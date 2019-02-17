import numpy as np
import os, time, glob
import utils
import matplotlib.pyplot as plt
from scipy.signal import convolve2d

def z_adjust(n, size, new):
	return ((n.astype(float) - 1.) * (float(size) / float(new)) - (float(size) / 2.))

class OpticalPerformance(object):

	def __init__(self, Defocus = 0, filename=None):
		# set the initial parameters of the system
		self.PARAMS = {'1': 400, # size of pupil aperture field in pixels (this defines the resolution of the calculation)
				       '5': .550, # imaging wavelength in microns
				       '6': 0, # number of pixels over which PSF is calculated
				       '7': 20 # increase to enhance the display of the wavefront (doesn't affect calculation)
				       }

		# --------Get the Zernike Coefficients -----------------------------------------------------------------------
		if filename is None:
			self.c = np.zeros(65, dtype=float)		
			self.PARAMS['2'] = 6. # size of pupil in mm for which PSF and MTF is to be calculated, 5mm default if no file chosen
			self.PARAMS['3'] = 6. # size of pupil in mm that Zernike coefficients define, 5mm default if no file chosen
		else:
			# read in the .zer file and extract information line by line
			fid = open(os.path.join('uploads', filename), 'r')
			version = fid.readline()
			instrument = fid.readline()
			manuf = fid.readline()
			oper = fid.readline()
			pupoff = fid.readline()
			geooff = fid.readline()
			datatype = fid.readline()
			Rfit = fid.readline()
			Rfit = float(Rfit[6:])
			Rmax = fid.readline()
			waverms = fid.readline()
			order = fid.readline()
			strehl = fid.readline()
			refent = fid.readline()
			refcor = fid.readline()
			resspec = fid.readline()
			data = fid.readline()
			# load in the data and reorder
			self.c = np.loadtxt(os.path.join('uploads', filename))[1:, 2] # ignore piston and extract just the coefficients

			# set the following parameters based on the .zer file 
			self.PARAMS['2'] = 2*Rfit # size of pupil in mm for which PSF and MTF is to be calculated
			self.PARAMS['3'] = 2*Rfit # size of pupil in mm that Zernike coefficients define
		
		self.PARAMS['4'] = 50; self.param4orig = self.PARAMS['4'] # automatically compute the field size

		# ----------- Calculate the RMS of the wave aberration --------------------------------------------------------

	def compute_intro(self, defocus=None, astigmatism=(None, None), coma=(None,None,None,None), high_order=True):

		# Remove the useless terms
		self.c[0] = 0; # this term is tilt and does not have any effect on image quality
		self.c[1] = 0; # this term is tilt and does not have any effect on image quality

		# remove or add other terms as desired
		if defocus is not None:
			self.c[3] =	defocus; # remove defocus
		if astigmatism[0] is not None:
			self.c[2] = astigmatism[0]
		if astigmatism[1] is not None:
			self.c[4] = astigmatism[1]
		if coma[0] is not None:
			self.c[5] = coma[0]
		if coma[1] is not None:
			self.c[6] = coma[1]
		if coma[2] is not None:
			self.c[7] = coma[2]
		if coma[3] is not None:
			self.c[8] = coma[3]
		if high_order is False:
			self.c[9:] = 0

		WaveAber_file = self.compWaveOSA()
		PSF_file = self.PSF()[0]

		return WaveAber_file, PSF_file

	def compute_defocus(self, start, stop, step, letter_size):
		letter_size = 5 * letter_size / 20.
		plotdimension = 60. * self.PARAMS['1'] * (180. * 60. / np.pi) * self.PARAMS['5'] *.001 / self.PARAMS['4']
		letter_size = self.PARAMS['1'] * letter_size * 60. / plotdimension 
		E = utils.make_E(letter_size, self.PARAMS['1'] ) # make an E of designated letter size

		self.c[:5] = 0 # set the correctable terms to zero
		defocus_range = np.arange(start, stop + step, step) # range of defoci in dioptres
		# fig, ax = plt.subplots(1, len(defocus_range), figsize=(6,2)) # intitialize figure

		if not os.path.isdir('static'):
			os.mkdir('static')
		else:
			# for filename in glob.glob(os.path.join('static', '*.png')):
			# 	os.remove(filename)
		
		fig, ax = plt.subplots(1,len(defocus_range))
		for i, defocus in enumerate(defocus_range):
			self.c[3] = utils.calculateDefocus(self.PARAMS['2'], defocus)
			PSF = self.PSF()[1]
			Conv = utils.convolve(E, PSF.T)
			ax[i].imshow(Conv, cmap='gray')
			ax[i].xaxis.set_visible(False)
			ax[i].yaxis.set_visible(False)
			ax[i].set_title(str(defocus) + " D")


		plotfile = os.path.join('static', 'DE_' + str(time.time()) + '.png')
		fig.savefig(plotfile)

		return plotfile




	def defocus_e(self):
		pass

	def defocus_streh_rms(self):
		pass

	def defocus_mtf(self):
		pass

	def PSF(self):
		pupilfunc = self.Zphase_MahajanOSA()

		Hamp = np.fft.fft2(pupilfunc)
		Hint = Hamp * np.conj(Hamp)


		plotdimension = 60 * self.PARAMS['1'] * (180 *60 / np.pi) * self.PARAMS['5'] *.001 / self.PARAMS['4'] # divided by the size of the pupil field.
		PSF = np.real(np.fft.fftshift(Hint)) # reorients the PSF so the origin is at the center of the image
		PSF = PSF / (self.PARAMS['6']**2) # scale the PSF so that peak represents the Strehl ratio

		axisPSF = np.arange(-plotdimension/2, ((plotdimension/2)-(plotdimension/self.PARAMS['1'])) + (plotdimension/self.PARAMS['1']), plotdimension/self.PARAMS['1'])

		# plot psf
		plt.figure(figsize=(8,4))
		plt.imshow(PSF, cmap='gray', extent=[-plotdimension/2, plotdimension/2, -plotdimension/2, plotdimension/2], aspect='auto')
		plt.title('Point Spread Function')
		plt.xlabel('arcsec')
		plt.ylabel('arcsec')
		plt.xlim([-plotdimension/2, plotdimension/2])
		plt.ylim([-plotdimension/2, plotdimension/2])

		# if not os.path.isdir('static'):
		# 	os.mkdir('static')
		# # else:
		# # 	for filename in glob.glob(os.path.join('static', '*.png')):
		# # 		os.remove(filename)
		plotfile = os.path.join('static', 'PSF_' + str(time.time()) + '.png')
		plt.savefig(plotfile)
		#plt.savefig("static/last/psf.png")
		return plotfile, PSF


	def compWaveOSA(self):
		waveabermap = self.Zwave_MahajanOSA()
		
		# waveabermap = waveabermap.T --> transpose for Austin's method
		n = waveabermap.shape

		# Define the axes of the wavefront plot in cyc/degree
		axispupil = np.arange(-self.PARAMS['2']/2, (self.PARAMS['2']/2)-(self.PARAMS['2']/n[0]) + self.PARAMS['2']/n[0], self.PARAMS['2']/n[0])

		# set parameters and display the wavefront aberration as a mesh plot
		v = np.arange(np.floor(np.min(waveabermap)) - 0.25, np.ceil(np.max(waveabermap)) + 0.25, 0.25)

		# plot wave aberration contour plot
		plt.figure(figsize=(8,4))
		plt.contourf(axispupil, axispupil, waveabermap)
		plt.colorbar()
		plt.title('Wavefront Aberration')
		plt.xlabel('mm (right-left)')
		plt.ylabel('mm (superior-inferior)')

		if not os.path.isdir('static'):
			os.mkdir('static')
		else:
			for filename in glob.glob(os.path.join('static', '*.png')):
				os.remove(filename)
		plotfile = os.path.join('static', 'WAM_' + str(time.time()) + '.png')
		plt.savefig(plotfile)
		#plt.savefig("static/last/wavemap.png")
		return plotfile

	def Zwave_MahajanOSA(self):
		newsize = np.ceil(self.PARAMS['7']*self.PARAMS['2']).astype(int)
		waveabermap = np.zeros((newsize, newsize))
		sizeoffield = self.PARAMS['2']

		self.PARAMS['6'] = 0

		#------ Joel's Method -------------------------------------------------------------
		nx = z_adjust(np.arange(1, newsize+1), sizeoffield, newsize)
		ny = z_adjust(np.arange(1, newsize+1), sizeoffield, newsize)

		xx, yy = np.meshgrid(nx, ny)

		angle, norm_radius = utils.cart2pol(xx,yy)
		norm_radius = norm_radius / (self.PARAMS['3'] / 2)
		r = norm_radius

		phase = utils.computePhase(self.c, angle, r)
		waveabermap[np.where(norm_radius > self.PARAMS['2']/self.PARAMS['3'])] = float('nan')
		waveabermap[np.where(norm_radius <= self.PARAMS['2']/self.PARAMS['3'])] = phase[np.where(norm_radius <= self.PARAMS['2']/self.PARAMS['3'])]


		# #---- Austin's Method --------------------------------------------------------------
		# for ny in range(newsize):
		# 	for nx in range(newsize):
		# 		xpos = ((nx - 1) * (sizeoffield / newsize) - (sizeoffield / 2))
		# 		print(xpos)
		# 		ypos = ((ny - 1) * (sizeoffield / newsize) - (sizeoffield / 2))

		# 		angle, norm_radius = utils.cart2pol(xpos,ypos)
		# 		norm_radius = norm_radius / (self.PARAMS['3'] / 2)
		# 		r = norm_radius

		# 		if norm_radius > self.PARAMS['2']/self.PARAMS['3']:
		# 			waveabermap[nx,ny] = float('nan')
		# 		else:
		# 			phase = utils.computePhase(self.c, angle, r)
		# 			waveabermap[nx,ny] = phase
		# 			self.PARAMS['6'] += 1

		return waveabermap

	def Zphase_MahajanOSA(self):
		phasemap = np.zeros((self.PARAMS['1'],self.PARAMS['1']))
		phasemap = phasemap + phasemap*1j 
		sizeoffield = self.PARAMS['4']

		B = 0.2 # fraction if the diffuse component;
		A = 1 - B # raction of directed component;
		peakx = 0 # normalized location of peak reflectance in x-direction
		peaky = 0 # normalized location of peak reflectance in x-direction
		p = 0 # 0.047 # set to zero to turn off the amplitude factor (Values from Burns et al in mm^(-2)

		self.PARAMS['6'] = 0

		# ----------Joel's Method --------------------------------
		nx = z_adjust(np.arange(1, self.PARAMS['1'] + 1), sizeoffield, self.PARAMS['1'])
		ny = z_adjust(np.arange(1, self.PARAMS['1'] + 1), sizeoffield, self.PARAMS['1'])
		xx, yy = np.meshgrid(nx, ny)

		angle, norm_radius = utils.cart2pol(xx,yy)
		norm_radius = norm_radius / (float(self.PARAMS['2']) / 2.)
		r = norm_radius

		phase = utils.computePhase(self.c, angle, r)
		d = np.sqrt((peakx - xx)**2 + (peaky - yy)**2);
		SCfactor = B + A * 10**(-p * ((self.PARAMS['2'])**2) * d**2)

		phase_result = SCfactor * np.exp(-1j * 2 * np.pi * phase / self.PARAMS['5'])

		phasemap[np.where(norm_radius <= 1)] = phase_result[np.where(norm_radius <= 1)]

		self.PARAMS['6'] = np.sum((norm_radius <= 1).astype(int))


		# #-------- Austin's Method ---------------------------------
		# newsize = self.PARAMS['1']
		# for ny in np.arange(1, self.PARAMS['1'] + 1):
		# 	for nx in np.arange(1, self.PARAMS['1'] + 1):
		# 		xpos = (((float(nx) - 1.) * (float(sizeoffield) / float(newsize))) - (sizeoffield / 2.))
		# 		ypos = (((float(ny) - 1.) * (float(sizeoffield) / float(newsize))) - (sizeoffield / 2.))

		# 		angle, norm_radius = utils.cart2pol(xpos,ypos)	
		# 		norm_radius = float(norm_radius) / (float(self.PARAMS['2']) / 2.)
		# 		r = norm_radius

				
		# 		if norm_radius > 1:
		# 			pass
		# 		else:
		# 			phase = utils.computePhase(self.c, angle, r)
		# 			d = np.sqrt((peakx - xpos)**2 + (peaky - ypos)**2);
		# 			SCfactor = B + A * 10**(-p * ((self.PARAMS['2'])**2) * d**2)
		# 			phasemap[nx,ny] = SCfactor * np.exp(-1j * 2 * np.pi * phase / self.PARAMS['5'])
		# 			self.PARAMS['6'] += 1

		return phasemap



if __name__ == '__main__':
	OP = OpticalPerformance(filename="Subject5_600.zer")
	# OP.compWaveOSA()
	# OP.PSF()
	OP.compute_defocus(-1,1,.5,20)
	plt.show()


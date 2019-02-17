import numpy as np
import matplotlib.pyplot as plt


def convolve(PSF, image):
	image = image + 0j
	PSF = PSF + 0j
	product = np.fft.fftshift(np.fft.fft2(image)) * np.fft.fftshift(np.fft.fft2(PSF))
	convolution = np.abs(np.fft.fftshift(np.fft.ifft2(product)))

	return convolution / np.max(convolution)

def make_E(letter_size, box_size):
	TumbE= np.ones((box_size,box_size))
	imcenter = np.array(TumbE.shape)/2;

	gapsize = round(letter_size/5.)

	TumbE[int(imcenter[0]-2*gapsize):int(imcenter[0]+3*gapsize-1), int(imcenter[1]-2*gapsize):int(imcenter[1]+3*gapsize-1)] = 0
	TumbE[int(imcenter[0]-1*gapsize):int(imcenter[0]-0*gapsize-1), int(imcenter[1]-1*gapsize):int(imcenter[1]+3*gapsize-1)] = 1
	TumbE[int(imcenter[0]+1*gapsize):int(imcenter[0]+2*gapsize-1), int(imcenter[1]-1*gapsize):int(imcenter[1]+3*gapsize-1)] = 1

	return TumbE


def calculateDefocus(pupilSize, Defocus):
	return (1e6 / (4 * np.sqrt(3))) * Defocus * ((pupilSize / 2000)**2) 

def RMS(coeff):
	return np.sqrt(np.sum(coeff**2))

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return theta, rho

def computePhase(c, angle, r):
	phase = 0
	phase = \
	c[0]*np.sqrt(4)*((1)*r**1)*np.sin(1*angle) + \
	c[1]*np.sqrt(4)*((1)*r**1)*np.cos(1*angle) + \
	c[2]*np.sqrt(6)*((1)*r**2)*np.sin(2*angle) + \
	c[3]*np.sqrt(3)*((2)*r**2+(-1)*r**0) + \
	c[4]*np.sqrt(6)*((1)*r**2)*np.cos(2*angle) + \
	c[5]*np.sqrt(8)*((1)*r**3)*np.sin(3*angle) + \
	c[6]*np.sqrt(8)*((3)*r**3+(-2)*r**1)*np.sin(1*angle) + \
	c[7]*np.sqrt(8)*((3)*r**3+(-2)*r**1)*np.cos(1*angle) + \
	c[8]*np.sqrt(8)*((1)*r**3)*np.cos(3*angle) + \
	c[9]*np.sqrt(10)*((1)*r**4)*np.sin(4*angle) + \
	c[10]*np.sqrt(10)*((4)*r**4+(-3)*r**2)*np.sin(2*angle) + \
	c[11]*np.sqrt(5)*((6)*r**4+(-6)*r**2+(1)*r**0) + \
	c[12]*np.sqrt(10)*((4)*r**4+(-3)*r**2)*np.cos(2*angle) + \
	c[13]*np.sqrt(10)*((1)*r**4)*np.cos(4*angle) + \
	c[14]*np.sqrt(12)*((1)*r**5)*np.sin(5*angle) + \
	c[15]*np.sqrt(12)*((5)*r**5+(-4)*r**3)*np.sin(3*angle) + \
	c[16]*np.sqrt(12)*((10)*r**5+(-12)*r**3+(3)*r**1)*np.sin(1*angle) + \
	c[17]*np.sqrt(12)*((10)*r**5+(-12)*r**3+(3)*r**1)*np.cos(1*angle) + \
	c[18]*np.sqrt(12)*((5)*r**5+(-4)*r**3)*np.cos(3*angle) + \
	c[19]*np.sqrt(12)*((1)*r**5)*np.cos(5*angle) + \
	c[20]*np.sqrt(14)*((1)*r**6)*np.sin(6*angle) + \
	c[21]*np.sqrt(14)*((6)*r**6+(-5)*r**4)*np.sin(4*angle) + \
	c[22]*np.sqrt(14)*((15)*r**6+(-20)*r**4+(6)*r**2)*np.sin(2*angle) + \
	c[23]*np.sqrt(7)*((20)*r**6+(-30)*r**4+(12)*r**2+(-1)*r**0) + \
	c[24]*np.sqrt(14)*((15)*r**6+(-20)*r**4+(6)*r**2)*np.cos(2*angle) + \
	c[25]*np.sqrt(14)*((6)*r**6+(-5)*r**4)*np.cos(4*angle) + \
	c[26]*np.sqrt(14)*((1)*r**6)*np.cos(6*angle) + \
	c[27]*np.sqrt(16)*((1)*r**7)*np.sin(7*angle) + \
	c[28]*np.sqrt(16)*((7)*r**7+(-6)*r**5)*np.sin(5*angle) + \
	c[29]*np.sqrt(16)*((21)*r**7+(-30)*r**5+(10)*r**3)*np.sin(3*angle) + \
	c[30]*np.sqrt(16)*((35)*r**7+(-60)*r**5+(30)*r**3+(-4)*r**1)*np.sin(1*angle) + \
	c[31]*np.sqrt(16)*((35)*r**7+(-60)*r**5+(30)*r**3+(-4)*r**1)*np.cos(1*angle) + \
	c[32]*np.sqrt(16)*((21)*r**7+(-30)*r**5+(10)*r**3)*np.cos(3*angle) + \
	c[33]*np.sqrt(16)*((7)*r**7+(-6)*r**5)*np.cos(5*angle) + \
	c[34]*np.sqrt(16)*((1)*r**7)*np.cos(7*angle) + \
	c[35]*np.sqrt(18)*((1)*r**8)*np.sin(8*angle) + \
	c[36]*np.sqrt(18)*((8)*r**8+(-7)*r**6)*np.sin(6*angle) + \
	c[37]*np.sqrt(18)*((28)*r**8+(-42)*r**6+(15)*r**4)*np.sin(4*angle) + \
	c[38]*np.sqrt(18)*((56)*r**8+(-105)*r**6+(60)*r**4+(-10)*r**2)*np.sin(2*angle) + \
	c[39]*np.sqrt(9)*((70)*r**8+(-140)*r**6+(90)*r**4+(-20)*r**2+(1)*r**0) + \
	c[40]*np.sqrt(18)*((56)*r**8+(-105)*r**6+(60)*r**4+(-10)*r**2)*np.cos(2*angle) + \
	c[41]*np.sqrt(18)*((28)*r**8+(-42)*r**6+(15)*r**4)*np.cos(4*angle) + \
	c[42]*np.sqrt(18)*((8)*r**8+(-7)*r**6)*np.cos(6*angle) + \
	c[43]*np.sqrt(18)*((1)*r**8)*np.cos(8*angle) + \
	c[44]*np.sqrt(20)*((1)*r**9)*np.sin(9*angle) + \
	c[45]*np.sqrt(20)*((9)*r**9+(-8)*r**7)*np.sin(7*angle) + \
	c[46]*np.sqrt(20)*((36)*r**9+(-56)*r**7+(21)*r**5)*np.sin(5*angle) + \
	c[47]*np.sqrt(20)*((84)*r**9+(-168)*r**7+(105)*r**5+(-20)*r**3)*np.sin(3*angle) + \
	c[48]*np.sqrt(20)*((126)*r**9+(-280)*r**7+(210)*r**5+(-60)*r**3+(5)*r**1)*np.sin(1*angle) + \
	c[49]*np.sqrt(20)*((126)*r**9+(-280)*r**7+(210)*r**5+(-60)*r**3+(5)*r**1)*np.cos(1*angle) + \
	c[50]*np.sqrt(20)*((84)*r**9+(-168)*r**7+(105)*r**5+(-20)*r**3)*np.cos(3*angle) + \
	c[51]*np.sqrt(20)*((36)*r**9+(-56)*r**7+(21)*r**5)*np.cos(5*angle) + \
	c[52]*np.sqrt(20)*((9)*r**9+(-8)*r**7)*np.cos(7*angle) + \
	c[53]*np.sqrt(20)*((1)*r**9)*np.cos(9*angle) + \
	c[54]*np.sqrt(22)*((1)*r**10)*np.sin(10*angle) + \
	c[55]*np.sqrt(22)*((10)*r**10+(-9)*r**8)*np.sin(8*angle) + \
	c[56]*np.sqrt(22)*((45)*r**10+(-72)*r**8+(28)*r**6)*np.sin(6*angle) + \
	c[57]*np.sqrt(22)*((120)*r**10+(-252)*r**8+(168)*r**6+(-35)*r**4)*np.sin(4*angle) + \
	c[58]*np.sqrt(22)*((210)*r**10+(-504)*r**8+(420)*r**6+(-140)*r**4+(15)*r**2)*np.sin(2*angle) + \
	c[59]*np.sqrt(11)*((252)*r**10+(-630)*r**8+(560)*r**6+(-210)*r**4+(30)*r**2+(-1)*r**0) + \
	c[60]*np.sqrt(22)*((210)*r**10+(-504)*r**8+(420)*r**6+(-140)*r**4+(15)*r**2)*np.cos(2*angle) + \
	c[61]*np.sqrt(22)*((120)*r**10+(-252)*r**8+(168)*r**6+(-35)*r**4)*np.cos(4*angle) + \
	c[62]*np.sqrt(22)*((45)*r**10+(-72)*r**8+(28)*r**6)*np.cos(6*angle) + \
	c[63]*np.sqrt(22)*((10)*r**10+(-9)*r**8)*np.cos(8*angle) + \
	c[64]*np.sqrt(22)*((1)*r**10)*np.cos(10*angle)

	return phase

import wtforms as wtf 

class Zernike(wtf.Form):
	filename = wtf.FileField(validators=[wtf.validators.InputRequired()])
	defocus = wtf.FloatField(label="Defocus", default=None)
	astigmatism_1 = wtf.FloatField(label="Astigmatism", default=None)
	astigmatism_2 = wtf.FloatField(label="Astigmatism", default = None) 
	coma_1 = wtf.FloatField(label="Coma", default=None)
	coma_2 = wtf.FloatField(label="Coma", default=None)
	coma_3 = wtf.FloatField(label="Coma", default=None)
	coma_4 = wtf.FloatField(label="Coma", default=None)
	higher_order = wtf.BooleanField(label="Higher-Order Aberrations", default=True)

class Defocus(wtf.Form):
	filename = wtf.FileField(validators=[wtf.validators.InputRequired()])
	start = wtf.FloatField(label="Start", default=-1)
	stop = wtf.FloatField(label="Stop", default=1)
	step = wtf.FloatField(label="Step", default=1)
	letter_size = wtf.FloatField(label="Letter Size", default=20)




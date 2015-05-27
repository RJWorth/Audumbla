###############################################################################
def mAnnulus(ri,ro,alpha,sigma0cgs):
	'''Calculate the mass of a dust disk with slope alpha, density sigma0.cgs
	in g/cm^2 at 1 AU, inner edge ri, and outer/truncation radius ro.'''
### Make sure the inputs are written as floats, not integers

### Two alpha value cases allowed
	assert ((alpha==1.5) | (alpha==1))
	assert (ro>= ri)

	import numpy as np
	from numpy import pi
	from cgs_constants import mEarth, AU

### Convert sigma to mSun/AU^2
	sigma0 = sigma0cgs*AU**2./mEarth
	if alpha==1:
		sigma0=0.303*sigma0	# normalization based on equal-mass 36-AU disks -- 5/4/15

### Calculate mTot in annulus with specified density profile
	if alpha==1.5:
		mTot = 4.*pi*sigma0* (1)**1.5 *(ro**0.5 - ri**0.5)
	elif alpha==1:
		mTot = 2.*pi*sigma0* (1)      *(ro      - ri)

	return mTot
###############################################################################


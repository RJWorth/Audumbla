###############################################################################

import random, copy
import DiskOOP as d
import numpy as np
from numpy import pi,sqrt
from cgs_constants import mSun,mEarth,AU,G
import multiprocessing as mp

###############################################################################
def AnnulusBody(disk, r0,r1):
	"""Return a body with the aggregate properties of a disk annulus."""

#		r0 = self.r_i + da*i	# inner edge of annulus
#		r1 = min(r0 + da,self.r_out)			# outer edge of annulus
	Annulus       = copy.deepcopy(disk)
	Annulus.r_i   = r0
	Annulus.r_out = r1
	m = Annulus.mDust()
	return d.Body(a=r0, m=m, M=disk.M)

###############################################################################
def mAnnulus(ri,ro,alpha,sigma0cgs):
	'''Calculate the mass in mSun of a dust disk with slope alpha, 
	density sigma0.cgs in g/cm^2 at 1 AU, inner edge ri, 
	and outer/truncation radius ro.'''
### Make sure the inputs are written as floats, not integers

### Two alpha value cases allowed
	assert ((alpha==1.5) | (alpha==1))
	assert (ro>= ri)

	import numpy as np
	from numpy import pi
	from cgs_constants import mSun, AU

### Convert sigma to mSun/AU^2
	sigma0 = sigma0cgs*AU**2./mSun
	if alpha==1:
		sigma0=0.303*sigma0	# normalization based on equal-mass 36-AU disks -- 5/4/15

### Calculate mTot in annulus with specified density profile
	if alpha==1.5:
		mTot = 4.*pi*sigma0* (1)**1.5 *(ro**0.5 - ri**0.5)
	elif alpha==1:
		mTot = 2.*pi*sigma0* (1)      *(ro      - ri)

	return mTot
###############################################################################
def rAnnulus(ri,m,alpha,sigma0cgs):
	'''Calculate the outer radius of a disk annulus with dust mass m, 
	slope alpha, density sigma0, and inner radius ri.'''
### Make sure the inputs are written as floats, not integers

### Two alpha value cases allowed
	assert ((alpha==1.5) | (alpha==1))

#	import numpy as np
#	from numpy import pi
#	from cgs_constants import mSun, AU

### Convert sigma to mSun/AU^2
	sigma0 = sigma0cgs*AU**2./mSun
	if alpha==1:
		sigma0=0.303*sigma0	# normalization based on equal-mass 36-AU disks -- 5/4/15

### Calculate mTot in annulus with specified density profile
	if alpha==1.5:
		ro = (m/(4.*pi*sigma0* (1)**1.5) + ri**0.5)**2.
	elif alpha==1:
		ro =  m/(2.*pi*sigma0* (1)     ) + ri

	return ro
###############################################################################
class Body(object):
	"""An orbiting body in a debris disk."""

#-----------------------------------------------------------------------------#
	def __init__(self, m, a, e=0, i=0, M=1.0):
		"""Instantialize object of Body class."""
			
		self.m = m # body's mass in mSun
		self.a = a # body's semimajor axis in AU
		self.e = e # body's eccentricity
		self.i = i # body's inclination
		self.M = M # mass of central object of body's orbit, in mSun
		self.mu = (self.m*self.M)/(self.m+self.M) # reduced mass

#-----------------------------------------------------------------------------#
	def eps(self):
		"""Specific orbital energy of Body around Central object"""
		eps = -self.mu/(2.*self.a)
		return eps

#-----------------------------------------------------------------------------#
	def P(self):
		"""Orbital period in days"""
		P = 2*pi * sqrt( (self.a*AU)**3. / (G*self.M*m) )/(24*3600.)
		return P

#-----------------------------------------------------------------------------#
#	def KE(self):
#		"""Kinetic energy of Body"""
#		KE = 0.
#		return KE		

#-----------------------------------------------------------------------------#
#	def PE(self):
#		"""Potential energy of Body"""
##		PE = -self.mu()/self.r
#		return PE		

#-----------------------------------------------------------------------------#
#	def E(self):
#		"""Total energy of Body"""
#		E = self.PE() + self.KE()
#		return E		

#-----------------------------------------------------------------------------#
	def RH(self):
		"""Hill radius of Body"""
		RH = self.a * (1-self.e) * (self.m/(3.*self.M))**(1./3.)
		return RH		

#-----------------------------------------------------------------------------#
	def RH2(self, other):
		"""Mutual Hill radius of two Bodies"""
  
		assert (self.M == other.M), 'Bodies do not have same central object'
		RH2 = ((self.m+other.m)/(3.*self.M))**(1./3.) * ((self.a+other.a)/2.)
		return RH2		
#-----------------------------------------------------------------------------#
	def IsClose(self,other,rh=10.):
		'''Check if two adjacent objects are within their mutual Hill radius
		of each other.'''


		da = abs(self.a - other.a)
		RH = self.RH2(other)

		if (rh*RH >= da):
			coll = True
		else:
			coll = False

		return coll


###############################################################################
#class DebrisDisk(list):
#	"""A list of Bodies in a debris disk."""

class DebrisDisk(list):

	def __init__(self, data, rh=10.):
		list.__init__(self, data)
		self.rh = rh
	
#	def times_two(self):
#		return self * 2 
#-----------------------------------------------------------------------------#
#	def __init__(self):
#		"""Instantialize object of DebrisDisk class."""
#		list.__init__(self)
		
#-----------------------------------------------------------------------------#
	def ListParams(self):
		'''Return list of masses/positions of debris objects.'''
		alist = np.array([i.a for i in self])
		mlist = np.array([i.m for i in self])

		return alist,mlist

##-----------------------------------------------------------------------------#
#	def IsClose(self,i1,out=True):
#		'''Check if two adjacent objects are within their mutual Hill radius
#		of each other.'''

#		if (out==True):
#			i2 = i1+1
#		else:
#			i2 = i1-1

#		da = abs(self[i1].a - self[i2].a)
#		RH = self[i1].RH2(self[i2])

#		if (RH >= da):
#			coll = True
#		else:
#			coll = False

#		return coll

#-----------------------------------------------------------------------------#
	def Merge(self,i,out=True):
		'''Merge debris object self.debris[i] with the outward or inward object'''

		if (out == True):
			i1 = i
			i2 = i+1
		else:
			i1 = i-1
			i2 = i

		deblist = copy.deepcopy(self)
		newobj  = copy.deepcopy(self[i1])
		oldobj1 = deblist[i1]
		oldobj2 = deblist[i2]

		newobj.m =  oldobj1.m + oldobj2.m
		newobj.a = np.average([oldobj1.a,oldobj2.a],
					weights = [oldobj1.m,oldobj2.m])

		deblist.remove(oldobj1)
		deblist.remove(oldobj2)
		deblist.insert(i1,newobj)
#		if (out == True):
#			deblist = d.DebrisDisk(deblist[:i]     + [newobj] + deblist[(i+2):])
#		else:
#			deblist = d.DebrisDisk(deblist[:(i-1)] + [newobj] + deblist[(i+1):])
		
		assert(len(deblist) == len(self)-1)

		return deblist

#-----------------------------------------------------------------------------#
	def PltFormRand(self):
		'''Collide random objects in the disk until they are all isolated, then
		output the resulting list of Bodies to Disk.planets.'''


		planets = [copy.deepcopy(self)]
		counter = 0
		N = len(planets[0])
		ind = range(N-1)

		while ((len(ind)>0) & (len(planets[-1])>1)):
			N = len(planets[-1])
			ind = range(N-1)
			merge = False
			while (merge == False):
				i = random.choice(ind)
#				print(counter, i,len(planets),len(planets[-1]))
				counter += 1
				this = planets[-1][i]
				that = planets[-1][i+1]
				if (this.IsClose(that, rh=self.rh)):
					planets.append(planets[-1].Merge(i))
#					ind = ind[0:i]+[j-1 for j in ind[(i+2):]]
					merge = True
				elif (i > 0):
					that = planets[-1][i-1]
					if (this.IsClose(that, rh=self.rh)):
						planets.append(planets[-1].Merge(i,out=False))
#						ind = ind[0:(i-1)]+[j-1 for j in ind[(i+1):]]
						merge = True
					else:
						ind.remove(i)
				else:
					ind.remove(i)
#				print(ind)

		return planets
#-----------------------------------------------------------------------------#
	def PltFormTilClear(self,i='default'):
		'''Pick one object, collide with nearby objects in the disk,
		alternating outward and inward, until it is isolated, then
		output the resulting list of Bodies.'''

		planets = copy.deepcopy(self)
		param_a = [planets.ListParams()[0]]
		param_m = [planets.ListParams()[1]]
		counter = 0
		N = len(planets[0])
		ind = range(N-1)
		if (i=='default'):
			i = random.choice(ind)
		out = False	# is planets[ind] isolated from next object outward?
		inw = False	# is planets[ind] isolated from next object inward?

		while ((out == False) & (inw == False)):
			print('counter={0} index={1} collisions={2} bodies remaining={3}'.format(
								counter, i,len(param_a),len(planets) ) )
			counter += 1
			if (i >= len(planets)-1):	
				out = True
			else:
				this = planets[i]
				that = planets[i+1]
				if (this.IsClose(that, rh=self.rh)):
					planets = planets.Merge(i)
					param_a.append(planets.ListParams()[0])
					param_m.append(planets.ListParams()[1])
				else:
					out = True
			if (i == 0):
				inw = True
			else:
				this = planets[i]
				that = planets[i-1]
				if (this.IsClose(that, rh=self.rh)):
					planets = planets.Merge(i,out=False)
					param_a.append(planets.ListParams()[0])
					param_m.append(planets.ListParams()[1])
					i += -1
				else:
					inw = True

		return planets, param_a, param_m
			
#-----------------------------------------------------------------------------#
	def PltFormSuccessiveClearings(self,i='default'):
		'''Pick one object, collide till clear, pick another and repeat
		until all are isolated. Then output the resulting list of Bodies.'''

		planets = copy.deepcopy(self)
		param_a = [planets.ListParams()[0]]
		param_m = [planets.ListParams()[1]]
		counter = 0
		N = len(planets)
		ind = range(N-1)

		while (len(ind)>0):
			if ((counter>0) | (i=='default')):
				i = random.choice(ind)
			BatchPlts, BatchPrm_a,BatchPrm_m = planets.PltFormTilClear(i)
			planets = BatchPlts
			param_a.extend(BatchPrm_a)
			param_m.extend(BatchPrm_m)

			ind  = range(len(planets)-1)
			ind2 = range(len(planets)-1)
			for i in ind2:
				this = planets[i]
				that = planets[i+1]
				if (not this.IsClose(that, rh=self.rh)):
					ind.remove(i)
			counter += 1

		return planets, param_a, param_m
		
###############################################################################
class Disk(object):
	"""The entire disk, with arguments holding disk properties."""

#-----------------------------------------------------------------------------#
	def __init__(self, unit='cgs', M=1., alpha=1.5, rh=10., 
		sigC=1.,sigma_g=1700., sigma_i=30., sigma_r=7.1,
		r_i=0.35, r_ice=2.7, r_out=36.,
		debris=[]):
		"""Instantialize object of Disk class."""
	
		assert ((alpha == 1.5) | (alpha==1.))
			
		self.unit    = unit		# over unit scheme used -- not using this yet
		self.M       = M       	# mass of central body (mSun)
		self.alpha   = alpha	# disk profile power law slope
		self.rh      = rh      	# separation w/in which bodies merge, in RH2
		self.sigC    = sigC    	# multiplier for surface density profile
		self.sigma_g = sigma_g	# gas surface density normalization (g/cm^2)
		self.sigma_i = sigma_i	# ice surface density normalization
		self.sigma_r = sigma_r	# rock surface density normalization
		self.sigma   = [sigma_r, sigma_i, sigma_g] # surface densities (rock, ice, gas)
		self.r_i     = r_i		# inner edge (AU)
		self.r_ice   = r_ice	# ice line
		self.r_out   = r_out	# outer edge
		self.r       = [r_i, r_ice, r_out] # radii of inner edge, ice line, outer edge
		self.debris  = d.DebrisDisk( [], rh) # list of objects in debris disk

#-----------------------------------------------------------------------------#
	def DebrisGenM(self,m,override=False):
		'''Populate the debris disk with Bodies with uniform mass m,
		based on the properties of the Disk.'''

		if override != True:
			assert(len(self.debris)==0),'Debris disk already exists!'
		else:
			print('Manually overwriting any existing debris disk')

	### Start at inner edge of disk, add objects incrementing outward until
	### ice line
	### with adjustments made crossing ice line and at outer edge
	### ri/ro = inner/outer radii of this annulus, m0 = mass of this annulus
		ri = self.r_i
		objlist = []
		while ((ri<=self.r_out) & (ri<=self.r_ice)):
			ro = d.rAnnulus(ri, m, self.alpha,self.sigC*self.sigma_r)
#			print(ri,ro)
			if (ro > self.r_ice):
				m1 = d.mAnnulus(ri, self.r_ice,  self.alpha, self.sigC*self.sigma_r)
				ro = d.rAnnulus(self.r_ice,m-m1,self.alpha, self.sigC*self.sigma_i)
				m0 = m
			if (ro > self.r_out):
				if (( ri < self.r_ice) & (self.r_ice < ro)):
					m1 = d.mAnnulus(ri,        self.r_ice,self.alpha,self.sigC*self.sigma_r)
					m2 = d.mAnnulus(self.r_ice,self.r_out,self.alpha,self.sigC*self.sigma_i)
					m0 = m1+m2
					assert ( (m0) <= m )
				else:
					m0 = d.mAnnulus(ri,self.r_out,self.alpha,self.sigC*self.sigma_r)
			else:
				m0 = m
			objlist.append(d.Body(a=ri, m = m0, M = self.M))
			ri = ro
	### Continue disk beyond ice line
		while ((ri<=self.r_out) & (ri>=self.r_ice)):
			ro = d.rAnnulus(ri, m, self.alpha,self.sigC*self.sigma_i)
#			print(ri,ro)
			if (ro <= self.r_out):
				m0 = m
			else:
				m0 = d.mAnnulus(ri,self.r_out,self.alpha,self.sigC*self.sigma_i)
			objlist.append(d.Body(a=ri, m = m0, M = self.M))
			ri = ro

#		print(len(objlist),m0,ri,ro)
		self.debris = d.DebrisDisk( objlist, self.rh )

#-----------------------------------------------------------------------------#
	def DebrisGenA(self,da,override=False):
		'''Populate the debris disk with Bodies with uniform spacing da,
		based on the properties of the Disk.'''

		if override != True:
			assert(len(self.debris)==0),'Debris disk already exists!'
		else:
			print('Manually overwriting any existing debris disk')

		expanse = self.r_out - self.r_i
		N = int(np.ceil(expanse/da))
#		self.debris    = d.DebrisDisk([d.Body(0,0) for i in range(N)])
#		for i in range(N):
#			self.debris[i] = self.AnnulusBody(self.r_i + da*i, 
#										  min(self.r_i + da*i + da,self.r_out)) 

		pool = mp.Pool(4)
		objlist = [pool.apply(d.AnnulusBody, args=(self, self.r_i + da*i,
						   min(self.r_i + da*i + da,self.r_out) ) )
														 for i in range(N)]
#		objlist = [d.AnnulusBody(self, self.r_i + da*i,
#						   min(self.r_i + da*i + da,self.r_out) )
#														 for i in range(N)]
		self.debris = d.DebrisDisk( objlist, self.rh )

#-----------------------------------------------------------------------------#
	def printargs(self):
		"""Print the arguments of Disk object."""
		print 'alpha  = {0}'.format(self.alpha)
		print 'sigma coefficient = {0}'.format(self.sigC)

#-----------------------------------------------------------------------------#
	def mGas(self):
		"""Gas mass of disk object."""

		s = self.sigC*self.sigma_g	
		mGas = d.mAnnulus(self.r_i,self.r_out,self.alpha,s)
		return mGas

#-----------------------------------------------------------------------------#
	def mDust(self):
		"""Dust mass of disk object."""

		sr = self.sigC*self.sigma_r	
		si = self.sigC*self.sigma_i	

	### Inside ice line
		if (self.r_i < self.r_ice):
			if (self.r_out >= self.r_ice):
				mDust1 = d.mAnnulus(self.r_i,self.r_ice,self.alpha,sr)
			elif (self.r_out < self.r_ice):
				mDust1 = d.mAnnulus(self.r_i,self.r_out,self.alpha,sr)
		else:
			mDust1 = 0.
#			print('no inner disk')

	### Outside ice line
		if (self.r_out > self.r_ice): 
			if (self.r_i <= self.r_ice):
				mDust2 = d.mAnnulus(self.r_ice,self.r_out,self.alpha,si)
			elif (self.r_i > self.r_ice):
				mDust2 = d.mAnnulus(self.r_i,self.r_out,self.alpha,si)
		else:
			mDust2 = 0.
#			print('no outer disk')
	
		return mDust1+mDust2


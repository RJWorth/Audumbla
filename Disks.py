###############################################################################

import random, copy
import Disks as d
import numpy as np
from numpy import pi,sqrt
from cgs_constants import mSun,mEarth,mMoon,mMars,AU,G
import multiprocessing as mp
tol = 1e-6	# tolerance in probability summing to 1.

###############################################################################
def Test():
	print('success')

###############################################################################
def AnnulusBody(disk, r0,r1):
	"""Return a body with the aggregate properties of a disk annulus."""

#		r0 = self.r_i + da*i	# inner edge of annulus
#		r1 = min(r0 + da,self.r_out)			# outer edge of annulus
	Annulus       = copy.deepcopy(disk)
	Annulus.r_i   = r0
	Annulus.r_out = r1
	m = Annulus.mDust
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
def PickFromList(x,prob):
	'''Pick an item from list x, where each item has the corresponding 
	probability in prob'''
	
	assert len(x)==len(prob)
	assert abs(1.-sum(prob)) < tol

	seed=random.random()
	for i in range(len(x)):
		if (seed <= sum(prob[0:(i+1)])):
			value = x[i]
			break
	assert value in x, 'invalid seed: {0}. Reduce tolerance?'.format(seed)

	return value

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
	def Test(self):
		print('success')

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
	def Test(self):
		print('success')

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
	def Merge(self,i,out=True, EjProb=0.,CorS=[0,0]):
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

		### Roll to see if objects merge or scatter
		if (random.random() >= EjProb):
			newobj.m =  oldobj1.m + oldobj2.m
			newobj.a = np.average([oldobj1.a,oldobj2.a],
						weights = [oldobj1.m,oldobj2.m])

			deblist.remove(oldobj1)
			deblist.remove(oldobj2)
			deblist.insert(i1,newobj)
			CorS[0] = CorS[0]+1
		else:							# if scattering, eject the smaller one
			if (oldobj1.m < oldobj2.m):
				deblist.remove(oldobj1)
			else:
				deblist.remove(oldobj2)	
			CorS[1] = CorS[1]+1

		assert(len(deblist) == len(self)-1)

		return deblist, CorS

#-----------------------------------------------------------------------------#
	def PltFormRand(self, EjProb=0.):
		'''Collide random objects in the disk until they are all isolated, then
		output the resulting list of Bodies to Disk.planets.'''


		planets = [copy.deepcopy(self)]
		counter = 0
		N = len(planets[0])
		ind = range(N-1)
		CorS = [0,0]

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
					deblist, CorS = planets[-1].Merge(i, EjProb=EjProb,CorS=CorS)
					planets.append(deblist)
					merge = True
					print(CorS)
				elif (i > 0):
					that = planets[-1][i-1]
					if (this.IsClose(that, rh=self.rh)):
						deblist, CorS = planets[-1].Merge(i,out=False, EjProb=EjProb,CorS=CorS)
						planets.append(deblist)
						print(CorS)
#						ind = ind[0:(i-1)]+[j-1 for j in ind[(i+1):]]
						merge = True
					else:
						ind.remove(i)
				else:
					ind.remove(i)
#				print(ind)

		return planets
#-----------------------------------------------------------------------------#
	def PltFormTilClear(self, i='default', EjProb=0.,CorS=[0,0]):
		'''Pick one object, collide with nearby objects in the disk,
		alternating outward and inward, until it is isolated, then
		output the resulting list of Bodies.'''

		planets = copy.deepcopy(self)
		param_a = [planets.ListParams()[0]]
		param_m = [planets.ListParams()[1]]
		counter = 0
		N = len(planets)
		ind = range(N-1)
		if (i=='default'):
			i = random.choice(ind)
		out = False	# is planets[ind] isolated from next object outward?
		inw = False	# is planets[ind] isolated from next object inward?

		while ((out == False) & (inw == False)):
			print('counter={0} index={1} coll={2} sctr={3} bodies remaining={4}'.format(
							counter, i,CorS[0],CorS[1],len(planets) ) )
			counter += 1
			if (i >= len(planets)-1):	
				out = True
			else:
				this = planets[i]
				that = planets[i+1]
				if (this.IsClose(that, rh=self.rh)):
					planets, CorS = planets.Merge(i, EjProb=EjProb,CorS=CorS)
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
#					planets = planets.Merge(i,out=False, EjProb=EjProb)
					planets, CorS = planets.Merge(i,out=False, EjProb=EjProb,CorS=CorS)
					param_a.append(planets.ListParams()[0])
					param_m.append(planets.ListParams()[1])
					i += -1
				else:
					inw = True

		return planets, param_a, param_m
			
#-----------------------------------------------------------------------------#
	def PltFormSuccessiveClearings(self,i='default', EjProb=0.):
		'''Pick one object, collide till clear, pick another and repeat
		until all are isolated. Then output the resulting list of Bodies.'''

		planets = copy.deepcopy(self)
		param_a = [planets.ListParams()[0]]
		param_m = [planets.ListParams()[1]]
		counter = 0
		N = len(planets)
		ind = range(N-1)
		CorS = [0,0]

		while (len(ind)>0):
			if ((counter>0) | (i=='default')):
				i = random.choice(ind)
			BatchPlts, BatchPrm_a,BatchPrm_m = planets.PltFormTilClear(i, EjProb=EjProb,CorS=CorS)
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

		return planets, param_a, param_m, CorS
		
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

		self.mGas = d.mAnnulus(self.r_i,self.r_out,self.alpha,self.sigC*self.sigma_g)
#-----------------------------------------------------------------------------#
### Calculate mDust
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
	### Outside ice line
		if (self.r_out > self.r_ice): 
			if (self.r_i <= self.r_ice):
				mDust2 = d.mAnnulus(self.r_ice,self.r_out,self.alpha,si)
			elif (self.r_i > self.r_ice):
				mDust2 = d.mAnnulus(self.r_i,self.r_out,self.alpha,si)
		else:
			mDust2 = 0.
	
		self.mDust = mDust1+mDust2
#-----------------------------------------------------------------------------#
	def Test(self):
		print('success')

#-----------------------------------------------------------------------------#
	def DebrisGenM(self,m=np.array([mMoon,mMars])/mSun,f=[.5,.5],override=False):
		'''Populate a debris disk based on the properties of the Disk 
		with Bodies of masses in m, with total mass fraction f in each.'''

		if override != True:
			assert(len(self.debris)==0),'Debris disk already exists!'
		else:
			print('Manually overwriting any existing debris disk')

		assert len(m) == len(f)
		assert abs(1.-sum(f)) <= tol

	### Calculate number of each size of planetesimal 
	#   and probability of a given planetesimal being each size,
	#   based on the total mass fraction that will be in each size
		N = [f[i]*self.mDust/m[i] for i in range(len(m))]
		prob = [n/sum(N) for n in N]

	### Start at inner edge of disk, add objects incrementing outward until
	#   ice line, with adjustments made crossing ice line and at outer edge
	### ri/ro = inner/outer radii of this annulus, m0 = mass of this annulus
		ri = self.r_i
		objlist = []
		while ((ri<=self.r_out) & (ri<=self.r_ice)):
			dm = d.PickFromList(m,prob)
			ro = d.rAnnulus(ri, dm, self.alpha,self.sigC*self.sigma_r)
			# if annulus crosses ice line, find the mass w/in ice line m1,
			# then annulus starting at ice line that contains remaining mass
			if (ro > self.r_ice):
				m1 = d.mAnnulus(ri, self.r_ice,  self.alpha, self.sigC*self.sigma_r)
				ro = d.rAnnulus(self.r_ice,dm-m1,self.alpha, self.sigC*self.sigma_i)
				m0 = dm
			# if outer radius is past edge of disk, if it also crosses ice line
			# integrate from ri to r_ice, then from r_ice to outer edge,
			# and this object has their combined mass. If not crossing ice line,
			# just integrate what's left of disk for this object
			if (ro > self.r_out):
				if (( ri < self.r_ice) & (self.r_ice < ro)):
					m1 = d.mAnnulus(ri,        self.r_ice,self.alpha,self.sigC*self.sigma_r)
					m2 = d.mAnnulus(self.r_ice,self.r_out,self.alpha,self.sigC*self.sigma_i)
					m0 = m1+m2
					assert ( (m0) <= dm )
				else:
					m0 = d.mAnnulus(ri,self.r_out,self.alpha,self.sigC*self.sigma_r)
			# if no complications, mass is as expected
			else:
				m0 = dm
			# add this object to the disk and begin next annulus at its outer edge
			objlist.append(d.Body(a=(ri+ro)/2., m = m0, M = self.M))
			ri = ro
	### Continue disk beyond ice line
		while ((ri<=self.r_out) & (ri>=self.r_ice)):
			dm = d.PickFromList(m,prob)
			ro = d.rAnnulus(ri, dm, self.alpha,self.sigC*self.sigma_i)
			if (ro <= self.r_out):
				m0 = dm
			else:
				m0 = d.mAnnulus(ri,self.r_out,self.alpha,self.sigC*self.sigma_i)
			objlist.append(d.Body(a=(ri+ro)/2., m = m0, M = self.M))
			ri = ro

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



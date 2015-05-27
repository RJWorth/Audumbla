###############################################################################

import random, copy
import DiskOOP as d
import numpy as np
from mDust import mAnnulus
from numpy import pi,sqrt
from cgs_constants import mSun,mEarth,AU,G

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
	def IsClose(self,other):
		'''Check if two adjacent objects are within their mutual Hill radius
		of each other.'''


		da = abs(self.a - other.a)
		RH = self.RH2(other)

		if (10.*RH >= da):
			coll = True
		else:
			coll = False

		return coll


###############################################################################
class DebrisDisk(list):
	"""A list of Bodies in a debris disk."""

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
				if (this.IsClose(that)):
					planets.append(planets[-1].Merge(i))
#					ind = ind[0:i]+[j-1 for j in ind[(i+2):]]
					merge = True
				elif (i > 0):
					that = planets[-1][i-1]
					if (this.IsClose(that)):
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

		planets = [copy.deepcopy(self)]
		counter = 0
		N = len(planets[0])
		ind = range(N-1)
		if (i=='default'):
			i = random.choice(ind)
		out = False
		inw = False

		while ((out == False) & (inw == False)):
			print('counter={0} index={1} collisions={2} bodies remaining={3}'.format(
								counter, i,len(planets),len(planets[-1]) ) )
			counter += 1
			if (i >= len(planets[-1])-1):	
				out = True
			else:
				this = planets[-1][i]
				that = planets[-1][i+1]
				if (this.IsClose(that)):
					planets.append(planets[-1].Merge(i))
				else:
					out = True
			if (i == 0):
				inw = True
			else:
				this = planets[-1][i]
				that = planets[-1][i-1]
				if (this.IsClose(that)):
					planets.append(planets[-1].Merge(i,out=False))
					i += -1
				else:
					inw = True

		return planets
			
#-----------------------------------------------------------------------------#
	def PltFormSuccessiveClearings(self,i='default'):
		'''Pick one object, collide till clear, pick another and repeat
		until all are isolated. Then output the resulting list of Bodies.'''

		planets = [copy.deepcopy(self)]
		counter = 0
		N = len(planets[0])
		ind = range(N-1)

		while (len(ind)>0):
			if ((counter>0) | (i=='default')):
				i = random.choice(ind)
			planets.extend(planets[-1].PltFormTilClear(i)[1:])

			ind  = range(len(planets[-1])-1)
			ind2 = range(len(planets[-1])-1)
			for i in ind2:
				this = planets[-1][i]
				that = planets[-1][i+1]
				if (not this.IsClose(that)):
					ind.remove(i)
			counter += 1

		return planets
		
###############################################################################
class Disk(object):
	"""The entire disk, with arguments holding disk properties."""

#-----------------------------------------------------------------------------#
	def __init__(self, unit='cgs', M=1., alpha=1.5, sig=1.,
		sigma_g=1700., sigma_i=30., sigma_r=7.1,
		r_i=0.35, r_ice=2.7, r_out=36.,
		debris=DebrisDisk()):
		"""Instantialize object of Disk class."""
	
		assert ((alpha == 1.5) | (alpha==1.))
			
		self.unit    = unit		# over unit scheme used -- not using this yet
		self.M       = M       	# mass of central body (mSun)
		self.alpha   = alpha	# disk profile power law slope
		self.sig     = sig     	# multiplier for surface density profile
		self.sigma_g = sigma_g	# gas surface density normalization (g/cm^2)
		self.sigma_i = sigma_i	# ice surface density normalization
		self.sigma_r = sigma_r	# rock surface density normalization
		self.r_i     = r_i		# inner edge (AU)
		self.r_ice   = r_ice	# ice line
		self.r_out   = r_out	# outer edge
		self.debris  = debris	# list of objects in debris disk

#-----------------------------------------------------------------------------#
	def DebrisGen(self,da,override=False):
		'''Populate the debris disk with Bodies with uniform spacing da,
		based on the properties of the Disk.'''

		if override != True:
			assert(len(self.debris)==0),'Debris disk already exists!'
		else:
			print('Manually overwriting any existing debris disk')

		expanse = self.r_out - self.r_i
		N = int(np.ceil(expanse/da))
#		self.debris = ['' for i in range(N)]
		self.debris = d.DebrisDisk()

		for i in range(N):
			r0 = self.r_i + da*i	# inner edge of annulus
			r1 = min(r0 + da,self.r_out)			# outer edge of annulus
			Annulus = copy.deepcopy(self)
			Annulus.r_i = r0
			Annulus.r_out = r1
			m = Annulus.mDust()*mEarth/mSun
			self.debris.append(Body(a=r0, m=m, M=self.M))

#-----------------------------------------------------------------------------#
	def printargs(self):
		"""Print the arguments of Disk object."""
		print 'alpha  = {0}'.format(self.alpha)
		print 'sigma coefficient = {0}'.format(self.sig)

#-----------------------------------------------------------------------------#
	def mGas(self):
		"""Gas mass of disk object."""

		s = self.sig*self.sigma_g	#*AU**2./mEarth
		mGas = mAnnulus(self.r_i,self.r_out,self.alpha,s)
		return mGas

#-----------------------------------------------------------------------------#
	def mDust(self):
		"""Dust mass of disk object."""

		sr = self.sig*self.sigma_r	#*AU**2./mEarth
		si = self.sig*self.sigma_i	#*AU**2./mEarth

	### Inside ice line
		if (self.r_i < self.r_ice):
			if (self.r_out >= self.r_ice):
				mDust1 = mAnnulus(self.r_i,self.r_ice,self.alpha,sr)
			elif (self.r_out < self.r_ice):
				mDust1 = mAnnulus(self.r_i,self.r_out,self.alpha,sr)
		else:
			mDust1 = 0.
#			print('no inner disk')

	### Outside ice line
		if (self.r_out > self.r_ice): 
			if (self.r_i <= self.r_ice):
				mDust2 = mAnnulus(self.r_ice,self.r_out,self.alpha,si)
			elif (self.r_i > self.r_ice):
				mDust2 = mAnnulus(self.r_i,self.r_out,self.alpha,si)
		else:
			mDust2 = 0.
#			print('no outer disk')
	
		return mDust1+mDust2


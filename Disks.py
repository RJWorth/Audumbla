###############################################################################

import random, copy
import Disks as d
import Merc as M
import numpy as np
from numpy import pi,sqrt
from cgs_constants import mSun,mEarth,mMoon,mMars,AU,G
import multiprocessing as mp
import matplotlib.pyplot as plt
import subprocess
tol = 1e-6	# tolerance in probability summing to 1.

###############################################################################
def Sep2D(planets):
	'''Takes a list of N DiskObjs. Returns an NxN array of the separation
	of each pair of objects, in mutual hill radii.'''

	N = len(planets)
	sep = np.zeros((N,N))

	for i in range(N):
		for j in range(N):
			if i != j:
				sep[i,j] = planets[i].sep(planets[j])
	return(sep)

###############################################################################
def A2D(planets):
	'''Takes a list of N DiskObjs. Returns an NxN array of the mass-weighted
	mean a of each pair of objects.'''

	N = len(planets)
	a2d = np.zeros((N,N))


	for i in range(N):
		for j in range(N):
			if i != j:
				a1, a2 = planets[i].a, planets[j].a
				m1, m2 = planets[j].m, planets[j].m
				a2d[i,j] = (a1*m1 + a2*m2)/(m1 + m2)
			else:
				a2d[i,j] = planets[i].a				

	return a2d

###############################################################################
def Prob2D_new(sep,a2):
	'''Takes a list of N DiskObjs. Returns an NxN array of the probabilities
	of each pair of objects colliding (unnormalized).'''

	assert sep.shape == a2.shape, "Prob2d error: Array shapes don't match"
	prob = np.zeros_like(a2)
	N1 = a2.shape[0]
	N2 = a2.shape[1]

	for i in range(N1):
		for j in range(N2):
			if i != j:
				prob[i,j] = 1./(sep[i,j]*a2[i,j])

	return(prob)

###############################################################################
def Prob2D(planets):
	'''Takes a list of N DiskObjs. Returns an NxN array of the probabilities
	of each pair of objects colliding (unnormalized).'''

	N = len(planets)
	prob = np.zeros((N,N))

	for i in range(N):
		for j in range(N):
			if i != j:
				prob[i,j] = planets[i].prob(planets[j])

	return(prob)

###############################################################################
def Roll(prob):
	'''Takes an NxN array of probabilities and a random float 
	between 0 and 1. Walks through the array in rows, subtracting
	each probability, until crossing zero, then returns those indices.'''

	x = random.uniform(0.,np.sum(prob))

	Ni,Nj = prob.shape
	assert Ni==Nj, 'array not square!'
	N = Ni

	# Step through array until reaching the random value
	for ind in range(N**2):
		i = (ind-ind%N)/N
		j = ind%N
		x = x-prob[ i, j]
		if x < 0.:
			return 	(ind-ind%N)/N, ind%N 
	raise ValueError, "Error in Roll: no match"

###############################################################################
def AnnulusBody(disk, r0,r1):
	"""Return a body with the aggregate properties of a disk annulus."""

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

############################################################################
def Comparison(WhichDir):
	'''Organize data from MERCURY run into list of interactions to 
	compare with my algorithm's output.'''

	Mcent = M.ReadMCent(WhichDir+'/')
	CorS = [0,0]

	### Get initial objects
	objs = M.ReadInObjList(WhichDir+'/In/','big.in')
	namelist = [objs[0].name]
	a_thistime, m_thistime = [objs[0].a()], [objs[0].mass]
	for o in objs[1:]:
		namelist.append(o.name)
		a_thistime.append(o.a())
		e_thistime.append(o.a())
		i_thistime.append(o.mass)
		m_thistime.append(o.mass)
#		for n in namelist:
#			a, m, tf = M.ReadAeiLine(WhichDir,n,tf)
	a_list = [a_thistime]
	e_list = [a_thistime]
	i_list = [m_thistime]
	m_list = [m_thistime]
	da = np.array([0.])
	Rh = np.array([0.])
	acoll = np.array([0.])
	mnsep = [np.mean([a_thistime[i+1]-a_thistime[i] for i in range(len(a_thistime)-1)])]
	mdsep = [np.median([a_thistime[i+1]-a_thistime[i] for i in range(len(a_thistime)-1)])]

	### Get interaction lists
	name, dest, time = M.ReadInfo(WhichDir)
#	namelist.remove(name[0])	# remove first object

	### Step through collisions in chronological order
#	for i, t in enumerate(time[1:]):
	for i, t in enumerate(time):
		a_other = 'reset'
		if dest[i] in ['ejected','Center']:
			CorS[1] =+ 1
			if dest[i] == 'ejected':
				a_other = 100.
			else:
				a_other = 0.
		else:
			CorS[0] =+ 1
	### Get last timestep of destroyed object
		a, m, tf = M.ReadAeiLine(WhichDir,name[i],t,iscoll=True)
		a_self, e_self, i_self, m_self = a, e, i, m
		a_thistime, m_thistime = [], []
		namelist.remove(name[i])
	### Get all other objects at the next timestep and record their positions
		for n in namelist:
			a, m, dummy = M.ReadAeiLine(WhichDir,n,tf)
			if not a  == None:
				a_thistime.append(a)
				m_thistime.append(m)
			if n == dest[i]:
				a_other = a
				m_other = m
		assert a_other != 'reset', 'a_other not matched!'
	### Put a, m into lists like IceCow
		a_thistime.sort
		m_thistime.sort
		a_list.append(a_thistime)
		m_list.append(m_thistime)
		acoll = np.append(acoll,np.mean([a_self,a_other]))
		da = np.append(da,a_self-a_other)
		rh2 = ((m_self+m_other)/(3.*Mcent))**(1./3.) * ((a_self+a_other)/2.)
		Rh = np.append(Rh, rh2)
		mnsep.append(np.mean([a_thistime[i+1]-a_thistime[i] for i in range(len(a_thistime)-1)]))
		mdsep.append(np.median([a_thistime[i+1]-a_thistime[i] for i in range(len(a_thistime)-1)]))
	acoll = np.append(acoll, 0.)
	da = np.append(da, 0.)
	Rh = np.append(Rh, 0.)

	print(namelist)
	### Get final masses and positions
	line = subprocess.check_output(['tail', '-1', WhichDir+'/Aei/'+namelist[0]+'.aei'])
	tMax = float(line.split()[0])
	print(tMax)
	a_thistime, m_thistime = [], []
	for n in namelist:
		a, m, tf = M.ReadAeiLine(WhichDir,n,tMax)
		a_thistime.append(a)
		m_thistime.append(m)
	a_thistime.sort
	m_thistime.sort
	a_list.append(a_thistime)
	m_list.append(m_thistime)
	mnsep.append(np.mean([a_thistime[i+1]-a_thistime[i] for i in range(len(a_thistime)-1)]))
	mdsep.append(np.median([a_thistime[i+1]-a_thistime[i] for i in range(len(a_thistime)-1)]))
		
	# insert 0 and tMax at beginning and end of timeline
	tlist = np.insert(np.append(time,tMax),0,0.)

	return np.array(a_list), np.array(m_list), tlist, acoll, da, Rh, CorS, mnsep, mdsep

###############################################################################
def PlotDiskMerge(a, m, tlist='default',
		acoll='default',da='default',Rh='default',
		mnsep='default',mdsep='default',
		fname='IceCow_000',ftype='png',ylim='default'):
	'''Plot the merge history of a list of semimajor axes at each step'''

	plt.close('all')

	### Set up plot objects
	if tlist == 'default':
		f, ax1 = plt.subplots(1)
	else:
		f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

	# Make x values array and plot each timestep
	x = np.array(len(a))
	for i in range(len(a)):
		if tlist == 'default':
			x[i] = [   i    for n in range(len(a[i]))]
			xlab = '# of Collisions'
		else:
			x[i] = [tlist[i] for n in range(len(a[i]))]
			xlab = 'Time'
#		ax1.scatter(x, a[i], s=1.,c=np.log(m[i]), lw=0.)
	ax1.scatter(x, a, s=1.,c=np.log(m[i]), lw=0.)

	(x1,x2),(y1,y2) = ax1.get_xlim(), ax1.get_ylim()
#	x1=0.1
	if not ylim == 'default':
		ax1.set_ylim(ylim)
		ax1.set_xscale('log')
	ax1.set_xlabel(xlab)
	ax1.set_ylabel('Semimajor Axis (AU)')
#	plt.title('Planetesimal Merging Heuristic')

	### If time provided, make t vs. n plot
	if tlist != 'default':
		dt = np.array([tlist[i+1]-tlist[i] for i in range(len(tlist)-1)])
		dRh = np.array(da)/np.array(Rh)
#		clip_dt = np.array([(i < 5e6) for i in dt])
#		clip_da = np.array([(i < 90) for i in da])
#		clip = np.array([ clip_dt[i] & clip_da[i] for i in range(len(dt))])
		da_sum = [np.sum(da[:i]) for i in range(len(da))]

		da2, dt2 = da[1:(len(da)-1)], dt[:(len(dt)-1)]
		ac2, Rh2, dRh2 = acoll[1:(len(a)-1)], Rh[1:(len(Rh)-1)], dRh[1:(len(dRh)-1)]

		test1,test2,test3 = np.ones(len(ac2)),np.ones(len(ac2)),np.ones(len(ac2))

		ax2.set_xscale('log')
		ax2.set_yscale('log')
		ax2.scatter(da2,dt2, c=np.log(ac2))
		ax2.set_xlabel('da')
		ax2.set_ylabel('dt')

		ax3.set_xscale('log')
		ax3.set_yscale('log')
		ax3.scatter(Rh2,dt2, c=np.log(ac2))
		ax3.set_xlabel('Rh')
		ax3.set_ylabel('dt')

		ax4.set_xscale('log')
		ax4.set_yscale('log')
		ax4.scatter(dRh2,dt2, c=np.log(ac2))
		ax4.set_xlabel('dRh')
		ax4.set_ylabel('dt')

#		plt.subplot(222)
#		plt.xscale('log')
#		plt.yscale('log')
##		plt.plot(da2,dt2,'ro')
#		plt.scatter(da2,dt2, c=test3)
#		plt.xlabel('da')
#		plt.ylabel('dt')

#		plt.subplot(223)
#		plt.plot(da2,dt2,'ro')
##		plt.scatter(Rh[1:(len(Rh)-1)],dt[:(len(dt)-1)], c=a[1:(len(da)-1)])
#		plt.xscale('log')
#		plt.yscale('log')
#		plt.xlabel('Rh')
#		plt.ylabel('dt')

#		plt.subplot(224)
#		plt.scatter(dRh[1:(len(dRh)-1)],dt[:(len(dt)-1)], c=a[1:(len(da)-1)])
#		plt.xscale('log')
#		plt.yscale('log')
#		plt.xlabel('dRh')
#		plt.ylabel('dt')

#		plt.subplot(222)
#		plt.plot(tlist, range(len(a)), 'bo')
#		plt.xscale('log')
#		plt.yscale('log')
#		plt.xlabel('# of Collisions')
#		plt.ylabel('Time')
#		plt.title('')
#		
#		plt.subplot(223)
#		plt.plot(da2[clip],dt[clip],'ro')
#		plt.xscale('log')
#		plt.yscale('log')
#		plt.xlabel('Object Separation')
#		plt.ylabel('Time to Collision')
#		plt.title('')

#		plt.subplot(224)
#		plt.plot(da_sum,tlist)
#		plt.xscale('log')
#		plt.yscale('log')
#		plt.xlabel('Cumulative Separations')
#		plt.ylabel('Time')
#		plt.title('')

#		plt.subplot(222)
#		plt.plot(tlist, mdsep, 'bo',
#				 tlist, mnsep, 'ro')
#		plt.xscale('log')
#		plt.yscale('log')
#		plt.xlabel('Time')
#		plt.ylabel('Average Separation')
#		plt.title('')

#		plt.subplot(223)
#		plt.plot(dt, mnsep[1:], 'ro')
#		plt.xscale('log')
#		plt.yscale('log')
#		plt.xlabel('dt')
#		plt.ylabel('Mean Separation')
#		plt.title('')


#		plt.subplot(224)
#		plt.plot(dt, mdsep[1:], 'bo')
#		plt.xscale('log')
#		plt.yscale('log')
#		plt.xlabel('dt')
#		plt.ylabel('Median Separation')
#		plt.title('')


	### Save
	plt.savefig(fname+'.'+ftype)
	plt.clf()
	plt.close(f)

###############################################################################
###############################################################################
#------------------------------- Objects --------------------------------------
###############################################################################
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
		P = 2*pi * sqrt( (self.a*AU)**3. / (G*(self.M+self.m)*mSun) )/(24*3600.)
		return P

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
	def sep(self, other):
		"""Separation of two Bodies, in mutual Hill radii"""

		da = self.a-other.a
		rh2 = ((self.m+other.m)/(3.*self.M))**(1./3.) * ((self.a+other.a)/2.)
		if da>0.:
			sep =  da/rh2
		else:
			sep = -da/rh2

		return sep

#-----------------------------------------------------------------------------#
	def prob(self, other):
		"""Probability of collision of two Bodies, defined as inverse
		of their separation in mutual Hill radii"""

		da = self.a-other.a
		rh2 = ((self.m+other.m)/(3.*self.M))**(1./3.) * ((self.a+other.a)/2.)
		if da>0.:
			prob =  rh2/da
		else:
			prob = -rh2/da

		return prob

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
		self.rh = rh	# # of hill radii at which objects won't collide
	
#-----------------------------------------------------------------------------#
	def ListParams(self):
		'''Return list of masses/positions of debris objects.'''
		alist = np.array([i.a for i in self])
		mlist = np.array([i.m for i in self])

		return alist,mlist

#-----------------------------------------------------------------------------#
	def SortObjs(self):
		'''Return list reordered from lowest to highest a.'''

		alist = self.ListParams()[0]
		order = np.argsort(alist)

		newlist = d.DebrisDisk([self[i] for i in order])

		return(newlist)

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
	def MergeAny(self, i1, i2, EjProb=0.,CorS=[0,0]):
		'''Merge debris object self.debris[i] with any other object'''

		deblist = copy.deepcopy(self)
		newobj  = copy.deepcopy(self[i1])
		oldobj1 = deblist[i1]
		oldobj2 = deblist[i2]

		### Roll to see if objects merge or scatter
		if (random.random() >= EjProb):	# if merging
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

	### Put the objects in order by a
		deblist = deblist.SortObjs()

		return deblist, CorS

#-----------------------------------------------------------------------------#
	def IceCow(self, EjProb=0.,CorS=[0,0], nRh=10.,drift=[0.,0.],vers='new'):
		'''Collide objects in the disk based on probability for each pair, 
		based on how close they are. Continue until all are isolated, then
		output the resulting list of Bodies to Disk.planets. Return lists
		of each a and e over the evolution, and list with counts of 
		collisions vs. scatters.'''

		planets = copy.deepcopy(self)
		param_a = [planets.ListParams()[0]]
		param_m = [planets.ListParams()[1]]
		counter = 0
		N = len(planets)
		ind = range(N-1)
		CorS = [0,0]
		t = [0.]

		# 2D array of probabilities ( = RH2/da)
		if vers=='old':
			prob = d.Prob2D(planets)
			minsep = 1./np.max(prob)
		elif vers=='new':
			sep2d = d.Sep2D(planets)
			a2d   = d.A2D(planets)
			prob  = d.Prob2D_new(sep2d, a2d)
			flat  = np.sort(sep2d.flatten())
			nzeros = len(flat[flat == 0.])
			print(flat[(nzeros-1):(nzeros+4)])
			minsep = flat[nzeros]

		print(minsep)
		# merge closest objects until all are spaced by nRh*RH2
		while (minsep <= nRh):
			print(t)
			i,j = d.Roll(prob)
			# merge i and j
			planets, CorS = planets.MergeAny(i,j, EjProb=EjProb,CorS=CorS)
			# add random drift in positions
			if drift != [0.,0.]:
				for k in range(len(planets)):
					planets[k].a = planets[k].a + random.uniform(
													-drift[0]/planets[k].m,
													 drift[1]/planets[k].m)
			# add a/e to record
			param_a.append(planets.ListParams()[0])
			param_m.append(planets.ListParams()[1])
			# advance timestep by minimum separation
			t.append(t[-1] + minsep)
			# recalculateprob
			print('{0: 4} t={1: .1e},{5: 3} remaining, ij=[{2: 4} {3: 4}], P={4: .1e}, minsep={6: .1e}'.format(
					counter, t[-1], i,j, prob[i,j], len(planets), minsep ) )
			if vers=='old':
				prob = d.Prob2D(planets)
				minsep = 1./np.max(prob)
			elif vers=='new':
				sep2d = d.Sep2D(planets)
				a2d   = d.A2D(planets)
				preva2d = a2d
				prob  = d.Prob2D_new(sep2d, a2d)
				flat  = np.sort(sep2d.flatten())
				nzeros = len(flat[flat == 0.])
				print(flat[(nzeros-1):(nzeros+4)])
				minsep = flat[nzeros]

				index = min(i,j)
				print(np.array(preva2d[:index,:index]) == np.array(a2d[:index,:index]))
				print(np.array(preva2d[:index,:index]))
				print(np.array(a2d[:index,:index]))
			counter += 1

		return planets, param_a, param_m, t, CorS
	

#-----------------------------------------------------------------------------#
	def PltFormTimed(self, EjProb=0.,CorS=[0,0], nRh=10.):
		'''Collide objects in the disk based on a 'time' parameter 
		based on how close they are. Continue until all are isolated, then
		output the resulting list of Bodies to Disk.planets.'''

		planets = copy.deepcopy(self)
		param_a = [planets.ListParams()[0]]
		param_m = [planets.ListParams()[1]]
		counter = 0
		N = len(planets)
		ind = range(N-1)
		CorS = [0,0]

		# 'merge time' (proportional to distance in rHs)
#		rh = np.array([planets[j].RH2(planets[j+1]) for j in range(len(planets)-1)])
#		da = np.array([abs(planets[j].a-planets[j+1].a) for j in range(len(planets)-1)])
		t = np.array([planets[j].sep(planets[j+1]) for j in range(len(planets)-1)])
		# merge closest objects until all are spaced by nRh*RH2
		now=0.
		while (min(t)-now <= nRh):
			now = min(t)
			i = np.where(t==min(t))[0][0]
			# merge i and i+1
			planets, CorS = planets.Merge(i, EjProb=EjProb,CorS=CorS)
			param_a.append(planets.ListParams()[0])
			param_m.append(planets.ListParams()[1])
			print('counter={0: 4} index={1: 4} coll={2: 4} sctr={3: 4} bodies remaining={4: 3}'.format(
							counter, i,CorS[0],CorS[1],len(planets) ) )
			# recalculate t
#			rh = np.array([planets[j].RH2(planets[j+1]) for j in range(len(planets)-1)])
#			da = np.array([abs(planets[j].a-planets[j+1].a) for j in range(len(planets)-1)])
###			t = [planets[j].sep(planets[j+1]) for j in range(len(planets)-1)]
			mask = np.ones(len(t), dtype=bool)
			if (i < len(t)-1):
				mask[i+1] = False
				t2 = t[mask]
				t2[i-1] = now + planets[i-1].sep(planets[i  ])
				t2[i]   = now + planets[i  ].sep(planets[i+1])
			else:
				mask[i] = False
				t2 = t[mask]
				t2[i-1] = now + planets[i-1].sep(planets[i  ])
#				t2[i]   = now + planets[i  ].sep(planets[i+1])
	
			t = t2
			print(i,len(t))

			counter += 1

		return planets, param_a, param_m, CorS
	

#-----------------------------------------------------------------------------#
	def PltFormRand(self, EjProb=0.):
		'''Collide random objects in the disk until they are all isolated, then
		output the resulting list of Bodies to Disk.planets.'''

		counter = 0
		planets = [copy.deepcopy(self)]
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


###############################################################################


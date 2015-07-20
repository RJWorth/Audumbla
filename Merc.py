from cgs_constants import mSun,mEarth,AU
from cgs_constants import G as G_cgs
from mks_constants import deg2rad,rad2deg
from cgs_constants import G as G_mks
import Disks as D
import Merc as M
import random as R
import numpy as np
from numpy import pi, sqrt, sin, cos

###############################################################################
class Obj(object):
	'''Object in a mercury simulation, input in asteroidal coords. Cartesian 
	coords are calculated based on them.'''

#-----------------------------------------------------------------------------#
	def __init__(self, units='Asteroidal', name='Earth', 
		mass=mEarth/mSun, mCent = mSun, density=5.51, 
		a=1., e=0., i=0., g=0., n=0., m=0., s=[0., 0., 0.]):

		self.name    = name
		self.mass    = mass
		self.mCent   = mCent
		self.density = density
		self.a       = a
		self.e       = e
		self.i       = i
		self.g       = g
		self.n       = n
		self.m       = m
		self.s       = s

		self.x, self.y, self.z, self.vx, self.vy, self.vz = M.Merc_El2X(
			[a,e,i, g*pi/180.,n*pi/180.,m*pi/180.], [mCent, mass])
		self.pos = [self.x,  self.y,  self.z]
		self.vel = [self.vx, self.vy, self.vz]

#-----------------------------------------------------------------------------#
	def RecalcCartesian(self):
		'''Run this after changing any Asteroidal parameters to make the 
		Cartesian ones consistent.'''

		self.x, self.y, self.z, self.vx, self.vy, self.vz = M.Merc_El2X(
			[self.a,self.e,self.i, 
				self.g*pi/180.,self.n*pi/180.,self.m*pi/180.],
			[self.mCent, self.mass])
		self.pos = [self.x,  self.y,  self.z]
		self.vel = [self.vx, self.vy, self.vz]

###############################################################################
def GetObjList(rtr = 5., sigC = 10., rh = 10., dm = 1.e-7, alpha = 1.5):
	'''Generate list of objects based on the Disks semi-analytic model'''

	disk = D.Disk(r_out=rtr, alpha=alpha, sigC=sigC, rh=rh)
	disk.DebrisGenM(dm)
	a, mass = disk.debris.ListParams()

#	digits = str(int(np.floor(np.log10(len(a))+1)))
	digits = str(len(str(len(a))))
	fmt = 'P{0:0'+digits+'}'

	objlist = []
	for i in range(len(a)):
		objlist.append( M.Obj(name=fmt.format(i), mass=mass[i], density=3.,
		a=a[i], 
		g=R.uniform(0.,360.), n=R.uniform(0.,360.), m=R.uniform(0.,360.)) )

	return objlist

###############################################################################
def WriteObjInFile(objlist='default', loc = 'Merc95/In/',infile='big', epoch=0.):
	'''Write a big.in or small.in file for mercury'''

### Currently non-variable variable
	style = 'Asteroidal'

### Make list of object parameters, unless provided
	if (objlist=='default'):
		objlist = M.GetObjList()

### Process big/small differences
	assert ((infile == 'big') | (infile=='small')), 'invalid infile: must be "big" or "small"'
	fname = loc+infile+'.in'
	if (infile == 'big'):
		vers = 0
	elif (infile == 'small'):
		vers = 1

	header = ')O+_06 '+['Big','Small'][vers]+ \
"-body initial data  (WARNING: Do not delete this line!!)\n) Lines beginning with `)' are ignored.\n)---------------------------------------------------------------------\n style (Cartesian, Asteroidal, Cometary) = {style}"\
+["\n epoch (in days) = {epoch}",""][vers]+ \
"\n)---------------------------------------------------------------------\n"

	objstr = '''  {0.name:16}  m={0.mass}  d={0.density}
    {0.a: .18e} {0.e: .18e} {0.i: .18e}
    {0.g: .18e} {0.n: .18e} {0.m: .18e}
    {0.s[0]: 19.18e} {0.s[1]: 19.18e} {0.s[2]: 19.18e}
'''

	with open(fname, 'w') as f:
		f.write(header.format(style=style,epoch=epoch))
		for i in range(len(objlist)):
			f.write(objstr.format(objlist[i]))

###############################################################################
def WriteParamInFile(loc = 'Merc95/In/', f = 'in', alg='hybrid', 
	ti=0., tf=365.25e3, tOut = 'default', dt=1., acc=1.e-12, 
	CEstop='no',CE='yes',CEfrag='no',tUnit='years',tRel='yes',prec='medium',
	rel='no',user='no',
	rEj=100, rStar=0.005, mStar=1.0, NrH=3.,dtDump='default',dtPer=100):
	'''Write a param.in file for mercury'''
	
	if (tOut == 'default'):
		tOut = tf/1.e3
	if (dtDump == 'default'):
		dtDump = int(tf/dt/100.)

	text = ''')O+_06 Integration parameters  (WARNING: Do not delete this line!!)
) Lines beginning with `)' are ignored.
)---------------------------------------------------------------------
) Important integration parameters:
)---------------------------------------------------------------------
 algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = {alg}
 start time (days)= {ti}
 stop time (days) = {tf}
 output interval (days) = {tOut}
 timestep (days) = {dt}
 accuracy parameter={acc}
)---------------------------------------------------------------------
) Integration options:
)---------------------------------------------------------------------
 stop integration after a close encounter = {CEstop}
 allow collisions to occur = {CE}
 include collisional fragmentation = {CEfrag}
 express time in days or years = {tUnit}
 express time relative to integration start time = {tRel}
 output precision = {prec}
 < not used at present >
 include relativity in integration= {rel}
 include user-defined force = {user}
)---------------------------------------------------------------------
) These parameters do not need to be adjusted often:
)---------------------------------------------------------------------
 ejection distance (AU)= {rEj}
 radius of central body (AU) = {rStar}
 central mass (solar) = {mStar}
 central J2 = 0
 central J4 = 0
 central J6 = 0
 < not used at present >
 < not used at present >
 Hybrid integrator changeover (Hill radii) = {NrH}
 number of timesteps between data dumps = {dtDump}
 number of timesteps between periodic effects = {dtPer}'''.format(
	alg=alg, ti=ti, tf=tf, tOut=tOut, dt=dt, acc=acc,
	CEstop=CEstop, CE=CE, CEfrag=CEfrag, tUnit=tUnit, tRel=tRel,prec=prec,
	rel=rel, user=user,
	rEj=rEj, rStar=rStar,mStar=mStar, NrH=NrH, dtDump=dtDump, dtPer=dtPer)

	fname=loc+'param.'+f
	with open(fname, 'w') as f:
		f.write(text)


###########################################################################
### Convert from cartesian (xyz uvw) to orbital elements (aei gnM)
#   BUT WHAT FRAME ARE THEY IN???
def Merc_El2X(el, mass):
	'''Convert orbital elements to cartesian for an ellipse (e < 1). 
Based on MCO_EL2X.FOR from J. Chambers' mercury6_2.for:
	gm = grav const * (central + secondary mass)
	q = perihelion distance
	e = eccentricity
	i = inclination                 )
	p = longitude of perihelion !!! )   in
	n = longitude of ascending node ) radians
	l = mean anomaly                )

	x,y,z = Cartesian positions  ( units the same as a )
	u,v,w =     "     velocities ( units the same as sqrt(gm/a) )'''
# Big.in format:
# Asteroidal = Keplerian orbital elements, in an `asteroidal' format.
#              i.e.  a e I g n M, where
#               a = semi-major axis (in AU)
#               e = eccentricity
#               I = inclination (degrees)
#               g = argument of pericentre (degrees)
#               n = longitude of the ascending node (degrees)
#               m = mean anomaly (degrees)

### Extract needed parameters from input list
	a,e,i,g,n,m = [float(i) for i in el]
	gm = G_mks*sum(mass)

### Convert input degrees to radians
#	i, g, n, m = deg2rad*i, deg2rad*g, deg2rad*n, deg2rad*m

### Pericenter
	q = (1.-e)*a

### Get p (longitude of pericenter) from g (argument of pericenter)
	p = g + n

### Rotation factors
	si, ci = sin(i), cos(i)
	sg, cg = sin(g), cos(g)
	sn, cn = sin(n), cos(n)
	z1 = cg * cn
	z2 = cg * sn
	z3 = sg * cn
	z4 = sg * sn
	d11 =  z1 - z4*ci
	d12 =  z2 + z3*ci
	d13 =       sg*si
	d21 = -z3 - z2*ci
	d22 = -z4 + z1*ci
	d23 =       cg*si

	assert (e < 1.0)	# haven't finished the other cases yet
### Calculate ellipse
	if (e < 1.0):
	### Ellipse
		romes = sqrt(1.0 - e*e)
		temp = M.Merc_KeplerEllipse(e,m)
		se, ce = sin(temp), cos(temp)
		z1 = a * (ce - e)
		z2 = a * romes * se
		temp = sqrt(gm/a) / (1.0 - e*ce)
		z3 = -se * temp
		z4 = romes * ce * temp
#	elif (e == 1.0) then
	### Parabola
#		ce = orbel_zget(l)
#		z1 = q * (1.0 - ce*ce)
#		z2 = 2.0 * q * ce
#		z4 = sqrt(2.0*gm/q) / (1.0 + ce*ce)
#	 	z3 = -ce * z4
#	elif (e > 1.):
	### Hyperbola
#		romes = sqrt(e*e - 1.0)
#		temp = orbel_fhybrid(e,l)
#		se, ce = sin(temp), cos(temp)
#		z1 = a * (ce - e)
#		z2 = -a * romes * se
#		temp = sqrt(gm/abs(a)) / (e*ce - 1.0)
#		z3 = -se * temp
#		z4 = romes * ce * temp
	
###
	x = d11 * z1  +  d21 * z2
	y = d12 * z1  +  d22 * z2
	z = d13 * z1  +  d23 * z2
	u = d11 * z3  +  d21 * z4
	v = d12 * z3  +  d22 * z4
	w = d13 * z3  +  d23 * z4

	return(x,y,z,u,v,w)

###############################################################################
def Merc_KeplerEllipse(e,oldl):
	'''Solves Kepler's equation for eccentricities less than one.
 Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.

  e = eccentricity
  l = mean anomaly      (radians)
  u = eccentric anomaly (   "   )'''
#------------------------------------------------------------------------------

	from math import pi, exp, log, sin, cos, sqrt
	twopi = 2.*pi
	piby2 = pi/2.

### Reduce mean anomaly to lie in the range 0 < l < pi
	if (oldl >= 0):
		l = oldl % twopi
	else:
		l = oldl % twopi + twopi
	sign = 1.0
	if (l > pi):
		l = twopi - l
		sign = -1.0

	ome = 1.0 - e
#-----------------------------
	if ((l >= .45) | (e < .55)):
### Regions A,B or C in Nijenhuis

### Rough starting value for eccentric anomaly
		if (l < ome):
			u1 = ome
		elif (l > (pi-1.0-e)):
			u1 = (l+e*pi)/(1.0+e)
		else:
			u1 = l + e

### Improved value using Halley's method
		flag = u1 > piby2
		if (flag):
			x = pi - u1
		else:
			x = u1
		x2 = x*x
		sn = x*(1.0 + x2*(-.16605 + x2*.00761) )
		dsn = 1.0 + x2*(-.49815 + x2*.03805)
		if (flag):
			dsn = -dsn
		f2 = e*sn
		f0 = u1 - f2 - l
		f1 = 1.0 - e*dsn
		u2 = u1 - f0/(f1 - .5*f0*f2/f1)

#---------------------
	else:
### Region D in Nijenhuis
### Rough starting value for eccentric anomaly
		z1 = 4.0*e + .50
		p = ome / z1
		q = .50 * l / z1
		p2 = p*p
		z2 = exp( log( sqrt( p2*p + q*q ) + q )/1.5 )
		u1 = 2.0*q / ( z2 + p + p2/z2 )

### Improved value using Newton's method
		z2 = u1*u1
		z3 = z2*z2
		u2 = u1 - .075*u1*z3 / (ome + z1*z2 + .375*z3)
		u2 = l + e*u2*( 3.0 - 4.0*u2*u2 )

### Accurate value using 3rd-order version of Newton's method
### N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!

### First get accurate values for u2 - sin(u2) and 1 - cos(u2)
	bigg = (u2 > piby2)
	if (bigg):
		z3 = pi - u2
	else:
		z3 = u2

	big = (z3 > (.5*piby2))
	if (big):
		x = piby2 - z3
	else:
		x = z3

	x2 = x*x
	ss = 1.0
	cc = 1.0

	ss = x*x2/6.*(1. - x2/20. *(1. - x2/42. *(1. - x2/72. *(1. -
         x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
	cc =   x2/2.*(1. - x2/12. *(1. - x2/30. *(1. - x2/56. *(1. -
         x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. -
         x2/306.))))))))

	if (big):
		z1 = cc + z3 - 1.0
		z2 = ss + z3 + 1.0 - piby2
	else:
		z1 = ss
		z2 = cc

	if (bigg):
		z1 = 2.0*u2 + z1 - pi
		z2 = 2.0 - z2

	f0 = l - u2*ome - e*z1
	f1 = ome + e*z2
	f2 = .5*e*(u2-z1)
	f3 = e/6.0*(1.0-z2)
	z1 = f0/f1
	z2 = f0/(f2*z1+f1)
	mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )

	return(mco_kep)

###########################################################################
### Convert from cartesian (xyz uvw) to orbital elements (aei gnM)
#   for orbit with no inclination
def El2X(el, mass):
	'''Convert orbital elements to cartesian for an ellipse (e < 1) with zero inclination
	mu = grav const * sum of masses
	q = perihelion distance
	e = eccentricity
	i = inclination                 )
	p = longitude of pericenter     )   in
	la = longitude of ascending node) radians
	f = true anomaly                )
	E = eccentric anomaly           )
	m = mean anomaly                )


	x,y,z = Cartesian positions  ( units the same as a )
	u,v,w =     "     velocities ( units the same as sqrt(mu/a) )'''


### Extract needed parameters from input list
	a,e,i,g,la,f = [float(i) for i in el]
	mu = G_mks*sum(mass)

### Convert input degrees to radians
#	i, g, la = deg2rad*i, deg2rad*g, deg2rad*la
	f = deg2rad*f

### Calculate radius
	r = a*(1-e**2)/(1+e*cos(f))

### Position coords
	x = r*cos(f)
	y = r*sin(f)
	z = 0.

### Period
	T = ( (4*pi**2/mu) * a**3 )**0.5

### Mean motion
	n = 2*pi/T

### Velocity coords
	u = -    sin(f)  * n*a/(1-e**2)**0.5
	v =   (e+cos(f)) * n*a/(1-e**2)**0.5
	w = 0.

	return(x,y,z,u,v,w)
	
###########################################################################


from cgs_constants import mSun,mEarth,AU,G
import DiskOOP as D
import Merc as M
import random as R
import numpy as np

###############################################################################
def GetObjList(rtr = 5., sigC = 10., rh = 10., dm = 1.e-7, alpha = 1.5):
	'''Generate list of objects based on the DiskOOP Heuristic model'''


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
def WriteObjInFile(objlist='default', loc = 'Merc95/In/',vers=0):
	'''Write a big.in or small.in file for mercury'''

	header = ')O+_06 '+['Big','Small'][vers]+ \
"-body initial data  (WARNING: Do not delete this line!!)\n) Lines beginning with `)' are ignored.\n)---------------------------------------------------------------------\n style (Cartesian, Asteroidal, Cometary) = {style}"\
+["\n epoch (in days) = {epoch}",""][vers]+ \
"\n)---------------------------------------------------------------------\n"

	objstr = '''  {0.name:16}  m={0.mass}  d={0.density}
    {0.a: .18e} {0.e: .18e} {0.i: .18e}
    {0.g: .18e} {0.n: .18e} {0.m: .18e}
    {0.s[0]: 19.18e} {0.s[1]: 19.18e} {0.s[2]: 19.18e}
'''

	style = 'Asteroidal'
	epoch = 0.0

	if (objlist=='default'):
		objlist = M.GetObjList()

	if (vers == 0):
		fname = loc+'big.in'
	elif (vers == 1):
		fname = loc+'small.in'
	else:
		print('invalid version number: must be 0 for big.in or 1 for small.in')

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


###############################################################################
class Obj(object):
	'''Object in a mercury simulation, currently in asteroidal coords.'''

#-----------------------------------------------------------------------------#
	def __init__(self, name='Earth', mass=mEarth/mSun, density=5.51,
		a=1., e=0., i=0., g=0., n=0., m=0., s=[0., 0., 0.]):

		self.name    = name
		self.mass    = mass
		self.density = density
		self.a       = a
		self.e       = e
		self.i       = i
		self.g       = g
		self.n       = n
		self.m       = m
		self.s       = s
#-----------------------------------------------------------------------------#

###############################################################################


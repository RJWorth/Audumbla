### Try machine learning techniques on Comparison mercury data
%matplotlib
import matplotlib.pyplot as plt
import seaborn as sb
import Disks
import Merc
from mks_constants import mEarth, mSun
import numpy as np
import pandas as pd
#from patsy import standardize
#from sklearn import metrics
#from sklearn import tree
#from sklearn import neighbors
#from sklearn import linear_model
#from sklearn.cross_validation import train_test_split, KFold
#from sklearn.grid_search import GridSearchCV

##############################################################################
WhichDir='Comparison2'
### Get mass of central object (Solar masses)
mCent = Merc.ReadMCent(WhichDir+'/')
### Get interaction lists
collisions = {}
collisions['name'], collisions['dest'], collisions['time'] = Merc.ReadInfo(WhichDir)
collisions = pd.DataFrame.from_dict(collisions)

### Get initial objects
objs = Merc.ReadInObjList(WhichDir+'/In/','big.in')
namelist = [o.name for o in objs]
nobjs = len(objs)
del objs

### Read in aei for each object to a dict
### too many rows to fit in memory, get multiple chunks and analyze separately?
### first chunk:
chunksize = 100
cube = {}
for n in namelist:
	aei = Merc.ReadAEI(WhichDir, n, nrows=chunksize)
	aei['m'] = aei.mass*mSun/mEarth # convert to Earth masses
	aei['p']  = aei.a * (1-aei.e)   # pericenter
	aei['ap'] = aei.a * (1+aei.e)   # apocenter
	# leave out density, Cartesian coordinates
	cube[n] = aei[['t','m','a','e','i','p','ap']]

### Turn dict into panel 
cube = pd.Panel.from_dict(cube)
nt = cube.shape[1]

### Get full list of timesteps from a complete object
tlist = cube[pd.notnull(cube[:,nt-1,'t']),:,'t'].iloc[:,0]

#items = namelist
#major_axis = times
#minor_axis = ['t','a','e','i','m']
# cube[obj].iloc[t, property#]
# cube[obj][property]
# cube.loc['P000',0:10,'m']

### Interaction terms between each pair of objects, at each timestep
da      = {}
RH2     = {}
overlap = {}

for i,obj1 in enumerate(namelist):
	# separation in a
	da[obj1]      = cube[:,:,'a'].apply(lambda x: x - cube[obj1,:,'a'])
	# sums of a, m
	suma          = cube[:,:,'a'].apply(lambda x: x + cube[obj1,:,'a'])
	summ          = cube[:,:,'m'].apply(lambda x: x + cube[obj1,:,'m'])
	# mutual hill radius
	RH2[obj1]     = (summ/(3.*mCent))**(1./3.) * (suma/2.)
	# overlap in orbital distance (between peri- and apocenters)
	minap         = cube[:,:,'ap'].apply(lambda x: np.min([x, cube[obj1,:,'ap']], axis=0 ))
	maxp          = cube[:,:, 'p'].apply(lambda x: np.max([x, cube[obj1,:, 'p']], axis=0 ))
	overlap[obj1] = minap - maxp

da      = pd.Panel.from_dict(da)
RH2     = pd.Panel.from_dict(RH2)
overlap = pd.Panel.from_dict(overlap)
del suma, summ, minap, maxp

### The above doesn't complete -- too much memory. 
### Reformulate to iterate over timesteps, analyzing each individually
### rather than saving all the info and analyzing all at the end
#for i,t in enumerate(tlist):
i = 1
t = tlist[i]
	# collisions this step
	CollThisStep = collisions.loc[(collisions.time>t) & (collisions.time<=tlist[i+1])]
	# separation in a
	a = cube[:,i,'a']
	m = cube[:,i,'m']
	da      = np.repeat([list(a)],len(a),axis=0) - np.repeat([list(a)],len(a),axis=0).transpose()
	# sums of a, m
	suma    = np.repeat([list(a)],len(a),axis=0) + np.repeat([list(a)],len(a),axis=0).transpose()
	summ    = np.repeat([list(m)],len(m),axis=0) + np.repeat([list(m)],len(m),axis=0).transpose()
	# mutual hill radius
	RH2     = (summ/(3.*mCent))**(1./3.) * (suma/2.)
	# overlap in orbital distance (between peri- and apocenters)
	minap = np.array([min(ap1, ap2) for ap1 in cube[:,i,'ap'] 
                                    for ap2 in cube[:,i,'ap']]).reshape((len(a),len(a)))
	maxp  = np.array([max( p1,  p2) for  p1 in cube[:,i, 'p'] 
                                    for  p2 in cube[:,i, 'p']]).reshape((len(a),len(a)))

	overlap = pd.DataFrame(minap - maxp, index=namelist, columns=namelist)
	ejoverlap = cube[:,i,'ap']-100.
	overlap['ejected'] = ejoverlap
	ejoverlap.loc['ejected'] = 0.
	overlap.loc['ejected'] = ejoverlap

### Compare da, RH2, overlap for name-dest collision pairs
	colloverlap = np.array([overlap.loc[collisions['name'].iloc[j],collisions['dest'].iloc[j]] for j in range(collisions.shape[0])])

	### future collision times vs current orbit overlap
	plt.scatter(collisions['time'], colloverlap)

	### All overlap values for this timestep in red, and ones that collide in blue
	x = np.array(overlap).flatten()
	x = x[~np.isnan(x) & (x>-80)]
	x.sort()

	CollOverlapThisStep = np.array([overlap.loc[ CollThisStep.dest.iloc[j], CollThisStep.name.iloc[j] ] for j in range(CollThisStep.shape[0])])
	CollTimeThisStep = CollThisStep.time
	plt.plot(CollOverlapThisStep, CollTimeThisStep, 'ro')

	### This shows that collisions did in fact happen mostly on the high end
	### of the distribution of overlap values => it is predictive of collisions
	### at least in the near term
	sb.distplot(x, bins=20, hist=True, kde=False, norm_hist=True)
	sb.distplot(y, hist=True, kde=False, norm_hist=True)
	for j in range(len(y)):
		plt.plot((y[j], y[j]), (0, 25000), 'k-')

### Look at each particle's drift to next step, vs neighbors?
	### Overtime, objects rearrange, have different ordering... how to use???
	alist = np.array(cube[:,i,'a'])
	order = np.argsort(alist)
	## indices of objects which are not nan, from smallest to greatest
	order = order[:(len(alist)-np.sum(np.isnan(alist)))]

	Step = {}
	Step['drift']  = cube[:,i+1,'a'] - cube[:,i,'a']
	Step = pd.DataFrame.from_dict(Step)
	Step['a']      = cube[:,i,'a']
	Step['masses'] = cube[:,i,'m']

	### RH2 btwn this object and outer/inner neighbor
	Step['RHout'] = np.array(     [RH2[k,k+1] for k in range(nobjs-1)]+[0.])
	Step['RHin']  = np.array([0.]+[RH2[k-1,k] for k in range(nobjs-1)]     )
	### Orbital overlap with outer/inner neighbor
	# have to add one entry to end/beginning respectively, for lack of neighbor
	Step['OOout'] = np.array([overlap.loc[namelist[k],namelist[k+1]] for k in range(nobjs-1)]+[0.])
	Step['OOin']  = np.array([0.]+
                     [overlap.loc[namelist[k-1],namelist[k]] for k in range(1,nobjs)] )
	### Orbital overlap in terms of RH
	Step['ORout'] = Step['OOout']/Step['RHout']
	Step['ORin' ] = Step['OOin' ]/Step['RHin' ]

	### Bin drift, color-code it, and plot; black = positive, white=negative
	driftbins = [-1, -.01, -.002, 0., 0.002, .01, 1]
	driftbinned = pd.cut( Step.drift, bins = driftbins, labels=False)
	Step['driftBinned'] = driftbinned
#	Step['driftBinned'] = Step['driftBinned']+1
#	Step['driftBinned'] = Step['driftBinned'].astype('category')

	### Visualizations
	# looking for trend with color (drift) and... anything else, but not seeing it
	sb.pairplot(vars=['drift','a','masses','RHout','RHin','OOout','OOin','ORout','ORin'],
		data=Step.dropna(), hue='driftBinned', 
		palette = sb.diverging_palette(10, 250, s=85, l=50, n=len(driftbins)-1, center="dark") )

	### See how drift magnitude is distributed over inner/outer spacing
	plt.scatter( Step.OOin, Step.OOout, c=driftbinned)


	plt.scatter( range(nobjs), RHout, c='blue')
	plt.scatter( range(nobjs), RHin , c='red')


	plt.scatter( range(nobjs), OOout, c='blue')
	plt.scatter( range(nobjs), OOin , c='red')
	plt.scatter( range(nobjs), drift, c='black')



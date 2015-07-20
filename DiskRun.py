
import Disks as d
import numpy as np
import matplotlib.pyplot as plt
from cgs_constants import mSun,mEarth,mMoon,mMars,AU,G
#-----------------------------------------------------------------------------#
rtr = 5.    # truncation radius
s   =10.    # coefficient of surface density (sigma = MMSN*s)
rh  =10.    # number of hill radii counting as interaction
m = [mMoon/mSun,mMars/mSun] # possible planetesimal masses
f = [0.5, 0.5]    # fraction of total mass in each size
a   =1.5    # alpha, slope of surface density profile
ej  = .5    # fraction of interactions that should lead to ejection rather than collision
#-----------------------------------------------------------------------------#
reload(d)
disk = d.Disk(r_out=rtr, alpha=a, sigC=s, rh=rh)
disk.DebrisGenM(m,f)
p, p_a, p_m, CorS = disk.debris.PltFormSuccessiveClearings(EjProb=ej)

### Plot
for i in range(len(p_a)):
	plt.plot([i for n in range(len(p_a[i]))], p_a[i], 'ro')

plt.xlabel('# of Collisions')
plt.ylabel('Semimajor Axis (AU)')
plt.title('Planetesimal Merging Heuristic')
plt.savefig('s{0}_r{1}_a{2}_ej{3}.png'.format(s,rh,a,ej))
plt.clf()

### To save and reload arrays
#import cPickle as pickle
#pickle.dump(   p1, open(   "p1.p", "wb" ) )
#pickle.dump( p_a1, open( "p_a1.p", "wb" ) )
#pickle.dump( p_m1, open( "p_m1.p", "wb" ) )
#p1 = pickle.load( open( "p1.p", "rb" ) )



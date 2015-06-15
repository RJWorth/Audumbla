
import Disks as d
import numpy as np
import matplotlib.pyplot as plt
from cgs_constants import mSun,mEarth,AU,G

rtr = 5.
s   =10.
rh  =10.
dm  =1.e-7
a   =1.5

reload(d)
disk = d.Disk(r_out=rtr, alpha=a, sigC=s, rh=rh)
disk.DebrisGenM(dm)
len(disk.debris)
#disk.DebrisGenA(da)

### RH/da ratio at center
#p3 = [disk.debris]
#mid = (len(p3[0])/2)
#print( abs(p3[0][mid].a-p3[0][mid+1].a),p3[0][mid].RH2(p3[0][mid+1]), 
#       abs(p3[0][mid].a-p3[0][mid+1].a)/p3[0][mid].RH2(p3[0][mid+1]) )

pl ,pa,pm = disk.debris.PltFormTilClear(1000)
pl1      = disk1.debris.PltFormTilClear(1000)

p, p_a, p_m = disk.debris.PltFormSuccessiveClearings()
#p1, p_a1, p_m1 = disk.debris.PltFormSuccessiveClearings()

import cPickle as pickle
pickle.dump(   p1, open(   "p1.p", "wb" ) )
pickle.dump( p_a1, open( "p_a1.p", "wb" ) )
pickle.dump( p_m1, open( "p_m1.p", "wb" ) )
p1 = pickle.load( open( "p1.p", "rb" ) )


for i in range(len(p_a)):
	plt.plot([i for n in range(len(p_a[i]))], p_a[i], 'ro')

plt.xlabel('# of Collisions')
plt.ylabel('Semimajor Axis (AU)')
plt.title('Planetesimal Merging Heuristic')
plt.savefig('s{0}_r{1}_m{2}_a{3}-1.png'.format(s,rh,dm,a))
plt.clf()

#plt.savefig('s{0}da{1}r10.png'.format(s,da))

### How does RH2 change with a?
r = [disk.debris[i].RH2(disk.debris[i+1]) for i in range(len(disk.debris)-1)]
a = [disk.debris[i].a                     for i in range(len(disk.debris)-1)]
plt.plot(a,r)
plt.show()


#reload(d)
#disk = d.Disk(r_out=5.)
#disk.DebrisGen(.1)
#p1 = disk.debris.PltFormRand()

#reload(d)
#disk = d.Disk(r_out=5.)
#disk.DebrisGen(.1)
#p2 = disk.debris.PltFormTilClear()


#for i in range(len(p)):
#	for j in range(len(p[i].ListParams()[0])):
#		plt.plot(i, p[i].ListParams()[0][j], 'ro')

#plt.show()

#for i in range(len(p2)):
#	for j in range(len(p2[i].ListParams()[0])):
#		plt.plot(i, p2[i].ListParams()[0][j], 'ro')

#plt.show()



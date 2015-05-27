
import DiskOOP as d
import numpy as np
import matplotlib.pyplot as plt
from cgs_constants import mSun,mEarth,AU,G

s=1.
da=0.005	# can't run plt form on this size -- python gets killed (memory probably)

reload(d)
disk = d.Disk(r_out=5.,sig=1.)
disk.DebrisGen(da)
### RH/da ratio at center
p3 = [disk.debris]
mid = (len(p3[0])/2)
print( abs(p3[0][mid].a-p3[0][mid+1].a),p3[0][mid].RH2(p3[0][mid+1]), 
       abs(p3[0][mid].a-p3[0][mid+1].a)/p3[0][mid].RH2(p3[0][mid+1]) )

p, p_a, p_m = disk.debris.PltFormSuccessiveClearings(10)



for i in range(len(p_a)):
	plt.plot([i for n in range(len(p_a[i]))], p_a[i], 'ro')

plt.savefig('s{0}da{0}r10.png'.format(s,da))
plt.clf()

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



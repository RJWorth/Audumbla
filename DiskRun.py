
import Disks as d
import Merc as M
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages
from cgs_constants import mSun,mEarth,mMoon,mMars,AU,G
import cProfile
import re
import pstats
#-----------------------------------------------------------------------------#
rtr = 1.0     # truncation radius
s   = 3.      # coefficient of surface density (sigma = MMSN*s)
rh  =10.      # number of hill radii counting as interaction
m = [mMoon/mSun,mMars/mSun] # possible planetesimal masses
f = [0.5, 0.5] # fraction of total mass in each size
a   =1.5      # alpha, slope of surface density profile
ej  = .2      # fraction of interactions that should lead to ejection rather than collision
dr  = [0.,0.] #[4.e-10,5.e-10]  # a drifts randomly by up to this much (in AU) each step, in/out respectively
#-----------------------------------------------------------------------------#
reload(d)
disk = d.Disk(r_out=rtr, alpha=a, sigC=s, rh=rh)
disk.DebrisGenM(m,f)
pl_n, a_n, m_n, t, CorS = disk.debris.IceCow(EjProb=ej,nRh=rh,drift=dr,vers="new")
d.PlotDiskMerge(a_n,fname='newprob_n')
d.PlotDiskMerge(a_n,fname='newprob_t',t=t)

cProfile.run('pl_n, a_n, m_n, t, CorS = disk.debris.IceCow(EjProb=ej,nrh=rh,drift=dr,vers="new")','new')
p_new = pstats.Stats('new')
d.PlotDiskMerge(a_n,fname='testnew')

cProfile.run('pl_o, a_o, m_o, t, CorS = disk.debris.IceCow(EjProb=ej, nrh=rh, drift=dr, vers="old")','old')
p_old = pstats.Stats('old')
d.PlotDiskMerge(a_o,fname='testold')

p_old.sort_stats('time').print_stats(10)
p_new.sort_stats('time').print_stats(10)

planets, a_icecow, m_icecow, CorS_icecow = disk.debris.IceCow(EjProb=ej,nrh=rh,drift=dr)
d.PlotDiskMerge(a_icecow,ylim=[0.,4.])

### Read in and plot merge history from MERCURY sim
a_merc, m_merc, t_merc, acoll, da, Rh, CorS_merc, mnsep, mdsep = d.Comparison('Comparison')
d.PlotDiskMerge(a=a_merc,m=m_merc,fname='Merc_n',ylim=[0.,4.])
d.PlotDiskMerge(a=a_merc,m=m_merc,tlist=t_merc,fname='Merc_t',ylim=[0.,4.])
d.PlotDiskMerge(a=a_merc,m=m_merc,tlist=t_merc,da=da,Rh=Rh,fname='Merc_t',ylim=[0.,4.])
d.PlotDiskMerge(a=a_merc,m=m_merc,
tlist=t_merc,acoll=acoll,da=da,Rh=Rh,mnsep=mnsep,mdsep=mdsep,fname='Merc_t',ylim=[0.,4.])

### To time and compare modifications:
#cProfile.run('p, p_a, p_m, CorS = disk.debris.PltFormProb(EjProb=ej,nrh=rh,drift=dr,vers="old")','old')
#p_old = pstats.Stats('old')
#cProfile.run('p, p_a, p_m, CorS = disk.debris.PltFormProb(EjProb=ej,nrh=rh,drift=dr,vers="new")','new')
#p_new = pstats.Stats('new')
#p_old.sort_stats('time').print_stats(10)
#p_new.sort_stats('time').print_stats(10)

### Notes on differences:
# mercury systems clear sooner in the inner part of the system => add an inverse-a term to probability?
# drift only occurs when objects are not yet well spaced? or when they are still small? it stops in inner system once they are mostly well separated
# drift happens preferentially outward, and is faster at later steps. Reintroduce a 'time' parameter and make drift proportional to it? also proportional to probabilities of interaction?




### To save and reload arrays
#import cPickle as pickle
#pickle.dump(   p1, open(   "p1.p", "wb" ) )
#pickle.dump( p_a1, open( "p_a1.p", "wb" ) )
#pickle.dump( p_m1, open( "p_m1.p", "wb" ) )
#p1 = pickle.load( open( "p1.p", "rb" ) )



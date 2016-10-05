
### Get packages (Disks, Merc, and cgs_constants are made by me)
import Disks as d
import Merc as M
import numpy as np
import matplotlib.pyplot as plt
from cgs_constants import mSun,mEarth,mMoon,mMars,AU,G
import cProfile
import re
import pstats
#-----------------------------------------------------------------------------#
### Input parameters
rtr = 3.08     # truncation radius
s   = 3.0      # coefficient of surface density (sigma = MMSN*s)
rh  = 9.       # number of hill radii counting as interaction
m = [mMoon/mSun,mMars/mSun] # possible planetesimal masses
f = [0.5, 0.5] # fraction of total mass in each size
a   =1.5       # alpha, slope of surface density profile
ej  = .2       # fraction of interactions that should lead to ejection rather than collision
dr  = [4.e-10,4.5e-10]  # a drifts randomly by up to this much (in AU) each step, in/out respectively
standardylim = [0., 4.5] # upper and lower values of a to use in plots
#-----------------------------------------------------------------------------#
# if tinkering with the code, reload and run from here
#reload(d)
### Set up disk model
disk = d.Disk(r_out=rtr, alpha=a, sigC=s, rh=rh)
### Generate debris disk objects
disk.DebrisGenM(m,f)
### Run the merging simulation
pl_n, a_n, m_n, t, CorS = disk.debris.IceCow(EjProb=ej,nRh=rh,drift=dr,vers="new")
### Plot
d.PlotDiskMerge(a=a_n, m=m_n, fname='newprob_n',ylim=standardylim)

#-----------------------------------------------------------------------------#
### For comparison, plot results of a MERCURY simulation the same way
### Read in and plot merge history from MERCURY sim
a_merc, e_merc, i_merc, m_merc, \
	t_merc, acoll, ecoll, icoll, mcoll, da, Rh, CorS_merc, mnsep, mdsep \
	= d.Comparison('Comparison')
### Make the plot
d.PlotDiskMerge(a=a_merc,m=m_merc,fname='Merc_n',ylim=standardylim)
#-----------------------------------------------------------------------------#


import matplotlib.pylab        as plt
import matplotlib.pyplot
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker       import ScalarFormatter
from   matplotlib.ticker       import FixedFormatter
import pylab                   as pl
import numpy                   as np
import math, os
import glob, pickle

formatter = ScalarFormatter(useMathText=True)
#formatter.set_scientific(True)
#formatter.set_powerlimits((-3,3))

fig_width_pt = 246.0*2 # Get this from LaTex using \the\columnwidth
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width  = fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize':10,
          'text.fontsize':8,
          'legend.fontsize':8,
          'xtick.labelsize':8.5,
          'ytick.labelsize':5.5,
          'figure.figsize':fig_size,
          'font.family': 'serif'}

pl.rcParams.update(params)
pl.clf()
pl.figure()
fig = pl.figure()
axes = pl.Axes(fig, [.2, .2, .7, .7])
fig.add_axes(axes)
axes.yaxis.set_major_formatter(formatter)

for i in xrange(1, 10, 1):
  data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_xvel_clipThreshold_1.0e+03_fullCube_BootStrap_kbin_0.010_00'+str(i)+'.dat')
  pl.loglog(data[:,0], data[:,1])

for i in xrange(10, 40, 1):
  data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_xvel_clipThreshold_1.0e+03_fullCube_BootStrap_kbin_0.010_0'+str(i)+'.dat')
  pl.loglog(data[:,0], data[:,1])
  
for i in xrange(1, 10, 1):
  data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_xvel_clipThreshold_1.0e+03_fullCube_BootStrap_kbin_0.010_00'+str(i)+'.dat')
  pl.loglog(data[:,0], data[:,2], '--')

for i in xrange(10, 40, 1):
  data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_xvel_clipThreshold_1.0e+03_fullCube_BootStrap_kbin_0.010_0'+str(i)+'.dat')
  pl.loglog(data[:,0], data[:,2], '--')

pl.xlabel(r'$k [h^{-1} Mpc]$')
pl.ylabel('P(k)')

pl.yscale('log')
pl.xscale('log')

pl.xlim(10**-2, 1.0)
pl.ylim(100, 5*10**4)

#pl.clf()

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Pk_BootStrap.pdf')
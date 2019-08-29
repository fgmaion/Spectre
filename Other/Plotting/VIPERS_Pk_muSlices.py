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
          'legend.fontsize':6,
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

HOD     = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/Pk_hod_20.0.dat')
linear  = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/camb_matterPk.dat')

beta    = 0.45
linbias = 1.495903**2

fig = pl.figure()

f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/muSlicedPk/zCube_xvel_clipThreshold_1.0e+03_fullCube_muSlicedPk_0.0mu0.2.dat')
ax1.loglog(data[:,0], data[:,1], 'r', label=r'$0.0 < \mu < 0.2$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/muSlicedPk/zCube_xvel_clipThreshold_5.0e+00_fullCube_muSlicedPk_0.0mu0.2.dat')
ax1.loglog(data[:,0], data[:,1], 'r--', label=r'$0.0 < \mu < 0.2$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/muSlicedPk/zCube_xvel_clipThreshold_1.5e+00_fullCube_muSlicedPk_0.0mu0.2.dat')
ax1.loglog(data[:,0], data[:,1], 'r', label=r'$0.0 < \mu < 0.2$')

ax1.loglog(HOD[:,0],               HOD[:,1]*(1. + beta*0.1*0.1)**2, 'k')
ax1.loglog(linear[:,0], linbias*linear[:,1]*(1. + beta*0.1*0.1)**2, 'k--')
ax1.loglog(linear[:,0],     0.45*linear[:,1]*(1. + beta*0.1*0.1)**2, 'k--')



data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/muSlicedPk/zCube_xvel_clipThreshold_1.0e+03_fullCube_muSlicedPk_0.3mu0.5.dat')
ax2.loglog(data[:,0], data[:,1], 'g', label=r'$0.3 < \mu < 0.5$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/muSlicedPk/zCube_xvel_clipThreshold_5.0e+00_fullCube_muSlicedPk_0.3mu0.5.dat')
ax2.loglog(data[:,0], data[:,1], 'g--', label=r'$0.3 < \mu < 0.5$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/muSlicedPk/zCube_xvel_clipThreshold_1.5e+00_fullCube_muSlicedPk_0.3mu0.5.dat')
ax2.loglog(data[:,0],   data[:,1], 'g-.', label=r'$0.3 < \mu < 0.5$')

ax2.loglog(HOD[:,0],               HOD[:,1]*(1. + beta*0.4*0.4)**2, 'k')
ax2.loglog(linear[:,0], linbias*linear[:,1]*(1. + beta*0.4*0.4)**2, 'k--')
ax2.loglog(linear[:,0],    0.45*linear[:,1]*(1. + beta*0.4*0.4)**2, 'k--')



data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/muSlicedPk/zCube_xvel_clipThreshold_1.0e+03_fullCube_muSlicedPk_0.6mu0.8.dat')
ax3.loglog(data[:,0], data[:,1], 'b', label=r'$\delta_0 = 10^{3}$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/muSlicedPk/zCube_xvel_clipThreshold_5.0e+00_fullCube_muSlicedPk_0.6mu0.8.dat')
ax3.loglog(data[:,0], data[:,1], 'b--', label=r'$\delta_0 = 5.$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/muSlicedPk/zCube_xvel_clipThreshold_1.5e+00_fullCube_muSlicedPk_0.6mu0.8.dat')
ax3.loglog(data[:,0], data[:,1], 'b:', label=r'$\delta_0 = 1.5$')

# data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/muSlicedPk/zCube_xvel_clipThreshold_2.0e+00_fullCube_muSlicedPk_0.6mu0.8.dat')
# ax3.loglog(data[:,0], data[:,1], 'b-.', label=r'$0.6 < \mu < 0.8$')

ax3.loglog(HOD[:,0],               HOD[:,1]*(1. + beta*0.7*0.7)**2, 'k')
ax3.loglog(linear[:,0], linbias*linear[:,1]*(1. + beta*0.7*0.7)**2, 'k--')
ax3.loglog(linear[:,0],    0.45*linear[:,1]*(1. + beta*0.7*0.7)**2, 'k--')

pl.yscale('log')
pl.xscale('log')

pl.xlim(0.01, 1.0)
pl.ylim(100, 9*10**4)

#pl.clf()

ax1.set_title('$0.0 < \mu < 0.2$', fontsize=8)
ax2.set_title('$0.3 < \mu < 0.5$', fontsize=8)
ax3.set_title('$0.6 < \mu < 0.8$', fontsize=8)

ax3.legend(loc=3, handlelength=2)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/HODCube_muSlicedPk.pdf')

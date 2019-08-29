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

fig_width_pt = 246.0*6 # Get this from LaTex using \the\columnwidth
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width  = fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize':14,
          'text.fontsize':8,
          'legend.fontsize':14,
          'xtick.labelsize':12.0,
          'ytick.labelsize':12.0,
          'figure.figsize':fig_size,
          'font.family': 'serif'}

pl.rcParams.update(params)
pl.clf()
pl.figure()
fig = pl.figure()
axes = pl.Axes(fig, [.2, .2, .7, .7])
fig.add_axes(axes)
axes.yaxis.set_major_formatter(formatter)

# import copy
# import scipy.interpolate

from   matplotlib.colors import LogNorm

# cmap = copy.copy(matplotlib.cm.jet)

# cmap.set_bad('w',1.)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Covariance/zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_Covariance.dat')
# data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Covariance/GaussCube_Realisation_clipThreshold_1.0e+03_fullCube_Covariance.dat')
# data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Covariance/GaussCube_BootStrap_clipThreshold_1.0e+03_fullCube_Covariance.dat')


# plt.imshow(data, origin='lower', cmap='bone_r', interpolation='nearest', norm = LogNorm(vmin=10**0, vmax=5*10**8))

plt.imshow(data, origin='lower', cmap='RdGy', interpolation='nearest', vmin=-1, vmax=1)

# print data.max()
# plt.imshow(data, origin='lower', cmap='YlGnBu', interpolation='nearest')

pl.xticks([10, 20, 30, 40, 50, 60, 70, 80], [str(10*0.01), str(20*0.01), str(30*0.01),  str(40*0.01), str(50*0.01), str(60*0.01), str(70*0.01), str(80*0.01)])

pl.yticks([10, 20, 30, 40, 50, 60, 70, 80], [str(10*0.01), str(20*0.01), str(30*0.01),  str(40*0.01), str(50*0.01), str(60*0.01), str(70*0.01), str(80*0.01)])

plt.colorbar()

pl.xlabel('k, [$h^{-1}$ Mpc]')

# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/GaussCube_Realisation_clipThreshold_1.0e+03_fullCube_Covariance.pdf', bbox_inches='tight', pad_inches=0.5)
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/zCube_BootStrap_clipThreshold_1.0e+03_fullCube_Covariance.pdf', bbox_inches='tight', pad_inches=0.5)

pl.clf()


'''
data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Covariance/GaussCube_Realisation_clipThreshold_1.0e+03_fullCube_DiagCov.dat')

pl.loglog(data[:,0], data[:,1], label='measured')

pl.loglog(data[:,0], (data[:,2]/np.sqrt(data[:,3]))**2, label='predicted')


pl.xlabel(r'$k [h^{-1} Mpc]$')
pl.ylabel(r'$\sigma^2(k)$')

pl.xlim(10**-2., 0.8)

pl.yscale('log')
pl.xscale('log')

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/GaussianRandomField_DiagCov.pdf')


pl.clf()

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Covariance/GaussCube_Realisation_clipThreshold_1.0e+03_fullCube_ChiSq.dat')

pl.loglog(data[:,0], data[:,1])
pl.loglog(data[:,0], data[:,0])

pl.loglog(data[:,0], data[:,2])


pl.xlabel(r'$N_u$')


# pl.yscale('log')
# pl.xscale('log')

# pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/GaussianRandomField_ChiSq.pdf')
'''

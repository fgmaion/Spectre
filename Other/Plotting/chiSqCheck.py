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


pl.xlabel(r'$k [h^{-1} Mpc]$')

# pl.xlim(10**-2, 1.0)
# pl.ylim(100, 5*10**4)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/bestfit_zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_beta_0.51_sigma_3.40_A11Sq_1.00_kmax_0.30_LorentzianRSD.dat')
pl.loglog(data[:,0], data[:,1])
pl.loglog(data[:,0], data[:,2])

pl.loglog(data[:,0], data[:,3], '--')
pl.loglog(data[:,0], data[:,4], '--', label='Lorentzian, $k_{max} = 0.3$')


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/bestfit_zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_beta_0.48_sigma_3.20_A11Sq_1.00_kmax_0.40_LorentzianRSD.dat')
pl.loglog(data[:,0], data[:,1])
pl.loglog(data[:,0], data[:,2])

pl.loglog(data[:,0], data[:,3], '--')
pl.loglog(data[:,0], data[:,4], '--', label='Lorentzian, $k_{max} = 0.4$')


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/bestfit_zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_beta_0.51_sigma_2.20_A11Sq_1.00_kmax_0.30_GaussianRSD.dat')
# pl.loglog(data[:,0], data[:,1])
# pl.loglog(data[:,0], data[:,2])

pl.loglog(data[:,0], data[:,3], '--')
pl.loglog(data[:,0], data[:,4], '--', label='Gaussian, $k_{max} = 0.3$')


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/bestfit_zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_beta_0.46_sigma_2.00_A11Sq_1.00_kmax_0.40_GaussianRSD.dat')
# pl.loglog(data[:,0], data[:,1])
# pl.loglog(data[:,0], data[:,2])

pl.loglog(data[:,0], data[:,3], '--')
pl.loglog(data[:,0], data[:,4], '--', label='Gaussian, $k_{max} = 0.4$')

pl.yscale('log')
pl.xscale('log')

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ChiSq_Check.pdf')

pl.clf()

pl.xlabel(r'$k [h^{-1} Mpc]$')

pl.errorbar(data[:,0], data[:,5], data[:,9])
pl.errorbar(data[:,0], data[:,6], data[:,10])

pl.loglog(data[:,0], data[:,7], '--')
pl.loglog(data[:,0], data[:,8], '--')

pl.yscale('linear')
pl.xscale('log')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ChiSq_Check_y.pdf')

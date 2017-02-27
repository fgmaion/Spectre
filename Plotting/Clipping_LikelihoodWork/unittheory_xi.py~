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

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/fftlog_hodpk_NoRSD_4096_1.00e-10_1.00e+14_11Aug_xi.dat')
pl.loglog(data[:,0],  np.abs(data[:,1]), 'y', label=r'$\xi_0$(r)')
pl.loglog(data[:,0],  np.abs(data[:,2]), 'r', label=r'$\xi_2$(r)')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/corrfn_multipoles_4.00_11Aug.dat')
pl.loglog(data[:,0], np.abs(data[:,1]), 'k^', markersize='2')
pl.loglog(data[:,0], np.abs(data[:,2]), 'k^', markersize='2')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/unitTheory_xi_mono_quad.dat')
pl.loglog(data[:,0], np.abs(data[:,1]), 'b--', markersize='2')
pl.loglog(data[:,0], np.abs(data[:,2]), 'g--', markersize='2')

pl.xlabel(r'$r [h^{-1} Mpc]$')
pl.ylabel(r'|$\xi$(r)|')

pl.yscale('log')
pl.xscale('log')
 
pl.xlim(10.**1., 300.)
pl.ylim(10.**-4.,  1.)

pl.legend(loc=1)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/quadrupole_xi.pdf')


pl.clf()

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/Si_minus_Si_Convergence.dat')
pl.plot(data[:,0], data[:,1] , 'y', label=r'Si - Si')

pl.ylim(-0.5, 1.01)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/Si_minusSiConvergence.pdf')

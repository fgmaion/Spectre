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

fig_width_pt = 246.0*2 # Get this from LaTex using \the\columnwidth
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width  = fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize':10,
          'text.fontsize':6,
          'legend.fontsize':6,
          'xtick.labelsize':11.0,
          'ytick.labelsize':11.0,
          'figure.figsize':fig_size,
          'font.family': 'serif'}

pl.rcParams.update(params)
pl.clf()
pl.figure()
fig = pl.figure()
axes = pl.Axes(fig, [.2, .2, .7, .7])
fig.add_axes(axes)
axes.yaxis.set_major_formatter(formatter)

# lines = {'linestyle': 'None'}
# plt.rc('lines', **lines)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/projected_xi_0.00_100.00.dat')
pl.loglog(data[:,0], 2.*np.abs(data[:,1]), 'k')
# pl.loglog(data[:,0], np.abs(data[:,2]),    'r--')

# data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/projected_xi_0.00_500.00.dat')
# pl.loglog(data[:,0], np.abs(data[:,1]), 'r')

'''
for i in xrange(2, 6, 1):
    Index = np.where(np.abs(data[:,0] - 10.**(-1. + 0.3*(i-2))) == np.abs(data[:,0] - 10.**(-1. + 0.3*(i-2))).min())
    
    pl.loglog(data[Index[0]:,0], np.abs(data[Index[0]:, i]), label=str(10.**(-1. + 0.3*(i-2))))
'''
pl.xscale('log')
pl.yscale('log')

pl.xlabel(r'$r_{\perp} [h^{-1} Mpc]$')
pl.ylabel(r'w$(r_{\perp})$')

# pl.xlabel(r'$r [h Mpc^{-1}]$')
# pl.ylabel(r'$\xi(r)$')

pl.xlim(10.**-1., 100.0)
pl.ylim(10.**-1.,  100.)

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/SSPOC/projected_xi.pdf')

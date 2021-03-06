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

#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/GalaxyCatalogue.dat')
#mod  = data[:,1]**2 + data[:,2]**2 + data[:,3]**2
#mod  = mod**0.5  

#a    = pl.hist(mod, bins=50)

#np.savetxt('/disk1/mjw/HOD_MockRun/Scripts/Plotting/SSPOC/height.dat', a[0])

#np.savetxt('/disk1/mjw/HOD_MockRun/Scripts/Plotting/SSPOC/mid.dat',    a[1])

height = np.loadtxt('/disk1/mjw/HOD_MockRun/Scripts/Plotting/SSPOC/height.dat')

mid    = np.loadtxt('/disk1/mjw/HOD_MockRun/Scripts/Plotting/SSPOC/mid.dat')

pl.loglog(mid[:-1], height/height.max())

# data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/GalaxyCatalogue_spec.dat')


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/fftlog_hodpk_NoRSD_4096_1.00e-10_1.00e+14_05Sep_xi.dat')
pl.loglog(data[:,0], data[:,0]**2*data[:,1]/100., 'r')

pl.xscale('log')
#pl.yscale('log')

# pl.xlabel(r'$r_{\perp} [h^{-1} Mpc]$')
# pl.ylabel(r'w$(r_{\perp})$')

# pl.xlabel(r'$r [h Mpc^{-1}]$')
# pl.ylabel(r'$\xi(r)$')

pl.xlim(10.**-1., 300.0)
pl.ylim(10.**-4.,  10.)

# pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/SSPOC/monoCorr_fromCat.pdf')

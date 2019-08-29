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

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/VipersMaskWf_rSpaceMonopole_newApproach.dat')
pl.loglog(data[:,0],  (7.5/10.)*10.**2.*data[:,1]*data[:,0], 'k-',         label=r'Monopole, new approach', markersize=2)

pl.loglog(data[:,0],  (7.5/10.)*10.**2.*np.abs(data[:,2])*data[:,0], 'r-', label=r'|Quadrupole|, new approach', markersize=2)

pl.loglog(data[:,0],  (7.5/10.)*10.**2.*data[:,3]*data[:,0], 'g-', label=r'Hexadecapole, new approach', markersize=2)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/VipersMaskWf_rSpaceMonopole.dat')
pl.loglog(data[:,0],  (7.5/10.)*10.**2.*data[:,1]*data[:,0], 'k^', label=r'Regression (Monopole order)', markersize=2)

pl.xlabel(r'$\mathbf{\Delta} [h^{-1} Mpc]$')
pl.ylabel(r'$\Delta W^2(\Delta)$')

pl.xscale('linear')
pl.yscale('linear')

pl.xlim(0.,  600.)
#pl.ylim(-0.005, 0.02)

pl.legend(loc=1)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/31JulyOnwards/windowfunc_rSpaceMultipoles.pdf')

pl.clf()

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/fftlog_hodpk_NoRSD_4096_1.00e-10_1.00e+14_04Sep_monoPkContributions_firstOrder_pk.dat')
pl.loglog(data[:,0],       np.abs(data[:,1]), 'k-', label=r'$|\xi_0(\Delta)W^{2}_{0}(\Delta)|$')
pl.loglog(data[:,0],       np.abs(data[:,2]), 'g-', label=r'$|\xi_2(\Delta)W^{2}_{2}(\Delta)|$')

pl.xlabel(r'$\mathbf{\Delta} [h^{-1} Mpc]$')
#pl.ylabel(r'$|\xi_{\ell}(\Delta)|$')

pl.xscale('log')
pl.yscale('log')

pl.xlim(10.,  400.)
pl.ylim(10.**-9., 10.**-2.)

pl.title(r'Components of cnvld. $\xi_{0}(\Delta)$')

pl.legend(loc=1)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/31JulyOnwards/corrfn_rSpaceMultipoles.pdf')

pl.clf()

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/fftlog_hodpk_NoRSD_4096_1.00e-10_1.00e+14_04Sep_quadPkContributions_firstOrder_pk.dat')
pl.loglog(data[:,0],       np.abs(data[:,1]), 'k--', label=r'$|\xi_0(\Delta)W^{2}_{0}(\Delta)|$')
pl.loglog(data[:,0],       np.abs(data[:,2]), 'r--', label=r'$|(1/5) \ \xi_2(\Delta)W^{2}_{0}(\Delta)|$')
pl.loglog(data[:,0],       np.abs(data[:,3]), 'g-', label=r'$|(2/7) \ \xi_2(\Delta)W^{2}_{2}(\Delta)|$')
pl.loglog(data[:,0],       np.abs(data[:,4]), 'y--', label=r'$|(18/35) \ \xi_2(\Delta)W^{2}_{4}(\Delta)|$')

pl.xlabel(r'$\mathbf{\Delta} [h^{-1} Mpc]$')
#pl.ylabel(r'$|\xi_{\ell}(\Delta)|$')

pl.xscale('log')
pl.yscale('log')

pl.xlim(10.,  400.)
pl.ylim(10.**-9., 10.**-3.)

pl.title(r'Components of cnvld. $\xi_{2}(\Delta)$')

pl.legend(loc=1)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/31JulyOnwards/corrfn_rSpaceMultipoles_quad.pdf')

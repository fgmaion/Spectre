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

TotalVolume = 800.**3.

# Stefano = np.loadtxt('/disk1/mjw/HOD_MockRun/Stefano/power/ps_mocks.txt')
# pl.errorbar(Stefano[:,0], Stefano[:,1]*np.exp(-0.5*(3.*Stefano[:,0])**2.), Stefano[:,2], c='k')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/ConvolvedPk_Multipoles_kbin_0.00500_kmax_0.80_CellSize_4.00_MonoInput.dat')
pl.loglog(data[:,0],  TotalVolume*data[:,1], 'k^', label=r'FFT estimated corr. fn., mono', markersize=2)
pl.loglog(data[:,0],  TotalVolume*data[:,2], 'k^', label=r'FFT estimated corr. fn., quad', markersize=2)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/mixingmatrix_convolvedPk_VipersMask.dat')
pl.loglog(data[:,0],  data[:,1], 'k',   label=r'Kaiser-Lorentz mono')
pl.loglog(data[:,0],  data[:,2], 'r',   label=r'Kaiser-Lorentz quad')
pl.loglog(data[:,0],  data[:,3], 'k--', label=r'conv. mono')
pl.loglog(data[:,0],  data[:,4], 'r--', label=r'conv. quad')
pl.loglog(data[:,0], -data[:,5], 'g--', label=r'|leak, mono to quad|')
pl.loglog(data[:,0],  data[:,4] -data[:,5], 'c--', label=r'leak corrected quad')

pl.xlabel(r'$k [h Mpc^{-1}]$')
pl.ylabel(r'P(k) e$^{-9k^2/2}$')

pl.xlim(10.**-2., 1.0)
pl.ylim(10.**0.,  10.**5.)

pl.legend(loc=1)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/31JulyOnwards/fullCube_pk.pdf')

pl.clf()

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/VipersMaskWf_Multipoles_VipersMaskCube_HOD_KaiserLorentzRSD_clipThreshold_2.0e-01_wfWork_fullCube_kbin_0.00500.dat')
pl.loglog(data[:,0],  data[:,1], 'k', label=r'VIPERS mask,   W${^2}_{0}$(k)')
pl.loglog(data[:,0], -data[:,2], 'r', label=r'VIPERS mask,  -W${^2}_{2}$(k)')

pl.xlabel(r'$k [h Mpc^{-1}]$')
pl.ylabel(r'W$^2$(k)')

pl.xlim(10.**-3., 1.0)
# pl.ylim(10.**1.,  5.*10.**4.)

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/31JulyOnwards/VIPERSmask_monopk.pdf')

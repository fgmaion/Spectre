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

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_GaussianCube_zvel_MonoQuadPk_powLawXi_KaiserRSD_clipThreshold_1.0e+03_fullCube_kbin_0.010_000.dat')
pl.loglog(data[:,0], data[:,1], 'k', label=r'$P_0(k)$, $\delta_0 = 10^3$')
pl.loglog(data[:,0], data[:,2], 'k', label=r'$P_2(k)$, $\delta_0 = 10^3$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_GaussianCube_zvel_MonoQuadPk_powLawXi_KaiserRSD_clipThreshold_1.5e+00_fullCube_kbin_0.010_000.dat')
pl.loglog(data[:,0], data[:,1], 'g', label=r'$P_0(k)$, $\delta_0 = 1.5$')
pl.loglog(data[:,0], data[:,2], 'g', label=r'$P_2(k)$, $\delta_0 = 1.5$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_clippedPrediction_GaussianCube_zvel_MonoQuadPk_powLawXi_KaiserRSD_clipThreshold_1.5e+00_fullCube_HOD_-20.0_beta_0.40_sigma_5.00_kbin_0.010.dat')
pl.loglog(data[:,0], data[:,1], 'k--', label=r'predicted $P_0$(k)')
pl.loglog(data[:,0], data[:,2], 'k--', label=r'predicted $P_2$(k)')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_clippedPrediction_SuppressionOnly_GaussianCube_zvel_MonoQuadPk_powLawXi_KaiserRSD_clipThreshold_1.5e+00_fullCube_HOD_-20.0_beta_0.40_sigma_5.00_kbin_0.010.dat')
pl.loglog(data[:,0], data[:,1], 'b-.', label=r'predicted $P_0$(k), suppressed')
pl.loglog(data[:,0], data[:,2], 'b-.', label=r'predicted $P_2$(k), suppressed')

pl.xlabel(r'$k [h^{-1}$ Mpc$]$')
pl.ylabel('P(k)')
pl.title(r'Kaiser model RSD, $\xi(r) = (5$h$^{-1}$Mpc$/r)^{1.8}$')

pl.yscale('log')
pl.xscale('log')

pl.xlim(10**-2, 1.0)
pl.ylim(10, 10**5)

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Clipping_LikelihoodWork/SpectralDistortion_Gaussian.pdf')

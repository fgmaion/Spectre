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

TotalVolume = 800.**3.

# Stefano = np.loadtxt('/disk1/mjw/HOD_MockRun/Stefano/power/ps_mocks.txt')
# pl.errorbar(Stefano[:,0], Stefano[:,1]*np.exp(-0.5*(3.*Stefano[:,0])**2.), Stefano[:,2], c='k')

MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_0_kbin_0.010_000.dat')
MockCat      = MockCat[:,1]

for i in xrange(1, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_'+str(i)+'_kbin_0.010_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))
   
for i in xrange(10, 100, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_'+str(i)+'_kbin_0.010_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))
   
for i in xrange(100, 500, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_'+str(i)+'_kbin_0.010_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))

kvals        = MockIn[:,0]

Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label='Ensemble Monopole', fmt='', color='k')

# Estimating mean using only cells in which window fn. in non-zero.
MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_IntegralConstraint_SubVolMean_0_kbin_0.010_000.dat')
MockCat      = MockCat[:,1]

for i in xrange(1, 7, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_IntegralConstraint_SubVolMean_'+str(i)+'_kbin_0.010_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))
   
#for i in xrange(10, 500, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_IntegralConstraint_SubVolMean_'+str(i)+'_kbin_0.010_0'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))

kvals        = MockIn[:,0]

Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals + 0.001, Mean, np.sqrt(Var)/np.sqrt(498), label='Ensemble Monopole, icc', fmt='', color='g')


MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_0_kbin_0.010_000.dat')
MockCat      = MockCat[:,2]

for i in xrange(1, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_'+str(i)+'_kbin_0.010_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))
   
for i in xrange(10, 100, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_'+str(i)+'_kbin_0.010_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))
   
for i in xrange(100, 500, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_'+str(i)+'_kbin_0.010_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))

kvals        = MockIn[:,0]

Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label='Ensemble Quadrupole', fmt='', color='r')

# Estimating mean using only cells in which window fn. in non-zero.
MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_IntegralConstraint_SubVolMean_0_kbin_0.010_000.dat')
MockCat      = MockCat[:,2]

for i in xrange(1, 7, 1):
  MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_IntegralConstraint_SubVolMean_'+str(i)+'_kbin_0.010_00'+str(i)+'.dat')
  MockCat  = np.vstack((MockIn[:,2], MockCat))
   
#for i in xrange(10, 500, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_MaskedMultiples_reducedMu_IntegralConstraint_SubVolMean_'+str(i)+'_kbin_0.010_0'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,2], MockCat))

kvals        = MockIn[:,0]

Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals + 0.001, np.abs(Mean), np.sqrt(Var)/np.sqrt(498), label='|Ensemble Quadrupole, icc|', fmt='', color='b')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/mixingmatrix_convolvedPk_VipersMask.dat')
pl.loglog(data[:,0],  data[:,3], 'k-', label=r'conv. mono')
pl.loglog(data[:,0],  data[:,4], 'r-', label=r'conv. quad')

pl.loglog(data[:,0],  data[:,1], 'k-',   label=r'Kaiser-Lorentz mono')
pl.loglog(data[:,0],  data[:,2], 'r-',   label=r'Kaiser-Lorentz quad')

pl.xlabel(r'$k [h Mpc^{-1}]$')
pl.ylabel(r'P(k) e$^{-9k^2/2}$')

pl.xlim(10.**-2., 0.6)
pl.ylim(10.**1.,  10.**6.)

pl.legend(loc=1)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/31JulyOnwards/maskedMultipoles.pdf')

pl.clf()


window = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/VipersMaskWf_Multipoles_VipersMask_HOD_KlRSD_clipThreshold_2.0e-01_IntConCor_apparentMean_kbin_0.00500.dat')
pl.loglog(window[:,0],  window[:,1]/window[0,1], 'k-', markersize=2, label=r'+ve monopole')
pl.loglog(window[:,0], -window[:,2]/window[0,1], 'r-', markersize=2, label=r'-ve quadrupole')

pl.xlabel(r'$k [h Mpc^{-1}]$')
pl.ylabel(r'W$^2$(k)')
pl.title('Vipers mask multipoles')

pl.xlim(5.*10.**-3., 0.6)

pl.legend(loc=1)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/31JulyOnwards/maskedMultipoles_window.pdf')

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

MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_fullCube_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_000.dat')
MockCat      = MockCat[:,1]

for i in xrange(1, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_fullCube_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))
   
for i in xrange(10, 100, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_fullCube_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
MonoVar          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(MonoVar), label='Ensemble Monopole', fmt='', color='g')

MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_fullCube_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_000.dat')
MockCat      = MockCat[:,2]

for i in xrange(1, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_fullCube_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))
   
for i in xrange(10, 100, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_fullCube_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
QuadVar      = np.var(MockCat,  axis=0)

#pl.errorbar(kvals - 0.002, Mean, np.sqrt(QuadVar), label='Ensemble Quadrupole', fmt='', color='b')

MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_fullCube_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_000.dat')
# pl.errorbar(MockCat[:,0] + 0.001, MockCat[:,1], MockCat[:,5], label='Monopole,   predicted errors', fmt='', color='k')
# pl.errorbar(MockCat[:,0] - 0.001, MockCat[:,2], MockCat[:,6], label='Quadrupole, predicted errors', fmt='', color='r')

pl.loglog(MockCat[:,0], MockCat[:,5], 'k-', label='predicted monopole   error')
pl.loglog(MockCat[:,0], MockCat[:,6], 'r-', label='predicted quadrupole error')

pl.loglog(kvals, np.sqrt(MonoVar), 'g-', label='monopole error from mocks')
pl.loglog(kvals, np.sqrt(QuadVar), 'b-', label='quadrupole error from mocks')

pl.xscale('log')
pl.yscale('log', nonposy='clip')

pl.xlabel(r'$k [h Mpc^{-1}]$')
# pl.ylabel(r'P(k) e$^{-9k^2/2}$')

pl.xlim(10.**-2., 0.4)
pl.ylim(10.**0.,  10.**5.)

pl.legend(loc=1)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/31JulyOnwards/quadrupoleError.pdf')

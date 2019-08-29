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

MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Today/Multipoles_250_fullCube_lnNorm_noRSD_CellSize_1.00_kbin_0.005_001.dat')
MockCat      = MockCat[:,1]

for i in xrange(1, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Today/Multipoles_250_fullCube_lnNorm_noRSD_CellSize_1.00_kbin_0.005_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))
   
for i in xrange(20, 50, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Today/Multipoles_250_fullCube_lnNorm_noRSD_CellSize_1.00_kbin_0.005_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label=r'Ensemble Monopole, 250 h$^{-1}$Mpc', fmt='', color='k')


MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Today/Multipoles_500_fullCube_lnNorm_noRSD_CellSize_2.00_kbin_0.005_001.dat')
MockCat      = MockCat[:,1]

for i in xrange(1, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Today/Multipoles_500_fullCube_lnNorm_noRSD_CellSize_2.00_kbin_0.005_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))
   
for i in xrange(20, 50, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Today/Multipoles_500_fullCube_lnNorm_noRSD_CellSize_2.00_kbin_0.005_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label=r'Ensemble Monopole, 500 h$^{-1}$Mpc', fmt='', color='g')


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Xi/fftlog_hodpk_4096_1.00e-10_1.00e+14_01Oct_lnnorm_pk.dat')
pl.loglog(data[:,0],  data[:,2],    'r-', label='ln normal')
pl.loglog(data[:,0],  data[:,1],    'k-', label='underlying Gaussian P(k)')

pl.xscale('log')
pl.yscale('log')

pl.xlabel(r'$k [h Mpc^{-1}]$')
pl.ylabel(r'P(k) e$^{-9k^2/2}$')

pl.xlim(4.*10.**-2.,  1.0)
pl.ylim(5.*10.**1.,  2.*10.**4.)

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/31JulyOnwards/lnNorm.pdf')

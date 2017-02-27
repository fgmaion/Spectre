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

ydata   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_yvel_clipThreshold_1.0e+03_fullCube_kbin_0.010_000.dat')
pl.loglog(ydata[:,0], ydata[:,1], 'g^', label='$P_0(k)$', markersize=4)

xdata   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_xvel_clipThreshold_1.0e+03_fullCube_kbin_0.010_000.dat')
pl.loglog(xdata[:,0], xdata[:,2], 'k^', label='$P_2(k)$, x vel. ', markersize=4)

pl.loglog(ydata[:,0], ydata[:,2], 'y^', label='$P_2(k)$, y vel.', markersize=4)

data   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_zvel_clipThreshold_1.0e+03_fullCube_kbin_0.010_000.dat')

pl.loglog(data[:,0], data[:,2], 'r^', label='$P_2(k)$, z vel.', markersize=4)


MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_clipThreshold_1.0e+03_subVol_0_kbin_0.010_000.dat')
MockCat      = MockCat[:,2]

for i in xrange(1, 8, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_clipThreshold_1.0e+03_subVol_'+str(i)+'_kbin_0.010_000.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(7), label='$P_2(k)$, x vel. sub volumes', fmt='k')


theory = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/Pk_hod_20.0.dat')
pl.loglog(theory[:,0], 1.37*theory[:,1], 'k')

theory = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_HOD_-20.0_beta_0.45_velDispersion_3.00.dat')
pl.loglog(theory[:,0], theory[:,2], 'k--', label=r'Kaiser-Lorentz, $\beta=0.45$, $\sigma=3.0$')

theory = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_HOD_-20.0_beta_0.45_velDispersion_3.00.dat')
pl.loglog(theory[:,0], theory[:,3], 'k--')

pl.xlabel(r'$k [h^{-1} Mpc]$')
pl.ylabel('P(k)')

pl.yscale('log')
pl.xscale('log')

pl.xlim(10**-2, 0.5)
pl.ylim(100, 10**5)

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Multipoles_Cube_clippingDraft.pdf')

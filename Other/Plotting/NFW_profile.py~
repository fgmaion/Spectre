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
fig_width  = 2.*fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize':14,
          'text.fontsize':8,
          'legend.fontsize':14,
          'xtick.labelsize':12.0,
          'ytick.labelsize':12.0,
          'figure.figsize':fig_size,
          'font.family': 'serif'}

pl.rcParams.update(params)
pl.clf()
pl.figure()
fig = pl.figure()
axes = pl.Axes(fig, [.2, .2, .7, .7])
fig.add_axes(axes)
axes.yaxis.set_major_formatter(formatter)


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_500_ukm_0.dat')
pl.semilogx(data[:,0], data[:,1]**0.5, 'k^', markersize=2)
pl.semilogx(data[:,0], data[:,3], 'r')

pl.xlabel(r'k [h$^{-1}$ Mpc]')
pl.ylabel('u(k|m)')

pl.yscale('linear')
pl.xscale('log')

pl.ylim(10.**-1., 1.)

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/NFW_window.pdf')

pl.clf()

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_fftlog_multipoles_urm.dat')
pl.loglog(data[:,0], data[:,1], 'k', label='u(r|m)')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_xi_multipoles_mean.dat')
pl.loglog(data[:,0], data[:,1]/(2.*10.**5.), 'r^', markersize=2, label=r'$\xi^{1-halo}$(r)')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_fftlog_multipoles_xi.dat')
pl.loglog(data[:,0], data[:,1], 'r', label=r'$\xi^{1-halo}$(r)')

pl.xlabel(r'r [h$^{-1}$ Mpc]')
#pl.ylabel(r'$\xi^{1-halo}$(r)')

pl.yscale('log')
pl.xscale('log')

pl.xlim(10.**-1., 100.0)
pl.ylim(10.**-7.,   0.2)

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/NFW_profile_xi0.pdf')

pl.clf()

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_500_pk_1.dat')
# pl.loglog(data[:,0],  data[:,1], 'k^', markersize=2)
pl.loglog(data[:,0],  data[:,3], '--k')
pl.loglog(data[:,0], -data[:,4], '--r')
pl.loglog(data[:,0],  data[:,5], 'k')
pl.loglog(data[:,0], -data[:,6], 'r')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_500_quad_1.dat')
# pl.loglog(data[:,0], -data[:,2], 'k^', markersize=2)

MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_500_pk_1.dat')
MockCat      = MockCat[:,1]

for i in xrange(2, 22, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_500_pk_'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(20), label='Ensemble Monopole', fmt='', color='k')


MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_500_quad_1.dat')
MockCat      = MockCat[:,2]

for i in xrange(2, 22, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_500_quad_'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, -Mean, np.sqrt(Var)/np.sqrt(20), label='Ensemble Quadrupole', fmt='', color='r')


pl.xlabel(r'k [h Mpc$^{-1}$]')
pl.ylabel('P(k)')

pl.yscale('log')
pl.xscale('log')

pl.xlim(10.**-1., 1.0)
pl.ylim(10.**3., 2.*10.**4.)

pl.legend(loc=1)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/NFW_halomodel_pk.pdf')

pl.clf()

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_monoxi.dat')
pl.loglog(data[:,0],         data[:,1],  'k-')
pl.loglog(data[:,0],  np.abs(data[:,2]), 'r-')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_xi_1_multipoles.dat')
#pl.loglog(data[:,0],         data[:,1],  'k^', markersize=6)
#pl.loglog(data[:,0],  np.abs(data[:,2]), 'r^', markersize=6)


MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_xi_1_multipoles.dat')
MockCat      = MockCat[:,1]

for i in xrange(2, 50, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_xi_'+str(i)+'_multipoles.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(48), label='Ensemble Monopole', fmt='', color='k')


MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_xi_1_multipoles.dat')
MockCat      = MockCat[:,2]

for i in xrange(2, 50, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/NFW_profile_xi_'+str(i)+'_multipoles.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, np.abs(Mean), np.sqrt(Var)/np.sqrt(48), label='Ensemble Monopole', fmt='', color='r')


pl.xlabel(r'r [h$^{-1}$ Mpc]')
pl.ylabel(r'$\xi$(r)')

pl.yscale('log')
pl.xscale('log')

pl.xlim(1., 50.)
pl.ylim(10.**-3., 120.)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/NFW_halomodel_xi.pdf')

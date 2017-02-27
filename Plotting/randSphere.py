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

'''
data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/randSphere_xi_multipoles.dat')
pl.plot(data[:,0],     data[:,1],    'k^', markersize=2)

pl.xlabel(r'$\Delta$ [h$^{-1}$ Mpc]')
# pl.ylabel('u(k|m)')

pl.yscale('log')
pl.xscale('log')

pl.xlim(0., 450.)
pl.ylim(10., 10.**6.)

pl.legend(loc=3)

nbar = 30000.*(500.)**-3.

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/spherical_RandDistribution_monoxi.dat')
pl.plot(data[:,0], 0.33*0.25*data[:,0]**3.*data[:,1], 'r')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/randSphere_xi.pdf')

pl.clf()
'''
'''
data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/spherical_RandDistribution_cnvldpk.dat')
pl.plot(data[:,0],     data[:,1],    'k')
pl.plot(data[:,0],     data[:,2],    'r')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/SphericalRands_mask_500_1.dat')
pl.plot(data[:,0],     data[:,1],    'k^', markersize=2)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/SphericalRands_mask_500_2.dat')
MockCat      = MockCat[:,1]

for i in xrange(3, 300, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/SphericalRands_mask_500_'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(298), label='cnvld. monopole', fmt='', color='r')

pl.xlabel(r'k [h Mpc$^{-1}$]')
pl.ylabel('P(k)')

pl.yscale('log')
pl.xscale('log')

pl.xlim(0.01,    0.3)
pl.ylim(10.**3., 4.*10.**4.)

pl.legend(loc=3)

#nbar = 30000.*(500.)**-3.

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/randSphere_cnvldpk.pdf')

pl.clf()
'''
'''
data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/spherical_tophat.dat')
pl.plot(data[:,0], data[:,1],    'k-')

pl.xlabel(r'k [h Mpc$^{-1}$]')

pl.yscale('log')
pl.xscale('log')

# pl.xlim(0., 450.)
# pl.ylim(10., 10.**6.)

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/spherical_tophat_pk.pdf')
'''
pl.clf()

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/randSphere_xi_multipoles.dat')
pl.plot(data[:,0],     data[:,1],    'k^', markersize=2)


data       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/spherical_RandDistribution_monoxi.dat')

data[:,1] /= data[0, 1]

Vsphere    = (4./3.)*pi*(200.**3.)

nbar       = 30000./Vsphere

print nbar

data[:,1] *= Vsphere

pl.plot(data[:,0], 0.5*nbar*nbar*4.*pi*data[:,0]**3.*np.log10(1.01)*data[:,1]*np.log(10.),    'k-')

pl.xlabel(r'$\Delta$ [h$^{-1}$ Mpc]')

pl.yscale('log')
pl.xscale('log')

# pl.xlim(1., 450.)
# pl.ylim(0.1, 1.5)

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/spherical_tophat_autocorrelation.pdf')

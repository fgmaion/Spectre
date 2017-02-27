import pylab as pl

from   scipy.special import eval_legendre as Leg

from   matplotlib.colors import LogNorm

import math,sys

kBinNumb  = 13
muBinNumb = 37

data      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/VipersMask_2Dpk_W2k.dat')

theta     = np.arccos(data[:muBinNumb*kBinNumb,0]).reshape(muBinNumb, kBinNumb)

kVals     = data[:muBinNumb*kBinNumb, 1].reshape(muBinNumb, kBinNumb)
    
PkArray   = data[:muBinNumb*kBinNumb, 2].reshape(muBinNumb, kBinNumb)


# plot the Temperature profile
fig, ax   = pl.subplots(subplot_kw=dict(projection='polar'))

#xVals    = np.random.multivariate_normal(np.array([2., 3.]), np.array([[2, 3], [4,5]]), 36*50)[:,0]

pax = ax.pcolormesh(theta, kVals, PkArray, norm=LogNorm(vmin=10.**0., vmax=data[:,2].max()))

ax.pcolormesh(np.pi*np.ones([muBinNumb, kBinNumb]) - theta, kVals, PkArray, norm=LogNorm(vmin=10.**0., vmax=data[:,2].max()))

ax.pcolormesh(np.pi*np.ones([muBinNumb, kBinNumb]) + theta, kVals, PkArray, norm=LogNorm(vmin=10.**0., vmax=data[:,2].max()))

ax.pcolormesh(-theta, kVals, PkArray, norm=LogNorm(vmin=10.**0., vmax=data[:,2].max()))


ax.set_theta_zero_location("N") # 'north' location for theta=0

ax.set_theta_direction(-1)      # angles increase clockwise

ax.set_xticklabels([str(np.cos(0.0)), str(np.cos(np.pi/4.)), str(np.cos(np.pi/2.)), str(np.cos(3.*np.pi/4.)), str(np.cos(np.pi)), '', '', ''])

fig.colorbar(pax)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/polar2Dpk_W2k.pdf')

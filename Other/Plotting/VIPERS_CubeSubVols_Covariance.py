# import copy
# import scipy.interpolate

from   matplotlib.colors import LogNorm

# cmap = copy.copy(matplotlib.cm.jet)

# cmap.set_bad('w',1.)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Covariance/zCube_clipThreshold_1.0e+03_subVols_Covariance.dat')

plt.imshow(data, origin='lower', cmap='YlGnBu', interpolation='nearest', norm = LogNorm(vmin=10**0, vmax=5*10**8))

# print data.max()
# plt.imshow(data, origin='lower', cmap='YlGnBu', interpolation='nearest')

pl.xticks([10, 20, 30, 40, 50, 60, 70], [str(10*0.01), str(20*0.01), str(30*0.01),  str(40*0.01), str(50*0.01), str(60*0.01), str(70*0.01)])

plt.colorbar()

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/zCube_subVols_clipThreshold_1.0e+03_Multipole_Covariance.pdf', bbox_inches='tight', pad_inches=0.5)

pl.clf()


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Covariance/zCube_clipThreshold_1.0e+03_subVols_InverseCovariance.dat')

plt.imshow(data, origin='lower', cmap='YlGnBu', interpolation='nearest', norm = LogNorm(vmin=10**-14, vmax=10**-8))

pl.xticks([10, 20, 30, 40, 50, 60, 70], [str(10*0.01), str(20*0.01), str(30*0.01),  str(40*0.01), str(50*0.01), str(60*0.01), str(70*0.01)])

print data.max(), data.min()
# plt.imshow(data, origin='lower', cmap='YlGnBu', interpolation='nearest')

plt.colorbar()

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/zCube_subVols_clipThreshold_1.0e+03_Multipole_InverseCovariance.pdf', bbox_inches='tight', pad_inches=0.5)

pl.clf()





data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Covariance/zCube_clipThreshold_1.0e+03_subVols_MeanMultipoles.dat')

# pl.loglog(data[:,0], data[:,1])
pl.loglog(data[:,0], data[:,2], 'k^')


MockCat     = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_clipThreshold_1.0e+03_subVol_0_kbin_0.010_000.dat')
MockCat     = MockCat[:,2]

for i in xrange(1, 8, 1):
    MockIn = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_clipThreshold_1.0e+03_subVol_'+str(i)+'_kbin_0.010_000.dat')
    # pl.loglog(single[:,0], single[:,1], 'm^')
    # pl.loglog(single[:,0], single[:,2])

    MockCat  = np.vstack((MockCat, MockIn[:,2]))
    
Mean         = np.mean(MockCat, axis=0)

MockCat     -= Mean

Var          = np.var(MockCat,  axis=0)

print Var

'''
single = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_zCube_clipThreshold_1.0e+03_subVol_1_kbin_0.010_000.dat')
# pl.loglog(single[:,0], single[:,1], 'g^')
pl.loglog(single[:,0], single[:,2] - data[:,2], 'b^')
'''


pl.savefig('/disk1/mjw/HOD_MockRun/Plots/zCube_subVols_clipThreshold_1.0e+03_Multipole_MeanMultipoles.pdf', bbox_inches='tight', pad_inches=0.5)

data = np.loadtxt('/disk1/mjw/VIPERS_ValueAddedHOD/mocks_W4_v1.2/mock_W4_001_ALLINFO.cat', usecols=[0,1,2])

ra   = data[:,1]
dec  = data[:,2]

print dec.min(), dec.max()
print  ra.min(),  ra.max()

print dec.mean(), ra.mean()

Area = (dec.max() - dec.min())*(ra.max() - ra.min())

print Area

pl.scatter(dec[::10], ra[::10])
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/AngularCoverageCheck_W4.pdf')

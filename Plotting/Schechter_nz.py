data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/Schechter_HODMocks_MockAvgNz_chiSliced_15.0_W1_22.00_50.0.dat')

pl.plot(data[:,0], data[:,3])
pl.plot(data[:,0], data[:,5])
pl.plot(data[:,0], data[:,6])

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Schechter_nz.pdf')

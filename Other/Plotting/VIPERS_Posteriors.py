data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_clipThreshold_1.0e+03_subVol_betaPosterior_LowRes_kmax_0.50.dat')
pl.plot(data[:,0], data[:,1], 'k^', label=r'$k_{max} = 0.5$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_clipThreshold_1.0e+03_subVol_betaPosterior_LowRes_kmax_0.10.dat')
pl.plot(data[:,0], data[:,1], 'y^', label=r'$k_{max} = 0.1$')

pl.axvline(x=0.542, ymin=0, ymax=1, c='r')

'''
data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_clipThreshold_1.0e+03_subVol_betaPosterior_LowRes_kmax_0.30.dat')
pl.plot(data[:,0], data[:,1], 'r')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_clipThreshold_1.0e+03_subVol_betaPosterior_LowRes_kmax_0.10.dat')
pl.plot(data[:,0], data[:,1], 'g')
'''

pl.xlabel(r'$\beta$')

pl.legend(loc=3)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/zCube_clipThreshold_1.0e+03_subVol_betaPosterior_LowRes.pdf')



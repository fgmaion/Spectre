'''

MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_VIPERSparent_StefanoBasis_Cic_GaussSmoothNz_50.0_SkinDepth_5.0_VolLim_22.00_Mesh_2.00_kInterval_0.01_001.dat')
MockCat      = MockCat[:,2]

for i in xrange(2, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_VIPERSparent_StefanoBasis_Cic_GaussSmoothNz_50.0_SkinDepth_5.0_VolLim_22.00_Mesh_2.00_kInterval_0.01_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))

for i in xrange(10, 27, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_VIPERSparent_StefanoBasis_Cic_GaussSmoothNz_50.0_SkinDepth_5.0_VolLim_22.00_Mesh_2.00_kInterval_0.01_0' +str(i)+ '.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))

kvals        = MockIn[:,0]

Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals[1:], Mean[1:], np.sqrt(Var[1:])/np.sqrt(25), c='g', fmt='--', label='dx=2., (12, 0.15, 0.1)')
'''
data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_VIPERSparent_StefanoBasis_Cic_GaussSmoothNz_50.0_SkinDepth_5.0_VolLim_22.00_Mesh_4.00_kInterval_0.01_001.dat')
pl.loglog(data[:,0], data[:,2], 'r^')

#Stefano      = np.loadtxt('/disk1/mjw/HOD_MockRun/Stefano/pk_06z09_H2/ps_1.txt') 
#pl.loglog(Stefano[:,0], Stefano[:,1], 'k^', label='Stefano')

#Stefano      = np.loadtxt('/disk1/mjw/HOD_MockRun/Stefano/pk_06z09_H4/ps_mocks.txt') 
#pl.errorbar(Stefano[:,0], Stefano[:,1], Stefano[:,2], c='k', fmt='--', label='dx=4., (12, 0.15, 0.1)')

pl.yscale('log')
pl.xscale('log')

pl.xlim([0.01, 0.8])
#pl.ylim([2*10**2, 5*10**4])

pl.ylabel(r'$P(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

#leg = pl.legend(loc=3, ncol=1, prop = FontProperties(size = '10'))
#leg.draw_frame(False)

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/CloudInCell.pdf')

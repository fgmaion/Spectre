xSlice       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/VIPERSparent_xSlice.dat')
ySlice       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/VIPERSparent_ySlice.dat')
zSlice       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/VIPERSparent_zSlice.dat')

abscissa     = np.arange(10**-3, 10**0, 10**-3)

SphereRadius = 250.

y            = abscissa*SphereRadius

analytic     = 3.*(y**-3)*(np.sin(y) - y*np.cos(y))

#pl.loglog(data[:,0], data[:,1],  'k')
#pl.loglog(abscissa, analytic**2, 'y')

pl.loglog(xSlice[:,1], xSlice[:,2], 'g')
pl.loglog(ySlice[:,1], ySlice[:,2], 'y')
pl.loglog(zSlice[:,1], zSlice[:,2], 'r^')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/Spherical/VIPERSparent_Windowfn.pdf')

##### 

data          = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_Spherical_Jenkins1.0.dat')
pl.loglog(data[:,0], data[:,2], 'b^', label='Measured')

theory = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/cambExtendedPk_hod_20.0.dat')
pl.loglog(theory[:,0], theory[:,1], 'y', label='HOD input')

pl.xlim([0.002, 0.4])
pl.ylim([10**3, 4*10**4])

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

leg = pl.legend(loc =2, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.ylabel(r'$P(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/Spherical/SphericalWf_Pk.pdf')


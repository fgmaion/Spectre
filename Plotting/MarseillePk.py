Mock001nz_03   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_VIPERSparent_Mock001nz_0.03.dat')
#pl.loglog(Mock001nz_03[:,0], Mock001nz_03[:,2], 'g^', label='Mock 001 nz, dz = 0.03')

Mock001nz_06   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_VIPERSparent_Mock001nz_0.06.dat')
#pl.loglog(Mock001nz_06[:,0], Mock001nz_06[:,2], 'b^', label='Mock 001 nz, dz = 0.06')

MockAvgnz_03   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_VIPERSparent_MockAvgnz_0.03.dat')
#pl.loglog(MockAvgnz_03[:,0], MockAvgnz_03[:,2], 'c^', label='Mock Avg nz, dz = 0.03')

MockAvgnz_06   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_VIPERSparent_MockAvgnz_0.06.dat')
#pl.loglog(MockAvgnz_06[:,0], MockAvgnz_06[:,2], 'm^', label='Mock Avg nz, dz = 0.06')

theory         = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/cambExtendedPk_hod_20.0.dat')
pl.loglog(theory[:,0], theory[:,1], 'y')

conv           = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Pk_Convolved_ICCcorrectedVIPERSparent.dat')
#pl.loglog(conv[:,0], conv[:,2], 'k^')

minChi2nz      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_VIPERSparent_Mock001_minChi2nz_0.03.dat')
pl.loglog(minChi2nz[:,0], minChi2nz[:,2], 'k^', label='min Chi^2 nz')

SphWf  = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSpherical/midK_W2q_VIPERSparent.dat')

ConvPkZeroPoint = 2.502092e+04
WfZeroPoint     = 1.836298*(10.**-5.)

IccCorrected    = conv[:,2] - (ConvPkZeroPoint/WfZeroPoint)*SphWf[:,1]
#pl.plot(conv[:,0], IccCorrected , 'r^')

pl.xlim([0.002, 0.4])
pl.ylim([10**3, 2*10**5])

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

leg = pl.legend(loc=3, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.ylabel(r'$P(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.title('Measured, parent 001, $0.5<z<0.9$')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/MarseillePk_minChi2nz.pdf')



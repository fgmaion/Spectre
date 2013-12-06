data   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_AnisoGauss_Jenkins1.0.dat')
pl.loglog(data[:,0], data[:,2], 'b^', label='Measured, volume limited sample.')

theory = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/cambExtendedPk_hod_20.0.dat')
pl.loglog(theory[:,0], theory[:,1], 'y', label='HOD input')

conv = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Pk_ConvolvedAnisoGauss_Jenkins1.0.dat')
pl.loglog(conv[:,0], conv[:,2], 'g^', label='convolved.')

pl.xlim([0.001, 0.3])
pl.ylim([10**3, 5*10**4])

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

leg = pl.legend(loc =2, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.ylabel(r'$P(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/MarseillePkAnisoGauss.pdf')



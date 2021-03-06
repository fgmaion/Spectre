data   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_FullCube.dat')
pl.loglog(data[:,0], data[:,2], 'r^', label='real space')

zdata   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_zFullCube.dat')
pl.loglog(zdata[:,0], zdata[:,2], 'g^', label='redshift space')

theory = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/cambExtendedPk_hod_20.0.dat')
pl.loglog(theory[:,0], theory[:,1], 'k')

matplotlib.pyplot.axhline(y=1./0.005481, c='y', xmin=0, xmax=1, label='expected shot noise')

leg = pl.legend(loc =1, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

pl.ylabel(r'$P(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/shotNoiseCheck.pdf')

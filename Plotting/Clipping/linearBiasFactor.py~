# del2k
#measured     = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_FullCube.dat')
#pl.loglog(measured[:, 0], measured[:, 1], 'b^', label = 'measured biased')

zspace        = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_zFullCube.dat')
pl.loglog(zspace[:, 0], zspace[:, 1], 'b^', label = 'z space')

# Pk
#linear        = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/LinearPk.dat')
#linear[:,1]  *= (0.5/np.pi**2)*linear[:,0]**3
#pl.loglog(linear[:,0],  2.9*linear[:,1], 'y', label='biased linear theory')

# clipping suppresed
suppressed     = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_zFullCube_clipped.dat')
pl.loglog(suppressed[:, 0], 2.6*suppressed[:, 1], 'g^', label = 'clipping suppressed, rescaled by suppression factor (2.6)')

leg = pl.legend(loc =2, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.xlim([0.01, 0.4])

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Clipping/linearBias.pdf')


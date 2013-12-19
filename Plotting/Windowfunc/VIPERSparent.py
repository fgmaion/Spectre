data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSpherical/midK_W2q_VIPERSparent.dat')
pl.loglog(data[:,0], data[:,1], 'g', label='spherical average')

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

leg = pl.legend(loc=3, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.ylabel(r'$W^2(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/MarseilleWfSphAvg.pdf')

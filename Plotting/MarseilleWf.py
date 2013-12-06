a = 300.
b = 400.
c = 500.

abscissa = np.arange(0.001, 0.3, 0.001)

xtheory  = np.exp(-(abscissa*a)**2.)
ytheory  = np.exp(-(abscissa*b)**2.)
ztheory  = np.exp(-(abscissa*c)**2.)

pl.loglog(abscissa, xtheory, 'b')
pl.loglog(abscissa, ytheory, 'g')
pl.loglog(abscissa, ztheory, 'y')

xdata    = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/AnisoGauss_Jenkins1.0_xSlice.dat')
pl.loglog(xdata[:,1], xdata[:,2], 'b^', label='x direction')

ydata    = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/AnisoGauss_Jenkins1.0_ySlice.dat')
pl.loglog(ydata[:,1], ydata[:,2], 'g^', label='y direction')

zdata    = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/AnisoGauss_Jenkins1.0_zSlice.dat')
pl.loglog(zdata[:,1], zdata[:,2], 'y^', label='z direction')

pl.xlim([0.001, 0.3])
pl.ylim([0.0001, 1.0])

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

leg = pl.legend(loc =2, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.ylabel(r'$W^2(k_i)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/MarseilleWfAnisoGauss.pdf')

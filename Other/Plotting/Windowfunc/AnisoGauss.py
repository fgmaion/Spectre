data          = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Pk_ConvolvedAnisoGauss_Jenkins1.0.dat')
pl.plot(data[:,0], data[:,2], 'b^', label='Measured')

theory = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/cambExtendedPk_hod_20.0.dat')
pl.plot(theory[:,0], theory[:,1], 'y', label='HOD input')

icc = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Pk_ConvolvedAnisoGauss_Jenkins1.0ICC.dat')
pl.plot(icc[:,0]+0.0001, icc[:,2], 'g^', label='icc corrected')

pl.xlim([0.002, 0.08])
pl.ylim([10**3, 4*10**4])

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

leg = pl.legend(loc=1, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.ylabel(r'$P(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ConvolvedPk/AnisoGauss_Pk.pdf')

pl.clf()

xSlice = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/AnisoGauss_Jenkins1.0_xSlice.dat')
ySlice = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/AnisoGauss_Jenkins1.0_ySlice.dat')
zSlice = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/AnisoGauss_Jenkins1.0_zSlice.dat')

pl.loglog(xSlice[:,1], xSlice[:,2], 'g^')
pl.loglog(ySlice[:,1], ySlice[:,2], 'r^')
pl.loglog(zSlice[:,1], zSlice[:,2], 'y^')

absicssa = np.arange(0.001, 1.0, 0.001)

a = 80.
b = 100.
c = 120.

xTheory  = np.exp(-a*a*absicssa**2.)
yTheory  = np.exp(-b*b*absicssa**2.)
zTheory  = np.exp(-c*c*absicssa**2.)

pl.loglog(absicssa, xTheory, 'g')
pl.loglog(absicssa, yTheory, 'r')
pl.loglog(absicssa, zTheory, 'y')

pl.ylim(10.**-24., 1.)

pl.xlabel('$k$')
pl.ylabel('$W^2(k_i)$')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/WindowfuncSlices/AnisoGauss.pdf')

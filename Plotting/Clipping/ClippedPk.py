NormalPk     = np.loadtxt('Scripts/Normal_PkTheory_KaiserLorentzian_0.6mu0.8.dat')
SuppressedPk = np.loadtxt('Scripts/Suppressed_PkTheory_KaiserLorentzian_0.6mu0.8.dat')
DistortedPk  = np.loadtxt('Scripts/Distorted_PkTheory_KaiserLorentzian_0.6mu0.8.dat')

pl.loglog(NormalPk[:, 1], NormalPk[:, 2], 'b', label = r'$ 0.6 < \mu < 0.8$')
pl.loglog(SuppressedPk[:,1], SuppressedPk[:,2], 'r', label='Clipping suppressed')
pl.loglog(DistortedPk[:,1], DistortedPk[:,2], 'g--', label='Clipping distorted')

leg = pl.legend(loc =1, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

#pl.xlim([0.00125, 0.154])
#pl.ylim([10**-8, 10**-5])

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

pl.ylabel(r'$P(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.title('Real space, linear P(k) clipped, u0 = 0.69')
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ClippedPk_0.6mu0.8.pdf')

########

pl.clf()

zNormalPk     = np.loadtxt('Scripts/zNormal_PkTheory_KaiserLorentzian_0.6mu0.8.dat')
zDistortedPk  = np.loadtxt('Scripts/zDistorted_PkTheory_KaiserLorentzian_0.6mu0.8.dat')

pl.loglog(zNormalPk[:, 1], zNormalPk[:, 2], 'r', label = r'$ 0.6 < \mu < 0.8$')
pl.loglog(zDistortedPk[:,1], zDistortedPk[:,2], 'g--', label='Clipping distorted')

leg = pl.legend(loc =1, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

#xx, locs = plt.xticks()
#ll = ['%.3f' % a for a in xx]
#plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

#pl.ylabel(r'$P(k)$', fontsize = '10')
#pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

#pl.title('Redshift space (Kaiser $[\Omega_M = 0.3, \ b=1]$ + Lorentzian $[300 kms^{-1}]$), linear P(k) clipped, u0 = 0.69')
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ClippedPk_KaiserLorentzian_0.6mu0.8.pdf')
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/zClippedPk_0.6mu0.8.pdf')

pl.clf()

pl.yscale('log')
pl.xscale('log')

# Real space. 
# Mono = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_Cube_Jenkins1.0_kInterval_0.01_000.dat')
# pl.loglog(Mono[:,0], Mono[:,2], 'g^')

# Multipoles = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Observed_Multipoles_zCube_Jenkins1.0_kbin_0.010_000.dat')
# pl.loglog(Multipoles[:,0], Multipoles[:,2], 'k^', markersize=4)

# MonoNorm   = Multipoles[:,4]

# Redshift space.
Mono = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_zCube_xvel_clipThreshold_1.0e+03_fullCube_kInterval_0.01_000.dat')
pl.loglog(Mono[:,0], Mono[:,2], 'k^', markersize=4, label='Unclipped Monopole')
pl.loglog(Mono[:,0], Mono[:,3], 'k^', markersize=4, label='Unclipped Monopole')

A11Sq = 1.55

Multipoles = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_fullCube_HOD_KlRSD_clipThreshold_2.0e-01_IntConCor_apparentMeanCorr_kbin_0.005_000.dat')
pl.loglog(Multipoles[:,0], A11Sq*Multipoles[:,1], 'k.', markersize=4, label='Clipped Monopole')
pl.loglog(Multipoles[:,0], A11Sq*Multipoles[:,2], 'k.', markersize=4)

# MB = -20.0
linearBias = 1.495903

beta       = 0.54  

KaiserMono =  1. + (2./3.)*beta + 0.2*beta*beta
KaiserQuad = (4./3.)*beta + (4./7.)*beta*beta

HOD = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/outdated_cosmology/Pk_hod_20.0.dat')

# pl.loglog(HOD[:,0], KaiserMono*HOD[:,1], 'y')
# pl.loglog(HOD[:,0], KaiserQuad*HOD[:,1], 'y')
# pl.loglog(HOD[:,0],            HOD[:,1], 'k')
'''
kaiserLorentzNorm = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_HOD_-20.0_beta_0.54_velDispersion_3.50.dat')
pl.loglog(kaiserLorentzNorm[:,0], kaiserLorentzNorm[:,2], 'r', label=r'HOD, $\sigma = 3.5$')

# kaiserLorentzNorm = kaiserLorentzNorm[:,2]

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_HOD_-20.0_beta_0.54_velDispersion_3.00.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,2], 'k', label=r'HOD, $\sigma = 3.0$')
#pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'k')

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_HOD_-20.0_beta_0.54_velDispersion_4.00.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,2], 'g', label=r'HOD, $\sigma = 4.0$')
#pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'g')

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_Linear_-20.0_beta_0.54_velDispersion_3.00.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,2], 'k--', label=r'linear, $\sigma = 3.0$')
#pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'k--')

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_Linear_-20.0_beta_0.54_velDispersion_2.75.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,2], 'r--', label=r'linear, $\sigma = 2.75$')
#pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'r--')

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_Linear_-20.0_beta_0.54_velDispersion_2.50.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,2], 'g--', label=r'linear, $\sigma = 2.5$')
#pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'g--')

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_Linear_-20.0_beta_0.54_velDispersion_2.00.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,2], 'y--', label=r'linear, $\sigma = 2.00$')
#pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'g--')
'''
'''
kaiserLorentzNorm = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_HOD_-20.0_beta_0.54_velDispersion_3.50.dat')
# kaiserLorentzNorm = kaiserLorentzNorm[:,2]

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_HOD_-20.0_beta_0.54_velDispersion_3.00.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'k', label=r'HOD, $\sigma = 3.0$')

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_HOD_-20.0_beta_0.54_velDispersion_3.50.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'r', label=r'HOD, $\sigma = 3.5$')

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_HOD_-20.0_beta_0.54_velDispersion_4.00.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'g', label=r'HOD, $\sigma = 4.0$')
'''
kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_Linear_-20.0_beta_0.54_velDispersion_3.00.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'k--', label=r'linear, $\sigma=3.00$')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,2], 'k--', label=r'linear, $\sigma=3.00$')

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_Linear_-20.0_beta_0.54_velDispersion_2.75.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,2], 'r--', label=r'linear, $\sigma=3.00$')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'r--', label=r'linear, $\sigma=2.75$')

kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_Linear_-20.0_beta_0.54_velDispersion_2.50.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,2], 'g--', label=r'linear, $\sigma=3.00$')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'g--', label=r'linear, $\sigma=2.50$')
'''
kaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_Pk_Linear_-20.0_beta_0.54_velDispersion_2.00.dat')
pl.loglog(kaiserLorentz[:,0], kaiserLorentz[:,3], 'y--', label=r'linear, $\sigma=2.00$')
'''
'''
Multipoles = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Observed_Multipoles_zCube_Jenkins1.0_kbin_0.010_000.dat')
pl.loglog(Multipoles[:,0], Multipoles[:,2], 'k^', markersize=6, label='Unclipped')

Multipoles = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Observed_Multipoles_zCube_Clipped_Jenkins1.0_kbin_0.010_000.dat')
pl.loglog(Multipoles[:,0], A11Sq*Multipoles[:,2], 'k.', markersize=6, label='Clipped')
'''
# pl.axhline(y=1.1, xmin=0, xmax=0.01, color='k')
# pl.axhline(y=0.9, xmin=0, xmax=0.01, color='k')

pl.xlim([0.01, 0.8])
# pl.ylim([10**1, 10**5])
pl.ylim(10.**2., 10.**5.)
# pl.title(r'Monopole, normed by HOD, $\sigma = 3.50$')
# pl.title(r'Quadrupole, normed by HOD, $\sigma = 3.50$')

pl.ylabel(r'$P_{\ell}(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h$ \ Mpc$^{-1}]$', fontsize = '10')

leg = pl.legend(loc=3, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/CubeMultipoles2.pdf')

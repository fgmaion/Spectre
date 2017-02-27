pl.clf()
'''
Pk = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_zPencilBeamCube_Jenkins1.0_xtrans_0.00_ytrans_0.00_kInterval_0.04_000.dat')

MockCat = Pk[:,2]

for i in xrange(0, 8, 1):
    for j in xrange(0, 8, 1):
        xtrans    = i*120
        ytrans    = j*120

        In        = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_zPencilBeamCube_Jenkins1.0_xtrans_'+str(xtrans)+'.00_ytrans_'+str(ytrans)+'.00_kInterval_0.04_000.dat')

        MockCat   = np.vstack((In[:,2], MockCat))

MeanPk = np.mean(MockCat, axis =0)
VarPk = np.var(MockCat, axis=0)
        
# pl.errorbar(In[1:, 0], MeanPk[1:], np.sqrt(VarPk[1:])/np.sqrt(48), fmt='^')
'''

pl.yscale('log')
pl.xscale('log')

theory         = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/cambExtendedPk_hod_20.00.dat')
# pl.loglog(theory[:,0], (6./5.)*theory[:,1], 'k', label='vol. limited, $M_B< -20.00$')
# pl.loglog(theory[:,0],         theory[:,1], 'k--', label='vol. limited, $M_B< -20.00$')

# conv           = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ConvolvedPk/IntegralConstraint_Uncorrected/PencilBeamCube_Jenkins1.0_kbinInterval_0.01.dat')
# pl.loglog(conv[:,0], conv[:,2], 'y^', label='Convolved.')

# ICC            = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ConvolvedPk/IntegralConstraint_Corrected/zPencilBeamCube_Jenkins1.0_kbinInterval_0.04.dat')
# pl.loglog(ICC[:,0], ICC[:,2], 'y', label='Conv. + ICC')

# first catalogue is being duplicated in the stack!
Multipoles     = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/Observed_Multipoles_zPencilBeamCube_Jenkins1.0_xtrans_0.00_ytrans_0.00_kbin_0.040_000.dat')
ClipMultipoles = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/Observed_Multipoles_Clipped_zPencilBeamCube_Jenkins1.0_xtrans_0.00_ytrans_0.00_kbin_0.040_000.dat')

Monopole   = Multipoles[:,1]
Quadrupole = Multipoles[:,2]

ClipMonopole   = ClipMultipoles[:,1]
ClipQuadrupole = ClipMultipoles[:,2]

for i in xrange(0, 8, 1):
    for j in xrange(0, 8, 1):
        xtrans    = i*120
        ytrans    = j*120

        In        = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/Observed_Multipoles_zPencilBeamCube_Jenkins1.0_xtrans_'+str(xtrans)+'.00_ytrans_'+str(ytrans)+'.00_kbin_0.040_000.dat')
        
        Monopole   = np.vstack((In[:,1], Monopole))
        Quadrupole = np.vstack((In[:,2], Quadrupole)) 

        In        = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/Observed_Multipoles_Clipped_zPencilBeamCube_Jenkins1.0_xtrans_'+str(xtrans)+'.00_ytrans_'+str(ytrans)+'.00_kbin_0.040_000.dat')

        ClipMonopole   = np.vstack((In[:,1], ClipMonopole))
        ClipQuadrupole = np.vstack((In[:,2], ClipQuadrupole))

MeanMono = np.mean(Monopole, axis =0)
MeanQuad = np.mean(Quadrupole, axis=0)

ClipMeanMono = np.mean(ClipMonopole, axis =0)
ClipMeanQuad = np.mean(ClipQuadrupole, axis=0)

print MeanMono
print ClipMeanMono

VarMono = np.var(Monopole, axis=0)
VarQuad = np.var(Quadrupole, axis=0)

ClipVarMono = np.var(ClipMonopole, axis=0)
ClipVarQuad = np.var(ClipQuadrupole, axis=0)

pl.errorbar(In[1:, 0] + 0.0005, MeanMono[1:], np.sqrt(VarMono[1:])/np.sqrt(42), label='Mono.', fmt='^', c='k')
pl.errorbar(In[1:, 0], MeanQuad[1:], np.sqrt(VarQuad[1:])/np.sqrt(42), fmt='^', label='Quad.', c='c')

A11sq = 2.58787686

pl.errorbar(In[1:, 0] + 0.0005, A11sq*ClipMeanMono[1:], np.sqrt(ClipVarMono[1:])/np.sqrt(48), label='Clipped Mono.', fmt='^')
pl.errorbar(In[1:, 0], A11sq*ClipMeanQuad[1:], np.sqrt(ClipVarQuad[1:])/np.sqrt(48), fmt='^', label='Clipped Quad.', c='r')


# MB = -20.0
linearBias = 1.495903
beta       = 0.54  

KaiserMono = 1. + (2./3.)*beta + 0.2*beta*beta

HOD = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/Pk_hod_20.0.dat')
pl.loglog(HOD[:,0], KaiserMono*HOD[:,1], 'y')

'''
linear     = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/camb_matterPk.dat')
# pl.loglog(linear[:,0], (6./5.)*linearBias*linearBias*linear[:,1], 'r')

# linearQuad = 0.23737*linearBias*linearBias*linear[:,1]
# pl.loglog(linear[:,0], linearQuad, 'r')

KaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_beta_0.35_velDispersion_1.20.dat')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,1], 'r--', label=r'Kaiser-Lorentz, $\beta=0.35$ & $\sigma=1.20$')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,2], 'r--')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,3], 'r--')

KaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_beta_0.35_velDispersion_1.34.dat')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,1], 'y--', label=r'Kaiser-Lorentz, $\beta=0.35$ & $\sigma=1.34$')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,2], 'y--')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,3], 'y--')

KaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_beta_0.35_velDispersion_1.06.dat')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,1], 'm--', label=r'Kaiser-Lorentz, $\beta=0.35$ & $\sigma=1.06$')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,2], 'm--')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,3], 'm--')

KaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_beta_0.35_velDispersion_1.45.dat')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,1], 'g--', label=r'Kaiser-Lorentz, $\beta=0.35$ & $\sigma=1.45$')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,2], 'g--')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,3], 'g--')

KaiserLorentz = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles_beta_0.35_velDispersion_1.75.dat')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,1], 'c--', label=r'Kaiser-Lorentz, $\beta=0.35$ & $\sigma=1.75$')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,2], 'c--')
pl.loglog(KaiserLorentz[:,0], KaiserLorentz[:,3], 'c--')

Kaiser        = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserMultipoles_beta_0.35.dat')
pl.loglog(KaiserLorentz[:,0], Kaiser[:,1], 'k', label='Real   Mono')
pl.loglog(KaiserLorentz[:,0], Kaiser[:,2], 'g', label='Kaiser Mono')
pl.loglog(KaiserLorentz[:,0], Kaiser[:,3], 'g', label='Kaiser Quad')
'''

pl.xlim([0.06, 0.8])
pl.ylim([10**1, 2*10**4])

pl.ylabel(r'$P(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

leg = pl.legend(loc=3, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/CubeSlices_AvgPk_beta_0.5.pdf')
'''
pl.clf()

Pk = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_PencilBeamCube_Jenkins1.0_xtrans_0.00_ytrans_0.00_kInterval_0.01_000.dat')

MockCat = Pk[:,2]

for i in xrange(0, 8, 1):
    for j in xrange(0, 8, 1):
        xtrans    = i*120
        ytrans    = j*120

        In        = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_PencilBeamCube_Jenkins1.0_xtrans_'+str(xtrans)+'.00_ytrans_'+str(ytrans)+'.00_kInterval_0.01_000.dat')

        MockCat   = np.vstack((In[:,2], MockCat))

MeanPk = np.mean(MockCat, axis =0)
VarPk = np.var(MockCat, axis=0)
        
pl.errorbar(In[1:, 0], MeanPk[1:], np.sqrt(VarPk[1:])/np.sqrt(48), fmt='^')

pl.yscale('log')
pl.xscale('log')

theory         = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/cambExtendedPk_hod_20.00.dat')
pl.loglog(theory[:,0], theory[:,1], 'k', label='vol. limited, $M_B< -20.00$')

ICC            = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ConvolvedPk/IntegralConstraint_Corrected/PencilBeamCube_Jenkins1.0_kbinInterval_0.01_29.dat')
pl.loglog(ICC[:,0], ICC[:,2], 'y', label='Conv. + ICC')

ICC            = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ConvolvedPk/IntegralConstraint_Corrected/PencilBeamCube_Jenkins1.0_kbinInterval_0.01.dat')
pl.loglog(ICC[:,0], ICC[:,2], 'g', label='Conv. + ICC')

FullCube       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_Cube_Jenkins1.0_kInterval_0.01_000.dat')
pl.loglog(FullCube[:,0], FullCube[:,2], 'k^')

linear = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/cambExtendedPk_hod_20.00.dat')
pl.loglog(linear[:,0], linear[:,1], 'y')

pl.xlim([0.01, 0.8])
pl.ylim([6*10**2, 5*10**4])

pl.ylabel(r'$P(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

leg = pl.legend(loc=3, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/CubeSlices_ConvolutionCheck.pdf')
'''

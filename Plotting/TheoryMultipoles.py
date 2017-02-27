theory         = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/cambExtendedPk_hod_20.00.dat')
#pl.loglog(theory[:,0], theory[:,1], 'k', label='Real space')

beta = 0.43

ztheory = (1. + (2./3.)*beta + 0.2*beta*beta)*theory[:,1]
#pl.loglog(theory[:,0], ztheory, 'k-.', label='Monopole')

expectedQuadrupole   = ((4./3.)*beta + (4./7.)*beta**2)*theory[:,1]
#pl.loglog(theory[:,0], expectedQuadrupole, 'k--', label='Quadrupole')

expectedHexadecapole = (8./35.)*beta*beta*theory[:,1]
#pl.loglog(theory[:,0], expectedHexadecapole, 'k-.', label='Hexadecapole')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Observed_Hexadecapole_VIPERSparent_zPad_12.5.5_GaussSmoothNz_50.0_SkinDepth_5.0_VolLim_22.00_Mesh_16.00_kbin_0.010_000.dat')
#pl.loglog(data[:,0], data[:,1], 'm^', label='recovered $P_0(k)$')
#pl.loglog(data[:,0], data[:,2], 'g^', label='recovered $P_2(k)$')
#pl.loglog(data[:,0], data[:,3], 'y^', label='recovered $P_4(k)$')

#kaiserGauss = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserLorentzMultipoles.dat')

#pl.loglog(kaiserGauss[:,0], kaiserGauss[:,1], 'c')         
#pl.loglog(kaiserGauss[:,0], kaiserGauss[:,2], 'm')        
#pl.loglog(kaiserGauss[:,0], kaiserGauss[:,3], 'g')        
#pl.loglog(kaiserGauss[:,0], kaiserGauss[:,4], 'y')       
# pl.loglog(kaiserGauss[:,0], kaiserGauss[:,5], 'r')       

#expansion = 1./3. - kaiserGauss[:,0]**2/5. + kaiserGauss[:,0]**4./14.
#pl.plot(kaiserGauss[:,0], 100000.*(expansion - kaiserGauss[:,2]))

muIntegrals = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/muOrder_nLorentzIntegrals.dat')

waveNumber = muIntegrals[:,0]
RealSpace  = muIntegrals[:,1]
OrderZero  = muIntegrals[:,2]
OrderTwo   = muIntegrals[:,3]
OrderFour  = muIntegrals[:,4]
OrderSix   = muIntegrals[:,5]
OrderEight = muIntegrals[:,6]

Mono = (OrderZero + 2.*beta*OrderTwo + beta*beta*OrderFour)

Quad = (5./2.)*(-OrderZero + OrderTwo*(3. -2.*beta) + OrderFour*(-beta*beta + 6.*beta) + 3.*beta*beta*OrderSix)

Hex  = (9./8.)*(35.*OrderEight*beta**2. + 3.*OrderFour*beta**2. + 70.*OrderSix*beta + 6.*OrderTwo*beta + 35.*OrderFour  + 3.*OrderZero - 30.*beta*beta*OrderSix - 60.*beta*OrderFour - 30.*OrderTwo)

# pl.loglog(waveNumber, Hex*RealSpace, 'y')

pl.clf()
'''
pl.semilogx(waveNumber,  (9./8.)*35.*beta**2.*OrderEight, 'b--')
pl.semilogx(waveNumber,  (9./8.)*30.*beta*beta*OrderSix, 'r--')
pl.semilogx(waveNumber,   (9./8.)*3.*beta**2.*OrderFour, 'b--')
pl.semilogx(waveNumber,  (9./8.)*70.*beta*OrderSix, 'b--')
pl.semilogx(waveNumber,  (9./8.)*60.*beta*OrderFour, 'r--')
pl.semilogx(waveNumber,   (9./8.)*6.*beta*OrderTwo, 'b--')
pl.semilogx(waveNumber,  (9./8.)*35.*OrderFour, 'b--')
pl.semilogx(waveNumber,  (9./8.)*30.*OrderTwo, 'r--')
pl.semilogx(waveNumber,   (9./8.)*3.*OrderZero, 'k--')

pl.loglog(waveNumber, Hex, 'y')

limit = (np.pi**0.5)*9.*3./(8.*2.*(waveNumber*velDispersion)**1.)
pl.loglog(waveNumber, limit, 'k')
'''

velDispersion = 2.*3./2.**0.5

pl.clf()

ks = waveNumber*velDispersion

pl.loglog(ks, OrderZero, 'y')
#pl.loglog(ks, OrderTwo, 'c')
#pl.loglog(ks, OrderFour, 'g')
#pl.loglog(ks, OrderSix, 'm')
#pl.loglog(ks, OrderEight, 'k')

# pl.axvline(x=0.015, ymin=0, ymax=1)

expansion = 1.0 - 0.166667*ks**2. + 0.05*ks**4. - 0.0178571*ks**6.
pl.loglog(ks, expansion, 'k--')

# pl.loglog(waveNumber, RealSpace, 'b')
#pl.loglog(waveNumber, 4.675*10.**6.*waveNumber**0.96, 'r')

measure =  10.**6.*abs(expansion - OrderZero)

pl.yscale('log')
pl.xscale('log')

#pl.xlim([0.01, 0.1.0])
# pl.ylim([0.31, 0.32])
# pl.ylim(0.1, 0.075)

# xx, locs = plt.xticks()
# ll = ['%.3f' % a for a in xx]
# plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

# leg = pl.legend(loc=1, ncol=1, prop = FontProperties(size = '10'))
# leg.draw_frame(False)

pl.ylabel(r'P$_n$(k $\sigma$)', fontsize = '10')
pl.xlabel(r'$k \sigma \ [h \ Mpc^{-1}]$', fontsize = '10')
pl.title(r'Kaiser-Lorentz model')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/kaiserLorentz_muIntegrals_0.01.pdf')
# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/kaiserLorentz_muIntegrals_0.01.pdf')
# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/TheoryMultipoles_0.01.pdf')

pl.clf()

kaiserGauss = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/KaiserGaussMultipoles_Beta_0.35_velDispersion_1.75.dat')
pl.loglog(kaiserGauss[:,0], kaiserGauss[:,1], 'c')         
pl.loglog(kaiserGauss[:,0], kaiserGauss[:,2], 'm')        
pl.loglog(kaiserGauss[:,0], kaiserGauss[:,3], 'g')        
pl.loglog(kaiserGauss[:,0], kaiserGauss[:,4], 'y')     

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Observed_Hexadecapole_VIPERSparent_zPad_12.0_8.0_8.0_GaussSmoothNz_50.0_SkinDepth_5.0_VolLim_22.00_Mesh_16.00_kbin_0.010_000.dat') 
pl.loglog(data[:,0], data[:,1], 'm^', label='Mono')
pl.loglog(data[:,0], data[:,2], 'g^', label='Quad')
pl.loglog(data[:,0], data[:,3], 'y^', label='Hex')

leg = pl.legend(loc=1, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.ylim(1., 10.**5.)

pl.ylabel(r'P$_n$(k)', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')
pl.title(r'Kaiser-Gauss model')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/KaiserGauss.pdf')

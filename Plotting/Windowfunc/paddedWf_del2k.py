def plotConvolvedPk():
    Datadir='/disk1/mjw/HOD_MockRun/Data/' 
    
    unpaddedConv = np.loadtxt(Datadir+'ConvolvedPk/IntegralConstraint_Uncorrected/midK_Del2k_Sphericalsplint_unpadded.dat')
    pl.loglog(unpaddedConv[:,0], unpaddedConv[:,1], 'k^', label='unpadded')
    
    paddedConv  = np.loadtxt(Datadir+'ConvolvedPk/IntegralConstraint_Uncorrected/midK_Del2k_Sphericalsplint_padded.dat')
    pl.loglog(paddedConv[:,0], paddedConv[:,1], 'r^', label='padded')
    
    analyticConv = np.loadtxt(Datadir+'ConvolvedPk/IntegralConstraint_Uncorrected/midK_Del2k_Sphericalsplint_analytic.dat')
    pl.loglog(analyticConv[:,0], analyticConv[:,1], label='analytic')
    
    regInterp = np.loadtxt(Datadir+'InterpTheoryPk/regInterpCambExtendedPk_0.0090_hod_20.0.dat')
    InputD2k  = (1./(2.*np.pi*np.pi))*regInterp[:,0]**3*regInterp[:,1]
    pl.loglog(regInterp[:,0], InputD2k, label='input P(k)')
    
    xx, locs = plt.xticks()
    ll = ['%.3f' % a for a in xx]
    plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

    leg = pl.legend(loc =2, ncol=1, prop = FontProperties(size = '10'))
    leg.draw_frame(False)
    
    pl.xlim([0.003, 0.3])
    pl.ylim([10**-4, 2.])
    
    pl.ylabel(r'$\Delta^2(k)$', fontsize = '10')
    pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')
    pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ConvolvedPk/Spherical_Del2k.eps')
    

plotConvolvedPk()

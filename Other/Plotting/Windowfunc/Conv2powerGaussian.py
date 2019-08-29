p0        = 150.
Rg        =  10.

convolved = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ConvolvedPk/IntegralConstraint_Uncorrected/midK_Pk_FullCube_Jenkins1.0_powerlaw2_Gauss.dat')

kvals     = convolved[:,0]

power2law = p0*(kvals**2) + p0*(3./2.)*(1./Rg)*(1./Rg)

residuals = (convolved[:,1] - power2law)/power2law

# pl.loglog(kvals, residuals, 'r^')

pl.plot(kvals + 0.001, convolved[:,1], 'r^')
pl.plot(kvals, power2law, 'g^')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ConvolvedPk/analytic2powerGauss.pdf')

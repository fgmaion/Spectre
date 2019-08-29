data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ConvolvedPk/IntegralConstraint_Uncorrected/midK_Pk_Spherical_Jenkins1.0_powerlaw2_Gauss.dat')

pl.loglog(data[:,0], data[:,1], 'b^')

p0 = 150.
R0 = 10.

abscissa = np.arange(0.001, 1.0, 0.001)

theory = p0*(abscissa**2. + (3./2.)/(R0*R0))

pl.loglog(abscissa, theory, 'g')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ConvolvedPk/SphConvAnalyticTestCase.pdf') 

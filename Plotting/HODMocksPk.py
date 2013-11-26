files     = glob.glob("/disk1/mjw/HOD_MockRun/*.txt")

binNumber = 362

Data      = np.zeros(binNumber*len(files)).reshape(binNumber, len(files))
midBin    = np.loadtxt(files[0])[:,0]

for i in xrange(0, len(files), 1):
  Data[:,i] = np.loadtxt(files[i])[:,1]

errors      = np.std(Data, axis  = 1)/np.sqrt(len(files) -1)
mean        = np.mean(Data, axis = 1) 

pl.errorbar(midBin, mean, yerr=errors, label = 'HOD mocks, VIPERS window fn.', fmt='o', markersize=2)

#matplotlib.pyplot.axhline(y=7.9*10**-6, c='y', label='Unweighted shot noise expectation.')
#matplotlib.pyplot.axhline(y=1.1*10**-6, c='g', label='FKP shot noise expectation.')

regInterpTheory  = np.loadtxt('/disk1/mjw/HOD_MockRun/regInterpCambExtendedPk_0.0050_hod_20.0.dat')
pl.loglog(regInterpTheory[:,0], regInterpTheory[:,1]*(regInterpTheory[:,0]**3)*4.*np.pi/((2.*np.pi)**3), 'r--', label='HOD theory.')

convolvedTheory  = np.loadtxt('/disk1/mjw/HOD_MockRun/ConvolvedPk.dat')

print regInterpTheory[100,1]/convolvedTheory[200,1]

convolvedTheory[:,1] /= convolvedTheory[200,1]
convolvedTheory[:,1] *= regInterpTheory[100,1]

pl.loglog(convolvedTheory[:,0], convolvedTheory[:,1]*(convolvedTheory[:,0]**3)*4.*np.pi/((2.*np.pi)**3), 'b--', label='convolved theory.')

#IntegralConstraintCorrected = np.loadtxt("/disk1/mjw/HOD_MockRun/ConvolvedPk_IntegralConstraintCorrected.dat")

#IntegralConstraintCorrected[:,1] /= IntegralConstraintCorrected[160,1]
#IntegralConstraintCorrected[:,1] *= regInterpTheory[160,1]

#pl.loglog(IntegralConstraintCorrected[:,0], IntegralConstraintCorrected[:,1]*(IntegralConstraintCorrected[:,0]**3)*4.*np.pi/((2.*np.pi)**3), 'b--', label='Integral constraint corrected.')

pl.xlim([0.001, 2.0])
#pl.ylim([10**-3, 2.5*10**2])

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

leg = pl.legend(loc =2, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/HODMocksW1_3hMpc_ConvolvedPk_RealSpace_0.7_0.9.eps')

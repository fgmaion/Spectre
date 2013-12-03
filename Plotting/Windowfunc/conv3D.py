data   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Pk/midK_Pk_ConvolvedAnsioGauss.dat')

a      = 1.0
b      = 0.5
c      = 0.2
          
theory = 150.*(data[:,0]**2 + 0.5*( a**-2. + b**-2. + c**-2. ))

pl.plot(data[:,0], data[:,1], 'r')
pl.plot(data[:,0], theory,    'b')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ConvolvedPk/3DconvolvedAnisoGauss.pdf')
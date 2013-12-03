data   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Pk/midK_Pk_ConvolvedAnisoGauss.dat')

P0     = 150.

a      = 1.0
b      = 0.5
c      = 0.2
          
inputPk  = P0*(data[:,0]**2)          
theory = P0*(data[:,0]**2 + 0.5*(a**-2. + b**-2. + c**-2.))

pl.plot(data[:,0], data[:,2], 'r')
pl.plot(data[:,0], theory,    'b')
pl.plot(data[:,0], inputPk,   'g')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ConvolvedPk/3DconvolvedAnisoGauss.pdf')
data       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Pk/midK_Pk_ConvolvedAnisoGauss.dat')

MeasuredWf = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Pk/midK_Pk_ConvolvedMeasuredAnisoGauss.dat')

P0         = 150.

a          =  80.0
b          =  100.0
c          =  120.0
          
inputPk    = P0*(data[:,0]**2)          
theory     = P0*(data[:,0]**2 + 0.5*(a**-2. + b**-2. + c**-2.))

pl.loglog(data[:,0], data[:,2], 'r^')
pl.loglog(MeasuredWf[:,0], MeasuredWf[:,2], 'y^')
pl.loglog(data[:,0], theory,    'b')
pl.loglog(data[:,0], inputPk,   'g')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ConvolvedPk/3DconvolvedAnisoGauss.pdf')

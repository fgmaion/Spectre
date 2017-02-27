data       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Pk_ConvolvedAnisoGaussTestCase2.dat')

# MeasuredWf = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Pk_ConvolvedAnisoGaussTestCase_MeasuredWf.dat')

P0         = 150.

a          =  80.0
b          =  100.0
c          =  120.0
          
inputPk    = P0*(data[:,0]**2)          
theory     = P0*(data[:,0]**2 + 0.5*(a**-2. + b**-2. + c**-2.))

pl.loglog(data[:,0], data[:,1], 'r^')

pl.loglog(data[:,0], theory,    'g')
pl.loglog(data[:,0], inputPk,   'k')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ConvolvedPk/3DWfConvolved_AnisoGaussTestCase2.pdf')

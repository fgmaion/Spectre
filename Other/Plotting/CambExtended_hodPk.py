one                 = np.loadtxt('/disk1/mjw/HOD_MockRun/Scripts/Pk_hod_20.0.dat')
linear              = np.loadtxt('/disk1/mjw/HOD_MockRun/Scripts/Camb_matterpower_100.dat')

biasIndexCamb       = 259
biasIndexHOD        = 1

linear[:,1]        *= one[biasIndexHOD,1]/linear[biasIndexCamb, 1] 

extendedK           = np.concatenate((linear[:259,0], one[1:,0]))
extendedPk          = np.concatenate((linear[:259,1], one[1:,1]))

Array               = np.zeros(len(extendedK)*2).reshape(len(extendedK), 2)
Array[:,0]          = extendedK
Array[:,1]          = extendedPk

np.savetxt('/disk1/mjw/HOD_MockRun/Scripts/cambExtendedPk_hod_20.0.dat', Array)

Del2L               = linear[:,1]*(linear[:,0]**3)*4*np.pi/((2*np.pi)**3)

pl.loglog(linear[:,0], Del2L, 'r', label = 'Camb extended HOD P(k), M_B < -20.0 magnitudes')
pl.loglog(one[:,0], 1./(2.*np.pi**2)*one[:,1]*one[:,0]**3, label='Original HOD P(k)')

#regInterpTheory  = np.loadtxt('/disk1/mjw/HOD_MockRun/regInterpCambExtendedPk_0.0025_hod_20.0.dat')
#pl.loglog(regInterpTheory[:,0], regInterpTheory[:,1]*(regInterpTheory[:,0]**3)*4.*np.pi/((2.*np.pi)**3), 'r--', label='HOD theory.')

leg = pl.legend(loc =2, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.xlim([0.001, 1.0])
#pl.ylim([10**-3, 2.5*10**2])

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

pl.ylabel(r'$\Delta^2(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/CambExtended_hodPk.eps')

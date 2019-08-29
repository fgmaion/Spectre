data             = np.loadtxt('/disk1/mjw/HOD_MockRun/midK_Del2k_HODPencilBeam.dat')
pl.loglog(data[:,0], data[:,1], 'r^', label='Measured - pencil beam survey.')

#data2             = np.loadtxt('/disk1/mjw/HOD_MockRun/midK_Del2k_HODCube.dat')
#pl.loglog(data2[:,0], data2[:,1], 'k^', label='Measured - 3D Cube.')

PencilBeam = np.loadtxt('/disk1/mjw/HOD_MockRun/midK_ConvDel2k_PencilBeam.dat')
pl.loglog(PencilBeam[:,0], PencilBeam[:,1], 'y--', label='Convolved theory, pencil beam survey.')

regInterpTheory  = np.loadtxt('/disk1/mjw/HOD_MockRun/regInterpCambExtendedPk_0.0090_hod_20.0.dat')
pl.loglog(regInterpTheory[:,0], regInterpTheory[:,1]*(regInterpTheory[:,0]**3)*4.*np.pi/((2.*np.pi)**3), 'g--', label='HOD theory.')

leg = pl.legend(loc =2, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.xlim([0.001, 5.0])

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

pl.ylabel(r'$\Delta^2(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/HODPencilBeam_del2k.eps')

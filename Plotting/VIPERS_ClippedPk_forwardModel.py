data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/Distorted/KaiserLorentzian_monopole_kInterval_0.01_beta_0.56_velDispersion_3.93_u0_-0.09.dat')

pl.loglog(data[:,0], data[:,5], 'y', label='HOD P(k)')
pl.loglog(data[:,0], data[:,6], 'k--', label='Suppressed')
pl.loglog(data[:,0], data[:,2], 'k', label='Distorted, u0 = -0.09')

'''
data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/Distorted/KaiserLorentzian_monopole_kInterval_0.01_beta_0.56_velDispersion_3.93_u0_6.28.dat')

pl.loglog(data[:,0], data[:,5], 'y', label='HOD P(k)')
pl.loglog(data[:,0], data[:,6], 'k--', label='Suppressed')
pl.loglog(data[:,0], data[:,2], 'k', label='Distorted, u0 = 6.28')
'''

pl.legend(loc=3)

pl.title('Forward modelling of Clipped P(k)')
pl.xlabel('k')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/ClippedPk_forwardModel_clippingDraft.pdf')

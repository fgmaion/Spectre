# Now using a single mock estimate of Nz, but a volume limited sample at -21.00;

fig  = plt.figure()

ax1  = fig.add_subplot(211)

surveyType = 'PencilBeamCube'

Data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_'+surveyType+'_Jenkins1.0_kInterval_0.01_000.dat')
ax1.loglog(Data[:,0], Data[:,2], 'b^', markersize='3', label='Observed')

# Data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_FullCube_Jenkins1.0_kInterval_0.01_000.dat')
# pl.loglog(Data[:,0], Data[:,2], 'g^', markersize='3')

pl.yscale('log')
pl.xscale('log')

theory         = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/HODTheoryPk/Pk_hod_20.15.dat')
ax1.loglog(theory[:,0], theory[:,1], 'c', label='vol. limited, $M_B< -20.15$')

# Cube      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Del2k/midK_Del2k_FullCube_Jenkins1.0.dat')
# pl.loglog(Cube[:,0], Cube[:,2], 'g^', label='Mock 001')

# conv           = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ConvolvedPk/IntegralConstraint_Uncorrected/FullCube_Jenkins1.0_kbinInterval_0.03.dat')
# pl.loglog(conv[:,0], conv[:,2], 'c^', label='theory conv.', markersize='3')

conv           = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ConvolvedPk/IntegralConstraint_Uncorrected/'+surveyType+'_Jenkins1.0_kbinInterval_0.01.dat')
ax1.loglog(conv[:,0], conv[:,2], 'c^', label='theory conv.', markersize='3')

leg = pl.legend(loc=1, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.title(surveyType+', ($z=0.8$)')

pl.xlim([0.01, 0.8])
pl.ylim([5.*10**2, 3.2*10**4])

pl.ylabel(r'$P(k)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

ax2  = fig.add_subplot(212)

pl.ylabel('%')

pl.xlim([0.01, 0.8])
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.xscale('log')
pl.ylim(-10,10)

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

theoryPrediction = conv[:,5]

theoryDiff = 100.*(conv[:,2] - conv[:,5])/conv[:,5]
dataDiff   = 100.*(Data[:,2] - Data[:,5])/Data[:,5]

ax2 = pl.plot(conv[:,0], theoryDiff, 'c^')
ax2 = pl.plot(conv[:,0], dataDiff, 'b^')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Pk_Windowfn_Applied2Cube_'+surveyType+'.pdf')

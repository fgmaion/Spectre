bestfit   = np.loadtxt("/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_MinChi2_MockAvgNz_0.03.dat")
pl.clf()

ax        = pl.subplot(111)

ax.bar(bestfit[:,0], bestfit[:, 1], color='k', width=0.03, label='mock avg.',  alpha=0.3, align='center')

pl.plot(bestfit[:, 0], bestfit[:,2], 'r')

pl.ylabel('N(z) $[deg^{-2} \ (\Delta z = 0.03)^{-1}]$')
pl.xlabel('z')

leg = pl.legend(loc =1, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.savefig("/disk1/mjw/HOD_MockRun/Plots/nz/MinChi2_Nz.pdf")

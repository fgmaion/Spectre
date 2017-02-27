mock001 = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_001_Nz_chiSliced_10.0_W1andW4_Gaussfiltered_100.0.dat')

ax      = pl.subplot(111)

ax.bar(mock001[:,1], mock001[:,3], color='g', width=10.0,  alpha=0.8)

pl.xlim(900, 2900)

pl.ylabel('N($\chi$)')
pl.xlabel(r'$ \chi $')
pl.title('Tophat filter check')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/nz/tophat_GaussianfilterCheck.pdf')

Mock001    = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_001_Nz_chiSliced_10.0_W1andW4_Gaussfiltered_100.0.dat')
MockVolLim = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_001_Nz_chiSliced_10.0_W1andW4_VolLim_-20.15_Gaussfiltered_100.0.dat')

ax         = pl.subplot(111)

# ax.bar(Mock001[:,1],  Mock001[:,2], color='k', width=10., label='mock 001.', alpha=0.3)
ax.bar(MockVolLim[:,1],  MockVolLim[:,2], color='k', width=10., label='mock 001 vol limited.', alpha=0.5)

ax.plot(MockVolLim[:,1], MockVolLim[:,3], color='r', label=r'Mock 001 vol. limited')
ax.plot(Mock001[:,1], Mock001[:,3], color='y', label=r'Mock 001 all galaxies')

# pl.axvline(x=, ymin=0, ymax=1)
# pl.axvline(x=, ymin=0, ymax=1)

# pl.ylim(0.,120.)

pl.xlabel('z')
pl.ylabel('N($\chi$) $[deg^{-2} \ (\Delta \chi = 10.0)^{-1}]$')
pl.title('Volume limited, $0.7<z<0.9$, $M_B<-20.15$')

leg = pl.legend(loc =1, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/nz/VolumeLimited_Nz.pdf')

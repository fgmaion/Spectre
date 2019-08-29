MockAvg  = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_MockAvgNz_chiSliced_15.0_W1andW4_VolLim_22.00_Gaussfiltered_50.0.dat')
Steffano = np.loadtxt('/disk1/mjw/HOD_MockRun/Steffano/nz.txt')

# Mock001 = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_001_SingleMockNz_chiSliced_10.0_W1andW4_VolLim_22.00_Gaussfiltered_50.0.dat')

ax      = pl.subplot(111)

#ax.plot(MockAvg[:,1], MockAvg[:,5], color='y', label=r'fitted func. form$')


W1area                =    10.875;
W4area                =       8.8;
TotalW1W4area         =    19.675; 

# ax.bar(Mock001[:,1],  Mock001[:,2], color='b', width=10., label='mock 001.', alpha=0.3)
ax.bar(MockAvg[:,1],  W1area*MockAvg[:,2], color='b', width=15., label='mock avg.', alpha=1.0)
ax.bar(Steffano[:,0], Steffano[:,1], color='y', width=15., alpha=0.4, label='Stefano')

pl.plot(Steffano[:,0], Steffano[:,2], 'y')
pl.plot(MockAvg[:,1],  W1area*MockAvg[:,3], c='b')

# ax.plot(MockAvg[:,1], MockAvg[:,3], color='r', label=r'Mock avg. gaussian filtered')
# ax.plot(Mock001[:,1], Mock001[:,3], color='c', label=r'Mock 001  gaussian filtered')
# ax.plot(Mock001[:,1], Mock001[:,5], color='g', label=r'fitting func. form')
# ax.plot(Mock001[:,1], Mock001[:,6], color='y', label=r'Mock 001, div func., ln etc.')

# Interim = np.log(Mock001[:,2]/Mock001[:,5])
# ax.plot(Mock001[:,1], Interim, color='g', label=r'Interim')

# pl.axvline(x=, ymin=0, ymax=1)
# pl.axvline(x=, ymin=0, ymax=1)

# pl.ylim(0.,120.)

pl.xlabel('$\chi$')
# pl.ylabel('N($\chi$) $[deg^{-2} \ (\Delta \chi = 10.0)^{-1}]$')

# pl.ylabel('smoothed ln(Observed/Theory)')

leg = pl.legend(loc =1, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/nz/MockAvg_Nz.pdf')

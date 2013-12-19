xdata    = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/VIPERSparent_xSlice.dat')
pl.loglog(xdata[:,1], xdata[:,2], 'b', label='x direction')

ydata    = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/VIPERSparent_ySlice.dat')
pl.loglog(ydata[:,1], ydata[:,2], 'g', label='y direction')

zdata    = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/VIPERSparent_zSlice.dat')
pl.loglog(zdata[:,1], zdata[:,2], 'y^', label='z direction')

#pl.xlim([0.001, 0.3])
#pl.ylim([0.0001, 1.0])

matplotlib.pyplot.axvline(x=0.0097, ymin=0, ymax=1, color='b')
matplotlib.pyplot.axvline(x=0.0180, ymin=0, ymax=1, color='g')
matplotlib.pyplot.axvline(x=0.0419, ymin=0, ymax=1, color='y')

xx, locs = plt.xticks()
ll = ['%.3f' % a for a in xx]
plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

leg = pl.legend(loc=3, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.ylabel(r'$W^2(k_i)$', fontsize = '10')
pl.xlabel(r'$k \ [h \ Mpc^{-1}]$', fontsize = '10')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/MarseilleWf.pdf')

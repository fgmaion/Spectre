xSlice       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/PencilBeamCube_Jenkins1.0_xSlice.dat')
ySlice       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/PencilBeamCube_Jenkins1.0_ySlice.dat')
zSlice       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/PencilBeamCube_Jenkins1.0_zSlice.dat')

pl.semilogx(xSlice[:,1], xSlice[:,2]/xSlice[0,2], 'g^', label='x')
pl.semilogx(ySlice[:,1] + 0.0005, ySlice[:,2]/ySlice[0,2], 'y^', label='y')
pl.semilogx(zSlice[:,1], zSlice[:,2]/zSlice[0,2], 'r^', label='z')

plt.axvline(x=0.0879646, ymin=0, ymax=1,  c='g')
plt.axvline(x=0.0879646, ymin=0, ymax=1,  c='y')
plt.axvline(x=0.00628319, ymin=0, ymax=1, c='r')

leg = pl.legend(loc=1, ncol=1, prop = FontProperties(size = '10'))
leg.draw_frame(False)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/PencilBeam_Windowfunc_Windowfn.pdf')

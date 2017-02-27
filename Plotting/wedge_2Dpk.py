import copy
import scipy.interpolate
from   matplotlib.colors import LogNorm

cmap = copy.copy(matplotlib.cm.jet)

cmap.set_bad('w',1.)

# plt.imshow(z2Pk, origin='lower', extent=[0.0, 1., 0.0, 1.], norm = LogNorm(vmin=10**2, vmax=3*10**4))

# plt.colorbar()

# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/linearPk.pdf', bbox_inches='tight', pad_inches=0.5)
# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/2D_kaiserPk.pdf', bbox_inches='tight', pad_inches=0.5)
# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/2d_kaiserLorentz_300.pdf', bbox_inches='tight', pad_inches=0.5)
# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/2D_kaiserLorentz_300Clipped.pdf', bbox_inches='tight', pad_inches=0.5)

# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/Observed2Dpk_FullCube_Jenkins2.0.pdf', bbox_inches='tight', pad_inches=0.5)
# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/Observed2Dpk_zFullCube_Jenkins4.0.pdf', bbox_inches='tight', pad_inches=0.5)
# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/Observed2Dpk_zFullCube_clipped.pdf', bbox_inches='tight', pad_inches=0.5)

# plt.scatter(data[:,0], data[:,1], c=np.log10(mean), vmin=1.5, vmax=5.5)
# plt.colorbar()

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Cartesian2Dpk_zCube_Clipped_Jenkins1.0_0.010_0.dat')
# plt.scatter(data[:,0], data[:,1], c=np.log10(data[:,2]), vmin=1.5, vmax=5.5)

plt.imshow(data[:,2].reshape(79, 79), origin='lower', interpolation='nearest', extent=[0.0, 0.8, 0.0, 0.8]) # norm = LogNorm(vmin=-10**3, vmax=3*10**4), 

plt.colorbar()

# pl.xlim(0., 0.6)
# pl.ylim(0., 0.6)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/zCube_Clipped_5.0Cartesian.pdf', bbox_inches='tight', pad_inches=0.5)

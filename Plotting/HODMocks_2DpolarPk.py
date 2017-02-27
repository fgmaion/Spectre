import copy
import scipy.interpolate

from   matplotlib.colors import LogNorm

cmap = copy.copy(matplotlib.cm.jet)

# cmap.set_bad('w',1.)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Polar2Dpk_VIPERSparent_zPad_12.0_0.15_0.1_GaussSmoothNz_50.0_SkinDepth_5.0_VolLim_22.00_Mesh_2.00_0.010_3.dat')

z2Pk = data[:,2].reshape(99, 99)

# z2Pk = np.rot90(z2Pk, 3)

plt.imshow(z2Pk, origin='lower', extent=[0.0, 1.0, 0.0, 1.0], norm = LogNorm(vmin=10**2, vmax=3*10**4))

# plt.yscale('log')

plt.colorbar()

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/HODMocks_2DpolarPk_VIPERSparent_003.pdf', bbox_inches='tight', pad_inches=0.5)


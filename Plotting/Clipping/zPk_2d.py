import scipy.interpolate
from   matplotlib.colors import LogNorm

#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/2D_linearPk.dat')
#data  = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/2D_kaiserPk.dat')
#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/2D_kaiserLorentz_300Pk.dat')
#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/2D_kaiserLorentz_300ClippedPk.dat')

#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_FullCube_Jenkins2.0.dat')
#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_zFullCube_Jenkins4.0.dat')
#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_zFullCube_clipped.dat')
data  = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/2D_legendre2weights.dat') 

perpk = data[:,0]
losk  = data[:,1] 
z2Pk  = data[:,2]

len   = len(data[:,2])**0.5

z2Pk  = z2Pk.reshape(len, len)

#plt.imshow(z2Pk, vmin=-0.5, vmax=10., origin='lower', extent=[perpk.min(), perpk.max(), losk.min(), losk.max()], norm=matplotlib.colors.LogNorm(), cmap='hsv')

plt.imshow(z2Pk, vmin=-0.5, vmax=1., origin='lower', extent=[perpk.min(), perpk.max(), losk.min(), losk.max()], cmap='hsv')

plt.colorbar()

#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/linearPk.pdf', bbox_inches='tight', pad_inches=0.5)
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/2D_kaiserPk.pdf', bbox_inches='tight', pad_inches=0.5)
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/2d_kaiserLorentz_300.pdf', bbox_inches='tight', pad_inches=0.5)
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/2D_kaiserLorentz_300Clipped.pdf', bbox_inches='tight', pad_inches=0.5)

#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/Observed2Dpk_FullCube_Jenkins2.0.pdf', bbox_inches='tight', pad_inches=0.5)
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/Observed2Dpk_zFullCube_Jenkins4.0.pdf', bbox_inches='tight', pad_inches=0.5)
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/Observed2Dpk_zFullCube_clipped.pdf', bbox_inches='tight', pad_inches=0.5)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/legendre2weights.pdf', bbox_inches='tight', pad_inches=0.5)

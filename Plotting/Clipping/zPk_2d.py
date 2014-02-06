import copy
import scipy.interpolate

from   matplotlib.colors import LogNorm

cmap = copy.copy(matplotlib.cm.jet)

cmap.set_bad('w',1.)

#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/2D_linearPk.dat')
#data  = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/2D_kaiserPk.dat')
#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/2D_kaiserLorentz_300Pk.dat')
#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/2D_kaiserLorentz_300ClippedPk.dat')

#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_FullCube_Jenkins2.0.dat')
#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_zFullCube_Jenkins4.0.dat')
#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_zFullCube_clipped.dat')
#data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/2D_legendre2weights.dat') 

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/MockAvgObserved2Dpk_VIPERSparent_Mock001_GaussSmoothNz_100.0_SkinDepth_5.0.dat')
# data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/ConvolvedTheory2Dpk_VIPERSparent_Mock001_GaussSmoothNz_100.0_SkinDepth_5.0.dat')
# data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/InputTheory2Dpk_VIPERSparent_Mock001_GaussSmoothNz_100.0_SkinDepth_5.0.dat')

# data= np.loadtxt('/disk1/mjw/HOD_MockRun/Data/ClippedPk/zSpace/2Dpk/IccConvolvedTheory2Dpk_VIPERSparent_Mock001_GaussSmoothNz_100.0_SkinDepth_5.0.dat')

#perpk    = data[:,0]
#losk     = data[:,1] 
z2Pk     = data
#ModeNumb = data[:,3]

z2Pk     = z2Pk.reshape(29, 29)

plt.imshow(z2Pk, origin='lower', extent=[0.0, 0.6, 0.0, 0.6], norm = LogNorm(vmin=10**2, vmax=3*10**4))

plt.colorbar()

#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/linearPk.pdf', bbox_inches='tight', pad_inches=0.5)
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/2D_kaiserPk.pdf', bbox_inches='tight', pad_inches=0.5)
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/2d_kaiserLorentz_300.pdf', bbox_inches='tight', pad_inches=0.5)
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/2D_kaiserLorentz_300Clipped.pdf', bbox_inches='tight', pad_inches=0.5)

#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/Observed2Dpk_FullCube_Jenkins2.0.pdf', bbox_inches='tight', pad_inches=0.5)
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/Observed2Dpk_zFullCube_Jenkins4.0.pdf', bbox_inches='tight', pad_inches=0.5)
#pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/Observed2Dpk_zFullCube_clipped.pdf', bbox_inches='tight', pad_inches=0.5)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/MockAvgObserved2Dpk_VIPERSparent_Mock001_GaussSmoothNz_100.0_SkinDepth_5.0.pdf', bbox_inches='tight', pad_inches=0.5)

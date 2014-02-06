import scipy.ndimage as ndi
import numpy as np
import pylab as pl

from matplotlib import pyplot as plt

data        = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_Nz_chiSliced_W1andW4_10.0.dat')

VIPERS_area = 21.47 # sq degs.
 
z           = data[:,0]
Chi         = data[:,1]
Nz          = data[:,2]

Nz         /= VIPERS_area

delta       = np.zeros(len(Chi))
delta[15]   = 1.0

filtNz      = ndi.gaussian_filter1d(Nz, sigma=16.0, output=np.float64, mode='nearest')

filtDel     = ndi.gaussian_filter1d(delta, sigma=1.0, output=np.float64, mode='nearest')

ax          = pl.subplot(111)

up          = np.roll(Chi, 1)
up[0]       = 0.0

ax.bar(Chi, Nz, color='b', width= 10, align='center') 
pl.plot(Chi, filtNz, 'r')
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/GaussianfilteredNz.jpg')

pl.xlim([Chi.min(), Chi.max()])

np.savetxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_Mock001_Gaussfilt16.0sig_Nz_10Mpc.dat', filtNz)

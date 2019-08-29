import scipy.ndimage as ndi
import numpy as np
import pylab as pl

from matplotlib import pyplot as plt

def Nz(sigma):
        pl.clf()
        
        ax          = pl.subplot(111)

        data        = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_001_Nz_chiSliced_10.0_W1andW4_Gaussfiltered_' + str(sigma) + '.0.dat')

        VIPERS_area = 21.47 # sq degs.
 
        z           = data[:,0]
        Chi         = data[:,1]
        Nz          = data[:,2]
        filtNz      = data[:,3]
        ComovDens   = data[:,4]
        Analytic    = data[:,5]
  
        divfn_ln    = data[:,6]
 
        Nz         /= VIPERS_area
        filtNz     /= VIPERS_area
  
        divfn_ln   /= VIPERS_area
        Analytic   /= VIPERS_area

        # 4*pi*Chi^2*dChi*               whole sky in sq degs. 
        Totalfilt = 4.*np.pi*Chi*Chi*10.*filtNz*(129600./np.pi)
        Total     = 4.*np.pi*Chi*Chi*10.*Nz*(129600./np.pi)

        ax.bar(Chi, Nz, color='k', width= 10, align='center', alpha=0.3) 
        
        pl.plot(Chi, filtNz, label = '$\sigma_{'+str(sigma)+'}$')
        pl.plot(Chi, divfn_ln, label = r'$ln(\frac{N(\chi)}{f(\chi)}), \ \sigma_{'+str(sigma)+'}$')

        pl.plot(Chi, Analytic, label='$f(\chi)$')

        leg = pl.legend(loc=1, ncol=1, prop = FontProperties(size = '10'))
        leg.draw_frame(False)

        pl.xlabel('$\chi$')
        pl.ylabel('$N(\chi) / [deg^{-2} \ (d\chi=10.0)^{-1}]$')

        pl.xlim([Chi.min(), Chi.max()])

        pl.savefig('/disk1/mjw/HOD_MockRun/Plots/GaussianfilteredNz' + str(sigma) + '.pdf')

Nz(100)
Nz(200)

# pl.clf()
# pl.hist(np.log(Nz/Analytic), bins=100, color='b', alpha=0.3, label=r'$ln(\frac{N(\chi)}{f(\chi)})$') 
# pl.hist(Nz/Analytic, bins=100, color='g', alpha=0.3, label=r'$\frac{N(\chi)}{f(\chi)}$') 

# leg = pl.legend(loc=1, ncol=1, prop = FontProperties(size = '10'))
# leg.draw_frame(False)

# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/NChi_logResiduals.pdf')

# pl.clf()
# pl.hist(np.log(divfn_ln/Analytic), bins=100, color='b', alpha=0.3, label=r'$ln(\frac{N_{filtered}(\chi)}{f(\chi)})$') 
# pl.xlim([-1., 2.0])

# leg = pl.legend(loc=1, ncol=1, prop = FontProperties(size = '10'))
# leg.draw_frame(False)

# pl.savefig('/disk1/mjw/HOD_MockRun/Plots/filteredNChi_logResiduals.pdf')

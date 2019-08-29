import numpy as np
import pylab as pl


vol = 1.e9

dat       = np.loadtxt('/global/homes/m/mjwilson/Spectrum-MC/output/output.txt')
dat[:,1] *= vol
dat[:,2] *= vol

pl.loglog(dat[:,0],        dat[:,1],  'k-', label=r'$P_0$')
pl.loglog(dat[:,0], np.abs(dat[:,2]), 'r-', label=r'$P_2$')

pl.legend(frameon=False)
pl.savefig('pk.pdf')

#!/home/mlam/anaconda/bin/
import sys, os;       sys.path.append(os.environ['PYTHONPATH'])
import matplotlib;    matplotlib.use('PDF') # matplotlib.use('Agg')

from   matplotlib     import style; style.use('classic')
from   emcee          import EnsembleSampler
from   emcee.utils    import MPIPool
from   scipy.optimize import minimize
from   ctypes         import cdll
from   ctypes         import c_double, c_int

import matplotlib     as mpl
import numpy          as np
import pylab          as pl
import mpi4py
import corner
import sys

'''
import getdist
from   getdist        import loadMCSamples, plots
'''

# load vipers likelihood as shared library.
likelihood         = cdll.LoadLibrary("/home/mjw/HOD_MockRun/likelihood.so")
calc_ChiSq         = likelihood.calc_ChiSq
calc_ChiSq.restype = c_double # return types are integer by default. 

def lnprior3D(x):
  if 0.0 < x[0] < 1.80 and 0.00 < x[1] < 3.5 and 0.0 < x[2] < 15.:
    return 0.0
  return -np.inf

def lnprior4D(x):
  if 0.0 < x[0] < 1.80 and 0.00 < x[1] < 3.5 and 0.0 < x[2] < 15. and -0.3 < x[3] < 0.8:
    return 0.0
  return -np.inf

def chiSq_3D(x):
    if not 0.0 < x[0] < 1.80 and not 0.00 < x[1] < 3.5 and not 0.0 < x[2] < 15.:
        return np.inf

    return -2.*calc_ChiSq(c_double(x[0]), c_double(x[1]), c_double(x[2]), c_double(0.0))
       
def chiSq_4D(x):
  if not 0.0 < x[0] < 1.80 and not 0.00 < x[1] < 3.5 and not 0.0 < x[2] < 15. and not -0.3 < x[3] < 0.8:
    return np.inf
  
  return -2.*calc_ChiSq(c_double(x[0]), c_double(x[1]), c_double(x[2]), c_double(x[3]))

def lnprob3D(x):
  lp = lnprior3D(x)

  if not np.isfinite(lp):
    return np.inf
  
  return lp + calc_ChiSq(c_double(x[0]), c_double(x[1]), c_double(x[2]), c_double(0.0))

def lnprob4D(x):
    lp = lnprior4D(x)

    if not np.isfinite(lp):
        return np.inf

    return lp + calc_ChiSq(c_double(x[0]), c_double(x[1]), c_double(x[2]), c_double(x[3]))


FIELD =    1
DZ    =  0.3
LOZ   =  0.6
HIZ   =  0.9
KMAX  =  0.8
D0    = 1000

for FIELD in [1, 4]:
    for LOZ in [0.6, 0.9]:
        HIZ = LOZ + DZ

        # all Python types except integers, strings and unicode strings have to be wrapped to corresponding ctypes type. 
        likelihood.get_main(D0, FIELD, c_double(LOZ), c_double(HIZ), c_double(KMAX))

        ## Reset k remapping between FFTlog and mock (to kmax).
        likelihood.prep_ctype_ChiSq();

        # print "Chi sq. result: %.6lf" % likelihood.calc_ChiSq(c_double(0.5), c_double(1.0), c_double(7.), c_double(0.01))

        results = []
        
        x03  = np.array([0.5, 0.85, 7.0]) 
        x04  = np.array([0.0, 1.00, 8.0, 0.4])

        for i in xrange(1, 10, 1):    
            likelihood.get_ydata(i)

            result = minimize(chiSq_3D, x03, method='Nelder-Mead', tol=1e-6)
            # result = minimize(chiSq_4D, x04, method='Nelder-Mead', tol=1e-6)
    
            results.append(result.x)

            print i, result.x 

        maxlike_points = np.array(results)

        np.savetxt('/home/mjw/HOD_MockRun/emcee_log/maxlike_points_W%d_%.1lf_%.1lf_d0_%d.dat' % (FIELD, LOZ, HIZ, D0), maxlike_points, fmt='%.6e')

        print "Field W%d, %.1lf < z < %.1lf, for k_max = %.3lf and d0 = %d" % (FIELD, LOZ, HIZ, KMAX, D0)
        print "\n\nMean values: (fs8, bs8, sp) = (%.4lf, %.4lf, %.4lf)"     % (maxlike_points[:,0].mean(), maxlike_points[:,1].mean(), maxlike_points[:,2].mean())

        likelihood.print_model(c_double(maxlike_points[:,0].mean()), c_double(maxlike_points[:,1].mean()), c_double(maxlike_points[:,2].mean()), c_double(0.00))

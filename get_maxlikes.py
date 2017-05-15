#!/home/mlam/anaconda/bin/
import sys, os;       sys.path.append(os.environ['PYTHONPATH'])
import matplotlib;    matplotlib.use('PDF') # matplotlib.use('Agg')

from   matplotlib     import style; style.use('classic')
import mpi4py
from   emcee          import EnsembleSampler
from   emcee.utils    import MPIPool
from   scipy.optimize import minimize

from   ctypes         import cdll
from   ctypes         import c_double, c_int

import matplotlib     as mpl
import numpy          as np
import pylab          as pl 
import corner
import sys

#import getdist
#from   getdist        import loadMCSamples, plots

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
          
outputdir =          os.environ['outputdir']
FIELD     =   np.int(os.environ['FIELDFLAG'])
D0        =   np.int(os.environ['d0'])
LOZ       = np.float(os.environ['LOZ'])
HIZ       = np.float(os.environ['HIZ'])
KMAX      = np.float(os.environ['KMAX'])

mockNum   = 153

print D0, FIELD, LOZ, HIZ, KMAX

# all Python types except integers, strings and unicode strings have to be wrapped to corresponding ctypes type. 
likelihood.get_main(D0, FIELD, c_double(LOZ), c_double(HIZ), c_double(KMAX))

## Reset k remapping between FFTlog and mock (to kmax).
likelihood.prep_ctype_ChiSq();

print "Chi sq. result: %.6lf" % likelihood.calc_ChiSq(c_double(0.5), c_double(0.8), c_double(8.), c_double(0.01))

results = []

x03     = np.array([0.5, 0.85, 7.0]) 
x04     = np.array([0.5, 0.85, 8.0, 0.01])

# for i in xrange(1, mockNum, 1):    
for i in xrange(1, 2, 1):  ## model that best matches the mean. 
    likelihood.get_ydata(i, 2)

    result = minimize(chiSq_3D, x03, method='Nelder-Mead', tol=1e-6)
    # result = minimize(chiSq_4D, x04, method='Nelder-Mead', tol=1e-6)
    
    results.append(result.x)

    print i, result.x 
    
maxlike  = np.array(results)
meanlike = [maxlike[:,i].mean() for i in range(3)]

print "\n\nMean values: (fs8, bs8, sp) = (%.4lf, %.4lf, %.4lf)" % (meanlike[0], meanlike[1], meanlike[2])

'''
header=""

for i in meanlike:
    header=header + r"%.6lf  " % i
    
np.savetxt(outputdir + "/maxlikes/maxlike_points_W%d_%.1lf_%.1lf_d0_%d.dat" % (FIELD, LOZ, HIZ, D0), maxlike, fmt='%.6e', header=header)
'''

## print models. 
likelihood.print_model(c_double(meanlike[0]), c_double(meanlike[1]), c_double(meanlike[2]), c_double(0.00))


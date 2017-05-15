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
    if -0.2 < x[0] < 1.80 and 0.00 < x[1] < 3.5 and 0.0 < x[2] < 15. and -0.3 < x[3] < 0.8:
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
D0        =   np.int(os.environ['D0'])
LOZ       = np.float(os.environ['LOZ'])
HIZ       = np.float(os.environ['HIZ'])
KMAX      = np.float(os.environ['KMAX'])

mockNum   = 152

# all Python types except integers, strings and unicode strings have to be wrapped to corresponding ctypes type. 
likelihood.get_main(D0, FIELD, c_double(LOZ), c_double(HIZ), c_double(KMAX))

## Reset k remapping between FFTlog and mock (to kmax).
likelihood.prep_ctype_ChiSq();

print "Chi sq. result: %.6lf" % likelihood.calc_ChiSq(c_double(0.5), c_double(1.0), c_double(7.), c_double(0.01))

results = []

x03     = np.array([0.5, 0.85, 7.0]) 
x04     = np.array([0.5, 0.85, 8.0, 0.01])

likelihood.get_ydata(0, 1)

# result = minimize(chiSq_3D, x03, method='Nelder-Mead', tol=1e-6)
result = minimize(chiSq_4D, x04, method='Nelder-Mead', tol=1e-6)
    
results.append(result.x)

print
print result.x 
    
## necessary?
np.seterr(invalid='warn')

ndim, nwalkers, nburn, nsteps = 4, 45, 10, 50

seed = 1
np.random.seed(seed)

# Choose an initial set of positions for the walkers.
# p0   = np.random.rand(nwalkers, ndim)
# p0   = result.x*(1. + 0.1*(p0 - 0.5))
# p0   = [np.array(p0.tolist()[i]) for i in range(nwalkers)] 

p0  = np.array([[j for j in np.random.rand(ndim)]  for i in xrange(nwalkers)])
p0  = result.x*(1. + 0.05*(p0 - 0.5))

# sampler = EnsembleSampler(nwalkers, ndim, lnprob4D)

pool = MPIPool() # loadbalance=True

if not pool.is_master():
    pool.wait()
    sys.exit(0)

sampler = EnsembleSampler(nwalkers, ndim, lnprob4D, pool=pool)

pos, lnprob, rand_state   = sampler.run_mcmc(p0, nburn) 

sampler.reset() # Reset the chain to remove the burn-in samples; keep walker positions.

pos, lnprob, rand_stateR  = sampler.run_mcmc(pos, nsteps, rstate0=rand_state) # rstate0=rstate

if pool.is_master():
    meanacceptance            = np.mean(sampler.acceptance_fraction)
    autocorrelationtimes      = sampler.get_autocorr_time()

    # Print out the mean acceptance fraction. In general, acceptance_fraction has an entry for each walker so, in this case, it is a
    # 50-dimensional vector.
    print "Mean acceptance fraction:" + str(meanacceptance)

    # Estimate the integrated autocorrelation time for the time series in each parameter.
    print "Autocorrelation time:" + str(autocorrelationtimes)

    samples        = sampler.flatchain
    lnprob         = sampler.flatlnprobability

    length =  len(samples[:,0])

    lnpriors = []

    for i in xrange(0, length, 1):    
        lnpriors.extend([lnprior4D(samples[i,:])])
 
    getdist_format = np.column_stack((np.ones(len(lnprob)), np.array(lnpriors), lnprob, samples))

    np.savetxt(outputdir+'/emcee_log/data_W%d_zlim_%.1lf_%.1lf_kmax_%.1lf_chain_%d.txt' % (FIELD, LOZ, HIZ, KMAX, seed), getdist_format, fmt='%.3lf', header="Mean acceptance: " + str(meanacceptance)+",  Autocorrelation times: " + str(autocorrelationtimes))
    
    fig = corner.corner(samples, labels=["$f \sigma_8$", "$b \sigma_8$", "$\sigma_p$", "$\epsilon$"])
    fig.savefig(outputdir + "/plots/data_W%d_zlim_%.1lf_%.1lf_kmax_%.1lf_triangle.pdf")

pool.close()

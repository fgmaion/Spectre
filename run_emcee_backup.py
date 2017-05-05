#!/usr/bin/env python
from   __future__   import print_function
import matplotlib;  matplotlib.use('PDF') # matplotlib.use('Agg')

from   matplotlib   import style; style.use('classic')
import mpi4py
from   emcee        import EnsembleSampler

from   ctypes       import cdll
from   ctypes       import c_double, c_int

import matplotlib   as mpl
import numpy        as np
import pylab        as pl
import sys 
import corner


def lnprior(x):
  if 0.0 < x[0] < 1.80 and 0.00 < x[1] < 3.5 and 0.0 < x[2] < 15.:
    return 0.0
  return -np.inf

# First, define the probability distribution that you would like to sample.
def lnprob_example(x, mu, icov):
    lp = lnprior(x)
    
    if not np.isfinite(lp):
      return -np.inf  
    
    diff = x - mu

    return lp - np.dot(diff, np.dot(icov,diff))/2.0

# load vipers likelihood as shared library. 
likelihood = cdll.LoadLibrary("/home/mjw/HOD_MockRun/likelihood.so")
calc_ChiSq = likelihood.calc_ChiSq

# return types are integer by default. 
calc_ChiSq.restype = c_double

def lnprob(x):
  lp = lnprior(x)

  if not np.isfinite(lp):
    return -np.inf
  
  return lp + calc_ChiSq(c_double(x[0]), c_double(x[1]), c_double(x[2]))

'''
# all Python types except integers, strings and unicode strings have to be wrapped to corresponding ctypes type. 
likelihood.get_main(1000, 1, c_double(0.6), c_double(0.9), c_double(0.8))

## Reset k remapping between FFTlog and mock (to kmax).
likelihood.prep_ctype_ChiSq();
'''
#### pool = em.utils.MPIPool()

# We'll sample im 3-dimensions...
ndim = 3

# We'll sample with N walkers.
nwalkers = 100

# Choose an initial set of positions for the walkers.
p0 = [(1.80*np.random.rand(), 3.5*np.random.rand(), 15.0*np.random.rand())  for i in xrange(nwalkers)]

#if not pool.is_master():
#    pool.wait()
#    sys.exit(0)


# Initialize the sampler with the chosen specs.
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
## sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=pool)
'''
# Run 5000 steps as a burn-in.
pos, prob, state = sampler.run_mcmc(p0, 1000)

# Reset the chain to remove the burn-in samples; keep walker positions. 
sampler.reset()

nsteps = 2000

for i, result in enumerate(sampler.sample(pos, iterations=nsteps, rstate0=state)):
  if (i+1) % 100 == 0:
    print("Run-in: {0:5.1%}".format(float(i)/nsteps))
                        
# Print out the mean acceptance fraction. In general, acceptance_fraction
# has an entry for each walker so, in this case, it is a 50-dimensional
# vector.
print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))

# Estimate the integrated autocorrelation time for the time series in each
# parameter.
print("Autocorrelation time:", sampler.get_autocorr_time())

## pool.close()

pl.clf()
pl.hist(sampler.flatchain[:,0], 100)
pl.savefig("/home/mjw/HOD_MockRun/W1_Spectro_V7_4/plots/emcee.pdf")

samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

pl.clf()
fig = corner.corner(samples, labels=["$f \sigma_8$", "$b \sigma_8$", "$\sigma_p$"])
fig.savefig("/home/mjw/HOD_MockRun/W1_Spectro_V7_4/plots/triangle.pdf")
'''

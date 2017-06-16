import os,        sys 
import numpy   as np
import pandas  as pd
from   scipy.optimize import minimize_scalar

root   = os.environ['outputdir']
root  += "/data_v1.7/"

args   = sys.argv

mocks  = 153

lo_zs  = np.array([np.float(args[1]), np.float(args[2])])
hi_zs  = np.array([np.float(args[3]), np.float(args[4])])

def simple_chi2(amp, unclipped, clipped, variance):
  return np.sum((unclipped - amp*clipped)**2./variance)

for a, field in enumerate(["W1", "W4"]):
  for b, lo_z in enumerate(lo_zs):
    hi_z     = hi_zs[b]

    ## mocks
    filename = os.environ['outputdir'] + "/mocks_v1.7/" + "pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (1000, field, lo_z, hi_z)
    shot     = np.loadtxt(filename)[:,1] # high k estimate.

    ## get unfolded.
    files    = [os.environ['outputdir']+"/mocks_v1.7/"+"pk/d0_%d/%s/mock_" % (1000, field) + "{0:03}".format(i+1) + "_zlim_%.1lf_%.1lf_Jf_0.dat" % (lo_z, hi_z) for i in xrange(0, mocks, 1)]
    dflist   = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N'])  for fname in files]

    for k in range(len(dflist)):
        dflist[k]['P0'] -= shot[k]

    noclipdf = pd.concat(dflist)

    indices  = dflist[0]['k'].between(0.08, 0.20) # fit at high k.

    mean     = noclipdf.groupby(noclipdf.index).mean()
    var      = noclipdf.groupby(noclipdf.index).var()['P0'][indices].values
    
    ## Data
    filename = root + "pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (1000, field, lo_z, hi_z)
    shot     = np.loadtxt(filename)[0] # high k estimate. 

    ## get unfolded.
    files    = [root + "pk/d0_%d/%s/data" % (1000, field) + "_zlim_%.1lf_%.1lf_Jf_0.dat" % (lo_z, hi_z)]
    dflist   = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N'])  for fname in files]

    for k in range(len(dflist)):
      dflist[k]['P0'] -= shot

    noclipdf = pd.concat(dflist)

    indices  = dflist[0]['k'].between(0.08, 0.20) # fit at high k. 
    
    print field, "%.1lf < z < %.1lf" % (lo_z, hi_z)

    ## save amplitude correction of unity for d0 = 1000.
    np.savetxt(root + "pk_derivedprops/d0_%d/%s/suppression_zlim_%.1lf_%.1lf.dat" % (1000, field, lo_z, hi_z), np.column_stack(np.ones(1)), fmt="%.6le")
   
    for c, d0 in enumerate([10, 6, 4]):                            
      filename = root + "pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z)
      shot     = np.loadtxt(filename)[0] # high k estimate. 
        
      ## get unfolded. 
      files   = [root + "pk/d0_%d/%s/data" % (d0, field) + "_zlim_%.1lf_%.1lf_Jf_0.dat" % (lo_z, hi_z)]
      mflist  = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N'])  for fname in files]
      
      for k in range(len(mflist)):
        mflist[k]['P0'] -= shot
    
      clipdf  = pd.concat(mflist)

      results = []
           
      answer = minimize_scalar(simple_chi2, bounds=(1., 100.), args=(dflist[0]['P0'][indices], mflist[0]['P0'][indices], var), method='brent')

      results.append(answer.x)
      
      np.savetxt(root + "pk_derivedprops/d0_%d/%s/suppression_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z), np.column_stack(np.array(results)), fmt="%.6le")
    

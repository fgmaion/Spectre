import os, sys 
import numpy   as np
import pandas  as pd
from   scipy.optimize import minimize_scalar

root   = os.environ['outputdir']
mocks  = 152
labels = np.array([i for i in xrange(1, mocks + 1)])

args   = sys.argv

zs     = [0.75, 1.05]
fsig8s = [0.4907, 0.4570]

def cost_function(amp, unclipped, clipped):
    return np.sum((unclipped - amp*clipped)**2.)

for a, field in enumerate(["W1", "W4"]):
    for b, lo_z in enumerate([np.float(args[1]), np.float(args[2])]):
        hi_z   = lo_z + np.float(args[3])

        filename = root + "/mocks_v1.7/pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (1000, field, lo_z, hi_z)
        shot     = np.loadtxt(filename)[:,1]
        
        ## get unfolded.
        files    = [root + "/mocks_v1.7/pk/d0_%d/%s/mock_" % (1000, field) + "{0:03}".format(i+1) + "_zlim_%.1lf_%.1lf_Jf_0.dat" % (lo_z, hi_z) for i in xrange(0, mocks, 1)]
        dflist   = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N'])  for fname in files]

        for k in range(len(dflist)):
          dflist[k]['P0'] -= shot[k]

        noclipdf = pd.concat(dflist)

        indices  = dflist[0]['k'].between(0.05, 0.1)

        print field, "%.1lf < z < %.1lf" % (lo_z, hi_z)

        np.savetxt(root + "/mocks_v1.7/pk_derivedprops/d0_%d/%s/suppression_zlim_%.1lf_%.1lf.dat" % (1000, field, lo_z, hi_z), np.column_stack((labels, np.ones(mocks))), fmt="%d \t %.6le")
        
        for c, d0 in enumerate([10, 6, 4]):                            
            filename = root + "/mocks_v1.7/pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z)
            shot     = np.loadtxt(filename)[:,1]

            ## get unfolded. 
            files  = [root + "/mocks_v1.7/pk/d0_%d/%s/mock_" % (d0, field) + "{0:03}".format(i+1) + "_zlim_%.1lf_%.1lf_Jf_0.dat" % (lo_z, hi_z) for i in xrange(0, mocks, 1)]
            mflist = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N'])  for fname in files]

            for k in range(len(mflist)):
              mflist[k]['P0'] -= shot[k]
    
            clipdf = pd.concat(mflist)

            results = []
            
            for i in xrange(1, mocks, 1):
                answer = minimize_scalar(cost_function, bounds=(1., 100.), args=(dflist[i]['P0'][indices], mflist[i]['P0'][indices]), method='brent')
                
                results.append(answer.x)

            np.savetxt(root + "/mocks_v1.7/pk_derivedprops/d0_%d/%s/suppression_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z), np.column_stack((labels, np.array(results))), fmt="%d \t %.6le")

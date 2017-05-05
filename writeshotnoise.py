import os 
import numpy   as np
import pandas  as pd

root   = os.environ['outputdir']
mocks  = 152

zs     = [0.75, 1.05]
fsig8s = [0.4907, 0.4570]

for c, d0 in enumerate([1000, 10, 6, 4]):
  for a, field in enumerate(["W1", "W4"]):
    for b, lo_z in enumerate([0.6, 0.9]):
      hi_z   = lo_z + 0.3

      ## mocks
      files  = [root + "/mocks_v1.7/pk/d0_" + str(d0) + "/" + field + "/" + "mock_" + "{0:03}".format(i+1) + "_zlim_" + str(lo_z) + "_" + str(hi_z) + "_Jf_4.dat" for i in xrange(0,mocks,1)]
      dflist = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N'])  for fname in files]
      catdf  =  pd.concat(dflist)
      shot   =  np.array([dflist[i][dflist[i]['k'].between(1.5, 3.0, inclusive=True)]['P0'].mean() for i in range(len(dflist))]) # shot noise for each mock.

      np.savetxt(root + "/mocks_v1.7/pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z), shot, fmt='%.6le')

      ## data
      files  = [root + "/data_v1.7/pk/d0_%d/%s/data_zlim_%.1lf_%.1lf_Jf_4.dat" % (d0, field, lo_z, hi_z)] 
      dflist = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N'])  for fname in files]
      catdf  =  pd.concat(dflist)
      shot   =  np.array([dflist[i][dflist[i]['k'].between(1.5, 3.0, inclusive=True)]['P0'].mean() for i in range(len(dflist))]) # shot noise for each mock.

      np.savetxt(root + "/data_v1.7/pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z), shot, fmt='%.6le')
      

import os, sys, glob 
import numpy   as np
import pandas  as pd

root   = os.environ['outputdir']
mocks  = 153

args   = sys.argv

zs     = [0.75, 1.05]
fsig8s = [0.4907, 0.4570]

for c, d0 in enumerate([1000]): # , 10, 6, 4]):
  for a, field in enumerate(["W1", "W4"]):
    for b, lo_z in enumerate([np.float(args[1]), np.float(args[2])]):
      hi_z   = lo_z + np.float(args[3])
      
      ## mocks
      files  = [root + "/mocks_v1.7/pk/d0_" + str(d0) + "/" + field + "/" + "mock_" + "{0:03}".format(i+1) + "_zlim_" + str(lo_z) + "_" + str(hi_z) + "_Jf_4.dat" for i in xrange(0,mocks,1)]
      dflist = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N'])  for fname in files]
      catdf  =  pd.concat(dflist)
      shot   =  np.array([dflist[i][dflist[i]['k'].between(1.5, 3.0, inclusive=True)]['P0'].mean() for i in range(len(dflist))]) # shot noise for each mock.

      labels = np.array([i for i in xrange(1, mocks + 1)])

      ## Get <n> shot noise estaimates. 
      files = glob.glob(root + "/mocks_v1.7/pk_derivedprops/d0_%d/%s/nbarshotnoise_mocks_*_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z))

      estimates = np.loadtxt(files[0])  
            
      for i, filename in enumerate(files[1:]):
        estimates = np.concatenate((estimates, np.loadtxt(files[i])))

      df = pd.DataFrame(estimates, columns=["num", "gal shot", "rand shot"])   
      df.drop_duplicates(keep="first", inplace=True)
      df["num"] = df["num"].astype(int)
      
      d = {key:value for key, value in zip(df["num"].values, df["gal shot"].values + df["rand shot"].values)} ## total shot noise

      np.savetxt(root + "/mocks_v1.7/pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z), np.column_stack((labels, shot, np.array([d[x] for x in labels]))), fmt='%d \t %.6le \t %.6le', header=" labels,    high k estimate,    <n> estimate")

      for filename in files:
        os.remove(filename) 

      ## data
      files  = [root + "/data_v1.7/pk/d0_%d/%s/data_zlim_%.1lf_%.1lf_Jf_4.dat" % (d0, field, lo_z, hi_z)] 
      dflist = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N'])  for fname in files]
      catdf  =  pd.concat(dflist)
      shot   =  np.array([dflist[i][dflist[i]['k'].between(1.5, 3.0, inclusive=True)]['P0'].mean() for i in range(len(dflist))]) # shot noise for each mock.

      np.savetxt(root + "/data_v1.7/pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z), shot, fmt='%.6le')
      

import os, sys, glob 

import numpy   as np
import pandas  as pd


root   = os.environ['outputdir']
mocks  = 305

args   = sys.argv         ## Redshift limits from command line.

lo_zs  = np.array([np.float(args[1]), np.float(args[2])])
hi_zs  = np.array([np.float(args[3]), np.float(args[4])])

for c, d0 in enumerate([1000]): #[1000 , 10, 6, 4]):
  for a, field in enumerate(["W1", "W4"]):
    for b, lo_z in enumerate(lo_zs):
      hi_z   = hi_zs[b]
      
      ## mocks
      files  = [root + "/mocks_v1.7/pk/d0_" + str(d0) + "/" + field + "/" + "mock_" + "{0:03}".format(i+1) + "_zlim_" + str(lo_z) + "_" + str(hi_z) + "_Jf_4.dat" for i in xrange(0,mocks,1)]

      ## Upper limit of 3.0 previously (15/08/2017); scales could be decided by minimising difference with <n> estimate. 
      dflist = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N'])  for fname in files]
      catdf  =  pd.concat(dflist)
      shot   =  np.array([dflist[i][dflist[i]['k'].between(1.2, 2.3, inclusive=True)]['P0'].mean() for i in range(len(dflist))]) # shot noise for each mock.

      labels =  np.array([i for i in xrange(1, mocks + 1)])

      ## Get <n> shot noise estimates. 
      files     = glob.glob(root + "/mocks_v1.7/pk_derivedprops/d0_%d/%s/shotnoise/nbarshotnoise_mocks_*_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z)) ## not ordered. 
      estimates = np.loadtxt(files[0])  
            
      for i, filename in enumerate(files[1:]):
        estimates = np.vstack((estimates, np.loadtxt(files[i+1]))) ## i in enumerate still starts at 0. 

      df = pd.DataFrame(estimates, columns=["num", "gal shot", "rand shot"])   

      df.drop_duplicates(keep="first", inplace=True)

      df["num"] = df["num"].astype(int)
      
      d = {key:value for key, value in zip(df["num"].values, df["gal shot"].values + df["rand shot"].values)} ## total shot noise
      
      np.savetxt(root + "/mocks_v1.7/pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z), np.column_stack((labels, shot, np.array([d[x] for x in labels]))), fmt='%d \t %.6le \t %.6le', header=" labels,    high k estimate,    <n> estimate")

      ## Clean up original nbar shotnoise estimate files. 
      #for filename in files:
      #  os.remove(filename) 
      
      ## data
      files  = [root + "/data_v1.7/pk/d0_%d/%s/data_zlim_%.1lf_%.1lf_Jf_4.dat" % (d0, field, lo_z, hi_z)] 

      dflist = [pd.read_csv(fname, sep='\t', header=None, names=['k', 'P0', 'P2', 'N']) for fname in files]

      catdf  =  pd.concat(dflist)

      shot   =  [dflist[i][dflist[i]['k'].between(1.2, 2.3, inclusive=True)]['P0'].mean() for i in range(len(dflist))] # shot noise for each mock.

      ## Get <n> shot noise estimates.
      filename  = root + "/data_v1.7/pk_derivedprops/d0_%d/%s/nbarshotnoise_data_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z)
      estimate  = np.loadtxt(filename)
            
      nbar_shot = estimate[0] + estimate[1]  

      output    = np.array([[shot[0], nbar_shot]])
               
      np.savetxt(root + "/data_v1.7/pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z), output, header=" high k estimate,    <n> estimate", delimiter='\t', fmt="%.6lf \t %.6lf")
      
      print root + "/data_v1.7/pk_derivedprops/d0_%d/%s/shotnoise_zlim_%.1lf_%.1lf.dat" % (d0, field, lo_z, hi_z), output
      

import matplotlib.pylab        as plt
import matplotlib.pyplot
from   matplotlib.font_manager import FontProperties
from   matplotlib.ticker       import ScalarFormatter
from   matplotlib.ticker       import FixedFormatter
import pylab                   as pl
import numpy                   as np
import math, os
import glob, pickle

formatter = ScalarFormatter(useMathText=True)

fig_width_pt = 246.0*2 # Get this from LaTex using \the\columnwidth
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width  = fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize':10,
          'text.fontsize':6,
          'legend.fontsize':6,
          'xtick.labelsize':11.0,
          'ytick.labelsize':11.0,
          'figure.figsize':fig_size,
          'font.family': 'serif'}

pl.rcParams.update(params)
pl.clf()
pl.figure()
fig = pl.figure()
axes = pl.Axes(fig, [.2, .2, .7, .7])
fig.add_axes(axes)
axes.yaxis.set_major_formatter(formatter)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_000.dat')
MockCat      = MockCat[:,1]

for i in xrange(1, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))
   
for i in xrange(10, 500, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label='Ensemble Monopole', fmt='', color='k')


MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_000.dat')
MockCat      = MockCat[:,2]

for i in xrange(1, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))
   
for i in xrange(10, 500, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_kbin_0.010_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label='Ensemble Quadrupole', fmt='', color='r')

'''
## Integral constraint corrected. 
MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_icc_kbin_0.010_000.dat')
MockCat      = MockCat[:,1]

for i in xrange(1, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_icc_kbin_0.010_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))
   
for i in xrange(10, 500, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_icc_kbin_0.010_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,1], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

# pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label='Ensemble Monopole, integral constraint', fmt='', color='g')
'''
'''
## Integral constraint corrected. 
MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_icc_kbin_0.010_000.dat')
MockCat      = MockCat[:,2]

for i in xrange(1, 10, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_icc_kbin_0.010_00'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))
   
for i in xrange(10, 500, 1):
   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VipersMask_HOD_KaiserLorentzRSD_clipThreshold_1.0e+03_MaskedMultiples_reducedMu_icc_kbin_0.010_0'+str(i)+'.dat')
   MockCat  = np.vstack((MockIn[:,2], MockCat))

kvals        = MockIn[:,0]
Mean         = np.mean(MockCat, axis=0)
Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label='Ensemble Quadrupole, integral constraint', fmt='', color='b')
'''
#data         = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/fftlog_hodpk_NoRSD_4096_1.00e-10_1.00e+14_18Sep_cnvld_firstOrder_pk.dat')
#pl.loglog(data[:,0],         data[:,1], 'k-', label='input mono')
#pl.loglog(data[:,0], (2./7.)*data[:,1], 'r-', label='input quad')

### these include window terms and hence have an arbitrary overall amplitude currently. 
#pl.loglog(data[:,0],  0.45*data[:,2], 'y-', label='cnvld mono')
#pl.loglog(data[:,0],  1.05*data[:,3], 'g-', label='cnvld quad')

#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_000.dat')
#MockCat      = MockCat[:,1]

#for i in xrange(1, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))
#   
##for i in xrange(10, 21, 1):
##   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_0'+str(i)+'.dat')
##   MockCat  = np.vstack((MockIn[:,1], MockCat))

#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(9), label='no mask, mono', fmt='', color='k')


#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_000.dat')
#MockCat      = MockCat[:,2]

#for i in xrange(1, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,2], MockCat))
#   
##for i in xrange(10, 21, 1):
##   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_0'+str(i)+'.dat')
##   MockCat  = np.vstack((MockIn[:,1], MockCat))

#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(9), label='no mask, quad', fmt='', color='y')


#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_NoRSD_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_000.dat')
#MockCat      = MockCat[:,1]

#for i in xrange(1, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_NoRSD_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))
#   
#for i in xrange(10, 177, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_NoRSD_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_0'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))

#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

## pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(175), label='Ensemble Monopole', fmt='', color='k')


#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_Binarymask_HOD_NoRSD_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_000.dat')
#MockCat      = MockCat[:,1]

#for i in xrange(1, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_Binarymask_HOD_NoRSD_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))
#   
#for i in xrange(10, 500, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_Binarymask_HOD_NoRSD_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_0'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))

#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label='Ensemble Monopole', fmt='', color='r')


#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_Binarymask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_000.dat')
#MockCat      = MockCat[:,1]

#for i in xrange(1, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_Binarymask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))
#   
#for i in xrange(10, 380, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_Binarymask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_0'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))

#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(378), label='mask, mono', fmt='', color='g')

#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_Binarymask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_000.dat')
#MockCat      = MockCat[:,2]

#for i in xrange(1, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_Binarymask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,2], MockCat))
#   
#for i in xrange(10, 5000, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_Binarymask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_0'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,2], MockCat))

#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(4998), label='mask, quad', fmt='', color='r')

#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_001.dat')
#MockCat      = MockCat[:,1]

#for i in xrange(2, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))
#   
#for i in xrange(10, 500, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_0'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))

#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label='Ensemble Monopole', fmt='', color='k')

#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_001.dat')
#MockCat      = MockCat[:,2]

#for i in xrange(2, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,2], MockCat))
#   
#for i in xrange(10, 500, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_kbin_0.010_0'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,2], MockCat))
#   
#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(498), label='Ensemble Quadrupole', fmt='', color='g')


#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_CellSize_2.00_kbin_0.010_001.dat')
#MockCat      = MockCat[:,1]

#for i in xrange(2, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_CellSize_2.00_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))
#   
#for i in xrange(10, 79, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_CellSize_2.00_kbin_0.010_0'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))
#   
#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals[70:], Mean[70:], np.sqrt(Var[70:])/np.sqrt(77), fmt='', color='k')


#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_CellSize_2.00_kbin_0.010_001.dat')
#MockCat      = MockCat[:,2]

#for i in xrange(2, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_CellSize_2.00_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,2], MockCat))
#   
#for i in xrange(10, 79, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_knownGRFmask_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_CellSize_2.00_kbin_0.010_0'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,2], MockCat))
#   
#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals[60:], Mean[60:], np.sqrt(Var[60:])/np.sqrt(77), fmt='', color='g')


#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_CellSize_4.00_kbin_0.010_001.dat')
#MockCat      = MockCat[:,1]

#for i in xrange(2, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_CellSize_4.00_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,1], MockCat))

#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(8), fmt='', color='k')


#MockCat      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_CellSize_4.00_kbin_0.010_001.dat')
#MockCat      = MockCat[:,2]

#for i in xrange(2, 10, 1):
#   MockIn   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Multipoles_VIPERS_fullCube_HOD_QuadrupoleOnly_clipThreshold_1.0e+03_MaskedMultiples_fullMu_CellSize_4.00_kbin_0.010_00'+str(i)+'.dat')
#   MockCat  = np.vstack((MockIn[:,2], MockCat))

#kvals        = MockIn[:,0]
#Mean         = np.mean(MockCat, axis=0)
#Var          = np.var(MockCat,  axis=0)

#pl.errorbar(kvals, Mean, np.sqrt(Var)/np.sqrt(8), fmt='', color='r')
'''
data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Xi/Multipoles_200_fullCube_mockGals_RSDPk_CellSize_2.00_kbin_0.005_001.dat')
pl.loglog(data[:,0],               data[:,1],    'r^', markersize=2)
pl.loglog(data[:,0],  0.25*(7./6.)*data[:,4],    'r-', markersize=2)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Xi/Multipoles_200_fullCube_mockGals_RSDPk_grid_CellSize_2.00_kbin_0.005_001.dat')
pl.loglog(data[:,0],               data[:,1],    'y^', markersize=2)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Xi/Multipoles_200_fullCube_mockGals_RSDPk_displaced_D.4_CellSize_2.00_kbin_0.005_001.dat')
pl.loglog(data[:,0],               data[:,1],    'k^', markersize=2)
pl.loglog(data[:,0],  0.16*(7./6.)*data[:,4],    'k-', markersize=2)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Xi/Multipoles_200_fullCube_mockGals_RSDPk_displaced_D.1_CellSize_2.00_kbin_0.005_001.dat')
pl.loglog(data[:,0],               data[:,1],    'g^', markersize=2)
pl.loglog(data[:,0],  0.01*(7./6.)*data[:,4],    'g-', markersize=2)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Xi/Multipoles_200_fullCube_mockGals_RSDPk_displaced_D.01_CellSize_2.00_kbin_0.005_001.dat')
pl.loglog(data[:,0],                 data[:,1],    'c^', markersize=2)
pl.loglog(data[:,0],  0.0001*(7./6.)*data[:,4],    'c-', markersize=2)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Xi/Multipoles_200_fullCube_mockGals_RSDPk_displaced_D.01_CellSize_1.00_kbin_0.005_001.dat')
pl.loglog(data[:,0],                 data[:,1],    'm^', markersize=2)
pl.loglog(data[:,0],  0.0001*(7./6.)*data[:,4],    'm--', markersize=2)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Xi/Multipoles_200_fullCube_mockGals_RSDPk_displaced_D.1_CellSize_1.00_kbin_0.005_001.dat')
pl.loglog(data[:,0],               data[:,1],    'b^', markersize=2)
pl.loglog(data[:,0],  0.01*(7./6.)*data[:,4],    'b--', markersize=2)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Xi/Multipoles_200_fullCube_mockGals_RSDPk_displaced_D.1_CellSize_1.00_onemillion_kbin_0.005_001.dat')
pl.loglog(data[:,0],               data[:,1],    'b^', markersize=2)
pl.loglog(data[:,0],  0.01*(7./6.)*data[:,4],    'b--', markersize=2)
'''
pl.xscale('log')
pl.yscale('log')

pl.xlabel(r'$k [h Mpc^{-1}]$')
pl.ylabel(r'P(k) e$^{-9k^2/2}$')

pl.xlim(10.**-2., 1.0)
pl.ylim(10.**-3.,  6.*10.**4.)

pl.legend(loc=1)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/31JulyOnwards/integralConstraint.pdf')

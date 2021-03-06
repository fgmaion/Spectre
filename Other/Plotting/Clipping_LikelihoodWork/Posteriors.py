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
#formatter.set_scientific(True)
#formatter.set_powerlimits((-3,3))

fig_width_pt = 240.0*3 # Get this from LaTex using \the\columnwidth
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width  = fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize':14,
          'text.fontsize':8,
          'legend.fontsize':14,
          'xtick.labelsize':12.0,
          'ytick.labelsize':12.0,
          'figure.figsize':fig_size,
          'font.family': 'serif'}

pl.rcParams.update(params)
pl.clf()
pl.figure()
fig = pl.figure()
axes = pl.Axes(fig, [.2, .2, .7, .7])
fig.add_axes(axes)
axes.yaxis.set_major_formatter(formatter)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_betaPosterior_Res_30_kmax_0.20_GaussianRSD.dat')
#pl.plot(data[:,0], data[:,1], label=r'Gaussian,   $k_{max} = 0.2$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_betaPosterior_Res_30_kmax_0.20_LorentzianRSD.dat')
plt.bar(data[:,0], data[:,1], 0.0333*10**-1., label=r'$k_{max} = 0.2$', color='g', log=1)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_betaPosterior_Res_30_kmax_0.30_GaussianRSD.dat')
#pl.plot(data[:,0], data[:,1], label=r'Gaussian,   $k_{max} = 0.3$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_betaPosterior_Res_30_kmax_0.30_LorentzianRSD.dat')
plt.bar(data[:,0], data[:,1], 0.0333*10**-1., label=r'$k_{max} = 0.3$', color='y', log=1, alpha=0.6)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_betaPosterior_Res_30_kmax_0.40_GaussianRSD.dat')
#pl.plot(data[:,0], data[:,1], label=r'Gaussian,   $k_{max} = 0.4$')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_xvel_BootStrap_clipThreshold_1.0e+03_fullCube_betaPosterior_Res_30_kmax_0.40_LorentzianRSD.dat')
plt.bar(data[:,0], data[:,1], 0.0333*10**-1., label=r'$k_{max} = 0.4$', color='r', log=1, alpha=1.0)

pl.axvline(x=0.542, ymin=0, ymax=1,c='k')

pl.xlabel(r'$\beta$')
pl.ylabel(r'PDF')

pl.ylim(10.**-6., 1.000)

lg = pl.legend(loc=3)
lg.draw_frame(False)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/zCube_xvel_BootStrap_clipThreshold_1.0e+03_betaPosterior_Res_30_kmax_0.40_betaPosterior_LowRes_clippingDraft.pdf')



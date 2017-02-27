from   matplotlib.colors import LogNorm

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

fig_width_pt = 246.0*2 # Get this from LaTex using \the\columnwidth
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width  = fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean # height in inches
fig_size = [fig_width, fig_height]
params = {'axes.labelsize':10,
          'text.fontsize':8,
          'legend.fontsize':8,
          'xtick.labelsize':8.5,
          'ytick.labelsize':5.5,
          'figure.figsize':fig_size,
          'font.family': 'serif'}

pl.rcParams.update(params)
pl.clf()
pl.figure()
fig = pl.figure()
axes = pl.Axes(fig, [.2, .2, .7, .7])
fig.add_axes(axes)
axes.yaxis.set_major_formatter(formatter)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/VipersMask_2Dpk.dat')

plt.imshow(data[:,2].reshape(15, 39), origin='lower', interpolation='nearest', extent=[0.0, 0.8, 0.0, 0.8],  norm = LogNorm(vmin=10.**2., vmax=10.**4.), cmap='PRGn')

plt.colorbar()

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/VipersMask_2Dpk.pdf', bbox_inches='tight', pad_inches=0.5)

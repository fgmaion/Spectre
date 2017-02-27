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


randoms = np.loadtxt('/disk1/mjw/HOD_MockRun/Stefano/xyzRan.dat')

Chi     = randoms[:,0]**2 + randoms[:,1]**2 + randoms[:,2]**2
Chi     = Chi**0.5

BINS = np.arange(574.28, 2800., 16.)


data   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_MockAvg_nz_chiSliced_16.0_W1_22.00_50.0.dat')

def plot(alpha):
  pl.hist(Chi, bins=BINS)
  
  pl.plot(data[:,1], alpha*data[:,4])

  pl.savefig('/disk1/mjw/HOD_MockRun/Plots/VIPERS_RandomsNz.pdf') 

plot(2.0)


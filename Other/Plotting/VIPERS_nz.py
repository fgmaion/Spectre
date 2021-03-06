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

data   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_001_SingleMockNz_chiSliced_15.0_W1andW4_VolLim_22.00_Gaussfiltered_50.0.dat')
#plt.bar(data[:,1], data[:,4], width=15.0, label='mock 001', alpha=0.4, color='r')
#plt.bar(data[:,1], data[:,4], width=15.0, label='mock 001', alpha=0.4, color='r')

data   = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_MockAvg_nz_chiSliced_16.0_W1_22.00_50.0.dat')
pl.plot(data[:,1], data[:,5], label='Gaussian filtered', color = 'k')
plt.bar(data[:,1], data[:,3], width=15.0, label='mock avg.', color='g')

pl.xlabel(r'$\chi [h^{-1} Mpc]$')
pl.ylabel(r'n$(\chi)$')

lg = pl.legend(loc=1)
lg.draw_frame(False)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/nz_clippingDraft.pdf')

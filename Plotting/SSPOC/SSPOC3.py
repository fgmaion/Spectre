from matplotlib.colors import LogNorm
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
inches_per_pt = 1.0/72.
golden_mean = (np.sqrt(5)-1.0)/2.0
#fig_width  = fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean # height in inches
fig_width = fig_height

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

bSpoc  = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/Spoc_angularcorrelation.dat')

bnoSpoc = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/noSpoc_angularcorrelation.dat')

bfiber  = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/SpectralDistortion/fiberCollision_angularcorrelation.dat') 

brandoms = np.loadtxt('/disk1/mjw/HOD_MockRun/randoms20_W1_Nagoya_angularCorrelation.dat')

Spoc   = bSpoc[0:50, 0:150]
noSpoc = bnoSpoc[0:50, 0:150]
fiber  = bfiber[0:50, 0:150]

randoms = brandoms[0:50, 0:150]

alpha = 529839/66548

randoms /= alpha**2

fig = plt.figure()

ax = fig.add_subplot(111, aspect='equal')

wCorrelation = noSpoc/randoms - 1.

cax = ax.matshow(np.transpose(wCorrelation), origin='lower')

ra_labels = ((1./(240.*50.))*60.*60.)*np.arange(0, 50, 10)

dec_labels = ((1./(20.*50.))*60.*60.)*20.*np.arange(0, 8, 1)

#for i in xrange(0, 5, 1):
#    ra_labels[i]  = "%0.2f" % ra_labels[i]

#    dec_labels[i] = "%0.2f" % dec_labels[i]

print ra_labels

strings = np.array(["%.1f\'\'" % number for number in ra_labels])

ra_labels = np.concatenate((np.array(['']), strings))

strings = np.array(["%.1f\'\'" % number for number in dec_labels])

dec_labels = np.concatenate((np.array(['']), strings))

ax.set_xticklabels(ra_labels)
ax.set_yticklabels(dec_labels)

#plt.imshow(1. + data[0:50,0:50], norm=LogNorm(vmin=.0, vmax=300.))

fig.colorbar(cax)

pl.savefig('wCorrelation_noSpoc.pdf')

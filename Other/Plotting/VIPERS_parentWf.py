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

xdata       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/NGP_Mesh_2.00_xSlice_NGPCorr.dat')
xdata[:,0] -= 3.0
xdata       = 10.**xdata
pl.loglog(xdata[:,0], xdata[:,1], 'r', label=r'declination')

ydata       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/NGP_Mesh_2.00_ySlice_NGPCorr.dat')
ydata[:,0] -= 1.
ydata       = 10.**ydata
pl.loglog(ydata[:,0], ydata[:,1], 'g', label='right ascension')

zdata       = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/NGP_Mesh_2.00_zSlice_NGPCorr.dat')
zdata       = 10.**(zdata)
pl.loglog(zdata[:,0], zdata[:,1], 'k', label='redshift')

xStef = np.loadtxt('/disk1/mjw/HOD_MockRun/Stefano/window/windowKx_Slice.txt')
#pl.plot(xStef[:,0], xStef[:,1]/(10.**6.7991541764558505), 'r^', label='x')

yStef = np.loadtxt('/disk1/mjw/HOD_MockRun/Stefano/window/windowKy_Slice.txt')
#pl.plot(yStef[:,0], yStef[:,1]/(10.**6.7991541764558505), 'g^', label='y')

zStef = np.loadtxt('/disk1/mjw/HOD_MockRun/Stefano/window/windowKz_Slice.txt')
#pl.plot(zStef[:,0], zStef[:,1]/(10.**6.7991541764558505), 'k^', label='z')

pl.xlim(10.**-5, 2.)
pl.ylim(10.**-7., 1.)

pl.xlabel(r'k, [$hMpc^{-1}$]')
pl.ylabel(r'$W^2$(k)')

lg = pl.legend(loc=3)
lg.draw_frame(False)

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/W2k_SlicesClippingDraft.pdf')

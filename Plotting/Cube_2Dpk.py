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


fig = pl.figure()

ax0 = fig.add_subplot(1,2,1)

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Multipoles/Cartesian2Dpk_zCube_xvel_clipThreshold_1.0e+03_fullCube_0.010_0.dat')

plt.imshow(data[:,4].reshape(49, 49), origin='lower', interpolation='nearest', extent=[0.0, 0.5, 0.0, 0.5], norm = LogNorm(vmin=10.**3., vmax=6.*10.**4.), cmap='PRGn')

plt.colorbar(shrink=0.57)

kperp         = np.arange(0.0, 0.49, 0.01)
kpara         = np.arange(0.0, 0.49, 0.01)

KPERP, KPARA  = np.meshgrid(kperp, kpara)

MODK          = np.sqrt(KPERP**2 + KPARA**2)

MU            = KPARA/np.sqrt(KPERP**2 + KPARA**2)

beta          = 0.542

SIG           = 3.2

KAISERFACTOR  = (1. + beta*MU**2)**2

KAISERLORENTZ = KAISERFACTOR*(1. + 0.5*MODK**2*MU**2*SIG**2)**-1

Pk            = data[:,1].reshape(49, 49)

Contours      = np.logspace(2.4, 3.9, num=8, base=10)

plt.contour(KPERP, KPARA, Pk*KAISERFACTOR,  Contours[3:], colors='k', norm = LogNorm(vmin=10.**2., vmax=10.**4.))

plt.contour(KPERP, KPARA, Pk*KAISERLORENTZ, Contours[2:], colors='k', norm = LogNorm(vmin=10.**2., vmax=10.**4.), linestyles='dashed')

ax0.set_aspect(1.)

pl.xlabel(r'$k_{\perp}$')
pl.ylabel(r'$k_{\parallel}$')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/z2Pk/zCube_Clipped_1000.0_Cartesian.pdf', bbox_inches='tight', pad_inches=0.5)

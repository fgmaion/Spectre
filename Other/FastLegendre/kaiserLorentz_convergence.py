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
formatter.set_scientific(True)
formatter.set_powerlimits((-3,3))

fig_width_pt = 240.0             # Get this from LaTex using \the\columnwidth
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width  = fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean # height in inches
fig_size = [fig_width, fig_height]
params = {'font.size':10,
          'axes.labelsize':10,
          'text.fontsize':10,
          'legend.fontsize':7.0,
          'xtick.labelsize':8.0,
          'ytick.labelsize':8.0,
          'figure.figsize':fig_size,
          'font.family': 'serif',
          # 'font.serif':'Computer Modern Roman',
          }

pl.rcParams.update(params)
pl.clf()
pl.figure()
fig = pl.figure()

axes = pl.Axes(fig, [0.125, 0.2, 0.95-0.125, 0.95-0.2])
fig.add_axes(axes)

axes.yaxis.set_major_formatter(formatter)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

oct_order          = np.loadtxt('/disk1/mjw/HOD_MockRun//Scripts/FastLegendre/kaiserLorentz_vipers_window_cnvldpk_beta_0.5_sigma_5.0_oct2.dat')
hex_order          = np.loadtxt('/disk1/mjw/HOD_MockRun//Scripts/FastLegendre/kaiserLorentz_vipers_window_cnvldpk_beta_0.5_sigma_5.0_hex2.dat')
quad_order         = np.loadtxt('/disk1/mjw/HOD_MockRun//Scripts/FastLegendre/kaiserLorentz_vipers_window_cnvldpk_beta_0.5_sigma_5.0_quad2.dat')

quad_fracdiff_mono = np.abs(quad_order[:,1] - oct_order[:,1])/np.abs(oct_order[:,1])
quad_fracdiff_quad = np.abs(quad_order[:,2] - oct_order[:,2])/np.abs(oct_order[:,2])

hex_fracdiff_mono  = np.abs( hex_order[:,1] - oct_order[:,1])/np.abs(oct_order[:,1])
hex_fracdiff_quad  = np.abs( hex_order[:,2] - oct_order[:,2])/np.abs(oct_order[:,2])

pl.plot(oct_order[:,0], quad_fracdiff_mono, 'k-', label=r"$|(M_2' - M_6')/M_6'|$")
pl.plot(oct_order[:,0], quad_fracdiff_quad, 'k--', label=r"$|(Q_2' - Q_6')/Q_6'|$", dashes=[1,1])

pl.plot(oct_order[:,0], hex_fracdiff_mono, 'r-', label=r"$|(M_4' - M_6')/M_6'|$")
pl.plot(oct_order[:,0], hex_fracdiff_quad, 'r--', label=r"$|(Q_4' - Q_6')/Q_6'|$", dashes=[1,1])

pl.legend(loc=2, ncol=2)

pl.xscale('log')

pl.xlim(0.01, 1.)
pl.ylim(-0.001, 0.06)

pl.xlabel(r'k [$h$ Mpc$^{-1}$]')

pl.savefig('/disk1/mjw/HOD_MockRun/Scripts/FastLegendre/kaiserLorentz_convergence.pdf')
pl.savefig('/disk1/mjw/Dropbox/maskedRSD/LaTeX_stuff/kaiserLorentz_convergence.pdf', bbox_inches='tight')

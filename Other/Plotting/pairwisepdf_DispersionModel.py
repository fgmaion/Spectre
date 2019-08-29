data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/pairwisepdf_DispersionModel.dat')

a    = plt.hist(data, bins=100, normed=True)

xbin  = -29.6271077 - -30.24586 

dispersion = 2.                   # sigma for one galaxy
pairwise   = 2.*dispersion**2.    # pair velocity pdf given by convolution, -> 'errors add in quadrature'

xvals =      np.arange(-40., 40., 0.1)
yvals = xbin*np.exp(-0.5*xvals*xvals/pairwise)*(2.*np.pi*pairwise)**-0.5

pl.plot(xvals, yvals, 'r')

pl.xlabel(r'v$_{\parallel}$ [h$^{-1}$Mpc]')
pl.ylabel(r'p(v$_{\parallel}$)')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/pairwisepdf_DispersionModel.pdf')

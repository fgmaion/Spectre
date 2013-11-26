data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_nz.dat')

# column 0: Lower redshift bound of bin. 
# column 1: Number of galaxies in redshift bin. 
# column 2: Comoving volume in the redshift interval of the bin.
# column 3: Comoving number density. 

pl.clf()
pl.plot(data[:,0], data[:,3], '^')
pl.xlabel(r'$z$',    fontsize = '10')
pl.ylabel(r'$\overline{n}(z)$', fontsize = '10')
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/FKPweights/HODMocks_W1_nz.pdf')

pl.clf()
pl.plot(data[:,0], data[:,1], '^')
pl.xlabel(r'$z$',    fontsize = '10')
pl.ylabel(r'$N(z)$', fontsize = '10')
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/FKPweights/HODMocks_W1_Nz.pdf')

pl.clf()
pl.semilogy(data[:,0], data[:,2])
pl.xlabel(r'$z$',    fontsize = '10')
pl.ylabel(r'$V(z)$', fontsize = '10')
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/FKPweights/HODMocks_W1_ComovVolume_z.pdf')


fkpPk = 2000.

pl.clf()
pl.plot(data[:,0], data[:,3]*fkpPk/(1. + data[:,3]*fkpPk))
pl.xlabel(r'$z$',    fontsize = '10')
pl.ylabel(r'$w(z)$', fontsize = '10')
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/FKPweights/HODMocks_W1_zFKPweights2000.pdf')

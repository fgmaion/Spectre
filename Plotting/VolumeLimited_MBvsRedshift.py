Catalogue = np.loadtxt('/disk1/mjw/VIPERS_ValueAddedHOD/mocks_W1_v1.2/mock_W1_001_ALLINFO.cat', usecols=[0,1,2,3,4,5,6,7,8])

pl.hexbin(Catalogue[:,3], Catalogue[:,7], bins='log', cmap='Oranges')

pl.xlabel('$z$')
pl.ylabel('$M_B(z)$')

pl.axvline(x=0.7, ymin=0, ymax=1, c='k')
pl.axvline(x=0.9, ymin=0, ymax=1, c='k')
pl.axhline(y=-20.15, xmin=0, xmax=1, c='k')

cb = plt.colorbar()
cb.set_label('log$_{10}$(Counts)')

plt.gca().invert_yaxis()
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/HOD_MocksVolLimSampleSelection.pdf')

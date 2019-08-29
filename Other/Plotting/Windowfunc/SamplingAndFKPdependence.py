data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/noFKP_Cic_Mesh_7.81_fkpPk_1000.00_sampl_1.00_xSlice_NGPuncorr.dat')
pl.plot(data[:,0], data[:,1], 'g', label='no FKP')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/noFKP_Cic_Mesh_7.81_fkpPk_1000.00_sampl_1.00_ySlice_NGPuncorr.dat')
pl.plot(data[:300,0], data[:300,1], 'g')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/noFKP_Cic_Mesh_7.81_fkpPk_1000.00_sampl_1.00_zSlice_NGPuncorr.dat')
pl.plot(data[:200,0] -1, data[:200,1], 'g')


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/FKP_Cic_Mesh_7.81_fkpPk_5000.00_sampl_1.00_xSlice_NGPuncorr.dat')
pl.plot(data[:,0], data[:,1], 'c', label='fkp Pk, 5000., 100% sampling')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/FKP_Cic_Mesh_7.81_fkpPk_5000.00_sampl_1.00_ySlice_NGPuncorr.dat')
pl.plot(data[:300,0], data[:300,1], 'c')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/FKP_Cic_Mesh_7.81_fkpPk_5000.00_sampl_1.00_zSlice_NGPuncorr.dat')
pl.plot(data[:200,0] -1, data[:200,1], 'c')


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/FKP_Cic_Mesh_7.81_fkpPk_1000.00_sampl_1.00_xSlice_NGPuncorr.dat')
pl.plot(data[:,0], data[:,1], 'm', label='fkp Pk, 1000., 100% sampling')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/FKP_Cic_Mesh_7.81_fkpPk_1000.00_sampl_1.00_ySlice_NGPuncorr.dat')
pl.plot(data[:300,0], data[:300,1], 'm')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/FKP_Cic_Mesh_7.81_fkpPk_1000.00_sampl_1.00_zSlice_NGPuncorr.dat')
pl.plot(data[:200,0] -1, data[:200,1], 'm')


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/FKP_Cic_Mesh_7.81_fkpPk_1000.00_sampl_0.40_xSlice_NGPuncorr.dat')
pl.plot(data[:,0], data[:,1], 'y', label='fkp Pk, 1000., 40% sampling')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/FKP_Cic_Mesh_7.81_fkpPk_1000.00_sampl_0.40_ySlice_NGPuncorr.dat')
pl.plot(data[:300,0], data[:300,1], 'y')

data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/FKP_Cic_Mesh_7.81_fkpPk_1000.00_sampl_0.40_zSlice_NGPuncorr.dat')
pl.plot(data[:200,0] -1, data[:200,1], 'y')

pl.legend(loc=3)
pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc_VIPERSparent_FKPweights_sampling_fkpPk_dependence.pdf')


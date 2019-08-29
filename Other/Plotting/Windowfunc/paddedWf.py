def plotWindowFunc(surveyType):
    plotSphericalAnalytic()

    #for i in ['x', 'y', 'z']:
    #    WindowFunc_Slice = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSlices/'+ surveyType +'_' + i + 'Slice.dat')
    #    pl.loglog(WindowFunc_Slice[:,1], WindowFunc_Slice[:,2], '^', label=i)
        
    WindowFunc  = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSpherical/midK_W2k_' + surveyType + '.dat')
    pl.loglog(WindowFunc[:,0], WindowFunc[:,1], 'r^', label='Spherical average')
            
    padWindowfn = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/WindowfuncSpherical/midK_pad3W2k_' + surveyType + '.dat')
    pl.loglog(padWindowfn[:,1], padWindowfn[:,2], 'g^', label='padded spherical average')
            
    xx, locs = plt.xticks()
    ll = ['%.3f' % a for a in xx]
    plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

    leg = pl.legend(loc =2, ncol=1, prop = FontProperties(size = '10'))
    leg.draw_frame(False)
    
    pl.xlim(0.0001, 0.1)
    pl.ylim([10**-6, 1.0])
    
    pl.xlabel('k')
    pl.ylabel(r'$W^2$')
    pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/' + surveyType + '/Windowfunc.eps')
    
def plotSphericalAnalytic():
    Radius              = 250.
    Volume              = (4./3.)*np.pi*250.**3
    
    kVals               = 0.0001 + 0.0001*np.arange(10000)
    y                   = kVals*Radius
    SphericalAnalytic   = (3./y**3)*(np.sin(y) - y*np.cos(y))
    
    SphericalAnalyticW2 = SphericalAnalytic**2 
    
    pl.clf()
    pl.loglog(kVals, SphericalAnalyticW2, 'y', label = 'analytic')
    
    xx, locs = plt.xticks()
    ll = ['%.3f' % a for a in xx]
    plt.gca().xaxis.set_major_formatter(FixedFormatter(ll))

    leg = pl.legend(loc =2, ncol=1, prop = FontProperties(size = '10'))
    leg.draw_frame(False)
    
    pl.xlabel('k')
    pl.ylabel(r'$W^2$')
    pl.savefig('/disk1/mjw/HOD_MockRun/Plots/Windowfunc/Spherical/Analytic.eps')
    
    
plotWindowFunc('Spherical')


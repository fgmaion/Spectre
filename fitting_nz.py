import numpy as np
import math
import scipy
import pylab as pl

from scipy import optimize

def CSR(z):
    return 0.5*(1. - math.erf(17.465*(0.424 - z)))

def N(params, z):
    A  = params[0]
    a  = params[1]
    b  = params[2]
    z0 = params[3]

    return A*((z/z0)**a)*math.exp(-1.*(z/z0)**b)*CSR(z)

def Ndz(params,z, dz):
    return N(params, z)*dz

def Chi2(params, zdata, ydata, errors, dz):
    Interim   = 0.0

    for i in xrange(lowLimit, HiLimit, 1):
        theory    = Ndz(params, zdata[i], dz)

        # Divides by zero.
        # Errors[i] = np.sqrt(theory)
        
        diff      = ydata[i] - theory

        Interim  += (diff**2)/Errors[i]**2

    return Interim 


def testfunc(params):
    return(params[0] - 2.)**2 + (params[1]-3.)**2 + (params[2]-5.)**2

    
data      = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_MockAvg_nz_0.03.dat')

zbinwidth = 0.03

# dz = 0.03 
lenData  = 27
lowLimit = 13
HiLimit  = 40

# dz = 0.06
# lenData  = 14
# lowLimit = 6
# HiLimit  = 20

midz   = data[:,0]
Counts = data[:,2]

Errors = np.ones(len(Counts))

for i in xrange(lowLimit, HiLimit, 1):
    Errors[i] = np.sqrt(Counts[i])
 
dof    = lenData - 4

def redChi2(params):
    return Chi2(params, midz, Counts, Errors, zbinwidth)/dof

# A, alpha, beta, z0

initParams = np.array([35.*3.103, 8.603, 1.448,0.191])
#initParams = np.array([3.103, 10.603, 5.448, 1.2])


minparams, fmin, xi, direc, iter, funcalls, allvecs = scipy.optimize.fmin_powell(redChi2, initParams, full_output=True, disp=True, retall=True, xtol=0.0001, ftol=0.0001)

allvecs = np.array(allvecs)

print minparams, fmin

# inittest   = np.array([1.0, 2.0, 3.0])

# minparams = scipy.optimize.fmin(testfunc, inittest)

abscissa = 0.2 + (zbinwidth/10.)*np.arange(400)
bestfit  = np.zeros(400)

for i in xrange(0, 400, 1):
    bestfit[i] = Ndz(minparams, abscissa[i], zbinwidth)

pl.clf()
ax        = pl.subplot(111)

ax.bar(midz, Counts, color='k', width=zbinwidth, label='mock avg.',  alpha=0.3, align='center')

pl.plot(abscissa, bestfit, 'b')
pl.xlim([0.2, 1.2])

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/MockAvgNz_0.03_PowellsPython.pdf')

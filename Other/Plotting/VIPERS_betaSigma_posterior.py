import copy
import scipy.interpolate

from  numpy.linalg import eig, inv

from   matplotlib.colors import LogNorm

cmap = copy.copy(matplotlib.cm.jet)

# cmap.set_bad('w',1.)


data = np.loadtxt('/disk1/mjw/HOD_MockRun//Data/Posteriors/Multipoles_zCube_clipThreshold_1.0e+03_subVol_betaSigmaPosterior_lowRes.dat')

plt.imshow(data, origin='lower', cmap='YlGnBu', extent=[0.3, 0.60, 2.7, 3.30], interpolation='nearest') # , norm = LogNorm(vmin=data.min(), vmax=data.max()))

# Maximum likelihood
# pl.plot(0.4633, 3.04, 'kx', label='Max. likelihood')

# True value
pl.plot(0.45, 3.00, 'rx', label='True')

# plt.axis('equal')

pl.xlabel(r'$\beta$')
pl.ylabel(r'$\sigma$')

pl.legend(loc=3, numpoints=1)

plt.colorbar()


# Fit and plot an ellipse to the projected chi sq. 

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a
    
    
def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])


def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))


def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])
    
    
arc = 2.
R   = np.arange(0, arc*np.pi, 0.01)

'''
x   = 1.5*np.cos(R) + 2 + 0.1*np.random.rand(len(R))
y   = np.sin(R) + 1. + 0.1*np.random.rand(len(R))


a = fitEllipse(x,y)
center = ellipse_center(a)
phi = ellipse_angle_of_rotation(a)
axes = ellipse_axis_length(a)
'''

confLimits = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/Posteriors/Multipoles_zCube_clipThreshold_1.0e+03_subVol_2DconfidenceLimits_dChiSq_2.71_lowRes.dat')
pl.plot(confLimits[:,0], confLimits[:,1], '^')

a = fitEllipse(confLimits[:,0], confLimits[:,1])
center = ellipse_center(a)
phi = ellipse_angle_of_rotation(a)
axes = ellipse_axis_length(a)

a, b = axes
xx = center[0] + a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi)
yy = center[1] + a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi)

# pl.plot(x, y)
pl.plot(xx, yy, c='r')

pl.savefig('/disk1/mjw/HOD_MockRun/Plots/zCube_clipThreshold_1.0e+03_subVol_betaSigmaPosterior_lowRes.pdf', bbox_inches='tight', pad_inches=0.5)

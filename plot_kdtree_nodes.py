import numpy as np
import matplotlib.pyplot as plt

from   itertools import product, combinations
from   mpl_toolkits.mplot3d import Axes3D


fig = plt.figure()

ax  = fig.add_subplot(111, projection='3d')

# ax.set_aspect("equal")


particles = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/Node_particle_positions.dat')

for k in xrange(0, len(particles[:,0]), 1):
    ax.scatter(particles[k,0], particles[k,1], particles[k,2], color="k", s=2)


data = np.loadtxt('/disk1/mjw/HOD_MockRun/Data/stacpolly/kdtree_nodeLimits_children.dat')

for j in xrange(0, len(data[:,0]), 1):
    r1 = [data[j,1], data[j, 4]]
    r2 = [data[j,2], data[j, 5]]
    r3 = [data[j,3], data[j, 6]]
    
    for s, e in combinations(np.array(list(product(r1,r2,r3))), 2):
        if np.sum(np.abs(s-e)) == r1[1]-r1[0]:
          ax.plot(*zip(s,e), color="b")

        if np.sum(np.abs(s-e)) == r2[1]-r2[0]:
          ax.plot(*zip(s,e), color="b")

        if np.sum(np.abs(s-e)) == r3[1]-r3[0]:
          ax.plot(*zip(s,e), color="b")


ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.set_xlim3d(0., 550.)
ax.set_ylim3d(0., 550.)
ax.set_zlim3d(0., 550.)

pl.show()

pl.savefig('kdtree_nodes.pdf')

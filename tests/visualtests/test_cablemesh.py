import sys
sys.path.append('../data')

import numpy as np
import matplotlib.pyplot as plt
from testcases import read_off
from orthogonalstreamlines import OrthogonalStreamlines, create_orthogonal_streamlines_mesh
from orthogonalstreamlines.intersection import pack_cables
import mtlheartsim as mtl
from mtlheartmesh.trisurf import TriSurf

ver, tri = read_off('../data/mesh2d.off')
# orient = np.fromfile('../data/mesh2d.orient', sep=' ').reshape((-1, 3))
# orient = 0.5*np.random.normal(size=orient.size).reshape((-1, 3))
# orient.tofile('orient')
orient = np.fromfile('../data/mesh2d.orient2', dtype=float).reshape((-1, 3))

mesh = create_orthogonal_streamlines_mesh(ver, tri, orient, 0.7, options={'random_seed': 1663374316})

x = mesh.vertices
for c in mesh.cables[:mesh.nlc]:
    plt.plot(x[c, 0], x[c, 1], c='C0')

for c in mesh.cables[mesh.nlc:]:
    plt.plot(x[c, 0], x[c, 1], c='C1')

plt.show()
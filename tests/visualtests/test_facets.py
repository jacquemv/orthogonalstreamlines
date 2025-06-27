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
orient = np.fromfile('../data/mesh2d.orient', sep=' ').reshape((-1, 3))
mesh = create_orthogonal_streamlines_mesh(ver, tri, orient, 0.7, options={'random_seed': 1663374316})

x = mesh.vertices
tri = mesh.triangles

print(mesh.info)

ngon = np.concatenate([np.full(f.shape[0], n) for n, f in enumerate(mesh.facets)])
start = np.cumsum([f.shape[0] for f in mesh.facets])

ngon = ngon[mesh.tri_to_facet]
facetid = mesh.tri_to_facet - start[ngon]

for i, t in enumerate(tri):
    f = mesh.facets[ngon[i]][facetid[i]]
    assert np.all(np.isin(t, f))


for n, f in enumerate(mesh.facets):
    for k in f:
        ke = np.concatenate((k, [k[0]]))
        plt.fill(x[ke, 0], x[ke, 1])
        plt.plot(x[ke, 0], x[ke, 1], 'w')

plt.show()
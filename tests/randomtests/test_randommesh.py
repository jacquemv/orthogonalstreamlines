import sys
sys.path.append('../data')

import numpy as np
from testcases import read_off, read_case, write_case, plot_case
from orthogonalstreamlines import create_orthogonal_streamlines_mesh
from orthogonalstreamlines.intersection import tri_normals
from orthogonalstreamlines.tessellation import facet_normals

def point_in_triangle(point, triangle):
    x, y = point
    (ax, ay), (bx, by), (cx, cy) = triangle
    side_1 = (x - bx) * (ay - by) - (ax - bx) * (y - by)
    side_2 = (x - cx) * (by - cy) - (bx - cx) * (y - cy)
    side_3 = (x - ax) * (cy - ay) - (cx - ax) * (y - ay)
    return (side_1 < 0.0) == (side_2 < 0.0) == (side_3 < 0.0)


ver, tri = read_off('../data/mesh2d.off')

niter = 0
while True:
    niter += 1
    print(f'ITERATION {niter}')
    orient = np.random.random((tri.shape[0], 3)) - 0.5
    radius = 0.2
    mesh = create_orthogonal_streamlines_mesh(ver, tri, orient, radius, 
                                              verbose=True, 
                                              add_ghost_nodes=True)

    # check types and shapes
    nv = mesh.vertices.shape[0]
    nt = mesh.triangles.shape[0]
    ne = sum(len(c)-1 for c in mesh.cables)
    assert isinstance(mesh.cables, list)
    assert isinstance(mesh.facets, list)
    assert isinstance(mesh.boundaries, list)
    assert mesh.vertices.shape[1] == 3
    assert mesh.triangles.shape[1] == 3
    assert mesh.neighbors.shape[1] == 4
    assert mesh.ver_to_orig_tri.shape == (nv,)
    assert mesh.tri_to_facet.shape == (nt,)
    assert mesh.sign.shape == (nv,)
    assert mesh.dx.shape == (ne,)

    # cables must be either longitudinal or transverse
    assert mesh.nlc + mesh.ntc == len(mesh.cables)

    # all cables must have a length >= 2
    assert min(len(c) for c in mesh.cables) >= 2

    # cables in the same direction must not intersect
    cables_long = np.concatenate(mesh.cables[:mesh.nlc])
    cables_trans = np.concatenate(mesh.cables[mesh.nlc:])
    assert np.bincount(cables_long).max() <= 1
    assert np.bincount(cables_trans).max() <= 1

    # cables must cover all the nodes
    cables = np.concatenate(mesh.cables)
    cnt = np.bincount(cables)
    assert cnt.min() >= 1
    assert cnt.size == nv

    # a vertex must not be part of more than 2 cables
    assert cnt.max() <= 2

    # the vertices must be inside the triangle they are associated with
    for i, x in enumerate(mesh.vertices):
        if mesh.is_node[i] == 0:
            continue
        j = mesh.ver_to_orig_tri[i]
        Tj = [ver[tri[j, k], :2] for k in range(3)]
        assert point_in_triangle(x[:2], Tj)

    # facets must have the appropriate number of nodes
    for n, f in enumerate(mesh.facets):
        assert f.shape[1] == n
    
    # each vertex must have at least 2 neighbors
    assert np.all(np.sum(mesh.neighbors != -1, 1) >= 2)

    # each edge creates a pair of neighbors
    assert np.sum(mesh.neighbors != -1) == 2*ne

    # check neighbors
    idx = np.full(nv, -1)
    idx[cables_long] = np.arange(cables_long.size)
    # column 0 of neighbors is the previous node along the longitudinal cable
    assert np.all((idx - idx[mesh.neighbors[:, 0]] == 1) 
                  | (mesh.neighbors[:, 0] == -1))
    # column 1 of neighbors is the next node alongn the longitudinal cable
    assert np.all((idx[mesh.neighbors[:, 1]] - idx == 1) 
                  | (mesh.neighbors[:, 1] == -1))
    idx[:] = -1
    idx[cables_trans] = np.arange(cables_trans.size)
    # column 2 of neighbors is the previous node along the transverse cable
    assert np.all((idx - idx[mesh.neighbors[:, 2]] == 1) 
                  | (mesh.neighbors[:, 2] == -1))
    # column 3 of neighbors is the next node along the transverse cable
    assert np.all((idx[mesh.neighbors[:, 3]] - idx == 1) 
                  | (mesh.neighbors[:, 3] == -1))

    # facets must be consistently oriented
    for n, norm in enumerate(facet_normals(mesh.vertices, mesh.facets)):
        if n < 3: continue
        assert np.all(norm[:, 2] > 0)
    
    # triangles must be associated with the correct facet
    nf = 0
    for n, f in enumerate(mesh.facets):
        for i, nodes in enumerate(f):
            # triangulation of the facet
            T = np.where(mesh.tri_to_facet == nf + i)[0]
            assert np.all(np.isin(mesh.triangles[T].ravel(), nodes))
        nf += f.shape[0]
    
    # triangles must be consistently oriented
    norm = tri_normals(mesh.vertices, mesh.triangles)
    assert np.all(norm[:, 2] >= 0)

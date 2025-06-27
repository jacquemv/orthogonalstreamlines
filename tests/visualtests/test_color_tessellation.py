import sys
sys.path.append('../data')

import numpy as np
from testcases import read_off, read_case, write_case, plot_case
from orthogonalstreamlines.intersection import create_cable_network, unpack_cables, tri_normals
from orthogonalstreamlines import tessellation, triangulation
from evenlyspacedstreamlines import evenly_spaced_streamlines
import matplotlib.pyplot as plt


ver, tri = read_off('../data/mesh2d.off')
orient = np.random.random((tri.shape[0], 3)) - 0.5
lines1, faces1, info1 = evenly_spaced_streamlines(ver, tri, orient, 0.2, orthogonal=False)
lines2, faces2, info2 = evenly_spaced_streamlines(ver, tri, orient, 0.2, orthogonal=True)
normal = tri_normals(ver, tri)
out = create_cable_network(normal, lines1, faces1, lines2, faces2, add_ghost_nodes=True)
facets, vneigh = tessellation.tessellate(out.cables, out.cables_len, out.nc_long, 
                                            out.sign, return_moves=True, return_neigh=True)
_, cutoff = tessellation.find_largest_facets(facets, 1)
triangles, facetid, triangulation_failures = \
    triangulation.triangulate_facets(out.vertices, facets, cutoff)
facets = [f[0] for f in facets if f[0].size < cutoff]


def adjacency_matrix(facets):
    adj = [set() for i in range(len(facets))]
    edges = {}
    for i, f in enumerate(facets):
        for j in range(len(f)):
            j1 = (j+1) % len(f)
            pair = f[j], f[j1]
            pair = min(pair), max(pair)
            if pair in edges:
                k = edges[pair]
                adj[i].add(k)
                adj[k].add(i)
            else:
                edges[pair] = i
    return adj

def greedy_coloring(adj):
    V = len(adj)
    result = np.full(V, -1, dtype=np.int32)
    result[0] = 0
    available = np.full(V, 0, dtype=np.uint8)
    for u in range(1, V):
        # Process all adjacent vertices and
        # flag their colors as unavailable
        for i in adj[u]:
            if (result[i] != -1):
                available[result[i]] = 1
        # Find the first available color
        cr = 0
        while cr < V:
            if (available[cr] == 0):
                break
            cr += 1
        # Assign the found color
        result[u] = cr 
        # Reset the values back to false 
        # for the next iteration
        for i in adj[u]:
            if (result[i] != -1):
                available[result[i]] = 0
    # Print the result
    return np.array(result)


print([f.shape for f in facets])
facet_color = tessellation.facet_coloring(facets)
tri_color = facet_color[facetid]

X = out.vertices
new_normals = tri_normals(X, triangles)
T = triangles
plt.tripcolor(X[:, 0], X[:, 1], facecolors=tri_color, triangles=T, cmap='hsv', clim=(0, tri_color.max()+1))
cables = unpack_cables(out.cables, out.cables_len)
for c in cables:
    plt.plot(X[c, 0], X[c, 1], '-', c='0.6', lw=2)
I = out.is_node == 1
plt.plot(X[I, 0], X[I, 1], 'o', c='0.4', ms=3)
#for i, f in enumerate(facets):
#    plt.text(X[f, 0].mean(), X[f, 1].mean(), str(i))
plt.show()

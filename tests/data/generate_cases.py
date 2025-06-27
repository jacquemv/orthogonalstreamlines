from itertools import chain
import numpy as np
from testcases import read_off, read_case, write_case, plot_case
from orthogonalstreamlines.intersection import create_cable_network, unpack_cables, tri_normals
from orthogonalstreamlines import tessellation, triangulation
from evenlyspacedstreamlines import evenly_spaced_streamlines
import matplotlib.pyplot as plt
import pickle
import surfacetopology


def create_case_simple1tri():
    ver = np.array([
        [0, 0, 0],
        [2, 1, 0],
        [1, 15, 0]
    ])
    tri = np.array([[0, 1, 2]])
    orient = np.array([[0.1, 1, 0]])
    lines1, faces1, info1 = evenly_spaced_streamlines(ver, tri, orient, 0.3, orthogonal=False)
    lines2, faces2, info2 = evenly_spaced_streamlines(ver, tri, orient, 0.3, orthogonal=True)
    write_case('simple1tri.case', ver, tri, orient, lines1, faces1, lines2, faces2)

def create_case_duplicates():
    ver = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 1, 0],
        [-1, 0, 0],
    ])
    tri = np.array([[0, 1, 2], [0, 2, 3]])
    orient = np.array([[1, 2, 0], [1, 2, 0]])
    lines1 = [np.array([[-0.25, -0.75, 0], [0, 0, 0], [0.25, 0.75, 0]]), 
              np.array([[-0.25, 0.75, 0], [-0.5, -0.5, 0]])]
    faces1 = [np.array([1, 0]), np.array([1])]
    lines2 = [np.array([[-0.75, 0.25, 0], [0, 0, 0], [0.75, -0.25, 0]]), 
              np.array([[0.5, 0.5, 0], [0, 0.5, 0], [-0.5, 0.5, 0]]),]
    faces2 = [np.array([1, 0]), np.array([0, 1])]
    write_case('duplicates.case', ver, tri, orient, lines1, faces1, lines2, faces2)

def create_mesh2d():
    ver, tri = read_off('mesh2d.off')
    orient = np.fromfile('mesh2d.orient', sep=' ').reshape((-1, 3))
    lines1, faces1, info1 = evenly_spaced_streamlines(ver, tri, orient, 0.3, orthogonal=False)
    lines2, faces2, info2 = evenly_spaced_streamlines(ver, tri, orient, 0.3, orthogonal=True)
    write_case('mesh2d.case', ver, tri, orient, lines1, faces1, lines2, faces2)

def check_nonconnected():
    ver, tri = read_off('mesh2d.off')
    orient = np.random.random((tri.shape[0], 3)) - 0.5
    lines1, faces1, info1 = evenly_spaced_streamlines(ver, tri, orient, 0.3, orthogonal=False)
    lines2, faces2, info2 = evenly_spaced_streamlines(ver, tri, orient, 0.3, orthogonal=True)
    normal = tri_normals(ver, tri)
    out = create_cable_network(normal, lines1, faces1, lines2, faces2)
    data = ver, tri, orient, lines1, faces1, lines2, faces2
    return out.nb_connected_components, data

def create_nonconnected():
    while True:
        cnt, data = check_nonconnected()
        if cnt > 3:
            plot_case(*data)
            write_case('nonconnected3.case', *data)
            plt.show()

def plot_nonconnected():
    data = read_case('nonconnected3.case')
    normal = tri_normals(*data[:2])
    out = create_cable_network(normal, *data[-4:], remove_isolated_regions=False)
    plot_case(*data)
    X = out.vertices
    plt.plot(X[:, 0], X[:, 1], '.', c='C1')
    cables = unpack_cables(out.cables, out.cables_len)
    for c in cables:
        plt.plot(X[c, 0], X[c, 1], c='C1')
    
    out = create_cable_network(normal, *data[-4:], remove_isolated_regions=True)
    cables = unpack_cables(out.cables, out.cables_len)
    X = out.vertices
    plt.plot(X[:, 0], X[:, 1], '.', c='C3')
    for c in cables:
        plt.plot(X[c, 0], X[c, 1], c='C3')
    plt.show()

def write_case_bin(fname, *args):
    with open(fname, 'wb') as file:
        pickle.dump(args, file)

def read_case_bin(fname):
    with open(fname, 'rb') as file:
        return pickle.load(file)

def check_tessellation_failure():
    ver, tri = read_off('mesh2d.off')
    orient = np.random.random((tri.shape[0], 3)) - 0.5
    print('streamlines')
    lines1, faces1, info1 = evenly_spaced_streamlines(ver, tri, orient, 0.2, orthogonal=False)
    lines2, faces2, info2 = evenly_spaced_streamlines(ver, tri, orient, 0.2, orthogonal=True)
    normal = tri_normals(ver, tri)
    #write_case('segfault.case', ver, tri, orient, lines1, faces1, lines2, faces2, fmt='%.16e')
    #write_case_bin('infiniteloop.bincase', ver, tri, orient, lines1, faces1, lines2, faces2)
    print('cables')
    out = create_cable_network(normal, lines1, faces1, lines2, faces2)
    print('tessellate')
    #print(out.cables, out.cables_len, out.nc_long)
    assert np.all(out.cables_len > 1)
    facets, vneigh = tessellation.tessellate(out.cables, out.cables_len, out.nc_long, 
                                             out.sign, return_moves=True, return_neigh=True)
    _, cutoff = tessellation.find_largest_facets(facets, 1)
    print('triangulate')
    triangles, facetid, triangulation_failures = \
        triangulation.triangulate_facets(out.vertices, facets, cutoff)
    if len(triangulation_failures) > 0:
        fail1, fail2 = 0, 0
        for idx in triangulation_failures:
            if np.unique(idx).size == idx.size:
                fail1 += 1
                if 0:
                    X = out.vertices
                    x, y = X[idx, 0], X[idx, 1]
                    plot_case(ver, tri, orient, lines1, faces1, lines2, faces2)
                    plt.plot(x, y, 'o')
                    x = np.concatenate((x, [x[0]]))
                    y = np.concatenate((y, [y[0]]))
                    plt.plot(x, y)
                    cables = unpack_cables(out.cables, out.cables_len)
                    for c in cables[:out.nc_long]:
                        plt.plot(X[c, 0], X[c, 1], c='C3', lw=2)
                    for c in cables[out.nc_long:]:
                        plt.plot(X[c, 0], X[c, 1], c='C2', lw=2)
                    plt.show()
            else:
                fail2 += 1
        return len(facets), fail1, fail2
    return len(facets), 0, 0

def debug_triangulation():
    data = read_case_bin('infiniteloop1.bincase')
    normal = tri_normals(*data[:2])
    out = create_cable_network(normal, *data[-4:], remove_duplicates=True)
    assert np.all(out.cables_len > 1)
    facets, vneigh = tessellation.tessellate(out.cables, out.cables_len, out.nc_long, 
                                             out.sign, return_moves=True, return_neigh=True)
    _, cutoff = tessellation.find_largest_facets(facets, 1)
    print('triangulate')
    triangles, facetid, triangulation_failures = \
        triangulation.triangulate_facets(out.vertices, facets, cutoff)
    triangles.tofile('triangles', sep=' ')

    X = out.vertices
    c = [148, 52, 53, 147]
    plot_case(*data)
    #plt.plot(X[c, 0], X[c, 1], c='C3', lw=2)
    for i in c:
        plt.text(X[i, 0], X[i, 1], str(i))
    i = 146
    plt.text(X[i, 0], X[i, 1], str(i))
    cables = unpack_cables(out.cables, out.cables_len)
    for c in cables:
        plt.plot(X[c, 0], X[c, 1], '.-', lw=2)
    plt.show()

def check_triangulation():
    ver, tri = read_off('mesh2d.off')
    orient = np.random.random((tri.shape[0], 3)) - 0.5
    print('streamlines')
    lines1, faces1, info1 = evenly_spaced_streamlines(ver, tri, orient, 0.2, orthogonal=False)
    lines2, faces2, info2 = evenly_spaced_streamlines(ver, tri, orient, 0.2, orthogonal=True)
    normal = tri_normals(ver, tri)
    print('cables')
    out = create_cable_network(normal, lines1, faces1, lines2, faces2)
    print('tessellate')
    assert np.all(out.cables_len > 1)
    facets, vneigh = tessellation.tessellate(out.cables, out.cables_len, out.nc_long, 
                                             out.sign, return_moves=True, return_neigh=True)
    _, cutoff = tessellation.find_largest_facets(facets, 1)
    print('triangulate')
    triangles, facetid, triangulation_failures = \
        triangulation.triangulate_facets(out.vertices, facets, cutoff)
    
    topo = surfacetopology.surface_topology(triangles)
    X = out.vertices
    new_normals = tri_normals(X, triangles)
    
    if not np.all(new_normals[:, 2] > -1e-8):
        print('normal min(z) =', new_normals[:, 2].min())
        print('orientable =', all(t.orientable for t in topo))
        plot_case(ver, tri, orient, lines1, faces1, lines2, faces2)
        T = triangles[new_normals[:, 2] < 0]
        plt.tripcolor(X[:, 0], X[:, 1], np.ones(T.shape[0]), triangles=T)
        cables = unpack_cables(out.cables, out.cables_len)
        for c in cables:
            plt.plot(X[c, 0], X[c, 1], '.-', c='C3', lw=2)
        plt.show()

def check_intermediate_nodes():
    ver, tri = read_off('mesh2d.off')
    orient = np.random.random((tri.shape[0], 3)) - 0.5
    print('streamlines')
    lines1, faces1, info1 = evenly_spaced_streamlines(ver, tri, orient, 0.2, orthogonal=False)
    lines2, faces2, info2 = evenly_spaced_streamlines(ver, tri, orient, 0.2, orthogonal=True)
    normal = tri_normals(ver, tri)
    print('cables')
    out = create_cable_network(normal, lines1, faces1, lines2, faces2, add_ghost_nodes=True)
    print('tessellate')
    assert np.all(out.cables_len > 1)
    facets, vneigh = tessellation.tessellate(out.cables, out.cables_len, out.nc_long, 
                                             out.sign, return_moves=True, return_neigh=True)
    _, cutoff = tessellation.find_largest_facets(facets, 1)
    print('triangulate')
    triangles, facetid, triangulation_failures = \
        triangulation.triangulate_facets(out.vertices, facets, cutoff)
        
    X = out.vertices
    new_normals = tri_normals(X, triangles)
    T = triangles
    #plot_case(ver, tri, orient, lines1, faces1, lines2, faces2)
    perm = np.random.permutation(len(facets))
    plt.tripcolor(X[:, 0], X[:, 1], perm[facetid], triangles=T, cmap='hsv', alpha=0.5)
    cables = unpack_cables(out.cables, out.cables_len)
    for c in cables:
        plt.plot(X[c, 0], X[c, 1], '-', c='0.6', lw=2)
    I = out.is_node == 1
    print(f'fraction additional vertices = {I.size/I.sum()-1:.3f}')
    plt.plot(X[I, 0], X[I, 1], 'o', c='0.4', ms=3)
    plt.show()



if __name__ == "__main__":

    check_intermediate_nodes()

    # RANDOM TESTING
    
    # i = 1
    # while True:
    #     print(f'ITERATION {i}')
    #     check_triangulation()
    #     i += 1

    # i = 1
    # fail1, fail2, total = 0, 0, 0
    # while True:
    #     print(f'ITERATION {i}')
    #     t1, f1, f2 = check_tessellation_failure()
    #     total += t1
    #     fail1 += f1
    #     fail2 += f2
    #     i += 1
    #     if fail1+fail2==0: continue
    #     print(f'ratio = {fail1/(fail1+fail2):.2f}    {fail1}  {fail2}  {total}')
    
    # data = read_case_bin('segfault.bincase')
    # #print(min(np.linalg.norm(line[1:] - line[:-1]).min() for line in data[3] + data[5]))
    
    # normal = tri_normals(*data[:2])
    # out = create_cable_network(normal, *data[-4:], remove_duplicates=True)
    # print(out.cables, out.cables_len, out.nc_long)
    #plot_case(*data)
    #line = data[3][2]
    #plt.plot(line[:, 0], line[:, 1], 'r')
    #print(line)
    #X = out.vertices
    #plt.plot(X[:, 0], X[:, 1], '.', c='C3')
    #plt.show()

    # while True:
    #     facets, vneigh = tessellation.tessellate(out.cables, out.cables_len, out.nc_long, 
    #                                             out.sign, return_moves=True, return_neigh=True)
    # plot_case(*data)
    # X = out.vertices
    # plt.plot(X[:, 0], X[:, 1], '.', c='C3')
    # cables = unpack_cables(out.cables, out.cables_len)
    # for c in cables[:out.nc_long]:
    #     plt.plot(X[c, 0], X[c, 1], c='C3', lw=2)
    # for c in cables[out.nc_long:]:
    #     plt.plot(X[c, 0], X[c, 1], c='C2', lw=2)
    # plt.show()
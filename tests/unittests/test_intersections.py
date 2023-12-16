import numpy as np
from orthogonalstreamlines.intersection import create_cable_network, tri_normals
from testcases import read_case, plot_case


def test_empty():
    out = create_cable_network(np.ones((1, 3)), [], [], [], [])
    assert out.vertices.shape[0] == 0
    assert out.cables.size == 0
    assert out.cables_len.size == 0
    out = create_cable_network(np.ones((0, 3)), [], [], [], [])
    assert out.vertices.shape[0] == 0

def test_simple1tri_intersections():
    data = read_case('../data/simple1tri.case')
    normal = tri_normals(*data[:2])
    X = create_cable_network(25*normal, *data[-4:], cut_loose_ends=True).vertices
    assert X.shape[0] == 107
    X = create_cable_network(normal, *data[-4:], cut_loose_ends=False).vertices
    assert X.shape[0] == 107+5

def test_simple1tri_robustness():
    data = read_case('../data/simple1tri.case')
    normal = tri_normals(*data[:2])
    for line in data[3] + data[5]: # add noise
        line += (np.random.random(line.shape)-0.5) * 1e-4
    X = create_cable_network(normal, *data[-4:]).vertices
    assert X.shape[0] == 107

def test_simple1tri_partial_segments():
    data = read_case('../data/simple1tri.case')
    normal = tri_normals(*data[:2])
    for line in data[3]:
        line[1] = (line[0] + line[1])/2
    out = create_cable_network(normal, *data[-4:], cut_loose_ends=True, remove_isolated_regions=True)
    X = out.vertices
    assert X.shape[0] == 50

def test_duplicates():
    data = read_case('../data/duplicates.case')
    normal = tri_normals(*data[:2])
    out = create_cable_network(normal, *data[-4:], cut_loose_ends=False, remove_duplicates=False)
    assert out.vertices.shape[0] == 5
    out = create_cable_network(normal, *data[-4:], cut_loose_ends=False, remove_duplicates=True)
    assert out.vertices.shape[0] == 4
    assert out.nc_long == 2 and out.nc_trans == 2
    assert np.all(out[1] == 2)
    assert out.cables.size == 8

def test_mesh2d():
    data = read_case('../data/mesh2d.case')
    normal = tri_normals(*data[:2])
    out = create_cable_network(normal, *data[-4:])
    assert out.vertices.shape[0] == 432

def test_nonconnected():
    for k in 1, 2, 3:
        data = read_case(f'../data/nonconnected{k}.case')
        normal = tri_normals(*data[:2])
        out = create_cable_network(normal, *data[-4:])
        assert out.nb_connected_components == k+1
        assert out.cnt_isolated_vertices == [10, 9, 16][k-1]

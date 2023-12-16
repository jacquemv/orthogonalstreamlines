from itertools import chain
import numpy as np
from orthogonalstreamlines.intersection import create_cable_network, tri_normals
from orthogonalstreamlines.tessellation import tessellate
from testcases import read_case

def _check_euler_characteristic(out, facets):
    euler = out.vertices.shape[0] - np.sum(out.cables_len-1) + len(facets)
    assert euler == 2

def test_empty():
    facets = tessellate(np.empty(0, np.int32), np.empty(0, np.int32), 0, 
                        np.empty(0, np.uint8))
    assert len(facets) == 0

def test_simple1tri():
    data = read_case('../data/simple1tri.case')
    normal = tri_normals(*data[:2])
    out = create_cable_network(normal, *data[-4:])
    facets = tessellate(out.cables, out.cables_len, out.nc_long, out.sign)
    _check_euler_characteristic(out, facets)
    
def test_duplicates():
    data = read_case('../data/duplicates.case')
    normal = tri_normals(*data[:2])
    out = create_cable_network(normal, *data[-4:])
    facets = tessellate(out.cables, out.cables_len, out.nc_long, out.sign)
    _check_euler_characteristic(out, facets)

def test_mesh2d():
    data = read_case('../data/mesh2d.case')
    normal = tri_normals(*data[:2])
    out = create_cable_network(normal, *data[-4:])
    facets = tessellate(out.cables, out.cables_len, out.nc_long, out.sign)
    _check_euler_characteristic(out, facets)

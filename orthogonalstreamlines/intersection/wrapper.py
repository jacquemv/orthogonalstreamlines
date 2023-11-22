from collections import namedtuple
import numpy as np
from .runengine import find_intersections

__all__ = ['create_cable_network', 'unpack_cables', 'pack_cables', 
           'edge_lengths', 'tri_normals']


CableNetworkOutput = namedtuple('CableNetworkOutput', ['cables', 'cables_len',
                                'nc_long', 'nc_trans', 'vertices', 
                                'indices_tri', 'sign', 
                                'cnt_loose_ends', 'cnt_empty_cables', 
                                'cnt_duplicates'])

#-----------------------------------------------------------------------------
def create_cable_network(face_normals, lines1, faces1, lines2, faces2,
                         cut_loose_ends=True, remove_empty_cables=True,
                         remove_duplicates=True, epsilon=1e-8):
    """Create a cable network from two sets of orthogonal streamlines built
    on a triangulated surface mesh. Each vertex of the cable network is the 
    intersection between two orthogonal streamlines.

    Args:
        face_normals (nt-by-3 float64 array): vector normal to each of the nt 
            triangles of the surface computed as the cross-product of edge 
            vectors (can be computed using the function 'tri_normals')
        lines1 (list of n-by-3 int32 arrays): for each streamline in the 
            longitudinal direction, the array represents a list of 
            n consecutive vertices, each lying on an edge of the triangulated 
            mesh
        faces1 (list of (n-1)-int32 vectors): for each streamline in the 
            longitudinal direction, the array represents a list of (n-1) 
            triangle indices in which the edge of the streamline lies
        lines2 (list of n-by-3 int32 arrays): same as lines1, but for the 
            transverse direction
        faces2 (list of (n-1)-int32 vectors): same as faces1, but for the 
            transverse direction
        cut_loose_ends (bool): remove nodes that only have one neighbor 
            (default: True)
        remove_empty_cables (bool): remove cables of length 0 or 1
            (default: True)
        remove_duplicates (bool): handles the case where the intersection is 
            exactly on an edge of a triangle (default: True)
        epsilon (float): tolerance for 'remove_duplicates' (default: 1e-8)
    
    Returns:
        namedtuple with the following fields:
        - cables (int array): concatenation of arrays of vertex indices, each
          describing a cable (all longitudinal cables are listed first, then 
          the transverse cables)
        - cables_len (int array): length of each cable, such that 
          sum(cables_len) == cables.size
        - nc_long (int): number of longitudinal cables
        - nc_trans (int): number of transverse cables
        - ver (nv-by-3 array): vertex 3D positions
        - idtri (int vector of size nv): triangle index in which each vertex 
          lies
        - sign (uint8 vector of size nv): gives 1 of the cross product of the 
          tangent vectors of the intersecting streamlines is in the same
          direction as the vector normal to the surface, and 0 otherwise
        - cnt_loose_ends (int): number of nodes removed
        - cnt_empty_cables (int): number ofcables removed
        - cnt_duplicates (int): number of duplicate nodes removed
    """
    face_normals = np.ascontiguousarray(face_normals, dtype=np.float64)
    out = find_intersections(face_normals, lines1, faces1, lines2, faces2,
                             cut_loose_ends=cut_loose_ends, 
                             remove_empty_cables=remove_empty_cables,
                             remove_duplicates=remove_duplicates,
                             epsilon=epsilon)
    return CableNetworkOutput(
        cables=out[0], cables_len=out[1], 
        nc_long=out[2][0], nc_trans=out[2][0], 
        vertices=out[3], indices_tri=out[4], sign=out[5],
        cnt_loose_ends=out[6][0], cnt_empty_cables=out[6][1], 
        cnt_duplicates=out[6][2]
    )

#-----------------------------------------------------------------------------
def tri_normals(vertices, triangles):
    """Compute normal vectors for 'create_cable_network'

    Args:
        vertices (nv-by-3 float array): vertex position
        triangles (nt-by-3 int array): vertex indices of triangles
    
    Returns:
        nt-bt-3 float array: normal vectors (non-normalized)
    """
    P = vertices[triangles[:, 0], :]
    u = vertices[triangles[:, 1], :] - P
    v = vertices[triangles[:, 2], :] - P
    return np.cross(u, v)

#-----------------------------------------------------------------------------
def unpack_cables(cables, cables_len):
    """Unpack cables encoded as a single array into a list of cables

    Args:
        cables (int array): concatenated array of vertex indices of all cables
        cables_len (int array): number of vertices in each cable
    
    Returns:
        list of int array: list of array of vertex indices
    """
    if cables_len.size == 0:
        return []
    return np.split(cables, np.cumsum(cables_len)[:-1])

#-----------------------------------------------------------------------------
def pack_cables(cables):
    """Pack cables encoded as a list of cables into a single array

    Args:
        list of int array: list of array of vertex indices
    
    Returns:
        cables (int array): concatenated array of vertex indices of all cables
        cables_len (int array): number of vertices in each cable
    """
    cables_len = np.array([c.size for c in cables], dtype=np.int32)
    return np.concatenate(cables), cables_len

#-----------------------------------------------------------------------------
def edge_lengths(cables, cables_len, ver):
    """Compute the distances between neighboring vertices

    Args:
        cables, cables_len, ver: outputs of 'create_cable_network'
    
    Returns:
        float64 array: all edge lengths (dx); the array has size 
        (cables_len-1).sum()
    """

    dx = np.linalg.norm(ver[cables[1:]] - ver[cables[:-1]], axis=1)
    dx[np.cumsum(cables_len)[:-1]-1] = np.nan
    return dx[np.isfinite(dx)]
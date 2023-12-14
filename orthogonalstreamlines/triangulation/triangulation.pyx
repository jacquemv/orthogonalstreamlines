import numpy as np

__all__ = ['triangulate_facets']

#-----------------------------------------------------------------------------
cdef extern from "triangulatefacets.h":
    cdef cppclass TriangulateFacets:
        TriangulateFacets()
        void initialize(int nv_, double* ver_, int nt_max_, int* tri_)
        double max_dihedral_thres, max_dihedral_incr
        int nt
        int triangulate_facet(int n, int* facet)

#-----------------------------------------------------------------------------
SUBDIVISION_TABLE = { # most common cases
    'FLLFLL':  (np.array([0, 1, 2, 3]),     np.array([0, 3, 4, 5])),
    'LFLLFL':  (np.array([1, 2, 3, 4]),     np.array([1, 4, 5, 0])),
    'LLFLLF':  (np.array([2, 3, 4, 5]),     np.array([2, 5, 0, 1])),
    'FLFLFLL': (np.array([0, 1, 2, 3, 4]),  np.array([0, 4, 5, 6])),
    'LFLFLFL': (np.array([1, 2, 3, 4, 5]),  np.array([1, 5, 6, 0])),
    'LLFLFLF': (np.array([2, 3, 4, 5, 6]),  np.array([2, 6, 0, 1])),
    'FLLFLFL': (np.array([3, 4, 5, 6, 0]),  np.array([3, 0, 1, 2])),
    'LFLLFLF': (np.array([4, 5, 6, 0, 1]),  np.array([4, 1, 2, 3])),
    'FLFLLFL': (np.array([5, 6, 0, 1, 2]),  np.array([5, 2, 3, 4])),
    'LFLFLLF': (np.array([6, 0, 1, 2, 3]),  np.array([6, 3, 4, 5])),
}

#-----------------------------------------------------------------------------
def triangulate_facets(double[:, ::1] vertices, list facets, 
                       int cutoff=1000000):
    cdef:
        TriangulateFacets engine
        int i, nf, k, nt_max, j, i_facet
        int is_tuple
        int[:, ::1] triangles_ptr
        int[::1] facet_ptr
        int nbtri
        list failures = []

    nf = len(facets)
    if nf == 0:
        return np.empty((0, 3), dtype=np.int32)
    is_tuple = isinstance(facets[0], tuple)

    # calculate the number of triangles
    nt_max = 0
    for i in range(nf):
        if is_tuple:
            k = facets[i][0].size
        else:
            k = facets[i].size
        if k >= cutoff:
            continue
        nt_max += max(k-2, 0)
    
    # initialize structures
    triangles = np.empty((nt_max, 3), dtype=np.int32)
    facetid = np.empty(nt_max, dtype=np.int32)
    facetid = -np.ones(nt_max, dtype=np.int32)
    triangles_ptr = triangles
    engine.initialize(vertices.shape[0], &vertices[0, 0], 
                      nt_max, &triangles_ptr[0, 0])

    # triangulate
    i_facet = 0
    for i in range(nf):
        if is_tuple:
            facet_ptr = facets[i][0]
        else:
            facet_ptr = facets[i]
        k = facet_ptr.size
        if k >= cutoff:
            continue
        if is_tuple and (k == 6 or k == 7):
            moves = facets[i][1]
            if moves in SUBDIVISION_TABLE:
                loop1, loop2 = SUBDIVISION_TABLE[moves]
                facet_ptr = facets[i][0][loop1]
                nbtri = engine.triangulate_facet(loop1.size, &facet_ptr[0])
                facet_ptr = facets[i][0][loop2]
                nbtri += engine.triangulate_facet(loop2.size, &facet_ptr[0])
                for j in range(nbtri):
                    facetid[engine.nt-j-1] = i_facet
                i_facet += 1
                continue
        nbtri = engine.triangulate_facet(k, &facet_ptr[0])
        if nbtri < 0:
            failures.append(np.array(facet_ptr))
        else:
            for j in range(nbtri):
                facetid[engine.nt-j-1] = i_facet
        i_facet += 1

    return triangles[:engine.nt], facetid[:engine.nt], failures
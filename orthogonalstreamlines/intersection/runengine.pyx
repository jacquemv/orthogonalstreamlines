import numpy as np

#-----------------------------------------------------------------------------
cdef extern from "intersection.h":
    cdef cppclass Intersection:
        Intersection()
        int add_ghost_nodes
        void set_normals(int nt_, double* face_normals_)
        void insert_streamlines(int orientation, int nb_curves, 
            int *nb_segments, double* vertices, int* triangle_idx)
        void allocate()
        void identify_intersections()
        int cut_loose_cable_ends()
        int remove_zero_length_cables()
        int remove_duplicates(double epsilon)
        int ncomp, niter
        int remove_isolated_regions()
        int remove_isolated_vertices()

        int get_number_of_vertices()
        void get_vertices(double* vertices)
        void get_triangle_id(int* idtri)
        void get_vertex_sign(char* sign)
        void get_vertex_is_node(char* sign)
        int get_number_of_cables(int dir)
        void get_cables_delimiters(int* cables_delimiters)
        void get_cables_length(int* cables_len, int* sum_of_len)
        void get_cables(int* cables)

#-----------------------------------------------------------------------------
def empty_output():
    return (np.empty(0, dtype=np.int32), 
            np.empty(0, dtype=np.int32),
            (0, 0),
            np.empty((0, 3), dtype=np.float64),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.uint8),
            np.empty(0, dtype=np.uint8),
            (0, 0, 0, 0, 0, 0, 0))

#-----------------------------------------------------------------------------
def find_intersections(double[:, ::1] face_normals, 
                       list lines1, list faces1, list lines2, list faces2,
                       int add_ghost_nodes=False,
                       int cut_loose_ends=True, int remove_empty_cables=True, 
                       int remove_duplicates=True, double epsilon=1e-8,
                       int remove_isolated_regions=True):
    cdef:
        int nt, nv, i, nc, nc1, nc2, size
        Intersection engine
        int[::1] nseg1_memview
        double[:, ::1] lines1_memview
        int[::1] faces1_memview
        int[::1] nseg2_memview
        double[:, ::1] lines2_memview
        int[::1] faces2_memview
        double[:, ::1] ver_memview
        int[::1] idtri_memview
        char[::1] sign_memview
        char[::1] is_node_memview
        int[::1] cables_memview
        int[::1] cables_len_memview

    # check input arguments
    if len(lines1) != len(faces1):
        raise ValueError('arguments lines1 and faces1 must be '
                         'two lists of same length')
    if len(lines2) != len(faces2):
        raise ValueError('arguments lines2 and faces2 must be '
                         'two lists of same length')

    for i in range(len(lines1)):
        if lines1[i].size != 3*faces1[i].size+3:
            raise ValueError(f'argument lines1[{i}] must have exactly one '
                             f'more row than faces1[{i}]')
    
    for i in range(len(lines2)):
        if lines2[i].size != 3*faces2[i].size+3:
            raise ValueError(f'argument lines2[{i}] must have exactly one '
                             f'more row than faces2[{i}]')
    
    # create the object 
    nt = <int>face_normals.shape[0]
    if nt == 0:
        return empty_output()
    engine.set_normals(nt, &face_normals[0, 0])

    # longitudinal streamlines
    nc = <int>len(lines1)
    if nc == 0:
        return empty_output()
    nseg1 = np.array([face.size for face in faces1], dtype=np.int32)
    nseg1_memview = nseg1
    lines1_stacked = np.vstack(lines1)
    lines1_memview = lines1_stacked
    faces1_stacked = np.concatenate(faces1)
    faces1_memview = faces1_stacked
    engine.insert_streamlines(1, nc, &nseg1_memview[0], 
                              &lines1_memview[0, 0], &faces1_memview[0])

    # transverse streamlines
    nc = <int>len(lines2)
    if nc == 0:
        return empty_output()
    nseg2 = np.array([face.size for face in faces2], dtype=np.int32)
    nseg2_memview = nseg2
    lines2_stacked = np.vstack(lines2)
    lines2_memview = lines2_stacked
    faces2_stacked = np.concatenate(faces2)
    faces2_memview = faces2_stacked
    engine.insert_streamlines(2, nc, &nseg2_memview[0], 
                              &lines2_memview[0, 0], &faces2_memview[0])

    # run the code
    engine.add_ghost_nodes = add_ghost_nodes
    engine.allocate()
    engine.identify_intersections()
    cdef:
        int cnt_loose_ends=0, cnt_empty_cables=0, cnt_duplicates=0
        int cnt_isolated_vertices=0, cnt_removed_vertices=0, n_comp=1, n_iter=0
    if remove_duplicates:
        cnt_duplicates = engine.remove_duplicates(epsilon)
    if cut_loose_ends:
        cnt_loose_ends = engine.cut_loose_cable_ends()
        n_iter = engine.niter
    if remove_isolated_regions:
        cnt_isolated_vertices = engine.remove_isolated_regions()
        n_comp = engine.ncomp
    if remove_empty_cables:
        cnt_empty_cables = engine.remove_zero_length_cables()
    cnt_removed_vertices = engine.remove_isolated_vertices()

    # return the output
    nv = engine.get_number_of_vertices()
    if nv == 0:
        return empty_output()
    
    ver = np.empty((nv, 3), dtype=np.float64)
    ver_memview = ver
    engine.get_vertices(&ver_memview[0, 0])

    idtri = np.empty(nv, dtype=np.int32)
    idtri_memview = idtri
    engine.get_triangle_id(&idtri_memview[0])

    sign = np.empty(nv, dtype=np.uint8)
    sign_memview = sign
    engine.get_vertex_sign(&sign_memview[0])

    is_node = np.empty(nv, dtype=np.uint8)
    is_node_memview = is_node
    engine.get_vertex_is_node(&is_node_memview[0])

    nc1 = engine.get_number_of_cables(1)
    nc2 = engine.get_number_of_cables(2)
    nc = nc1 + nc2

    cables_len = np.empty(nc, dtype=np.int32)
    cables_len_memview = cables_len
    engine.get_cables_length(&cables_len_memview[0], &size)

    cables = np.empty(size, dtype=np.int32)
    cables_memview = cables
    engine.get_cables(&cables_memview[0])

    return (cables, cables_len, (nc1, nc2), ver, idtri, sign, is_node,
            (cnt_loose_ends, n_iter, cnt_empty_cables, 
             cnt_duplicates, cnt_isolated_vertices, n_comp,
             cnt_removed_vertices))
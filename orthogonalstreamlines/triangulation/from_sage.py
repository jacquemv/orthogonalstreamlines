import numpy as np
from mtlheartmesh.trisurf import TriSurf

nmax = 10
nb_triangles = np.zeros(nmax+1, dtype=int)
nb_dihedra = np.zeros(nmax+1, dtype=int)

with open('sage_triangulations.txt', 'rt') as file:
    n = 3
    matrix_list = [0] * (nmax+1)
    while True:
        cnt = 0
        matrix = []
        for line in file:
            if not line.strip():
                break
            s = line.translate(line.maketrans("", "", "()<>[]\n")).strip()
            matrix.append(np.fromstring(s, sep=','))
            cnt += 1
        matrix_list[n] = np.row_stack(matrix).astype(np.int32)

        if n == nmax:
            break
        n += 1
        if not cnt:
            break

    print("int empty_tri[] = {0, 0, 0, 0};")
    for n in range(3, nmax+1):
        t = np.arange(n) * 2*np.pi/n
        ver = np.column_stack((np.cos(t), np.sin(t)))
        tri_list = matrix_list[n]
        print("int all_triangulations_tri_"+str(n)+"[] = {")
        for tri in tri_list:
            S = TriSurf((ver, tri.reshape((-1, 3))))
            assert S.is_oriented
            assert np.all(S.face_normal[:, 2] > 0
            idx = S.tri.ravel().astype(int)
            s = repr(idx)[7:-2]+','
            print(s.replace(" ","").replace("\n", ""))
        print("};")

    line = "int* all_triangulations_tri["+str(nmax+1)+"] = {"
    line += "empty_tri, " * 3
    for n in range(3, nmax+1):
        line += "all_triangulations_tri_"+str(n)+", "
    print(line[:-2]+"};")

    print('')

    nb_triangles[3] = 1
    for n in range(4, nmax+1):
        t = np.arange(n) * 2*np.pi/n
        ver = np.column_stack((np.cos(t), np.sin(t)))
        tri_list = matrix_list[n]
        nb_triangles[n] = tri_list.shape[0]
        print("int all_triangulations_edges_"+str(n)+"[] = {")
        for tri in tri_list:
            S = TriSurf((ver, tri.reshape((-1, 3))))
            dihedra = np.hstack((S.edges, S.vertices_opposite_to_edge))
            dihedra = dihedra[dihedra[:, -1] != -1].ravel()
            s = repr(dihedra)[7:-2]+','
            print(s.replace(" ","").replace("\n", ""))
            nb_dihedra[n] = dihedra.size // 4
        print("};")

    line = "int* all_triangulations_edges["+str(nmax+1)+"] = {"
    line += "empty_tri, " * 4
    for n in range(4, nmax+1):
        line += "all_triangulations_edges_"+str(n)+", "
    print(line[:-2]+"};")

    print('')

    print('int nb_triangulations['+str(nmax+1)+'] = {',repr(nb_triangles)[7:-2],'};')
    print('int nb_dihedra['+str(nmax+1)+']   = {',repr(nb_dihedra)[7:-2],'};')

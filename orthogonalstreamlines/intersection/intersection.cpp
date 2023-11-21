#include <stdlib.h>
#include <stdio.h>
#include <cassert>

#include "intersection.h"
#include "algebra.cpp"
#include "streamlinecollection.cpp"

//-----------------------------------------------------------------------------
Intersection::Intersection()
{
    nt = -1;
    face_normals = NULL;
    dict_size = NULL;
    _ver1_ptr = _ver2_ptr = NULL;
    _tri1_ptr = _tri2_ptr = NULL;
    ver = NULL;
    ver_idtri = NULL;
    ver_sign = NULL;
    cables = NULL;
    cables_split = NULL;
    buffer = NULL;
}

//-----------------------------------------------------------------------------
void Intersection::set_normals(int nt_, double* face_normals_)
{
    nt = nt_;
    face_normals = face_normals_;
}

//-----------------------------------------------------------------------------
void Intersection::insert_streamlines(int orientation, int nb_curves, 
                    int *nb_segments, double* vertices, int* triangle_idx)
{
    int j = 0;
    double** ver_ptr = new double* [nb_curves];
    int** tri_ptr = new int* [nb_curves];
    for (int i=0;i<nb_curves;i++) {
        ver_ptr[i] = vertices + 3*(i+j);
        tri_ptr[i] = triangle_idx + j;
        j += nb_segments[i];
    }

    if (orientation == 1) {
        set1.initialize(nb_curves, nb_segments, ver_ptr, tri_ptr);
        _ver1_ptr = ver_ptr; // save pointers for freeing memory at the end
        _tri1_ptr = tri_ptr;
    } else if (orientation == 2) {
        set2.initialize(nb_curves, nb_segments, ver_ptr, tri_ptr);
        _ver2_ptr = ver_ptr;
        _tri2_ptr = tri_ptr;
    }
}

//-----------------------------------------------------------------------------
void Intersection::insert_streamlines(int orientation, int nb_curves, 
                    int *nb_segments, double** vertices, int** triangle_idx)
{
    if (orientation == 1) {
        set1.initialize(nb_curves, nb_segments, vertices, triangle_idx);
    } else if (orientation == 2) {
        set2.initialize(nb_curves, nb_segments, vertices, triangle_idx);
    }
}

//-----------------------------------------------------------------------------
void Intersection::allocate()
{
    // number of triangles
    if (nt < 0) {
        nt = set1.nt;
        if (set2.nt > nt) nt = set2.nt;
    }

    // upper bound for the number of vertices
    nv_max = 0;
    int min_nt = set1.nt < set2.nt ? set1.nt : set2.nt;
    for (int i=0;i<min_nt;i++)
        nv_max += set1.bin_size[i] * set2.bin_size[i];
    ver = new double [nv_max*3];
    ver_idtri = new int [nv_max];
    ver_sign = new char [nv_max];
    nv = 0;

    // allocate cables
    nc = set1.nb_curves + set2.nb_curves;
    cables = new int [nv_max*2];
    cables_split = new int [nc + 1];
    buffer = new int [nv_max*2];

    // create dictionary
    dict_size = new int [nt];
    dict_keys = new int* [nt];
    dict_values = new int* [nt];
    for (int i=0;i<nt;i++) {
        dict_size[i] = 0;
        int max_size = 1;
        if (i < min_nt)
            max_size = set1.bin_size[i] * set2.bin_size[i];
        dict_keys[i] = new int [4*max_size];
        dict_values[i] = new int [max_size];
    }
}

//-----------------------------------------------------------------------------
Intersection::~Intersection()
{
    delete [] ver;
    delete [] ver_idtri;
    delete [] ver_sign;
    delete [] cables;
    delete [] cables_split;
    delete [] buffer;
    if (dict_size) {
        delete [] dict_size;
        for (int i=0;i<nt;i++) {
            delete [] dict_keys[i];
            delete [] dict_values[i];
        }
        delete [] dict_keys;
        delete [] dict_values;
    }
    delete [] _ver1_ptr;
    delete [] _ver2_ptr;
    delete [] _tri1_ptr;
    delete [] _tri2_ptr;
}

//-----------------------------------------------------------------------------
int Intersection::get_number_of_vertices()
{
    return nv;
}

//-----------------------------------------------------------------------------
void Intersection::get_vertices(double* vertices)
{
    for (int i=0;i<3*nv;i++) vertices[i] = ver[i];
}

//-----------------------------------------------------------------------------
void Intersection::get_triangle_id(int* idtri)
{
    for (int i=0;i<nv;i++) idtri[i] = ver_idtri[i];
}

//-----------------------------------------------------------------------------
void Intersection::get_vertex_sign(char* sign)
{
    for (int i=0;i<nv;i++) sign[i] = ver_sign[i];
}

//-----------------------------------------------------------------------------
int Intersection::get_number_of_cables(int orientation)
{
    if (orientation == 1) return nc1;
    if (orientation == 2) return nc2;
    return nc;
}

//-----------------------------------------------------------------------------
void Intersection::get_cables_delimiters(int* cables_delimiters)
{
    for (int i=0;i<nc+1;i++) cables_delimiters[i] = cables_split[i];
}

//-----------------------------------------------------------------------------
void Intersection::get_cables_length(int* cable_len, int* sum_of_len)
{
    for (int i=0;i<nc;i++) cable_len[i] = cables_split[i+1] - cables_split[i];
    *sum_of_len = cables_split[nc];
}

//-----------------------------------------------------------------------------
void Intersection::get_cables(int* cables_)
{
    for (int i=0;i<cables_split[nc];i++) cables_[i] = cables[i];
}

//-----------------------------------------------------------------------------
void Intersection::dict_add(int idtri, int a, int b, int c, int d, int value)
{
    int k = dict_size[idtri];
    int* p = dict_keys[idtri] + 4*k;
    p[0] = a; p[1] = b; p[2] = c; p[3] = d;
    dict_values[idtri][k] = value;
    dict_size[idtri]++;
}

//-----------------------------------------------------------------------------
int Intersection::dict_find(int idtri, int a, int b, int c, int d)
{
    for (int k=0;k<dict_size[idtri];k++) {
        int* p = dict_keys[idtri] + 4*k;
        if (p[0] == a && p[1] == b && p[2] == c && p[3] == d)
            return dict_values[idtri][k];
    }
    return -1;
}

//-----------------------------------------------------------------------------
int segments_intersect(double* A, double* B, double* C, double* D,
                       double* X, double* n)
// returns the number of intersections (infinity -> 2)
// n = triangle normal vector
{
    double r[3] = {B[0]-A[0], B[1]-A[1], B[2]-A[2]};
    double s[3] = {D[0]-C[0], D[1]-C[1], D[2]-C[2]};
    double pq[3] = {C[0]-A[0], C[1]-A[1], C[2]-A[2]};
    double y[3], n2, u, v;
    vcross(n, r, s);
    n2 = n[0]*n[0] + n[1]*n[1] + n[2]*n[2];
    if (n2 == 0)
        return 0;
    vcross(y, pq, s);
    u = (y[0]*n[0] + y[1]*n[1] + y[2]*n[2])/n2;
    if ((u < 0) || (u > 1))
        return 0;
    vcross(y, pq, r);
    v = (y[0]*n[0] + y[1]*n[1] + y[2]*n[2])/n2;
    if ((v < 0) || (v > 1))
        return 0;
    X[0] = A[0]+u*r[0];
    X[1] = A[1]+u*r[1];
    X[2] = A[2]+u*r[2];
    return 1;
}

//-----------------------------------------------------------------------------
void Intersection::identify_intersections()
{
    if (!dict_size)
        allocate();
    int ne = 0;
    cables_split[0] = 0;
    nc = 0;

    // first pass: longitudinal cables
    for (int i=0;i<set1.nb_curves;i++) {
        for (int j=0;j<set1.nb_segments[i];j++) {
            int idtri = set1.get_triangle_id(i, j);
            double* seg1 = set1.get_segment_position(i, j);
            int n, *idcurv, *idseg;
            set2.get_segment_ordered_list(idtri, seg1, n, idcurv, idseg);
            for (int k=0;k<n;k++) {
                double* seg2 = set2.get_segment_position(idcurv[k], idseg[k]);
                double cross[3], x[3];
                int ni;
                ni = segments_intersect(seg1, seg1+3, seg2, seg2+3, x, cross);
                if (ni == 1) {
                    ver[3*nv] = x[0];
                    ver[3*nv+1] = x[1];
                    ver[3*nv+2] = x[2];
                    ver_idtri[nv] = idtri;
                    ver_sign[nv] = vdot(cross, face_normals+3*idtri) > 0;
                    dict_add(idtri, i, j, idcurv[k], idseg[k], ne);
                    cables[ne] = ne;
                    ne++;
                    nv++;
                }
            }
        }
        cables_split[++nc] = ne;
    }
    nc1 = nc;

    // second pass: transverse cables
    for (int i=0;i<set2.nb_curves;i++) {
        for (int j=0;j<set2.nb_segments[i];j++) {
            int idtri = set2.get_triangle_id(i, j);
            double* seg2 = set2.get_segment_position(i, j);
            int n, *idcurv, *idseg;
            set1.get_segment_ordered_list(idtri, seg2, n, idcurv, idseg);
            for (int k=0;k<n;k++) {
                int idv = dict_find(idtri, idcurv[k], idseg[k], i, j);
                if (idv >= 0)
                    cables[ne++] = idv;
            }
        }
        cables_split[++nc] = ne;
    }
    nc2 = nc-nc1;
}

//-----------------------------------------------------------------------------
void Intersection::count_neighbors(int* nb_neigh)
{
    for (int i=0;i<nv;i++) nb_neigh[i] = 0;

    for (int n=0;n<nc;n++) {
        for (int i=cables_split[n];i<cables_split[n+1]-1;i++) {
            if ((cables[i] == -1) || (cables[i+1] == -1)) continue;
            nb_neigh[cables[i]]++;
            nb_neigh[cables[i+1]]++;
        }
    }
}

//-----------------------------------------------------------------------------
void Intersection::print_neighbors_stat(int* nb_neigh)
{
    int hist[6] = {0, 0, 0, 0, 0, 0};
    for (int i=0;i<nv;i++) {
        if (nb_neigh[i] < 5)
            hist[nb_neigh[i]]++;
        else
            hist[5]++;
    }
    printf("#neigh: {0: %d, 1: %d, 2: %d, 3: %d, 4: %d, 5+: %d}\n", 
            hist[0], hist[1], hist[2], hist[3], hist[4], hist[5]);
}

//-----------------------------------------------------------------------------
void Intersection::print_cable_length_stat()
{
    int hist[6] = {0, 0, 0, 0, 0, 0};
    for (int n=0;n<nc;n++) {
        int size = cables_split[n+1] - cables_split[n];
        if (size < 5)
            hist[size]++;
        else
            hist[5]++;
    }
    printf("cable len: {0: %d, 1: %d, 2: %d, 3: %d, 4: %d, 5+: %d}\n", 
            hist[0], hist[1], hist[2], hist[3], hist[4], hist[5]);
}

//-----------------------------------------------------------------------------
int Intersection::remove_tagged_cable_nodes()
{
    int total_size = cables_split[nc];
    int cum_size = 0;
    int k1 = 0, k2;
    // update cable lengths
    for (int n=0;n<nc;n++) {
        k2 = cables_split[n+1];
        for (int i=k1;i<k2;i++)
            if (cables[i] != -1) cum_size++;
        k1 = k2;
        cables_split[n+1] = cum_size;
    }
    // update cable indices
    int j = 0, cnt = 0;
    for (int i=0;i<total_size;i++) {
        if (cables[i] == -1)
            cnt++;
        else
            cables[j++] = cables[i];
    }
    assert( j == cum_size );
    return cnt;
}

//-----------------------------------------------------------------------------
int Intersection::remove_isolated_vertices()
{
    int* nb_neigh = buffer;
    int* new_idx = buffer; // Beware! new_idx overwrites nb_neigh (but it's OK)
    count_neighbors(nb_neigh);

    // remove vertices
    int j = 0;
    for (int i=0;i<nv;i++) {
        if (nb_neigh[i] == 0) {
            new_idx[i] = 999999999;
        } else {
            ver[3*j] = ver[3*i];
            ver[3*j+1] = ver[3*i+1];
            ver[3*j+2] = ver[3*i+2];
            ver_idtri[j] = ver_idtri[i];
            ver_sign[j] = ver_sign[i];
            new_idx[i] = j;
            j++;
        }
    }
    int cnt = nv - j;
    nv = j;

    // renumber vertices
    for (int i=0;i<cables_split[nc];i++) {
        cables[i] = new_idx[cables[i]];
    }
    return cnt;
}

//-----------------------------------------------------------------------------
int Intersection::cut_loose_cable_ends()
{
    int cnt = -1, total_cnt = 0;
    int* nb_neigh = buffer;

    while (cnt != 0) {
        count_neighbors(nb_neigh);

        // tag loose cable ends
        cnt = 0;
        for (int n=0;n<nc;n++) {
            int k1 = cables_split[n];
            int k2 = cables_split[n+1]-1;
            int size = k2-k1+1;
            if (size < 1) continue;
            // remove cables of size one
            if (size == 1) {
                cables[k1] = -1;
                cnt++;
                continue;
            }
            // cut the beginning of the cable
            if (nb_neigh[cables[k1]] == 1) {
                cables[k1] = -1;
                cnt++;
                k1++;
                while (k1 < k2 && nb_neigh[cables[k1]] == 2) {
                    cables[k1] = -1;
                    cnt++;
                    k1++;
                }
            }
            // cut the end of the cable
            k1 = cables_split[n];
            if (nb_neigh[cables[k2]] == 1) {
                cables[k2] = -1;
                cnt++;
                k2--;
                while (k1 < k2 && nb_neigh[cables[k2]] == 2) {
                    if (cables[k2] == -1) break;
                    cables[k2] = -1;
                    cnt++;
                    k2--;
                }
            }
        }
        remove_tagged_cable_nodes();
        total_cnt += cnt;
    }
    remove_isolated_vertices();
    return total_cnt;
}

//-----------------------------------------------------------------------------
int Intersection::remove_zero_length_cables()
{
    int m = 0, cnt1 = 0, cnt2 = 0;
    for (int n=0;n<nc;n++) {
        int size = cables_split[n+1] - cables_split[n];
        if (size == 0) {
            if (n < nc1) cnt1++; else cnt2++;
        } else {
            cables_split[m++] = cables_split[n];
        }
    }
    cables_split[m] = cables_split[nc];
    nc -= cnt1 + cnt2;
    nc1 -= cnt1;
    nc2 -= cnt2;
    return cnt1 + cnt2;
}

//-----------------------------------------------------------------------------
 int Intersection::remove_duplicates(double epsilon)
 {
    int cnt = 0;
    double eps2 = epsilon*epsilon;

    for (int n=0;n<nc;n++) {
        for (int i = cables_split[n]; i < cables_split[n+1]-1; i++) {
            double dx[3];
            vdiff(dx, ver+3*cables[i], ver+3*cables[i+1]);
            if (vnorm2(dx) < eps2) {
                if (cables[i] < cables[i+1]) // remove the largest index
                    cables[i+1] = -1;
                else
                    cables[i] = -1;
                cnt++;
            }
        }
    }
    if (cnt) {
        remove_tagged_cable_nodes();
        remove_isolated_vertices();
    }
    return cnt;
 }
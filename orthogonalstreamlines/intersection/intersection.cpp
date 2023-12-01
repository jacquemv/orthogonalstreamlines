#include <stdlib.h>
#include <stdio.h>
#include <cassert>

#include "intersection.h"
#include "algebra.cpp"
#include "streamlinecollection.cpp"
#include "fixedsizedict.cpp"
#include "vertices.cpp"
#include "cablenetwork.cpp"

#define MINIMUM_DETERMINANT 1e-12

//-----------------------------------------------------------------------------
Intersection::Intersection()
{
    nt = -1;
    add_ghost_nodes = 0;
    face_normals = nullptr;
}

//-----------------------------------------------------------------------------
Intersection::~Intersection()
{
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
    int min_nt = (set1.nt < set2.nt) ? set1.nt : set2.nt;

    // upper bound for the number of cables and vertices
    int nv_max = 0;
    for (int i=0;i<min_nt;i++)
        nv_max += set1.bin_size[i] * set2.bin_size[i];
    if (add_ghost_nodes) {
        for (int i=0;i<set1.nb_curves;i++)
            nv_max += set1.nb_segments[i]+1;
        for (int i=0;i<set2.nb_curves;i++)
            nv_max += set2.nb_segments[i]+1;
    }
    int nc_max = set1.nb_curves + set2.nb_curves;
    
    ver.allocate(nv_max);
    cnet.allocate(nc_max, nv_max*2);

    // create dictionary
    dict.allocate_bins(nt);
    for (int i=0;i<nt;i++) {
        int cap = i < min_nt ? set1.bin_size[i] * set2.bin_size[i] : 1;
        dict.set_bin_capacity(i, cap);
    }
}

//-----------------------------------------------------------------------------
int segments_intersect(double* A, double* B, double* C, double* D,
                       double* normal, double* X, double* cross)
// returns the number of intersections (0 or 1)
// X = intersection point
// cross = (B-A) x (D-C) normal to the triangle
{
    double r[3], s[3], t[3];
    vdiff(r, A, B); vdiff(s, C, D); vdiff(t, A, C);
    vrmcomp(r, normal); vrmcomp(s, normal); 
    vcross(cross, r, s);
    double r2 = vnorm2(r), s2 = vnorm2(s);
    double rs = vdot(r, s), rt = vdot(r, t), st = vdot(s, t);
    double u, v;
    double delta = solve2x2(r2, -rs, -rs, s2, rt, -st, u, v);
    if (delta < MINIMUM_DETERMINANT)
        return 0;
    if ((u < 0) || (u > 1) || (v < 0) || (v > 1))
        return 0;
    vcopyadd(X, A, u, r);
    return 1;
}

//-----------------------------------------------------------------------------
void Intersection::identify_intersections()
{
    cnet.clear();

    // first pass: longitudinal cables
    for (int i=0;i<set1.nb_curves;i++) {
        cnet.new_cable();
        for (int j=0;j<set1.nb_segments[i];j++) {
            int idtri = set1.get_triangle_id(i, j);
            double* seg1 = set1.get_segment_position(i, j);
            if (add_ghost_nodes) {
                if (j > 0) {
                    cnet.append(ver.size);
                    ver.append_ghost(seg1, idtri, 1);
                }
            }
            int n, *idcurv, *idseg;
            set2.get_segment_ordered_list(idtri, seg1, n, idcurv, idseg);
            double* normal = face_normals+3*idtri;
            for (int k=0;k<n;k++) {
                double* seg2 = set2.get_segment_position(idcurv[k], idseg[k]);
                double cross[3], x[3];
                int ni = segments_intersect(seg1, seg1+3, seg2, seg2+3, normal, x, cross);
                if (ni == 1) {
                    cnet.append(ver.size);
                    ver.append_node(x, idtri, vdot(cross, normal) > 0);
                    dict.add(idtri, i, j, idcurv[k], idseg[k], ver.size-1);
                }
            }
        }
    }
    cnet.new_group();

    // second pass: transverse cables
    for (int i=0;i<set2.nb_curves;i++) {
        cnet.new_cable();
        for (int j=0;j<set2.nb_segments[i];j++) {
            int idtri = set2.get_triangle_id(i, j);
            double* seg2 = set2.get_segment_position(i, j);
            if (add_ghost_nodes) {
                if (j > 0) {
                    cnet.append(ver.size);
                    ver.append_ghost(seg2, idtri, 1);
                }
            }
            int n, *idcurv, *idseg;
            set1.get_segment_ordered_list(idtri, seg2, n, idcurv, idseg);
            for (int k=0;k<n;k++) {
                int idv = dict.find(idtri, idcurv[k], idseg[k], i, j);
                if (idv >= 0)
                    cnet.append(idv);
            }
        }
    }
}

//-----------------------------------------------------------------------------
int Intersection::remove_isolated_vertices()
// isolated = in none of the cables or in one/two cable(s) with a single node
// there may still be cables of length 0 or 1 at the end
{
    int* nb_neigh = new int[ver.size];
    cnet.count_neighbors(ver.size, nb_neigh);
    int cnt = ver.remove_if_zero(nb_neigh);
    int invalid = ver.renumber_indices(cnet.size(), cnet.idx);
    if (invalid) printf("Error: %d invalid node indices", invalid);
    delete [] nb_neigh;
    return cnt;
}

//-----------------------------------------------------------------------------
int Intersection::cut_loose_cable_ends()
{
    return cnet.cut_loose_ends(niter);
}

//-----------------------------------------------------------------------------
int Intersection::remove_zero_length_cables()
{
    return cnet.squeeze();
}

//-----------------------------------------------------------------------------
 int Intersection::remove_duplicates(double epsilon)
 {
    int cnt = 0;
    double eps2 = epsilon*epsilon;

    for (int n=0;n<cnet.nc;n++) {
        for (int i = cnet.sep[n]; i < cnet.sep[n+1]-1; i++) {
            double dx[3];
            vdiff(dx, ver.pos+3*cnet.idx[i], ver.pos+3*cnet.idx[i+1]);
            if (vnorm2(dx) <= eps2) {
                if (cnet.idx[i] < cnet.idx[i+1]) // remove the largest index
                    cnet.idx[i+1] = -1;
                else
                    cnet.idx[i] = -1;
                cnt++;
            }
        }
    }
    if (cnt)
        cnet.remove_negative_indices();
    return cnt;
 }

//-----------------------------------------------------------------------------
int Intersection::remove_isolated_regions()
{
    return cnet.remove_isolated_regions(ncomp);
}

//-----------------------------------------------------------------------------
bool Intersection::check_cable_indices()
{
    for (int i=0;i<cnet.size();i++)
        if ((cnet.idx[i] < 0) || (cnet.idx[i] >= ver.size)) return 0;
    return 1;
}

//-----------------------------------------------------------------------------
int Intersection::get_number_of_vertices()
{
    return ver.size;
}

//-----------------------------------------------------------------------------
void Intersection::get_vertices(double* vertices)
{
    for (int i=0;i<3*ver.size;i++) vertices[i] = ver.pos[i];
}

//-----------------------------------------------------------------------------
void Intersection::get_triangle_id(int* idtri)
{
    for (int i=0;i<ver.size;i++) idtri[i] = ver.idtri[i];
}

//-----------------------------------------------------------------------------
void Intersection::get_vertex_sign(char* sign)
{
    for (int i=0;i<ver.size;i++) sign[i] = ver.sign[i];
}

//-----------------------------------------------------------------------------
void Intersection::get_vertex_is_node(char* is_node)
{
    for (int i=0;i<ver.size;i++) is_node[i] = ver.is_node[i];
}

//-----------------------------------------------------------------------------
int Intersection::get_number_of_cables(int orientation)
{
    if (orientation == 1) return cnet.group_sep;
    if (orientation == 2) return cnet.nc - cnet.group_sep;
    return cnet.nc;
}

//-----------------------------------------------------------------------------
void Intersection::get_cables_delimiters(int* cables_delimiters)
{
    for (int i=0;i<=cnet.nc;i++) cables_delimiters[i] = cnet.sep[i];
}

//-----------------------------------------------------------------------------
void Intersection::get_cables_length(int* cable_len, int* sum_of_len)
{
    for (int i=0;i<cnet.nc;i++) cable_len[i] = cnet.sep[i+1] - cnet.sep[i];
    *sum_of_len = cnet.size();
}

//-----------------------------------------------------------------------------
void Intersection::get_cables(int* cables)
{
    for (int i=0;i<cnet.size();i++) cables[i] = cnet.idx[i];
}
#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include "algebra.h"
#include "streamlinecollection.h"

class Intersection {
public:
    Intersection();
    ~Intersection();

    void set_normals(int nt_, double* face_normals_);
    // convenient for C
    void insert_streamlines(int orientation, int nb_curves, 
                            int *nb_segments, double** vertices, int** triangle_idx);
    // convenient for python
    void insert_streamlines(int orientation, int nb_curves, 
                            int *nb_segments, double* vertices, int* triangle_idx);

    void identify_intersections();
    int cut_loose_cable_ends();
    int remove_zero_length_cables();

    // export output data
    int get_number_of_vertices();
    void get_vertices(double* vertices);
    void get_triangle_id(int* idtri);
    void get_vertex_sign(char* sign);
    int get_number_of_cables(int orientation);
    void get_cables_length(int* cable_len, int* sum_of_len);
    void get_cables_delimiters(int* cables_delimiters);
    void get_cables(int* cables);

//private:
    StreamlineCollection set1, set2;
    double* ver;
    int* ver_idtri;
    char* ver_sign;
    double* face_normals;
    int nv, nv_max, nt;
    int nc, nc1, nc2;
    int* cables;
    int* cables_split;
    int** dict_keys;
    int** dict_values;
    int* dict_size;
    int* buffer;
    double **_ver1_ptr, **_ver2_ptr;
    int **_tri1_ptr, **_tri2_ptr;

    void allocate(); // called by identify_intersections

    void dict_add(int idtri, int a, int b, int c, int d, int value);
    int dict_find(int idtri, int a, int b, int c, int d);

    void count_neighbors(int* nb_neigh);
    int remove_tagged_cable_nodes();
    int remove_isolated_vertices();

    void print_neighbors_stat(int* nb_neigh);
    void print_cable_length_stat();
};

#endif
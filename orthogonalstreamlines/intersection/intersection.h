#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include "algebra.h"
#include "streamlinecollection.h"
#include "cablenetwork.h"
#include "fixedsizedict.h"
#include "vertices.h"

class Intersection {
public:
    Intersection();
    ~Intersection();
    int add_ghost_nodes;

    void set_normals(int nt_, double* face_normals_);
    void insert_streamlines(int orientation, int nb_curves, 
                            int *nb_segments, double* vertices, int* triangle_idx);

    void identify_intersections();
    int cut_loose_cable_ends();
    int remove_zero_length_cables();
    int remove_duplicates(double epsilon);
    int remove_isolated_regions();

    // export output data
    int get_number_of_vertices();
    void get_vertices(double* vertices);
    void get_triangle_id(int* idtri);
    void get_vertex_sign(char* sign);
    void get_vertex_is_node(char* is_node);
    int get_number_of_cables(int orientation);
    void get_cables_length(int* cable_len, int* sum_of_len);
    void get_cables_delimiters(int* cables_delimiters);
    void get_cables(int* cables);

//private:
    StreamlineCollection set1, set2;
    Vertices ver;
    double* face_normals;
    int nt;
    CableNetwork cnet;
    FixedSizeDict dict;
    int ncomp, niter;

    void allocate(); // called by identify_intersections
    int remove_isolated_vertices();

    bool check_cable_indices();
};

#endif
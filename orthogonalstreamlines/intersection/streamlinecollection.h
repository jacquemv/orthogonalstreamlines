#ifndef STREAMLINECOLLECTION_H_
#define STREAMLINECOLLECTION_H_

class StreamlineCollection {
public:
    int nb_curves;
    int* nb_segments;
    int nt;

    StreamlineCollection();
    StreamlineCollection(int nb_curves_, int *nb_segments_, double** vertices_, int** triangle_idx_);
    ~StreamlineCollection();

    void initialize(int nb_curves_, int *nb_segments_, double** vertices_, int** triangle_idx_);
    void initialize(int nb_curves_, int *nb_segments_, double* vertices_, int* triangle_id_);

    inline double* get_segment_position(int idcurv, int idseg);
    inline int get_triangle_id(int idcurv, int idseg);
    void get_segment_list(int idtri, int &n, int* &idcurv, int* &idseg);
    void get_segment_ordered_list(int idtri, double* direction, int &n, int* &idcurv, int* &idseg);

//private:
    int** triangle_id;
    double** ver;
    int* bin_size;
    int** bin_idcurv;
    int** bin_idseg;
    int max_bin_size;
    double* buffer;
    double** ver_ptr;
    int** tri_ptr;
};

#endif
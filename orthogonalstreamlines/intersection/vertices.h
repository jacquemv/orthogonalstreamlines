#ifndef VERTICES_H_
#define VERTICES_H_

#define INVALID_INDEX 999999999

class Vertices {
public:
    int size;
    double* pos;
    int* idtri;
    char* sign;
    char* is_node;

    Vertices();
    ~Vertices();
    void allocate(int capacity_);

    void append_node(double* pos_, int idtri_, char sign_);
    void append_ghost(double* pos_, int idtri_, char sign_);
    int remove_if_zero(int* flag); // flag must be of length at least 'size'
    int renumber_indices(int n, int* indices); // indices has size n

//private:
    int capacity;
    int* new_idx;
    int size_new_idx;
};

#endif
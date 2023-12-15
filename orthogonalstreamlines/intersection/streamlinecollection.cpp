#include "streamlinecollection.h"
#include "algebra.h"
#include <cfloat>

//-----------------------------------------------------------------------------
StreamlineCollection::StreamlineCollection(int nb_curves_, int *nb_segments_, 
                                    double** vertices_, int** triangle_id_)
{
    initialize(nb_curves_, nb_segments_, vertices_, triangle_id_);
}

//-----------------------------------------------------------------------------
StreamlineCollection::StreamlineCollection()
{
    nt = -1;
    bin_size = nullptr;
    bin_idcurv = nullptr;
    bin_idseg = nullptr;
    buffer = nullptr;
    ver_ptr = nullptr;
    tri_ptr = nullptr;
}

//-----------------------------------------------------------------------------
StreamlineCollection::~StreamlineCollection()
{
    delete [] bin_size;
    for (int i=0;i<nt;i++) {
        delete [] bin_idcurv[i];
        delete [] bin_idseg[i];
    }
    delete [] bin_idcurv;
    delete [] bin_idseg;
    delete [] buffer;
    delete [] ver_ptr;
    delete [] tri_ptr;
}

//-----------------------------------------------------------------------------
void StreamlineCollection::initialize(int nb_curves_, int *nb_segments_, 
                                    double** vertices_, int** triangle_id_)
{
    nb_curves = nb_curves_;
    nb_segments = nb_segments_;
    triangle_id = triangle_id_;
    ver = vertices_;

    // count triangles/bins
    nt = -1;
    for (int i=0;i<nb_curves;i++) for (int j=0;j<nb_segments[i];j++) {
        if (triangle_id[i][j] > nt)
            nt = triangle_id[i][j];
    }
    nt++;

    // determine max bin size
    bin_size = new int [nt];
    for (int i=0;i<nt;i++)
        bin_size[i] = 0;
    for (int i=0;i<nb_curves;i++) for (int j=0;j<nb_segments[i];j++)
        bin_size[triangle_id[i][j]]++;
    
    // create bins
    bin_idcurv = new int* [nt];
    bin_idseg = new int* [nt];
    for (int i=0;i<nt;i++) {
        bin_idcurv[i] = new int [bin_size[i]];
        bin_idseg[i] = new int [bin_size[i]];
    }

    // reset bin size
    for (int i=0;i<nt;i++)
        bin_size[i] = 0;

    // fill bins
    for (int i=0;i<nb_curves;i++) for (int j=0;j<nb_segments[i];j++) {
        int k = triangle_id[i][j];
        int n = bin_size[k];
        bin_idcurv[k][n] = i;
        bin_idseg[k][n] = j;
        bin_size[k]++;
    }

    // compute max bin size
    max_bin_size = 0;
    for (int i=0;i<nt;i++)
        if (bin_size[i] > max_bin_size)
            max_bin_size = bin_size[i];
    buffer = new double [max_bin_size];
}

//-----------------------------------------------------------------------------
void StreamlineCollection::initialize(int nb_curves_, int *nb_segments_, 
                                    double* vertices_, int* triangle_id_)
{
    int j = 0;
    ver_ptr = new double* [nb_curves_];
    tri_ptr = new int* [nb_curves_];
    for (int i=0;i<nb_curves_;i++) {
        ver_ptr[i] = vertices_ + 3*(i+j);
        tri_ptr[i] = triangle_id_ + j;
        j += nb_segments_[i];
    }
    initialize(nb_curves_, nb_segments_, ver_ptr, tri_ptr);
}

//-----------------------------------------------------------------------------
inline double* StreamlineCollection::get_segment_position(int idline, 
                                                          int idseg)
{
    return ver[idline] + 3*idseg;
}

//-----------------------------------------------------------------------------
inline int StreamlineCollection::get_triangle_id(int idline, int idseg)
{
    return triangle_id[idline][idseg];
}

//-----------------------------------------------------------------------------
void StreamlineCollection::get_segment_list(int idtri, int &n, int* &idcurv, 
                                            int* &idseg)
{
    if (idtri >= nt) {
        n = 0;
        return;
    }
    n = bin_size[idtri];
    idcurv = bin_idcurv[idtri];
    idseg = bin_idseg[idtri];
}

//-----------------------------------------------------------------------------
double vector_coordinate(double *A, double *B, double *C, double *D)
// compute the coordinate along (AB) of the intersection between the lines
// (AB) and (CD)
{
    double u[3], v[3], w[3];
    vdiff(u, A, B); vdiff(v, C, D); vdiff(w, A, C);
    double u2 = vnorm2(u), v2 = vnorm2(v), uv = vdot(u, v);
    double delta = u2*v2 - uv*uv;
    double coord = (vdot(u, w)*v2 - vdot(w, v)*uv) / delta;
    if (isfinite(coord)) return coord;
    return DBL_MAX;
}

//-----------------------------------------------------------------------------
void insertion_sort_multiple(double arr[], int arr2[], int arr3[], int n)
{
    int i, j, key2, key3;
    double key;
    for (i = 1; i < n; i++) {
        key = arr[i];
        key2 = arr2[i];
        key3 = arr3[i];
        j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            arr2[j + 1] = arr2[j];
            arr3[j + 1] = arr3[j];
            j = j - 1;
        }
        arr[j + 1] = key;
        arr2[j + 1] = key2;
        arr3[j + 1] = key3;
    }
}

//-----------------------------------------------------------------------------
void StreamlineCollection::get_segment_ordered_list(int idtri, 
                        double* direction, int &n, int* &idcurv, int* &idseg)
{
    get_segment_list(idtri, n, idcurv, idseg);
    for (int i=0;i<n;i++) {
        double* x = get_segment_position(idcurv[i], idseg[i]);
        buffer[i] = vector_coordinate(direction, direction+3, x, x+3);
    }
    insertion_sort_multiple(buffer, idcurv, idseg, n);
    while ((buffer[n-1] == DBL_MAX) && (n > 0)) n--;
}
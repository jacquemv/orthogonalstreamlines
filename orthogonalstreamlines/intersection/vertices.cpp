#include <cassert>
#include "vertices.h"

//-----------------------------------------------------------------------------
Vertices::Vertices()
{
    size = capacity = 0;
    pos = nullptr;
    idtri = nullptr;
    sign = nullptr;
    is_node = nullptr;
    new_idx = nullptr;
}

//-----------------------------------------------------------------------------
Vertices::~Vertices()
{
    delete [] pos;
    delete [] idtri;
    delete [] sign;
    delete [] is_node;
    delete [] new_idx;
}

//-----------------------------------------------------------------------------
void Vertices::allocate(int capacity_)
{
    capacity = capacity_;
    pos = new double [3*capacity];
    idtri = new int [capacity];
    sign = new char [capacity];
    is_node = new char [capacity];
    new_idx = new int [capacity];
    for (int i=0;i<capacity;i++) {
        new_idx[i] = INVALID_INDEX;
        is_node[i] = 1;
    }
    size = 0;
    size_new_idx = 0;
}

//-----------------------------------------------------------------------------
void Vertices::append_node(double* pos_, int idtri_, char sign_)
{
    assert(size < capacity);
    pos[3*size] = pos_[0];
    pos[3*size+1] = pos_[1];
    pos[3*size+2] = pos_[2];
    idtri[size] = idtri_;
    sign[size] = sign_;
    size++;
}

//-----------------------------------------------------------------------------
void Vertices::append_ghost(double* pos_, int idtri_, char sign_)
{
    assert(size < capacity);
    pos[3*size] = pos_[0];
    pos[3*size+1] = pos_[1];
    pos[3*size+2] = pos_[2];
    idtri[size] = idtri_;
    sign[size] = sign_;
    is_node[size] = 0;
    size++;
}

//-----------------------------------------------------------------------------
int Vertices::remove_if_zero(int* flag)
{
    int j = 0;
    size_new_idx = size;
    for (int i=0;i<size;i++) {
        if (flag[i] == 0) {
            new_idx[i] = INVALID_INDEX;
        } else {
            pos[3*j] = pos[3*i];
            pos[3*j+1] = pos[3*i+1];
            pos[3*j+2] = pos[3*i+2];
            idtri[j] = idtri[i];
            sign[j] = sign[i];
            is_node[j] = is_node[i];
            new_idx[i] = j;
            j++;
        }
    }
    size = j;
    return size_new_idx - size;
}

//-----------------------------------------------------------------------------
int Vertices::renumber_indices(int n, int* indices)
{
    int cnt_invalid = 0;
    for (int i=0;i<n;i++) {
        assert((indices[i] >= 0) && (indices[i] < size_new_idx));
        int idx = new_idx[indices[i]];
        cnt_invalid += (idx == INVALID_INDEX);
        indices[i] = idx;
    }
    return cnt_invalid;
}
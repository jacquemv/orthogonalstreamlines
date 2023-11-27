#include <cassert>
#include "fixedsizedict.h"

//-----------------------------------------------------------------------------
FixedSizeDict::FixedSizeDict()
{
    nbins = -1;
}

//-----------------------------------------------------------------------------
FixedSizeDict::~FixedSizeDict()
{
    if (nbins < 0) return;
    delete [] capacity;
    delete [] size;
    for (int i=0;i<nbins;i++) {
        delete [] keys[i];
        delete [] values[i];
    }
    delete [] keys;
    delete [] values;
}

//-----------------------------------------------------------------------------
void FixedSizeDict::allocate_bins(int nbins_)
{
    nbins= nbins_;
    size = new int [nbins];
    capacity = new int [nbins];
    keys = new int* [nbins];
    values = new int* [nbins];
    for (int i=0;i<nbins;i++) {
        size[i] = capacity[i] = 0;
        keys[i] = values[i] = nullptr;
    }
}

//-----------------------------------------------------------------------------
void FixedSizeDict::set_bin_capacity(int bin_id, int capacity_)
{
    assert((bin_id >= 0) && (bin_id < nbins));
    capacity[bin_id] = capacity_;
    keys[bin_id] = new int [4 * capacity_];
    values[bin_id] = new int [capacity_];
}

//-----------------------------------------------------------------------------
void FixedSizeDict::add(int bin_id, int a, int b, int c, int d, int value)
{
    int k = size[bin_id];
    int* p = keys[bin_id] + 4*k;
    p[0] = a; p[1] = b; p[2] = c; p[3] = d;
    values[bin_id][k] = value;
    size[bin_id]++;
}

//-----------------------------------------------------------------------------
int FixedSizeDict::find(int bin_id, int a, int b, int c, int d)
{
    for (int k=0;k<size[bin_id];k++) {
        int* p = keys[bin_id] + 4*k;
        if (p[0] == a && p[1] == b && p[2] == c && p[3] == d)
            return values[bin_id][k];
    }
    return -1;
}
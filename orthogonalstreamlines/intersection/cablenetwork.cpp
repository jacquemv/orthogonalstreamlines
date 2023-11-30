#include <cassert>
#include "cablenetwork.h"
#include "unionfind.cpp"

//-----------------------------------------------------------------------------
CableNetwork::CableNetwork()
{
    nc = 0;
    idx = nullptr;
    sep = nullptr;
}

//-----------------------------------------------------------------------------
CableNetwork::~CableNetwork()
{
    delete [] idx;
    delete [] sep;
}

//-----------------------------------------------------------------------------
void CableNetwork::allocate(int capacity_cables, int capacity_indices)
{
    nc_max = capacity_cables;
    nv_max = capacity_indices;
    idx = new int [nv_max];
    sep = new int [nc_max+1];
}

//-----------------------------------------------------------------------------
void CableNetwork::clear()
{
    nc = 0;
    sep[0] = 0;
    group_sep = 0;
}

//-----------------------------------------------------------------------------
void CableNetwork::append(int index)
{
    idx[sep[nc]++] = index;
}

//-----------------------------------------------------------------------------
void CableNetwork::new_cable()
{
    nc++;
    assert(nc < nc_max);
    sep[nc] = sep[nc-1];
}

//-----------------------------------------------------------------------------
void CableNetwork::new_group()
{
    group_sep = nc;
}

//-----------------------------------------------------------------------------
int CableNetwork::size()
{
    return sep[nc];
}

//-----------------------------------------------------------------------------
int CableNetwork::squeeze()
{
    int c = 0, i = 0, cnt = 0;
    for (int c0 = 0; c0 < nc; c0++) {
        if (c0 == group_sep)
            group_sep -= cnt;
        if (sep[c0+1]-sep[c0] > 1) {
            for (int i0 = sep[c0]; i0 < sep[c0+1]; i0++) {
                idx[i++] = idx[i0];
            }
            sep[++c] = i;
        } else {
            cnt++;
        }
    }
    sep[c] = i;
    nc = c;
    if (group_sep >= nc)
        group_sep = nc;
    return cnt;
}

//-----------------------------------------------------------------------------
int CableNetwork::remove_negative_indices()
{
    int i = 0, cnt = 0;
    for (int c = 0; c < nc; c++) {
        int start = i;
        for (int i0 = sep[c]; i0 < sep[c+1]; i0++) {
            if (idx[i0] >= 0)
                idx[i++] = idx[i0];
            else
                cnt++;
        }
        sep[c] = start;
    }
    sep[nc] = i;
    return cnt;
}

//-----------------------------------------------------------------------------
int CableNetwork::largest_index()
{
    int imax = 0;
    for (int i=0; i<sep[nc]; i++) if (idx[i] > imax) imax = idx[i];
    return imax;
}

//-----------------------------------------------------------------------------
void CableNetwork::count_neighbors(int nv, int* nb_neigh)
{
    for (int i=0;i<nv;i++)
        nb_neigh[i] = 0;

    for (int n=0;n<nc;n++) {
        for (int i=sep[n];i<sep[n+1]-1;i++) {
            if ((idx[i] == -1) || (idx[i+1] == -1)) continue;
            nb_neigh[idx[i]]++;
            nb_neigh[idx[i+1]]++;
        }
    }
}

//-----------------------------------------------------------------------------
int most_freq_val(int n, int *x, int bins, int &cnt_less_freq)
// negative values are discarded
{
    int* count = new int [bins];
    for (int i=0;i<bins;i++)
        count[i] = 0;
    for (int i=0;i<n;i++)
        if (x[i] >= 0) count[x[i]]++;
    
    int maxcount = 0, argmax = -1;
    for (int i=0;i<bins;i++)
        if (count[i] > maxcount) {
            maxcount = count[i];
            argmax = i;
        }

    cnt_less_freq = 0;
    for (int i=0;i<bins;i++)
        if (i != argmax) cnt_less_freq += count[i];

    delete [] count;
    return argmax;
}

//-----------------------------------------------------------------------------
 int CableNetwork::remove_isolated_regions(int &n_comp)
 {
    UnionFind UF;
    int nv = largest_index() + 1;
    UF.allocate(nv);
    int *label = new int [nv];
    UF.reset(label);

    for (int n=0;n<nc;n++) {
        if (sep[n+1] - sep[n] == 1) 
            UF.mark(idx[sep[n]]);
        else
            for (int i = sep[n]; i < sep[n+1]-1; i++)
                UF.unite(idx[i], idx[i+1]);
    }
    
    int cnt = 0;
    n_comp = UF.assign_label();

    if (n_comp > 1) {
        // find the biggest connection region
        int main_comp = most_freq_val(nv, label, n_comp, cnt);

        // remove these nodes
        for (int i=0;i<sep[nc];i++) {
            if (label[idx[i]] != main_comp) {
                idx[i] = -1;
            }
        }
        remove_negative_indices();
    }

    delete [] label;
    return cnt;
 }

//-----------------------------------------------------------------------------
int CableNetwork::tag_loose_ends(int* nb_neigh)
{
    int cnt = 0;
    for (int n=0;n<nc;n++) {
        int k1 = sep[n];
        int k2 = sep[n+1]-1;
        int size = k2-k1+1;
        if (size < 1) continue;
        // remove cables of size one
        if (size == 1) {
            idx[k1] = -1;
            cnt++;
            continue;
        }
        // cut the beginning of the cable
        if (nb_neigh[idx[k1]] == 1) {
            idx[k1] = -1;
            cnt++;
            k1++;
            while (k1 < k2 && nb_neigh[idx[k1]] == 2) {
                idx[k1] = -1;
                cnt++;
                k1++;
            }
        }
        // cut the end of the cable
        k1 = sep[n];
        if (nb_neigh[idx[k2]] == 1) {
            idx[k2] = -1;
            cnt++;
            k2--;
            while (k1 < k2 && nb_neigh[idx[k2]] == 2) {
                if (idx[k2] == -1) break;
                idx[k2] = -1;
                cnt++;
                k2--;
            }
        }
    }
    return cnt;
}

//-----------------------------------------------------------------------------
int CableNetwork::cut_loose_ends(int &niter)
{
    int nv = largest_index() + 1;
    int* nb_neigh = new int [nv];
    int cnt = -1, total_cnt = 0;
    niter = -1;
    while (cnt != 0) {
        niter ++;
        count_neighbors(nv, nb_neigh);
        cnt = tag_loose_ends(nb_neigh);
        remove_negative_indices();
        total_cnt += cnt;
    }

    delete [] nb_neigh;
    return total_cnt;
}
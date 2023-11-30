#ifndef CABLE_NETWORK_H_
#define CABLE_NETWORK_H_

class CableNetwork {
public:
    // the indices of cable c are:
    // idx[sep[c]], idx[sep[c]+1], ..., idx[sep[c+1]-1]
    int nc; // numbre of cables
    int* idx; // node indices of all cables (concatenated)
    int* sep; // size nc+1
    int group_sep; // number of cables in the first group

    CableNetwork();
    ~CableNetwork();

    void allocate(int capacity_cables, int capacity_indices);
    void clear();
    void append(int index);
    void new_cable();
    void new_group();

    int size(); // total number of indices
    int largest_index();
    void count_neighbors(int nv, int* nb_neigh);

    int squeeze(); // remove cables of length 0 or 1
    int remove_negative_indices();
    int remove_isolated_regions(int &n_comp);
    int tag_loose_ends(int* nb_neigh);
    int cut_loose_ends(int &niter);

//private:
    int nc_max, nv_max;
};

#endif
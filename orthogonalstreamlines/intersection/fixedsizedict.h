#ifndef FIXEDSIZEDICT_H_
#define FIXEDSIZEDICT_H_

class FixedSizeDict {
public:
    FixedSizeDict();
    ~FixedSizeDict();

    void allocate_bins(int nbins_);
    void set_bin_capacity(int bin_id, int capacity);

    void add(int bin_id, int a, int b, int c, int d, int value);
    int find(int bin_id, int a, int b, int c, int d);

//private:
    int nbins;
    int *capacity;
    int *size;
    int **keys;
    int **values;
};

#endif
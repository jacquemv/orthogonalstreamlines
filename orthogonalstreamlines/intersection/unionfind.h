#ifndef _Union_Find
#define _Union_Find

class UnionFind {
public:
	int n;
	int* parent; // parent[i] = parent of i
	int* size;	 // size[i] = number of sites in tree rooted at i
	int count;	 // number of components
	int* label; // component id

	UnionFind();
	~UnionFind();
	void allocate(int n_);
	void reset(int* label_);
	int find(int p);
	void mark(int p);
	void unite(int p, int q);
	bool connected(int p, int q);
	int assign_label(); // overwrite size
};
#endif
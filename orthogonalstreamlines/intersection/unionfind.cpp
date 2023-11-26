// Weighted Quick-Union with path compression
// based on code by Robert Sedgewick and Kevin Wayne
#include <stdlib.h>
#include <stdio.h>
#include "unionfind.h"

UnionFind::UnionFind() {
	n = count = 0;
	parent = size = label = NULL;
}

UnionFind::~UnionFind()
{
	delete [] parent;
	delete [] size;
}

void UnionFind::allocate(int n_) {
	n = n_;
	parent = new int [n];
	size = new int [n];
}

void UnionFind::reset(int* label_)
{
	label = label_;
	for (int i = 0; i < n; i++) {
		parent[i] = i;
		size[i] = 1;
		label[i] = -1;
	}
	count = n;
}

int UnionFind::find(int p) {
	int root = p;
	while (root != parent[root])
		root = parent[root];
	while (p != root) {
		int newp = parent[p];
		parent[p] = root;
		p = newp;
	}
	return root;
}

void UnionFind::mark(int p) {
	label[p] = 0;
}

void UnionFind::unite(int p, int q) {
	int rootP = find(p);
	int rootQ = find(q);
	label[p] = label[q] = 0;
	if (rootP == rootQ) return;

	// make smaller root point to larger one
	if (size[rootP] < size[rootQ]) {
		parent[rootP] = rootQ;
		size[rootQ] += size[rootP];
	}
	else {
		parent[rootQ] = rootP;
		size[rootP] += size[rootQ];
	}
	count--;
}

int UnionFind::assign_label()
{
	int i, max_label = 0;
	for (i = 0; i < n; i++) size[i] = -1;
	for (i = 0; i < n; i++) if (label[i] >= 0) {
		int p = find(i);
		if (size[p] < 0) size[p] = max_label++;
		label[i] = size[p];
	}
	return max_label;
}


	 

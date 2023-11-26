#ifndef TRIANGULATE_FACETS_H_
#define TRIANGULATE_FACETS_H_

class TriangulateFacets {
public:
	double* ver;
	int nv;
	int* tri;
	int nt, nt_max;
	
	// parameters
	double max_dihedral_thres, max_dihedral_incr;
	
	// triangulation
	TriangulateFacets() {}
	~TriangulateFacets();
	TriangulateFacets(int nv_, double* ver_, int nt_max_, int* tri_);
	void initialize(int nv_, double* ver_, int nt_max_, int* tri_);

	int triangulate_facet(int n, int* facet);

//private:
	// edge and triangle calculations
	void edge_vec(int i, int j, double* x);
	inline double distance2(int i, int j);
	void tri_normal(int i, int j, int k, double* n);
	void facet_normal(int n, int* facet, double* normal);
	double edge_angle(int i, int j, int k1, int k2);
	double min_angle(int i, int j, int k);
	double vertex_angle(int i, int j, int k);
	
	// triangulation
	void insert_triangle(int i, int j, int k);
	int triangulate_polygon2d(double* ver, int nv, int** tri, int* nt);
	void triangulate_small_facet(int n, int* facet, double thres);
	int triangulate_large_facet(int n, int* facet);
	
	// utilities
	double to_degree(double angle);
	double from_degree(double angle);
	
	bool _allocate;
};

#endif
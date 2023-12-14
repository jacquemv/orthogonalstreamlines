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
	int max_iter;
	
	// triangulation
	TriangulateFacets() {}
	~TriangulateFacets();
	TriangulateFacets(int nv_, double* ver_, int nt_max_, int* tri_);
	void initialize(int nv_, double* ver_, int nt_max_, int* tri_);

	int triangulate_facet(int n, int* facet); // returns #tri or -1

//private:
	// edge and triangle calculations
	void edge_vec(int i, int j, double* x);
	inline double distance2(int i, int j);
	void tri_normal(int i, int j, int k, double* n);
	void facet_normal(int n, int* facet, double* normal);
	double edge_angle(int i, int j, int k1, int k2);
	double min_angle(int i, int j, int k);
	double vertex_angle(int i, int j, int k);
	double min_of_min_angles(int n, int* facet, int* idx);
	double max_of_dihedral_angles(int n, int* facet, int* idx);

	// triangulation
	void insert_triangle(int i, int j, int k);
	void insert_facet_triangulation(int n, int* facet, int* idx);
	int triangulate_polygon2d_polypart(double* ver, int nv, int** tri, int* nt);
	int triangulate_polygon2d_earcut(double* ver, int nv, int** tri, int* nt);
	int triangulate_small_facet(int n, int* facet, double thres, int iter);
	int triangulate_large_facet(int n, int* facet);
	void fix_triangle_orientation(int n, int* facet, int* idx, double* normal);

	// utilities
	double to_degree(double angle);
	double from_degree(double angle);
	
	bool _allocate;
};

#endif
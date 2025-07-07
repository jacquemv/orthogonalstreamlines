#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cassert>
#include <array>
#include <cfloat>

#include "lookup_tables.h"
#include "triangulatefacets.h"
#include "algebra.cpp"

#include "polypartition.cpp"	
#include "earcut.cpp"
using Coord = double;
using N = int;
using Point = std::array<Coord, 2>;
#ifndef M_PI // for MSVC compiler
	#define M_PI 3.14159265358979323846264338327950288419716
#endif

//-----------------------------------------------------------------------------
TriangulateFacets::TriangulateFacets(int nv_, double* ver_, int nt_max_, 
									 int* tri_)
{
	initialize(nv_, ver_, nt_max_, tri_);
}

//-----------------------------------------------------------------------------
void TriangulateFacets::initialize(int nv_, double* ver_, int nt_max_, 
								   int* tri_)
{
	nv = nv_;
	nt = 0;
	nt_max = nt_max_;
	ver = ver_;
	_allocate = 0;
	if (tri_ == NULL) {
		tri = new int [3*nt_max];
		_allocate = 1;
	} else
		tri = tri_;
	max_dihedral_thres = 0;
	max_dihedral_incr = 0.1;
	max_iter = (int)(1.0/max_dihedral_incr+0.5) + 1;
}

//-----------------------------------------------------------------------------
TriangulateFacets::~TriangulateFacets()
{
	if (_allocate) delete [] tri;
}

//-----------------------------------------------------------------------------
inline double TriangulateFacets::distance2(int i, int j)
{
	double d2, dx;
	dx = ver[3*i]-ver[3*j];
	d2 = dx*dx;
	dx = ver[3*i+1]-ver[3*j+1];
	d2 += dx*dx;
	dx = ver[3*i+2]-ver[3*j+2];
	d2 += dx*dx;
	return d2;
}

//-----------------------------------------------------------------------------
void TriangulateFacets::edge_vec(int i, int j, double* x)
{
	x[0] = ver[3*j]   - ver[3*i];
	x[1] = ver[3*j+1] - ver[3*i+1];
	x[2] = ver[3*j+2] - ver[3*i+2];
}

//-----------------------------------------------------------------------------
void TriangulateFacets::tri_normal(int i, int j, int k, double* normal)
{
	double x1[3], x2[3];
	edge_vec(i, j, x1);
	edge_vec(i, k, x2);
	vcross(normal, x1, x2);
	vnormalize(normal);
}

//-----------------------------------------------------------------------------
void TriangulateFacets::facet_normal(int n, int* facet, double* normal)
{
	double L[3], dx[3];
	edge_vec(facet[n-1], facet[0], dx);
	vcross(normal, &ver[3*facet[n-1]], dx);
	for (int i=0;i<n-1;i++) {
		edge_vec(facet[i], facet[i+1], dx);
		vcross(L, &ver[3*facet[i]], dx);
		normal[0] += L[0]; normal[1] += L[1]; normal[2] += L[2];
	}
	vnormalize(normal);
}

//-----------------------------------------------------------------------------
double TriangulateFacets::edge_angle(int i, int j, int k1, int k2)
// actually computes -cos(angle), which is a monotonic function of the angle
// on the interval 0 to pi
{
	double n1[3], n2[3];
	tri_normal(i, j, k1, n1);
	tri_normal(j, i, k2, n2);
	double s = -vdot(n1, n2);
	if (isfinite(s)) return s;
	return 1;
}

//-----------------------------------------------------------------------------
double TriangulateFacets::to_degree(double angle)
{
	return acos(-angle)/M_PI*180; // note the minus sign
}

//-----------------------------------------------------------------------------
double TriangulateFacets::from_degree(double angle)
{
	return -cos(angle*M_PI/180); // note the minus sign
}

//-----------------------------------------------------------------------------
double TriangulateFacets::min_angle(int i, int j, int k)
// actually computes -cos(angle), which is a monotonic function of the angle
// on the interval 0 to pi
{
	double a2 = distance2(i, j);
	double b2 = distance2(j, k);
	double c2 = distance2(k, i);
	double a = sqrt(a2);
	double b = sqrt(b2);
	double c = sqrt(c2);
	double angc = -(a2+b2-c2)/a/b/2; // note the minus sign
	double angb = -(a2+c2-b2)/a/c/2;
	double anga = -(b2+c2-a2)/b/c/2;
	double minang = angc;
	if (angb < minang) minang = angb;
	if (anga < minang) minang = anga;
	if (isfinite(minang)) return minang; 
	return -1; // equivalent to an angle of zero
}

//-----------------------------------------------------------------------------
double TriangulateFacets::vertex_angle(int i, int j, int k)
{ // -cos(angle) at j in the triangle (i,j,k)
	double a2 = distance2(i, j);
	double b2 = distance2(j, k);
	double c2 = distance2(k, i);
	return -(a2+b2-c2)/sqrt(a2*b2)/2; // note the minus sign
}

//-----------------------------------------------------------------------------
void orthogonal_basis(double* normal, double* v1, double* v2)
{
    double v0[3] = {0, 0, 1};
    vcross(v1, v0, normal);
    double norm2 = vnorm2(v1);
    if (norm2 == 0) {
        v0[1] = 1; v0[2] = 0;
        vcross(v1, v0, normal);
        norm2 = vnorm2(v1);
    }
    vscale(v1, 1/sqrt(norm2));
    vcross(v2, normal, v1);
    norm2 = vnorm2(v2);
    vscale(v2, 1/sqrt(norm2));
}

//-----------------------------------------------------------------------------
void points3d_to_2d(double* ver, int nv, double* normal)
{
    //center_points(ver, nv);
    double v1[3], v2[3];
    orthogonal_basis(normal, v1, v2);
    for (int i=0;i<nv;i++) {
        double x = v1[0]*ver[3*i] + v1[1]*ver[3*i+1] + v1[2]*ver[3*i+2];
        double y = v2[0]*ver[3*i] + v2[1]*ver[3*i+1] + v2[2]*ver[3*i+2];
        ver[2*i] = x;
        ver[2*i+1] = y;
    }
}

//-----------------------------------------------------------------------------
void TriangulateFacets::insert_triangle(int i, int j, int k)
{
	assert(nt < nt_max);
	tri[3*nt] = i;
	tri[3*nt+1] = j;
	tri[3*nt+2] = k;
	nt++;
}

//-----------------------------------------------------------------------------
int TriangulateFacets::triangulate_polygon2d_polypart(double* ver, int nv, 
											 		  int** tri, int* nt)
{
	*tri = nullptr;
	*nt = 0;
	TPPLPoly tp;
	tp.Init(nv);
	for (int i=0;i<nv;i++) {
		tp[i].x = ver[2*i];
		tp[i].y = ver[2*i+1];
		tp[i].id = i;
	}
	tp.SetOrientation(TPPL_CCW);

	TPPLPolyList pl;
	TPPLPartition pp;
	int ok = pp.Triangulate_OPT(&tp, &pl);
	if (!ok) return -1;

	TPPLPolyList::iterator iter;
	*nt = pl.size();
	*tri = new int [3*(*nt)];
	int k = 0;
	for (iter=pl.begin(); iter!=pl.end(); iter++) {
		if (iter->GetNumPoints() != 3) return -3;
		for (int j=0;j<3;j++)
			(*tri)[k++] = iter->GetPoint(j).id;
	}
	return 0;
}

//-----------------------------------------------------------------------------
int TriangulateFacets::triangulate_polygon2d_earcut(double* ver, int nv, 
											 		int** tri, int* nt)
{
	*tri = nullptr;
	*nt = 0;
	std::vector<Point> vertices;
	std::vector<std::vector<Point>> polygon;
	for (int i=0;i<nv;i++) {
		Point P = {ver[2*i], ver[2*i+1]};
		vertices.push_back(P);
	}
	polygon.push_back(vertices);
	std::vector<N> indices = mapbox::earcut<N>(polygon);
	int n = (int)(indices.size());
	*nt = n / 3;
	*tri = new int [n];
	for (int k=0;k<n;k++)
		(*tri)[k] = indices[k];
	return 0;
}

//-----------------------------------------------------------------------------
void TriangulateFacets::fix_triangle_orientation(int n, int* facet, int* idx, 
												 double* normal)
{
	for (int j=0;j<n;j++) {
		double vec[3];
		tri_normal(facet[idx[0]], facet[idx[1]], facet[idx[2]], vec);
		if (vdot(vec, normal) < 0) {
			int tmp = idx[1];
			idx[1] = idx[2];
			idx[2] = tmp;
		}
		idx += 3;
	}
}

//-----------------------------------------------------------------------------
double TriangulateFacets::min_of_min_angles(int n, int* facet, int* idx)
// minimum over the triangles of the minimum angle of each triangle
{
	double angle = DBL_MAX;
	for (int j=0;j<n-2;j++) {
		double a = min_angle(facet[idx[0]], facet[idx[1]], facet[idx[2]]);
		if (a < angle) angle = a;
		idx += 3;
	}
	return angle;
}

//-----------------------------------------------------------------------------
double TriangulateFacets::max_of_dihedral_angles(int n, int* facet, int* idx)
{
	double dihedral_angle = -DBL_MAX;
	for (int j=0;j<n-3;j++) {
		double a = edge_angle(facet[idx[0]], facet[idx[1]], 
							  facet[idx[2]], facet[idx[3]]);
		if (a > dihedral_angle) dihedral_angle = a;
		idx += 4;
	}
	return dihedral_angle;
}

//-----------------------------------------------------------------------------
void TriangulateFacets::insert_facet_triangulation(int nt, int* facet, 
														  int* idx)
{
	for (int j=0;j<nt;j++) {
		insert_triangle(facet[idx[0]], facet[idx[1]], facet[idx[2]]);
		idx += 3;
	}
}

//-----------------------------------------------------------------------------
int TriangulateFacets::triangulate_small_facet(int n, int* facet, 
												double thres, int iter)
{
	double best_value = -DBL_MAX;
	int best_id = -1;
	int* table_tri = all_triangulations_tri[n];
	int* table_adj = all_triangulations_edges[n];

	for (int i=0;i<nb_triangulations[n];i++) {
		double inner_angle = min_of_min_angles(n, facet, &table_tri[3*(n-2)*i]);
		double dihedral_angle = max_of_dihedral_angles(n, facet, 
													   &table_adj[4*(n-3)*i]);
		// consider first only "flat" trangulations, then progressively allow
		// larger dihedral angles
		if (dihedral_angle > thres) continue;
		
		if (inner_angle > best_value) {
			best_id = i;
			best_value = inner_angle;
		}
	}

	if ((best_id < 0) && (iter < max_iter)) {
		// if it doesn't work, releave the constraint
		return triangulate_small_facet(n, facet, thres + max_dihedral_incr, 
									   iter+1);
	}

	int nt = n-2;
	int *idx = table_tri + 3*nt*best_id;
	if (iter <= max_iter / 2) {
		insert_facet_triangulation(nt, facet, idx);
	} else {
		int* idx_copy = new int [3*nt];
		for (int j=0;j<3*nt;j++) idx_copy[j] = idx[j];

		double normal[3];
		facet_normal(n, facet, normal);

		fix_triangle_orientation(nt, facet, idx_copy, normal);
		insert_facet_triangulation(nt, facet, idx_copy);
		delete [] idx_copy;
	}
	return nt;
}

//-----------------------------------------------------------------------------
int TriangulateFacets::triangulate_large_facet(int n, int* facet)
{
	double normal[3];
	facet_normal(n, facet, normal);
	
	double* polyver = new double [3*n];
	for (int i=0;i<n;i++) {
		int j = facet[i];
		polyver[3*i] = ver[3*j];
		polyver[3*i+1] = ver[3*j+1];
		polyver[3*i+2] = ver[3*j+2];
	}
	points3d_to_2d(polyver, n, normal);
	int* idx;
	int nt, err;

	err = triangulate_polygon2d_polypart(polyver, n, &idx, &nt);
	if (err < 0)
		err = triangulate_polygon2d_earcut(polyver, n, &idx, &nt);
	
	if (err == 0) {
		insert_facet_triangulation(nt, facet, idx);
		fix_triangle_orientation(nt, facet, idx, normal);
		delete [] idx; // allocated only if no error
	}
	delete [] polyver;
	return err < 0 ? err : nt;
}

//-----------------------------------------------------------------------------
int TriangulateFacets::triangulate_facet(int n, int* facet)
{
	if (n<3) return 0;
	if (n == 3) {
		insert_triangle(facet[0], facet[1], facet[2]);
		return 1;
	}
	if (n <= 10)
		return triangulate_small_facet(n, facet, max_dihedral_thres, 0);
	return triangulate_large_facet(n, facet);
}
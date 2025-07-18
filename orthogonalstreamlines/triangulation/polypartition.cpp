//Copyright (C) 2011 by Ivan Fratric
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in
//all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//THE SOFTWARE.

// https://github.com/ivanfratric/polypartition/tree/master

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <list>
#include <algorithm>
#include <set>
#include <vector>
#include <stdexcept>

using namespace std;

#include "polypartition.h"

#define TPPL_VERTEXTYPE_REGULAR 0
#define TPPL_VERTEXTYPE_START 1
#define TPPL_VERTEXTYPE_END 2
#define TPPL_VERTEXTYPE_SPLIT 3
#define TPPL_VERTEXTYPE_MERGE 4

TPPLPoly::TPPLPoly() {
	hole = false;
	numpoints = 0;
	points = NULL;
}

TPPLPoly::~TPPLPoly() {
	if(points) delete [] points;
}

void TPPLPoly::Clear() {
	if(points) delete [] points;
	hole = false;
	numpoints = 0;
	points = NULL;
}

void TPPLPoly::Init(long numpoints) {
	Clear();
	this->numpoints = numpoints;
	points = new TPPLPoint[numpoints];
}

void TPPLPoly::Triangle(TPPLPoint &p1, TPPLPoint &p2, TPPLPoint &p3) {
	Init(3);
	points[0] = p1;
	points[1] = p2;
	points[2] = p3;
}

TPPLPoly::TPPLPoly(const TPPLPoly &src) : TPPLPoly() {
	hole = src.hole;
	numpoints = src.numpoints;

	if(numpoints > 0) {
		points = new TPPLPoint[numpoints];
		memcpy(points, src.points, numpoints*sizeof(TPPLPoint));
	}
}

TPPLPoly& TPPLPoly::operator=(const TPPLPoly &src) {
	Clear();
	hole = src.hole;
	numpoints = src.numpoints;

	if(numpoints > 0) {
		points = new TPPLPoint[numpoints];
		memcpy(points, src.points, numpoints*sizeof(TPPLPoint));
	}

	return *this;
}

int TPPLPoly::GetOrientation() const {
	long i1,i2;
	tppl_float area = 0;
	for(i1=0; i1<numpoints; i1++) {
		i2 = i1+1;
		if(i2 == numpoints) i2 = 0;
		area += points[i1].x * points[i2].y - points[i1].y * points[i2].x;
	}
	if(area>0) return TPPL_CCW;
	if(area<0) return TPPL_CW;
	return 0;
}

void TPPLPoly::SetOrientation(int orientation) {
	int polyorientation = GetOrientation();
	if(polyorientation&&(polyorientation!=orientation)) {
		Invert();
	}
}

void TPPLPoly::Invert() {
	std::reverse(points, points + numpoints);
}

TPPLPartition::PartitionVertex::PartitionVertex() : previous(NULL), next(NULL) {

}

TPPLPoint TPPLPartition::Normalize(const TPPLPoint &p) {
	TPPLPoint r;
	tppl_float n = sqrt(p.x*p.x + p.y*p.y);
	if(n!=0) {
		r = p/n;
	} else {
		r.x = 0;
		r.y = 0;
	}
	return r;
}

tppl_float TPPLPartition::Distance(const TPPLPoint &p1, const TPPLPoint &p2) {
	tppl_float dx,dy;
	dx = p2.x - p1.x;
	dy = p2.y - p1.y;
	return(sqrt(dx*dx + dy*dy));
}

//checks if two lines intersect
int TPPLPartition::Intersects(TPPLPoint &p11, TPPLPoint &p12, TPPLPoint &p21, TPPLPoint &p22) {
	if((p11.x == p21.x)&&(p11.y == p21.y)) return 0;
	if((p11.x == p22.x)&&(p11.y == p22.y)) return 0;
	if((p12.x == p21.x)&&(p12.y == p21.y)) return 0;
	if((p12.x == p22.x)&&(p12.y == p22.y)) return 0;

	TPPLPoint v1ort,v2ort,v;
	tppl_float dot11,dot12,dot21,dot22;

	v1ort.x = p12.y-p11.y;
	v1ort.y = p11.x-p12.x;

	v2ort.x = p22.y-p21.y;
	v2ort.y = p21.x-p22.x;

	v = p21-p11;
	dot21 = v.x*v1ort.x + v.y*v1ort.y;
	v = p22-p11;
	dot22 = v.x*v1ort.x + v.y*v1ort.y;

	v = p11-p21;
	dot11 = v.x*v2ort.x + v.y*v2ort.y;
	v = p12-p21;
	dot12 = v.x*v2ort.x + v.y*v2ort.y;

	if(dot11*dot12>0) return 0;
	if(dot21*dot22>0) return 0;

	return 1;
}

//removes holes from inpolys by merging them with non-holes
int TPPLPartition::RemoveHoles(TPPLPolyList *inpolys, TPPLPolyList *outpolys) {
	TPPLPolyList polys;
	TPPLPolyList::iterator holeiter,polyiter,iter,iter2;
	long i,i2,holepointindex,polypointindex=0;
	TPPLPoint holepoint,polypoint,bestpolypoint={0.0, 0.0, 0};
	TPPLPoint linep1,linep2;
	TPPLPoint v1,v2;
	TPPLPoly newpoly;
	bool hasholes;
	bool pointvisible;
	bool pointfound;

	//check for trivial case (no holes)
	hasholes = false;
	for(iter = inpolys->begin(); iter!=inpolys->end(); iter++) {
		if(iter->IsHole()) {
			hasholes = true;
			break;
		}
	}
	if(!hasholes) {
		for(iter = inpolys->begin(); iter!=inpolys->end(); iter++) {
			outpolys->push_back(*iter);
		}
		return 1;
	}

	polys = *inpolys;

	while(1) {
		//find the hole point with the largest x
		hasholes = false;
		for(iter = polys.begin(); iter!=polys.end(); iter++) {
			if(!iter->IsHole()) continue;

			if(!hasholes) {
				hasholes = true;
				holeiter = iter;
				holepointindex = 0;
			}

			for(i=0; i < iter->GetNumPoints(); i++) {
				if(iter->GetPoint(i).x > holeiter->GetPoint(holepointindex).x) {
					holeiter = iter;
					holepointindex = i;
				}
			}
		}
		if(!hasholes) break;
		holepoint = holeiter->GetPoint(holepointindex);

		pointfound = false;
		for(iter = polys.begin(); iter!=polys.end(); iter++) {
			if(iter->IsHole()) continue;
			for(i=0; i < iter->GetNumPoints(); i++) {
				if(iter->GetPoint(i).x <= holepoint.x) continue;
				if(!InCone(iter->GetPoint((i+iter->GetNumPoints()-1)%(iter->GetNumPoints())),
					iter->GetPoint(i),
					iter->GetPoint((i+1)%(iter->GetNumPoints())),
					holepoint))
					continue;
				polypoint = iter->GetPoint(i);
				if(pointfound) {
					v1 = Normalize(polypoint-holepoint);
					v2 = Normalize(bestpolypoint-holepoint);
					if(v2.x > v1.x) continue;
				}
				pointvisible = true;
				for(iter2 = polys.begin(); iter2!=polys.end(); iter2++) {
					if(iter2->IsHole()) continue;
					for(i2=0; i2 < iter2->GetNumPoints(); i2++) {
						linep1 = iter2->GetPoint(i2);
						linep2 = iter2->GetPoint((i2+1)%(iter2->GetNumPoints()));
						if(Intersects(holepoint,polypoint,linep1,linep2)) {
							pointvisible = false;
							break;
						}
					}
					if(!pointvisible) break;
				}
				if(pointvisible) {
					pointfound = true;
					bestpolypoint = polypoint;
					polyiter = iter;
					polypointindex = i;
				}
			}
		}

		if(!pointfound) return 0;

		newpoly.Init(holeiter->GetNumPoints() + polyiter->GetNumPoints() + 2);
		i2 = 0;
		for(i=0;i<=polypointindex;i++) {
			newpoly[i2] = polyiter->GetPoint(i);
			i2++;
		}
		for(i=0;i<=holeiter->GetNumPoints();i++) {
			newpoly[i2] = holeiter->GetPoint((i+holepointindex)%holeiter->GetNumPoints());
			i2++;
		}
		for(i=polypointindex;i<polyiter->GetNumPoints();i++) {
			newpoly[i2] = polyiter->GetPoint(i);
			i2++;
		}

		polys.erase(holeiter);
		polys.erase(polyiter);
		polys.push_back(newpoly);
	}

	for(iter = polys.begin(); iter!=polys.end(); iter++) {
		outpolys->push_back(*iter);
	}

	return 1;
}

bool TPPLPartition::IsConvex(TPPLPoint& p1, TPPLPoint& p2, TPPLPoint& p3) {
	tppl_float tmp;
	tmp = (p3.y-p1.y)*(p2.x-p1.x)-(p3.x-p1.x)*(p2.y-p1.y);
	if(tmp>0) return 1;
	else return 0;
}

bool TPPLPartition::IsReflex(TPPLPoint& p1, TPPLPoint& p2, TPPLPoint& p3) {
	tppl_float tmp;
	tmp = (p3.y-p1.y)*(p2.x-p1.x)-(p3.x-p1.x)*(p2.y-p1.y);
	if(tmp<0) return 1;
	else return 0;
}

bool TPPLPartition::IsInside(TPPLPoint& p1, TPPLPoint& p2, TPPLPoint& p3, TPPLPoint &p) {
	if(IsConvex(p1,p,p2)) return false;
	if(IsConvex(p2,p,p3)) return false;
	if(IsConvex(p3,p,p1)) return false;
	return true;
}

bool TPPLPartition::InCone(TPPLPoint &p1, TPPLPoint &p2, TPPLPoint &p3, TPPLPoint &p) {
	bool convex;

	convex = IsConvex(p1,p2,p3);

	if(convex) {
		if(!IsConvex(p1,p2,p)) return false;
		if(!IsConvex(p2,p3,p)) return false;
		return true;
	} else {
		if(IsConvex(p1,p2,p)) return true;
		if(IsConvex(p2,p3,p)) return true;
		return false;
	}
}

bool TPPLPartition::InCone(PartitionVertex *v, TPPLPoint &p) {
	TPPLPoint p1,p2,p3;

	p1 = v->previous->p;
	p2 = v->p;
	p3 = v->next->p;

	return InCone(p1,p2,p3,p);
}

void TPPLPartition::UpdateVertexReflexity(PartitionVertex *v) {
	PartitionVertex *v1 = NULL,*v3 = NULL;
	v1 = v->previous;
	v3 = v->next;
	v->isConvex = !IsReflex(v1->p,v->p,v3->p);
}

void TPPLPartition::UpdateVertex(PartitionVertex *v, PartitionVertex *vertices, long numvertices) {
	long i;
	PartitionVertex *v1 = NULL,*v3 = NULL;
	TPPLPoint vec1,vec3;

	v1 = v->previous;
	v3 = v->next;

	v->isConvex = IsConvex(v1->p,v->p,v3->p);

	vec1 = Normalize(v1->p - v->p);
	vec3 = Normalize(v3->p - v->p);
	v->angle = vec1.x*vec3.x + vec1.y*vec3.y;

	if(v->isConvex) {
		v->isEar = true;
		for(i=0;i<numvertices;i++) {
			if((vertices[i].p.x==v->p.x)&&(vertices[i].p.y==v->p.y)) continue;
			if((vertices[i].p.x==v1->p.x)&&(vertices[i].p.y==v1->p.y)) continue;
			if((vertices[i].p.x==v3->p.x)&&(vertices[i].p.y==v3->p.y)) continue;
			if(IsInside(v1->p,v->p,v3->p,vertices[i].p)) {
				v->isEar = false;
				break;
			}
		}
	} else {
		v->isEar = false;
	}
}

//triangulation by ear removal
int TPPLPartition::Triangulate_EC(TPPLPoly *poly, TPPLPolyList *triangles) {
	if(!poly->Valid()) return 0;

	long numvertices;
	PartitionVertex *vertices = NULL;
	PartitionVertex *ear = NULL;
	TPPLPoly triangle;
	long i,j;
	bool earfound;

	if(poly->GetNumPoints() < 3) return 0;
	if(poly->GetNumPoints() == 3) {
		triangles->push_back(*poly);
		return 1;
	}

	numvertices = poly->GetNumPoints();

	vertices = new PartitionVertex[numvertices];
	for(i=0;i<numvertices;i++) {
		vertices[i].isActive = true;
		vertices[i].p = poly->GetPoint(i);
		if(i==(numvertices-1)) vertices[i].next=&(vertices[0]);
		else vertices[i].next=&(vertices[i+1]);
		if(i==0) vertices[i].previous = &(vertices[numvertices-1]);
		else vertices[i].previous = &(vertices[i-1]);
	}
	for(i=0;i<numvertices;i++) {
		UpdateVertex(&vertices[i],vertices,numvertices);
	}

	for(i=0;i<numvertices-3;i++) {
		earfound = false;
		//find the most extruded ear
		for(j=0;j<numvertices;j++) {
			if(!vertices[j].isActive) continue;
			if(!vertices[j].isEar) continue;
			if(!earfound) {
				earfound = true;
				ear = &(vertices[j]);
			} else {
				if(vertices[j].angle > ear->angle) {
					ear = &(vertices[j]);
				}
			}
		}
		if(!earfound) {
			delete [] vertices;
			return 0;
		}

		triangle.Triangle(ear->previous->p,ear->p,ear->next->p);
		triangles->push_back(triangle);

		ear->isActive = false;
		ear->previous->next = ear->next;
		ear->next->previous = ear->previous;

		if(i==numvertices-4) break;

		UpdateVertex(ear->previous,vertices,numvertices);
		UpdateVertex(ear->next,vertices,numvertices);
	}
	for(i=0;i<numvertices;i++) {
		if(vertices[i].isActive) {
			triangle.Triangle(vertices[i].previous->p,vertices[i].p,vertices[i].next->p);
			triangles->push_back(triangle);
			break;
		}
	}

	delete [] vertices;

	return 1;
}

int TPPLPartition::Triangulate_EC(TPPLPolyList *inpolys, TPPLPolyList *triangles) {
	TPPLPolyList outpolys;
	TPPLPolyList::iterator iter;

	if(!RemoveHoles(inpolys,&outpolys)) return 0;
	for(iter=outpolys.begin();iter!=outpolys.end();iter++) {
		if(!Triangulate_EC(&(*iter),triangles)) return 0;
	}
	return 1;
}

int TPPLPartition::ConvexPartition_HM(TPPLPoly *poly, TPPLPolyList *parts) {
	if(!poly->Valid()) return 0;

	TPPLPolyList triangles;
	TPPLPolyList::iterator iter1,iter2;
	TPPLPoly *poly1 = NULL,*poly2 = NULL;
	TPPLPoly newpoly;
	TPPLPoint d1,d2,p1,p2,p3;
	long i11,i12,i21,i22,i13,i23,j,k;
	bool isdiagonal;
	long numreflex;

	//check if the poly is already convex
	numreflex = 0;
	for(i11=0;i11<poly->GetNumPoints();i11++) {
		if(i11==0) i12 = poly->GetNumPoints()-1;
		else i12=i11-1;
		if(i11==(poly->GetNumPoints()-1)) i13=0;
		else i13=i11+1;
		if(IsReflex(poly->GetPoint(i12),poly->GetPoint(i11),poly->GetPoint(i13))) {
			numreflex = 1;
			break;
		}
	}
	if(numreflex == 0) {
		parts->push_back(*poly);
		return 1;
	}

	if(!Triangulate_EC(poly,&triangles)) return 0;

	for(iter1 = triangles.begin(); iter1 != triangles.end(); iter1++) {
		poly1 = &(*iter1);
		for(i11=0;i11<poly1->GetNumPoints();i11++) {
			d1 = poly1->GetPoint(i11);
			i12 = (i11+1)%(poly1->GetNumPoints());
			d2 = poly1->GetPoint(i12);

			isdiagonal = false;
			for(iter2 = iter1; iter2 != triangles.end(); iter2++) {
				if(iter1 == iter2) continue;
				poly2 = &(*iter2);

				for(i21=0;i21<poly2->GetNumPoints();i21++) {
					if((d2.x != poly2->GetPoint(i21).x)||(d2.y != poly2->GetPoint(i21).y)) continue;
					i22 = (i21+1)%(poly2->GetNumPoints());
					if((d1.x != poly2->GetPoint(i22).x)||(d1.y != poly2->GetPoint(i22).y)) continue;
					isdiagonal = true;
					break;
				}
				if(isdiagonal) break;
			}

			if(!isdiagonal) continue;

			p2 = poly1->GetPoint(i11);
			if(i11 == 0) i13 = poly1->GetNumPoints()-1;
			else i13 = i11-1;
			p1 = poly1->GetPoint(i13);
			if(i22 == (poly2->GetNumPoints()-1)) i23 = 0;
			else i23 = i22+1;
			p3 = poly2->GetPoint(i23);

			if(!IsConvex(p1,p2,p3)) continue;

			p2 = poly1->GetPoint(i12);
			if(i12 == (poly1->GetNumPoints()-1)) i13 = 0;
			else i13 = i12+1;
			p3 = poly1->GetPoint(i13);
			if(i21 == 0) i23 = poly2->GetNumPoints()-1;
			else i23 = i21-1;
			p1 = poly2->GetPoint(i23);

			if(!IsConvex(p1,p2,p3)) continue;

			newpoly.Init(poly1->GetNumPoints()+poly2->GetNumPoints()-2);
			k = 0;
			for(j=i12;j!=i11;j=(j+1)%(poly1->GetNumPoints())) {
				newpoly[k] = poly1->GetPoint(j);
				k++;
			}
			for(j=i22;j!=i21;j=(j+1)%(poly2->GetNumPoints())) {
				newpoly[k] = poly2->GetPoint(j);
				k++;
			}

			triangles.erase(iter2);
			*iter1 = newpoly;
			poly1 = &(*iter1);
			i11 = -1;

			continue;
		}
	}

	for(iter1 = triangles.begin(); iter1 != triangles.end(); iter1++) {
		parts->push_back(*iter1);
	}

	return 1;
}

int TPPLPartition::ConvexPartition_HM(TPPLPolyList *inpolys, TPPLPolyList *parts) {
	TPPLPolyList outpolys;
	TPPLPolyList::iterator iter;

	if(!RemoveHoles(inpolys,&outpolys)) return 0;
	for(iter=outpolys.begin();iter!=outpolys.end();iter++) {
		if(!ConvexPartition_HM(&(*iter),parts)) return 0;
	}
	return 1;
}

//minimum-weight polygon triangulation by dynamic programming
//O(n^3) time complexity
//O(n^2) space complexity
int TPPLPartition::Triangulate_OPT(TPPLPoly *poly, TPPLPolyList *triangles) {
	if(!poly->Valid()) return 0;

	long i,j,k,gap,n;
	DPState **dpstates = NULL;
	TPPLPoint p1,p2,p3,p4;
	long bestvertex;
	tppl_float weight,minweight=0,d1,d2;
	Diagonal diagonal,newdiagonal;
	DiagonalList diagonals;
	TPPLPoly triangle;
	int ret = 1;

	n = poly->GetNumPoints();
	dpstates = new DPState *[n];
	for(i=1;i<n;i++) {
		dpstates[i] = new DPState[i];
	}

	//init states and visibility
	for(i=0;i<(n-1);i++) {
		p1 = poly->GetPoint(i);
		for(j=i+1;j<n;j++) {
			dpstates[j][i].visible = true;
			dpstates[j][i].weight = 0;
			dpstates[j][i].bestvertex = -1;
			if(j!=(i+1)) {
				p2 = poly->GetPoint(j);

				//visibility check
				if(i==0) p3 = poly->GetPoint(n-1);
				else p3 = poly->GetPoint(i-1);
				if(i==(n-1)) p4 = poly->GetPoint(0);
				else p4 = poly->GetPoint(i+1);
				if(!InCone(p3,p1,p4,p2)) {
					dpstates[j][i].visible = false;
					continue;
				}

				if(j==0) p3 = poly->GetPoint(n-1);
				else p3 = poly->GetPoint(j-1);
				if(j==(n-1)) p4 = poly->GetPoint(0);
				else p4 = poly->GetPoint(j+1);
				if(!InCone(p3,p2,p4,p1)) {
					dpstates[j][i].visible = false;
					continue;
				}

				for(k=0;k<n;k++) {
					p3 = poly->GetPoint(k);
					if(k==(n-1)) p4 = poly->GetPoint(0);
					else p4 = poly->GetPoint(k+1);
					if(Intersects(p1,p2,p3,p4)) {
						dpstates[j][i].visible = false;
						break;
					}
				}
			}
		}
	}
	dpstates[n-1][0].visible = true;
	dpstates[n-1][0].weight = 0;
	dpstates[n-1][0].bestvertex = -1;

	for(gap = 2; gap<n; gap++) {
		for(i=0; i<(n-gap); i++) {
			j = i+gap;
			if(!dpstates[j][i].visible) continue;
			bestvertex = -1;
			for(k=(i+1);k<j;k++) {
				if(!dpstates[k][i].visible) continue;
				if(!dpstates[j][k].visible) continue;

				if(k<=(i+1)) d1=0;
				else d1 = Distance(poly->GetPoint(i),poly->GetPoint(k));
				if(j<=(k+1)) d2=0;
				else d2 = Distance(poly->GetPoint(k),poly->GetPoint(j));

				weight = dpstates[k][i].weight + dpstates[j][k].weight + d1 + d2;

				if((bestvertex == -1)||(weight<minweight)) {
					bestvertex = k;
					minweight = weight;
				}
			}
			if(bestvertex == -1) {
				for(i=1;i<n;i++) {
					delete [] dpstates[i];
				}
				delete [] dpstates;

				return 0;
			}

			dpstates[j][i].bestvertex = bestvertex;
			dpstates[j][i].weight = minweight;
		}
	}

	newdiagonal.index1 = 0;
	newdiagonal.index2 = n-1;
	diagonals.push_back(newdiagonal);
	while(!diagonals.empty()) {
		diagonal = *(diagonals.begin());
		diagonals.pop_front();
		bestvertex = dpstates[diagonal.index2][diagonal.index1].bestvertex;
		if(bestvertex == -1) {
			ret = 0;
			break;
		}
		triangle.Triangle(poly->GetPoint(diagonal.index1),poly->GetPoint(bestvertex),poly->GetPoint(diagonal.index2));
		triangles->push_back(triangle);
		if(bestvertex > (diagonal.index1+1)) {
			newdiagonal.index1 = diagonal.index1;
			newdiagonal.index2 = bestvertex;
			diagonals.push_back(newdiagonal);
		}
		if(diagonal.index2 > (bestvertex+1)) {
			newdiagonal.index1 = bestvertex;
			newdiagonal.index2 = diagonal.index2;
			diagonals.push_back(newdiagonal);
		}
	}

	for(i=1;i<n;i++) {
		delete [] dpstates[i];
	}
	delete [] dpstates;

	return ret;
}

void TPPLPartition::UpdateState(long a, long b, long w, long i, long j, DPState2 **dpstates) {
	Diagonal newdiagonal;
	DiagonalList *pairs = NULL;
	long w2;

	w2 = dpstates[a][b].weight;
	if(w>w2) return;

	pairs = &(dpstates[a][b].pairs);
	newdiagonal.index1 = i;
	newdiagonal.index2 = j;

	if(w<w2) {
		pairs->clear();
		pairs->push_front(newdiagonal);
		dpstates[a][b].weight = w;
	} else {
		if((!pairs->empty())&&(i <= pairs->begin()->index1)) return;
		while((!pairs->empty())&&(pairs->begin()->index2 >= j)) pairs->pop_front();
		pairs->push_front(newdiagonal);
	}
}

void TPPLPartition::TypeA(long i, long j, long k, PartitionVertex *vertices, DPState2 **dpstates) {
	DiagonalList *pairs = NULL;
	DiagonalList::iterator iter,lastiter;
	long top;
	long w;

	if(!dpstates[i][j].visible) return;
	top = j;
	w = dpstates[i][j].weight;
    if(k-j > 1) {
		if (!dpstates[j][k].visible) return;
		w += dpstates[j][k].weight + 1;
    }
	if(j-i > 1) {
		pairs = &(dpstates[i][j].pairs);
		iter = pairs->end();
		lastiter = pairs->end();
		while(iter!=pairs->begin()) {
			iter--;
			if(!IsReflex(vertices[iter->index2].p,vertices[j].p,vertices[k].p)) lastiter = iter;
			else break;
		}
		if(lastiter == pairs->end()) w++;
		else {
			if(IsReflex(vertices[k].p,vertices[i].p,vertices[lastiter->index1].p)) w++;
			else top = lastiter->index1;
		}
	}
	UpdateState(i,k,w,top,j,dpstates);
}

void TPPLPartition::TypeB(long i, long j, long k, PartitionVertex *vertices, DPState2 **dpstates) {
	DiagonalList *pairs = NULL;
	DiagonalList::iterator iter,lastiter;
	long top;
	long w;

	if(!dpstates[j][k].visible) return;
	top = j;
	w = dpstates[j][k].weight;

	if (j-i > 1) {
		if (!dpstates[i][j].visible) return;
		w += dpstates[i][j].weight + 1;
	}
	if (k-j > 1) {
		pairs = &(dpstates[j][k].pairs);

		iter = pairs->begin();
		if((!pairs->empty())&&(!IsReflex(vertices[i].p,vertices[j].p,vertices[iter->index1].p))) {
			lastiter = iter;
			while(iter!=pairs->end()) {
				if(!IsReflex(vertices[i].p,vertices[j].p,vertices[iter->index1].p)) {
					lastiter = iter;
					iter++;
				}
				else break;
			}
			if(IsReflex(vertices[lastiter->index2].p,vertices[k].p,vertices[i].p)) w++;
			else top = lastiter->index2;
		} else w++;
	}
	UpdateState(i,k,w,j,top,dpstates);
}

int TPPLPartition::ConvexPartition_OPT(TPPLPoly *poly, TPPLPolyList *parts) {
	if(!poly->Valid()) return 0;

	TPPLPoint p1,p2,p3,p4;
	PartitionVertex *vertices = NULL;
	DPState2 **dpstates = NULL;
	long i,j,k,n,gap;
	DiagonalList diagonals,diagonals2;
	Diagonal diagonal,newdiagonal;
	DiagonalList *pairs = NULL,*pairs2 = NULL;
	DiagonalList::iterator iter,iter2;
	int ret;
	TPPLPoly newpoly;
	vector<long> indices;
	vector<long>::iterator iiter;
	bool ijreal,jkreal;

	n = poly->GetNumPoints();
	vertices = new PartitionVertex[n];

	dpstates = new DPState2 *[n];
	for(i=0;i<n;i++) {
		dpstates[i] = new DPState2[n];
	}

	//init vertex information
	for(i=0;i<n;i++) {
		vertices[i].p = poly->GetPoint(i);
		vertices[i].isActive = true;
		if(i==0) vertices[i].previous = &(vertices[n-1]);
		else vertices[i].previous = &(vertices[i-1]);
		if(i==(poly->GetNumPoints()-1)) vertices[i].next = &(vertices[0]);
		else vertices[i].next = &(vertices[i+1]);
	}
	for(i=1;i<n;i++) {
		UpdateVertexReflexity(&(vertices[i]));
	}

	//init states and visibility
	for(i=0;i<(n-1);i++) {
		p1 = poly->GetPoint(i);
		for(j=i+1;j<n;j++) {
			dpstates[i][j].visible = true;
			if(j==i+1) {
				dpstates[i][j].weight = 0;
			} else {
				dpstates[i][j].weight = 2147483647;
			}
			if(j!=(i+1)) {
				p2 = poly->GetPoint(j);

				//visibility check
				if(!InCone(&vertices[i],p2)) {
					dpstates[i][j].visible = false;
					continue;
				}
				if(!InCone(&vertices[j],p1)) {
					dpstates[i][j].visible = false;
					continue;
				}

				for(k=0;k<n;k++) {
					p3 = poly->GetPoint(k);
					if(k==(n-1)) p4 = poly->GetPoint(0);
					else p4 = poly->GetPoint(k+1);
					if(Intersects(p1,p2,p3,p4)) {
						dpstates[i][j].visible = false;
						break;
					}
				}
			}
		}
	}
	for(i=0;i<(n-2);i++) {
		j = i+2;
		if(dpstates[i][j].visible) {
			dpstates[i][j].weight = 0;
			newdiagonal.index1 = i+1;
			newdiagonal.index2 = i+1;
			dpstates[i][j].pairs.push_back(newdiagonal);
		}
	}

	dpstates[0][n-1].visible = true;
	vertices[0].isConvex = false; //by convention

	for(gap=3; gap<n; gap++) {
		for(i=0;i<n-gap;i++) {
			if(vertices[i].isConvex) continue;
			k = i+gap;
			if(dpstates[i][k].visible) {
				if(!vertices[k].isConvex) {
					for(j=i+1;j<k;j++) TypeA(i,j,k,vertices,dpstates);
				} else {
					for(j=i+1;j<(k-1);j++) {
						if(vertices[j].isConvex) continue;
						TypeA(i,j,k,vertices,dpstates);
					}
					TypeA(i,k-1,k,vertices,dpstates);
				}
			}
		}
		for(k=gap;k<n;k++) {
			if(vertices[k].isConvex) continue;
			i = k-gap;
			if((vertices[i].isConvex)&&(dpstates[i][k].visible)) {
				TypeB(i,i+1,k,vertices,dpstates);
				for(j=i+2;j<k;j++) {
					if(vertices[j].isConvex) continue;
					TypeB(i,j,k,vertices,dpstates);
				}
			}
		}
	}


	//recover solution
	ret = 1;
	newdiagonal.index1 = 0;
	newdiagonal.index2 = n-1;
	diagonals.push_front(newdiagonal);
	while(!diagonals.empty()) {
		diagonal = *(diagonals.begin());
		diagonals.pop_front();
		if((diagonal.index2 - diagonal.index1) <=1) continue;
		pairs = &(dpstates[diagonal.index1][diagonal.index2].pairs);
		if(pairs->empty()) {
			ret = 0;
			break;
		}
		if(!vertices[diagonal.index1].isConvex) {
			iter = pairs->end();
			iter--;
			j = iter->index2;
			newdiagonal.index1 = j;
			newdiagonal.index2 = diagonal.index2;
			diagonals.push_front(newdiagonal);
			if((j - diagonal.index1)>1) {
				if(iter->index1 != iter->index2) {
					pairs2 = &(dpstates[diagonal.index1][j].pairs);
					while(1) {
						if(pairs2->empty()) {
							ret = 0;
							break;
						}
						iter2 = pairs2->end();
						iter2--;
						if(iter->index1 != iter2->index1) pairs2->pop_back();
						else break;
					}
					if(ret == 0) break;
				}
				newdiagonal.index1 = diagonal.index1;
				newdiagonal.index2 = j;
				diagonals.push_front(newdiagonal);
			}
		} else {
			iter = pairs->begin();
			j = iter->index1;
			newdiagonal.index1 = diagonal.index1;
			newdiagonal.index2 = j;
			diagonals.push_front(newdiagonal);
			if((diagonal.index2 - j) > 1) {
				if(iter->index1 != iter->index2) {
					pairs2 = &(dpstates[j][diagonal.index2].pairs);
					while(1) {
						if(pairs2->empty()) {
							ret = 0;
							break;
						}
						iter2 = pairs2->begin();
						if(iter->index2 != iter2->index2) pairs2->pop_front();
						else break;
					}
					if(ret == 0) break;
				}
				newdiagonal.index1 = j;
				newdiagonal.index2 = diagonal.index2;
				diagonals.push_front(newdiagonal);
			}
		}
	}

	if(ret == 0) {
		for(i=0;i<n;i++) {
			delete [] dpstates[i];
		}
		delete [] dpstates;
		delete [] vertices;

		return ret;
	}

	newdiagonal.index1 = 0;
	newdiagonal.index2 = n-1;
	diagonals.push_front(newdiagonal);
	while(!diagonals.empty()) {
		diagonal = *(diagonals.begin());
		diagonals.pop_front();
		if((diagonal.index2 - diagonal.index1) <= 1) continue;

		indices.clear();
		diagonals2.clear();
		indices.push_back(diagonal.index1);
		indices.push_back(diagonal.index2);
		diagonals2.push_front(diagonal);

		while(!diagonals2.empty()) {
			diagonal = *(diagonals2.begin());
			diagonals2.pop_front();
			if((diagonal.index2 - diagonal.index1) <= 1) continue;
			ijreal = true;
			jkreal = true;
			pairs = &(dpstates[diagonal.index1][diagonal.index2].pairs);
			if(!vertices[diagonal.index1].isConvex) {
				iter = pairs->end();
				iter--;
				j = iter->index2;
				if(iter->index1 != iter->index2) ijreal = false;
			} else {
				iter = pairs->begin();
				j = iter->index1;
				if(iter->index1 != iter->index2) jkreal = false;
			}

			newdiagonal.index1 = diagonal.index1;
			newdiagonal.index2 = j;
			if(ijreal) {
				diagonals.push_back(newdiagonal);
			} else {
				diagonals2.push_back(newdiagonal);
			}

			newdiagonal.index1 = j;
			newdiagonal.index2 = diagonal.index2;
			if(jkreal) {
				diagonals.push_back(newdiagonal);
			} else {
				diagonals2.push_back(newdiagonal);
			}

			indices.push_back(j);
		}

		std::sort(indices.begin(), indices.end());
		newpoly.Init((long)indices.size());
		k=0;
		for(iiter = indices.begin();iiter!=indices.end();iiter++) {
			newpoly[k] = vertices[*iiter].p;
			k++;
		}
		parts->push_back(newpoly);
	}

	for(i=0;i<n;i++) {
		delete [] dpstates[i];
	}
	delete [] dpstates;
	delete [] vertices;

	return ret;
}

//triangulates a set of polygons by first partitioning them into monotone polygons
//O(n*log(n)) time complexity, O(n) space complexity
//the algorithm used here is outlined in the book
//"Computational Geometry: Algorithms and Applications"
//by Mark de Berg, Otfried Cheong, Marc van Kreveld and Mark Overmars
int TPPLPartition::MonotonePartition(TPPLPolyList *inpolys, TPPLPolyList *monotonePolys) {
	TPPLPolyList::iterator iter;
	MonotoneVertex *vertices = NULL;
	long i,numvertices,vindex,vindex2,newnumvertices,maxnumvertices;
	long polystartindex, polyendindex;
	TPPLPoly *poly = NULL;
	MonotoneVertex *v = NULL,*v2 = NULL,*vprev = NULL,*vnext = NULL;
	ScanLineEdge newedge;
	bool error = false;

	numvertices = 0;
	for(iter = inpolys->begin(); iter != inpolys->end(); iter++) {
		if(!iter->Valid()) return 0;
		numvertices += iter->GetNumPoints();
	}

	maxnumvertices = numvertices*3;
	vertices = new MonotoneVertex[maxnumvertices];
	newnumvertices = numvertices;

	polystartindex = 0;
	for(iter = inpolys->begin(); iter != inpolys->end(); iter++) {
		poly = &(*iter);
		polyendindex = polystartindex + poly->GetNumPoints()-1;
		for(i=0;i<poly->GetNumPoints();i++) {
			vertices[i+polystartindex].p = poly->GetPoint(i);
			if(i==0) vertices[i+polystartindex].previous = polyendindex;
			else vertices[i+polystartindex].previous = i+polystartindex-1;
			if(i==(poly->GetNumPoints()-1)) vertices[i+polystartindex].next = polystartindex;
			else vertices[i+polystartindex].next = i+polystartindex+1;
		}
		polystartindex = polyendindex+1;
	}

	//construct the priority queue
	long *priority = new long [numvertices];
	for(i=0;i<numvertices;i++) priority[i] = i;
	std::sort(priority,&(priority[numvertices]),VertexSorter(vertices));

	//determine vertex types
	char *vertextypes = new char[maxnumvertices];
	for(i=0;i<numvertices;i++) {
		v = &(vertices[i]);
		vprev = &(vertices[v->previous]);
		vnext = &(vertices[v->next]);

		if(Below(vprev->p,v->p)&&Below(vnext->p,v->p)) {
			if(IsConvex(vnext->p,vprev->p,v->p)) {
				vertextypes[i] = TPPL_VERTEXTYPE_START;
			} else {
				vertextypes[i] = TPPL_VERTEXTYPE_SPLIT;
			}
		} else if(Below(v->p,vprev->p)&&Below(v->p,vnext->p)) {
			if(IsConvex(vnext->p,vprev->p,v->p))
			{
				vertextypes[i] = TPPL_VERTEXTYPE_END;
			} else {
				vertextypes[i] = TPPL_VERTEXTYPE_MERGE;
			}
		} else {
			vertextypes[i] = TPPL_VERTEXTYPE_REGULAR;
		}
	}

	//helpers
	long *helpers = new long[maxnumvertices];

	//binary search tree that holds edges intersecting the scanline
	//note that while set doesn't actually have to be implemented as a tree
	//complexity requirements for operations are the same as for the balanced binary search tree
	set<ScanLineEdge> edgeTree;
	//store iterators to the edge tree elements
	//this makes deleting existing edges much faster
	set<ScanLineEdge>::iterator *edgeTreeIterators,edgeIter;
	edgeTreeIterators = new set<ScanLineEdge>::iterator[maxnumvertices];
	pair<set<ScanLineEdge>::iterator,bool> edgeTreeRet;
	for(i = 0; i<numvertices; i++) edgeTreeIterators[i] = edgeTree.end();

	//for each vertex
	for(i=0;i<numvertices;i++) {
		vindex = priority[i];
		v = &(vertices[vindex]);
		vindex2 = vindex;
		v2 = v;

		//depending on the vertex type, do the appropriate action
		//comments in the following sections are copied from "Computational Geometry: Algorithms and Applications"
		switch(vertextypes[vindex]) {
			case TPPL_VERTEXTYPE_START:
				//Insert ei in T and set helper(ei) to vi.
				newedge.p1 = v->p;
				newedge.p2 = vertices[v->next].p;
				newedge.index = vindex;
				edgeTreeRet = edgeTree.insert(newedge);
				edgeTreeIterators[vindex] = edgeTreeRet.first;
				helpers[vindex] = vindex;
				break;

			case TPPL_VERTEXTYPE_END:
				if (edgeTreeIterators[v->previous] == edgeTree.end()) {
					error = true;
					break;
				}
				//if helper(ei-1) is a merge vertex
				if(vertextypes[helpers[v->previous]]==TPPL_VERTEXTYPE_MERGE) {
					//Insert the diagonal connecting vi to helper(ei-1) in D.
					AddDiagonal(vertices,&newnumvertices,vindex,helpers[v->previous],
						vertextypes, edgeTreeIterators, &edgeTree, helpers);
				}
				//Delete ei-1 from T
				edgeTree.erase(edgeTreeIterators[v->previous]);
				break;

			case TPPL_VERTEXTYPE_SPLIT:
				//Search in T to find the edge e j directly left of vi.
				newedge.p1 = v->p;
				newedge.p2 = v->p;
				edgeIter = edgeTree.lower_bound(newedge);
				if(edgeIter == edgeTree.begin()) {
					error = true;
					break;
				}
				edgeIter--;
				//Insert the diagonal connecting vi to helper(ej) in D.
				AddDiagonal(vertices,&newnumvertices,vindex,helpers[edgeIter->index],
					vertextypes, edgeTreeIterators, &edgeTree, helpers);
				vindex2 = newnumvertices-2;
				v2 = &(vertices[vindex2]);
				//helper(e j)�vi
				helpers[edgeIter->index] = vindex;
				//Insert ei in T and set helper(ei) to vi.
				newedge.p1 = v2->p;
				newedge.p2 = vertices[v2->next].p;
				newedge.index = vindex2;
				edgeTreeRet = edgeTree.insert(newedge);
				edgeTreeIterators[vindex2] = edgeTreeRet.first;
				helpers[vindex2] = vindex2;
				break;

			case TPPL_VERTEXTYPE_MERGE:
				if (edgeTreeIterators[v->previous] == edgeTree.end()) {
					error = true;
					break;
				}
				//if helper(ei-1) is a merge vertex
				if(vertextypes[helpers[v->previous]]==TPPL_VERTEXTYPE_MERGE) {
					//Insert the diagonal connecting vi to helper(ei-1) in D.
					AddDiagonal(vertices,&newnumvertices,vindex,helpers[v->previous],
						vertextypes, edgeTreeIterators, &edgeTree, helpers);
					vindex2 = newnumvertices-2;
					v2 = &(vertices[vindex2]);
				}
				//Delete ei-1 from T.
				edgeTree.erase(edgeTreeIterators[v->previous]);
				//Search in T to find the edge e j directly left of vi.
				newedge.p1 = v->p;
				newedge.p2 = v->p;
				edgeIter = edgeTree.lower_bound(newedge);
				if(edgeIter == edgeTree.begin()) {
					error = true;
					break;
				}
				edgeIter--;
				//if helper(ej) is a merge vertex
				if(vertextypes[helpers[edgeIter->index]]==TPPL_VERTEXTYPE_MERGE) {
					//Insert the diagonal connecting vi to helper(e j) in D.
					AddDiagonal(vertices,&newnumvertices,vindex2,helpers[edgeIter->index],
						vertextypes, edgeTreeIterators, &edgeTree, helpers);
				}
				//helper(e j)�vi
				helpers[edgeIter->index] = vindex2;
				break;

			case TPPL_VERTEXTYPE_REGULAR:
				//if the interior of P lies to the right of vi
				if(Below(v->p,vertices[v->previous].p)) {
					if (edgeTreeIterators[v->previous] == edgeTree.end()) {
						error = true;
						break;
					}
					//if helper(ei-1) is a merge vertex
					if(vertextypes[helpers[v->previous]]==TPPL_VERTEXTYPE_MERGE) {
						//Insert the diagonal connecting vi to helper(ei-1) in D.
						AddDiagonal(vertices,&newnumvertices,vindex,helpers[v->previous],
							vertextypes, edgeTreeIterators, &edgeTree, helpers);
						vindex2 = newnumvertices-2;
						v2 = &(vertices[vindex2]);
					}
					//Delete ei-1 from T.
					edgeTree.erase(edgeTreeIterators[v->previous]);
					//Insert ei in T and set helper(ei) to vi.
					newedge.p1 = v2->p;
					newedge.p2 = vertices[v2->next].p;
					newedge.index = vindex2;
					edgeTreeRet = edgeTree.insert(newedge);
					edgeTreeIterators[vindex2] = edgeTreeRet.first;
					helpers[vindex2] = vindex;
				} else {
					//Search in T to find the edge ej directly left of vi.
					newedge.p1 = v->p;
					newedge.p2 = v->p;
					edgeIter = edgeTree.lower_bound(newedge);
					if(edgeIter == edgeTree.begin()) {
						error = true;
						break;
					}
					edgeIter--;
					//if helper(ej) is a merge vertex
					if(vertextypes[helpers[edgeIter->index]]==TPPL_VERTEXTYPE_MERGE) {
						//Insert the diagonal connecting vi to helper(e j) in D.
						AddDiagonal(vertices,&newnumvertices,vindex,helpers[edgeIter->index],
							vertextypes, edgeTreeIterators, &edgeTree, helpers);
					}
					//helper(e j)�vi
					helpers[edgeIter->index] = vindex;
				}
				break;
		}

		if(error) break;
	}

	char *used = new char[newnumvertices];
	memset(used,0,newnumvertices*sizeof(char));

	if(!error) {
		//return result
		long size;
		TPPLPoly mpoly;
		for(i=0;i<newnumvertices;i++) {
			if(used[i]) continue;
			v = &(vertices[i]);
			vnext = &(vertices[v->next]);
			size = 1;
			while(vnext!=v) {
				vnext = &(vertices[vnext->next]);
				size++;
			}
			mpoly.Init(size);
			v = &(vertices[i]);
			mpoly[0] = v->p;
			vnext = &(vertices[v->next]);
			size = 1;
			used[i] = 1;
			used[v->next] = 1;
			while(vnext!=v) {
				mpoly[size] = vnext->p;
				used[vnext->next] = 1;
				vnext = &(vertices[vnext->next]);
				size++;
			}
			monotonePolys->push_back(mpoly);
		}
	}

	//cleanup
	delete [] vertices;
	delete [] priority;
	delete [] vertextypes;
	delete [] edgeTreeIterators;
	delete [] helpers;
	delete [] used;

	if(error) {
		return 0;
	} else {
		return 1;
	}
}

//adds a diagonal to the doubly-connected list of vertices
void TPPLPartition::AddDiagonal(MonotoneVertex *vertices, long *numvertices, long index1, long index2,
								char *vertextypes, set<ScanLineEdge>::iterator *edgeTreeIterators,
								set<ScanLineEdge> *edgeTree, long *helpers)
{
	long newindex1,newindex2;

	newindex1 = *numvertices;
	(*numvertices)++;
	newindex2 = *numvertices;
	(*numvertices)++;

	vertices[newindex1].p = vertices[index1].p;
	vertices[newindex2].p = vertices[index2].p;

	vertices[newindex2].next = vertices[index2].next;
	vertices[newindex1].next = vertices[index1].next;

	vertices[vertices[index2].next].previous = newindex2;
	vertices[vertices[index1].next].previous = newindex1;

	vertices[index1].next = newindex2;
	vertices[newindex2].previous = index1;

	vertices[index2].next = newindex1;
	vertices[newindex1].previous = index2;

	//update all relevant structures
	vertextypes[newindex1] = vertextypes[index1];
	edgeTreeIterators[newindex1] = edgeTreeIterators[index1];
	helpers[newindex1] = helpers[index1];
	if(edgeTreeIterators[newindex1] != edgeTree->end())
		edgeTreeIterators[newindex1]->index = newindex1;
	vertextypes[newindex2] = vertextypes[index2];
	edgeTreeIterators[newindex2] = edgeTreeIterators[index2];
	helpers[newindex2] = helpers[index2];
	if(edgeTreeIterators[newindex2] != edgeTree->end())
		edgeTreeIterators[newindex2]->index = newindex2;
}

bool TPPLPartition::Below(TPPLPoint &p1, TPPLPoint &p2) {
	if(p1.y < p2.y) return true;
	else if(p1.y == p2.y) {
		if(p1.x < p2.x) return true;
	}
	return false;
}

//sorts in the falling order of y values, if y is equal, x is used instead
bool TPPLPartition::VertexSorter::operator() (long index1, long index2) {
	if(vertices[index1].p.y > vertices[index2].p.y) return true;
	else if(vertices[index1].p.y == vertices[index2].p.y) {
		if(vertices[index1].p.x > vertices[index2].p.x) return true;
	}
	return false;
}

bool TPPLPartition::ScanLineEdge::IsConvex(const TPPLPoint& p1, const TPPLPoint& p2, const TPPLPoint& p3) const {
	tppl_float tmp;
	tmp = (p3.y-p1.y)*(p2.x-p1.x)-(p3.x-p1.x)*(p2.y-p1.y);
	if(tmp>0) return 1;
	else return 0;
}

bool TPPLPartition::ScanLineEdge::operator < (const ScanLineEdge & other) const {
	if(other.p1.y == other.p2.y) {
		if(p1.y == p2.y) {
			if(p1.y < other.p1.y) return true;
			else return false;
		}
		if(IsConvex(p1,p2,other.p1)) return true;
		else return false;
	} else if(p1.y == p2.y) {
		if(IsConvex(other.p1,other.p2,p1)) return false;
		else return true;
	} else if(p1.y < other.p1.y) {
		if(IsConvex(other.p1,other.p2,p1)) return false;
		else return true;
	} else {
		if(IsConvex(p1,p2,other.p1)) return true;
		else return false;
	}
}

//triangulates monotone polygon
//O(n) time, O(n) space complexity
int TPPLPartition::TriangulateMonotone(TPPLPoly *inPoly, TPPLPolyList *triangles) {
	if(!inPoly->Valid()) return 0;

	long i,i2,j,topindex,bottomindex,leftindex,rightindex,vindex;
	TPPLPoint *points = NULL;
	long numpoints;
	TPPLPoly triangle;

	numpoints = inPoly->GetNumPoints();
	points = inPoly->GetPoints();

	//trivial case
	if(numpoints == 3) {
		triangles->push_back(*inPoly);
		return 1;
	}

	topindex = 0; bottomindex=0;
	for(i=1;i<numpoints;i++) {
		if(Below(points[i],points[bottomindex])) bottomindex = i;
		if(Below(points[topindex],points[i])) topindex = i;
	}

	//check if the poly is really monotone
	i = topindex;
	while(i!=bottomindex) {
		i2 = i+1; if(i2>=numpoints) i2 = 0;
		if(!Below(points[i2],points[i])) return 0;
		i = i2;
	}
	i = bottomindex;
	while(i!=topindex) {
		i2 = i+1; if(i2>=numpoints) i2 = 0;
		if(!Below(points[i],points[i2])) return 0;
		i = i2;
	}

	char *vertextypes = new char[numpoints];
	long *priority = new long[numpoints];

	//merge left and right vertex chains
	priority[0] = topindex;
	vertextypes[topindex] = 0;
	leftindex = topindex+1; if(leftindex>=numpoints) leftindex = 0;
	rightindex = topindex-1; if(rightindex<0) rightindex = numpoints-1;
	for(i=1;i<(numpoints-1);i++) {
		if(leftindex==bottomindex) {
			priority[i] = rightindex;
			rightindex--; if(rightindex<0) rightindex = numpoints-1;
			vertextypes[priority[i]] = -1;
		} else if(rightindex==bottomindex) {
			priority[i] = leftindex;
			leftindex++;  if(leftindex>=numpoints) leftindex = 0;
			vertextypes[priority[i]] = 1;
		} else {
			if(Below(points[leftindex],points[rightindex])) {
				priority[i] = rightindex;
				rightindex--; if(rightindex<0) rightindex = numpoints-1;
				vertextypes[priority[i]] = -1;
			} else {
				priority[i] = leftindex;
				leftindex++;  if(leftindex>=numpoints) leftindex = 0;
				vertextypes[priority[i]] = 1;
			}
		}
	}
	priority[i] = bottomindex;
	vertextypes[bottomindex] = 0;

	long *stack = new long[numpoints];
	long stackptr = 0;

	stack[0] = priority[0];
	stack[1] = priority[1];
	stackptr = 2;

	//for each vertex from top to bottom trim as many triangles as possible
	for(i=2;i<(numpoints-1);i++) {
		vindex = priority[i];
		if(vertextypes[vindex]!=vertextypes[stack[stackptr-1]]) {
			for(j=0;j<(stackptr-1);j++) {
				if(vertextypes[vindex]==1) {
					triangle.Triangle(points[stack[j+1]],points[stack[j]],points[vindex]);
				} else {
					triangle.Triangle(points[stack[j]],points[stack[j+1]],points[vindex]);
				}
				triangles->push_back(triangle);
			}
			stack[0] = priority[i-1];
			stack[1] = priority[i];
			stackptr = 2;
		} else {
			stackptr--;
			while(stackptr>0) {
				if(vertextypes[vindex]==1) {
					if(IsConvex(points[vindex],points[stack[stackptr-1]],points[stack[stackptr]])) {
						triangle.Triangle(points[vindex],points[stack[stackptr-1]],points[stack[stackptr]]);
						triangles->push_back(triangle);
						stackptr--;
					} else {
						break;
					}
				} else {
					if(IsConvex(points[vindex],points[stack[stackptr]],points[stack[stackptr-1]])) {
						triangle.Triangle(points[vindex],points[stack[stackptr]],points[stack[stackptr-1]]);
						triangles->push_back(triangle);
						stackptr--;
					} else {
						break;
					}
				}
			}
			stackptr++;
			stack[stackptr] = vindex;
			stackptr++;
		}
	}
	vindex = priority[i];
	for(j=0;j<(stackptr-1);j++) {
		if(vertextypes[stack[j+1]]==1) {
			triangle.Triangle(points[stack[j]],points[stack[j+1]],points[vindex]);
		} else {
			triangle.Triangle(points[stack[j+1]],points[stack[j]],points[vindex]);
		}
		triangles->push_back(triangle);
	}

	delete [] priority;
	delete [] vertextypes;
	delete [] stack;

	return 1;
}

int TPPLPartition::Triangulate_MONO(TPPLPolyList *inpolys, TPPLPolyList *triangles) {
	TPPLPolyList monotone;
	TPPLPolyList::iterator iter;

	if(!MonotonePartition(inpolys,&monotone)) return 0;
	for(iter = monotone.begin(); iter!=monotone.end();iter++) {
		if(!TriangulateMonotone(&(*iter),triangles)) return 0;
	}
	return 1;
}

int TPPLPartition::Triangulate_MONO(TPPLPoly *poly, TPPLPolyList *triangles) {
	TPPLPolyList polys;
	polys.push_back(*poly);

	return Triangulate_MONO(&polys, triangles);
}

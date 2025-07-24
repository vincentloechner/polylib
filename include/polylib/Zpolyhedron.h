#ifndef _Zpolyhedron_h_
#define _Zpolyhedron_h_

#if defined(__cplusplus)
extern "C" {
#endif

// DEFINITIONS USED IN THOSE FUNCTIONS:
//
// - a ZPolyhedron is a single lattice associated to a polyhedral domain
//   (a polyhedral domain can be a union of polyhedra)
//
// - a ZDomain is a chained list of ZPolyhedra (with possibly multiple lattices;
//   using the ->next structure element).

extern void CanonicalForm(ZPolyhedron *Zpol, ZPolyhedron **Result,
                          Matrix **Basis);
extern ZPolyhedron *EmptyZPolyhedron(int dimension);
extern ZPolyhedron *IntegraliseLattice(ZPolyhedron *A);
extern Bool isEmptyZPolyhedron(ZPolyhedron *Zpol);
extern ZPolyhedron *ZDomainDifference(ZPolyhedron *A, ZPolyhedron *B);
extern ZPolyhedron *ZDomainImage(ZPolyhedron *A, Matrix *Func);
extern Bool ZDomainIncludes(ZPolyhedron *A, ZPolyhedron *B);
extern ZPolyhedron *ZDomainIntersection(ZPolyhedron *A, ZPolyhedron *B);
extern ZPolyhedron *ZDomainPreimage(ZPolyhedron *A, Matrix *Func);
extern void ZDomainPrint(FILE *fp, const char *format, ZPolyhedron *A);
extern ZPolyhedron *ZDomainSimplify(ZPolyhedron *ZDom);
extern ZPolyhedron *ZDomainUnion(ZPolyhedron *A, ZPolyhedron *B);
extern ZPolyhedron *ZDomain_Copy(ZPolyhedron *Head);
extern void ZDomain_Free(ZPolyhedron *Head);
extern ZPolyhedron *ZPolyhedronAlloc(Lattice *Lat, Polyhedron *Poly);
extern ZPolyhedron *SplitZpolyhedron(ZPolyhedron *ZPol, Lattice *B);
extern void Matrix_Move_Homogeneous_Dim_First(Matrix* A);
extern void Matrix_Move_Homogeneous_Dim_Last(Matrix *A);
extern void Canonical_ZDomain(ZPolyhedron* A);

#if defined(__cplusplus)
}
#endif

#endif /* _Zpolyhedron_h_ */

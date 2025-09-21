#ifndef _LBL_h_
#define _LBL_h_

#if defined(__cplusplus)
extern "C" {
#endif

// DEFINITIONS:
// - a 'single LBL' is a single affine integer function (as the matrix Lat)
//   associated to a polyhedral domain P (a polyhedral domain can be a union
//   of polyhedra). As such, it represents the set of all points:
//      {x = Lat z with z \in P and z \in Z^d}
//   Lat is in homogeneous HNF and P does not contain equalities (rational).
//
// - the name 'LBL' is used to define a union of LBLs, as a chained list of
//   LBLs (with possibly multiple _different_ lattices).
// *************************************************************************
// All user exposed functions manipulate canonical unions of LBLs by default
// *************************************************************************
//
// - a 'Z-polyhedron' is the intersection of an integer lattice and an
//   integer polyhedron. A Z-polyhedron *is* also a specific single LBL and
//   can be represented as such. The difference between an LBL and a
//   Z-polyhedron is that matrix Lat does not contain any column of zeros in
//   a Z-polyhedron (no dimension is eliminated by the Lat function).
//
// - a 'Z-domain' is a union of Z-polyhedra.
//
// All those objects are represented using the same data structure (LBL *),
// so the functions have explicit prefix/suffix names depending on what they
// handle: sLBL, LBL, Zpolyhedron, ZDomain.

extern LBL *LBLAlloc(Matrix *Lat, Polyhedron *Poly);
extern void LBLFree(LBL *Head);
extern void LBLPrint(FILE *fp, const char *format, LBL *A);
extern LBL *LBLCopy(LBL *Head);
extern LBL *EmptyLBL(int dimension);
extern Bool isEmptyLBL(LBL *Zpol);

extern LBL *LBLUnion(LBL *A, LBL *B);
extern Bool LBLIncludes(LBL *A, LBL *B);
extern LBL *LBLIntersection(LBL *A, LBL *B);
extern LBL *LBLDifference(LBL *A, LBL *B);
extern LBL *LBLImage(LBL *A, Matrix *Func);
extern LBL *LBLPreimage(LBL *A, Matrix *Func);
extern LBL *LBLComplement(LBL *A);

extern LBL *LBLSimplify(LBL *A);

// All above functions always return a canonical LBL.
// The CanonicalLBL() function is exposed to the user to enable transforming
// a self-built non-canonical union of LBLs into a canonical LBL:
extern void CanonicalLBL(LBL* A);

// This function transforms a union of LBLs into a union of Z-domains:
extern LBL *LBL2ZDomain(LBL *A);


// removed:
// extern LBL *LBLSimplify(LBL *ZDom);
// extern LBL *SplitLBL(LBL *ZPol, Matrix *B);
// extern LBL *IntegraliseMatrix(LBL *A);

#if defined(__cplusplus)
}
#endif

#endif /* _LBL_h_ */

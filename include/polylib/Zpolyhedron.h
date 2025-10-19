#ifndef _LBL_h_
#define _LBL_h_

#if defined(__cplusplus)
extern "C" {
#endif

// DEFINITIONS:
// - a 'single LBL' is a single affine integer function (as PolyLib matrix Lat)
//   associated to a polyhedral domain P (a union of polyhedra).
//   As such, it represents the set of all points:
//      {x = Lat z | z \in P and z \in Z^d}
//   In its canonical form, Lat is in homogeneous HNF and P does not contain
//   equalities (rational).
//
// - the name 'LBL' is used to describe a union of canonical LBLs, as a linked
//   list of single LBLs with possibly multiple _different_ lattices
//   ************************************************************
//    All user exposed functions manipulate canonical LBL unions
//   ************************************************************
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
// but all the user-exposed functions handle LBL unions.

extern LBL *LBLAlloc(Matrix *Lat, Polyhedron *Domain);
extern void LBLFree(LBL *A);
extern void LBLPrint(FILE *fp, const char *format, LBL *A);
extern LBL *LBLCopy(LBL *A);
extern LBL *EmptyLBL(int dimension);
extern LBL *UniverseLBL(int dimension);
extern Bool isEmptyLBL(LBL *A);

extern LBL *LBLUnion(LBL *A, LBL *B);
extern Bool LBLIncluded(LBL *A, LBL *B);       // True if A \in B
extern LBL *LBLIntersection(LBL *A, LBL *B);
extern LBL *LBLDifference(LBL *A, LBL *B);     // A - B
extern LBL *LBLImage(LBL *A, Matrix *Func);
extern LBL *LBLPreimage(LBL *A, Matrix *Func);
extern LBL *LBLComplement(LBL *A);

extern void LBLSimplify(LBL *A);  // in place

// All above functions always manipulate canonical LBLs.
// The CanonicalLBL() function is exposed to the user to enable transforming
// a self-built non-canonical union of LBLs into a canonical LBL:
extern void CanonicalLBL(LBL* A); // in place

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

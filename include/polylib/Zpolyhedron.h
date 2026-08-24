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

// allocates the single LBL (Lat, Domain).
extern LBL *LBLAlloc(Matrix *Lat, Polyhedron *Domain);

// free memory
extern void LBLFree(LBL *A);

// print out, format depends on the Value type ("%d", "ld", etc.):
extern void LBLPrint(FILE *fp, const char *format, LBL *A);

// copy of A
extern LBL *LBLCopy(LBL *A);

// return the empty/universe LBL of given dimension:
extern LBL *EmptyLBL(int dimension);
extern LBL *UniverseLBL(int dimension);

// check for emptiness (if the LBL is the result of a complex operation, it is
// advised to call LBLSimplifyEmpty(A) before):
extern Bool isEmptyLBL(LBL *A);

// union A U B
extern LBL *LBLUnion(LBL *A, LBL *B);

// True if A \in B
extern Bool LBLIncluded(LBL *A, LBL *B);

// True if pt \in A (pt an integer vector of same dimension as A)
extern Bool LBLContainsPoint(LBL *A, Value *pt);

// intersection A \cap B
extern LBL *LBLIntersection(LBL *A, LBL *B);

// difference A - B
extern LBL *LBLDifference(LBL *A, LBL *B);

// image by affine function
extern LBL *LBLImage(LBL *A, Matrix *Func);

// preimage by affine function
extern LBL *LBLPreimage(LBL *A, Matrix *Func);

// complement = Universe - A
extern LBL *LBLComplement(LBL *A);

// check for empty coordinate polyhedra and remove them (in place)
extern void LBLSimplifyEmpty(LBL *A);

// tries to simplify the union of polyhedra representing A (in place)
extern void LBLSimplify(LBL *A);


// All above functions always manipulate canonical LBLs.
// The CanonicalLBL() function is only exposed to the user to enable
// normalization of a self-built non-canonical union of LBLs (in place)
extern void CanonicalLBL(LBL* A);

// This function transforms a union of LBLs into a union of Z-domains:
extern LBL *LBL2ZDomain(LBL *A);

// This function transforms a union of LBLs into a disjoint union of LBLs:
extern LBL *LBLDisjointUnion(LBL *A);

// removed:
// extern LBL *SplitLBL(LBL *ZPol, Matrix *B);
// extern LBL *IntegraliseMatrix(LBL *A);

#if defined(__cplusplus)
}
#endif

#endif /* _LBL_h_ */

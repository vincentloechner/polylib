#ifndef _LBL_h_
#define _LBL_h_

#if defined(__cplusplus)
extern "C" {
#endif

// DEFINITIONS USED IN THOSE FUNCTIONS:
//
// - a LBL is a single lattice associated to a polyhedral domain
//   (a polyhedral domain can be a union of polyhedra)
//
// - a LBL is a chained list of ZPolyhedra
//   (with possibly multiple lattices; using the ->next structure element)
//   * It is an LBL according to *
//
// A Z-polyhedron is the intersection of an integer lattice and an integer
// polyhedron:
// L = { z ∈ Qd | ∃z′, z = Lz′ + l, Cz + c ≥ 0 }.
// An LBL, or linearly bounded lattice, is the affine integer image (L, l) of an
// integer polyhedron, called the coordinate polyhedron:
// Z = { z = Lx + l | Cx + c ≥ 0, x ∈ Z^d }.


extern LBL *EmptyLBL(int dimension);
extern Bool isEmptyLBL(LBL *Zpol);
extern LBL *LBLDifference(LBL *A, LBL *B);
extern LBL *LBLImage(LBL *A, Matrix *Func);
extern Bool LBLIncludes(LBL *A, LBL *B);
extern LBL *LBLIntersection(LBL *A, LBL *B);
extern LBL *LBLPreimage(LBL *A, Matrix *Func);
extern void LBLPrint(FILE *fp, const char *format, LBL *A);
extern LBL *LBLUnion(LBL *A, LBL *B);
extern LBL *LBLCopy(LBL *Head);
extern void LBLFree(LBL *Head);
extern LBL *LBLAlloc(Lattice *Lat, Polyhedron *Poly);
extern void CanonicalLBL(LBL* A);

// removed:
// extern LBL *LBLSimplify(LBL *ZDom);
// extern LBL *SplitLBL(LBL *ZPol, Lattice *B);
// extern LBL *IntegraliseLattice(LBL *A);

#if defined(__cplusplus)
}
#endif

#endif /* _LBL_h_ */

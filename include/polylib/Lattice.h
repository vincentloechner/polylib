#ifndef _Lattice_h_
#define _Lattice_h_

#if defined(__cplusplus)
extern "C" {
#endif

extern void AffineHermite(Lattice *A, Lattice **H, Matrix **U);
// extern void AffineSmith(Lattice *A, Lattice **U, Lattice **V, Lattice **Diag);
extern int intcompare(const void *a, const void *b);
extern Bool isEmptyLattice(Lattice *A);
extern Bool isLinear(Lattice *A);
extern LatticeUnion *LatticeDifference(Lattice *A, Lattice *B);
extern Bool LatticeIncludes(Lattice *A, Lattice *B);
extern Lattice *LatticeIntersection(Lattice *X, Lattice *Y);
// extern LatticeUnion *LatticeSimplify(LatticeUnion *latlist);
extern LatticeUnion *LatticeUnion_Alloc(void);
extern void LatticeUnion_Free(LatticeUnion *Head);
extern void PrintLatticeUnion(FILE *fp, char *format, LatticeUnion *Head);
extern Bool sameLattice(Lattice *A, Lattice *B);
extern void Matrix_Move_Homogeneous_Dim_First(Matrix *A);
extern void Matrix_Move_Homogeneous_Dim_Last(Matrix *A);
extern Vector* get_pivots(Matrix* A);
#if defined(__cplusplus)
}
#endif

#endif /* _Lattice_h_ */

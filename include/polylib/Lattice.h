#ifndef _Lattice_h_
#define _Lattice_h_

#if defined(__cplusplus)
extern "C" {
#endif

extern void AffineHermite(Matrix *A, Matrix **H, Matrix **U);
// extern void AffineSmith(Matrix *A, Matrix **U, Matrix **V, Matrix **Diag);
// extern LatticeUnion *LatticeSimplify(LatticeUnion *latlist);
extern int intcompare(const void *a, const void *b);
extern Bool isEmptyLattice(Matrix *A);
extern LatticeUnion *LatticeDifference(Matrix *A, Matrix *B);
extern Bool LatticeIncludes(Matrix *A, Matrix *B);
extern Matrix *LatticeIntersection(Matrix *X, Matrix *Y);
extern LatticeUnion *LatticeUnion_Alloc(void);
extern void LatticeUnion_Free(LatticeUnion *Head);
extern void PrintLatticeUnion(FILE *fp, char *format, LatticeUnion *Head);
extern Bool sameLattice(Matrix *A, Matrix *B);
extern void Matrix_Move_Homogeneous_Dim_First(Matrix *A);
extern void Matrix_Move_Homogeneous_Dim_Last(Matrix *A);
extern Vector* get_pivots(Matrix* A);
#if defined(__cplusplus)
}
#endif

#endif /* _Lattice_h_ */

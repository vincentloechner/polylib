#ifndef _Lattice_h_
#define _Lattice_h_

#if defined(__cplusplus)
extern "C" {
#endif

extern void AffineHermite(Matrix *A, Matrix **H, Matrix **U);
extern Bool isEmptyLattice(Matrix *A);
extern LatticeUnion *LatticeUnion_Alloc(void);
extern void LatticeUnion_Free(LatticeUnion *Head);
extern LatticeUnion *LatticeDifference(Matrix *A, Matrix *B);
extern Bool LatticeIncluded(Matrix *A, Matrix *B);    // True if A \in B
extern Matrix *LatticeIntersection(Matrix *X, Matrix *Y);
extern void PrintLatticeUnion(FILE *fp, char *format, LatticeUnion *Head);
extern int LatCountZeroCols(Matrix* M);
extern Bool isEqualLattice(Matrix *A, Matrix *B);     // exact equality
extern Bool isSameLatticeSpace(Matrix *A, Matrix *B); // spread the same points
extern void Matrix_Move_Homogeneous_Dim_First(Matrix *A);
extern void Matrix_Move_Homogeneous_Dim_Last(Matrix *A);
extern Bool isNormalLattice(Matrix *A);

// removed:
// extern void AffineSmith(Matrix *A, Matrix **U, Matrix **V, Matrix **Diag);
// extern LatticeUnion *LatticeSimplify(LatticeUnion *latlist);

#if defined(__cplusplus)
}
#endif

#endif /* _Lattice_h_ */

/*
 *   This file was automatically generated by version 1.7 of cextract.
 *   Manual editing not recommended.
 *
 *   Created: Mon Mar 30 11:33:19 1998
 */
#ifndef _matrix_H_
#define _matrix_H_
#if __STDC__

extern Matrix *Matrix_Alloc(unsigned NbRows, unsigned NbColumns);
extern void Matrix_Free(Matrix *Mat);
extern void Matrix_Print(FILE * Dst,char *Format,Matrix *Mat);
extern void Matrix_Read_Input(Matrix *Mat);
extern Matrix *Matrix_Read(void);
extern void right_hermite(Matrix *A,Matrix **Hp,Matrix **Up,Matrix
			  **Qp);
extern void left_hermite(Matrix *A,Matrix **Hp,Matrix **Qp,Matrix
			 **Up);
extern int MatInverse(Matrix *M,Matrix *MInv);
extern void rat_prodmat(Matrix *S,Matrix *X,Matrix *P);
extern void Matrix_Vector_Product(Matrix *mat,Value *p1,Value *p2);
extern void Vector_Matrix_Product(Value *p1,Matrix *mat,Value *p2);
extern void Matrix_Product(Matrix *mat1,Matrix *mat2,Matrix *mat3);
extern int Matrix_Inverse(Matrix *Mat,Matrix *MatInv);

#else /* __STDC__ */

extern Matrix *Matrix_Alloc(/* unsigned NbRows, unsigned NbColumns */);
extern void Matrix_Free(/* Matrix *Mat */);
extern void Matrix_Print(/* FILE * Dst,char *Format,Matrix *Mat */);
extern void Matrix_Read_Input(/* Matrix *Mat */);
extern Matrix *Matrix_Read(/* void */);
extern void right_hermite(/* Matrix *A,Matrix **Hp,Matrix **Up,Matrix
			  **Qp */);
extern void left_hermite(/* Matrix *A,Matrix **Hp,Matrix **Qp,Matrix
			 **Up */);
extern int MatInverse(/* Matrix *M,Matrix *MInv */);
extern void rat_prodmat(/* Matrix *S,Matrix *X,Matrix *P */);
extern void Matrix_Vector_Product(/* Matrix *mat,Value *p1,Value *p2 */);
extern void Vector_Matrix_Product(/* Value *p1,Matrix *mat,Value *p2 */);
extern void Matrix_Product(/* Matrix *mat1,Matrix *mat2,Matrix *mat3 */);
extern int Matrix_Inverse(/* Matrix *Mat,Matrix *MatInv */);

#endif /* __STDC__ */
#endif /* _matrix_H_ */
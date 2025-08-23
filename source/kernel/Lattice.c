#include <polylib/polylib.h>
#include <stdlib.h>

// #define LATINTER_DEBUG 1
#define LATDIF_DEBUG 1

typedef struct {
  int count;
  int *fac;
} factor;


static factor allfactors(int num);
static LatticeUnion *generate_lattice_union_line(int line_nb, int *pivots_columns,
    Lattice *A, Lattice *Intersection, Lattice *L, LatticeUnion *Result);
static void get_pivots_columns(Matrix* A, int *columns);

/*
 * Print the contents of a list of Lattices 'Head'
 */
void PrintLatticeUnion(FILE *fp, char *format, LatticeUnion *Head) {

  LatticeUnion *temp;

  for (temp = Head; temp != NULL; temp = temp->next)
    Matrix_Print(fp, format, (Matrix *)temp->M);
  return;
} /* PrintLatticeUnion */

/*
 * Free the memory allocated to a list of lattices 'Head'
 */
void LatticeUnion_Free(LatticeUnion *Head) {

  LatticeUnion *temp;

  while (Head != NULL) {
    temp = Head;
    Head = temp->next;
    Matrix_Free(temp->M);
    free(temp);
  }
  return;
} /* LatticeUnion_Free */

/*
 * Allocate a heads for a list of Lattices
 */
LatticeUnion *LatticeUnion_Alloc(void) {

  LatticeUnion *temp;

  temp = (LatticeUnion *)malloc(sizeof(LatticeUnion));
  temp->M = NULL;
  temp->next = NULL;
  return temp;
} /* LatticeUnion_Alloc */

/*
 * Given two Lattices 'A' and 'B', return True if they have the same affine
 * part (the last column) otherwise return 'False'.
 */
Bool sameAffinepart(Lattice *A, Lattice *B) {

  int i;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered SAMEAFFINEPART \n");
  fclose(fp);
#endif

  for (i = 0; i < A->NbRows; i++)
    if (value_ne(A->p[i][A->NbColumns - 1], B->p[i][B->NbColumns - 1]))
      return False;
  return True;
} /* sameAffinepart */

/*
 * Return an empty lattice of dimension 'dimension-1'. An empty lattice is
 * represented as [[0 0 ... 0] .... [0 ... 0][0 0.....0 1]].
 */
Lattice *EmptyLattice(int dimension) {

  Lattice *result;
  int j;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered NULLATTICE \n");
  fclose(fp);
#endif

  return NULL;
  result = Matrix_Alloc(1, dimension);
  for (j = 0; j < dimension; j++)
    value_set_si(result->p[0][j], 0);
  value_set_si(result->p[0][dimension-1], 1);
  return result;
} /* EmptyLattice */

/*
 * Return True if lattice 'A' is empty, otherwise return False.
 * A has been normalized.
 */
Bool isEmptyLattice(Lattice *A) {

  if(A == NULL || A->NbColumns == 0) {
    return True;
  }
  // A is empty if there are only zero's in the first column
  for (int j = 0; j < A->NbRows-1; j++) {
    if (value_notzero_p(A->p[0][j])) {
      return False;
    }
  }
  return True;
  
} /* isEmptyLattice */

/*
 * Given a Lattice 'A', check whether it is linear or not, i.e. whether the
 * affine part is NULL or not. If affine part is empty, it returns True other-
 * wise it returns False.
 */
Bool isLinear(Lattice *A) {

  int i;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ISLINEAR \n");
  fclose(fp);
#endif

  for (i = 0; i < A->NbRows - 1; i++)
    if (value_notzero_p(A->p[i][A->NbColumns - 1])) {
      return False;
    }
  return True;
} /* isLinear */

/*
 * Return the affine Hermite normal form of the affine lattice 'A'. The unique
 * affine Hermite form if a lattice is stored in 'H' and the unimodular matrix
 * corresponding to 'A = H . U' is stored in the matrix 'U' (if not NULL).
 * Algorithm:
 *     1) move the homogeneous dimensions first (on top-left)
 *     2) compute left_hermite
 *     3) move back the homogeneous dimensions (bottom-right)
 * -> works also on non square matrices (lattices having less row than columns)
 */
void AffineHermite(Lattice *A, Lattice **H, Matrix **U) {

  // DEBUG
  // printf("Entering AffineHermite: A= ");
  // Matrix_Print(stdout, P_VALUE_FMT, A);

  // for left hermite to include the constant, move it on top-left:
  Matrix_Move_Homogeneous_Dim_First(A);
  left_hermite(A, H, U, NULL);
  Matrix_Move_Homogeneous_Dim_Last(*H);
  if(U)
    Matrix_Move_Homogeneous_Dim_Last(*U);
  Matrix_Move_Homogeneous_Dim_Last(A); // restore A as it was

  // OLD VERSION, working fine on square matrices, but not on non-square ones...
  // Lattice *temp;
  // Bool flag = True;

  // #ifdef DOMDEBUG
  //   FILE *fp;
  //   fp = fopen("_debug", "a");
  //   fprintf(fp, "\nEntered AFFINEHERMITE \n");
  //   fclose(fp);
  // #endif

  // if (isLinear(A) == False)
  //   temp = Homogenise(A, True);
  // else {
  //   flag = False;
  //   temp = Matrix_Copy(A);
  // }
  // Hermite(temp, H, U);
  // if (flag == True) {
  //   Matrix_Free(temp);
  //   temp = Homogenise(H[0], False);
  //   Matrix_Free(H[0]);
  //   H[0] = Matrix_Copy(temp);
  //   Matrix_Free(temp);
  //   temp = Homogenise(U[0], False);
  //   Matrix_Free(U[0]);
  //   U[0] = Matrix_Copy(temp);
  // }
  // Matrix_Free(temp);


  // DEBUG
  // printf("Exit AffineHermite: H = ");
  // Matrix_Print(stdout, P_VALUE_FMT, *H);
  // printf("                    U = ");
  // Matrix_Print(stdout, P_VALUE_FMT, *U);

  return;
} /* AffineHermite */

// /*
//  * Given a Polylib matrix 'A' that represents an affine function, return the
//  * affine Smith normal form 'Delta' of 'A' and unimodular matrices 'U' and 'V'
//  * such that 'A = U*Delta*V'.
//  * Algorithm:
//  *           (1) Homogenize the Lattice.
//  *           (2) Call Smith
//  *           (3) The Smith Normal Form Delta must be Dehomogenized and also
//  *               corresponding changes must be made to the Unimodular Matrices
//  *               U and V.
//  *           4) Bring Delta into AffineSmith Form.
//  */
// void AffineSmith(Lattice *A, Lattice **U, Lattice **V, Lattice **Diag) {

//   Lattice *temp;
//   Lattice *Uinv;
//   int i, j;
//   Value sum, quo, rem;

// #ifdef DOMDEBUG
//   FILE *fp;
//   fp = fopen("_debug", "a");
//   fprintf(fp, "\nEntered AFFINESMITH \n");
//   fclose(fp);
// #endif

//   value_init(sum);
//   value_init(quo);
//   value_init(rem);
//   temp = Homogenise(A, True);
//   Smith(temp, U, V, Diag);
//   Matrix_Free(temp);

//   temp = Homogenise(*U, False);
//   Matrix_Free(*U);
//   *U = temp;

//   temp = Homogenise(*V, False);
//   Matrix_Free(*V);
//   *V = temp;

//   temp = Homogenise(*Diag, False);
//   Matrix_Free(*Diag);
//   *Diag = temp;

//   temp = Matrix_Copy(*U);
//   Uinv = Matrix_Alloc((*U)->NbColumns, (*U)->NbRows);
//   Matrix_Inverse(temp, Uinv);
//   Matrix_Free(temp);

//   for (i = 0; i < U[0]->NbRows - 1; i++) {
//     value_set_si(sum, 0);
//     for (j = 0; j < U[0]->NbColumns - 1; j++) {
//       value_addmul(sum, Uinv->p[i][j], U[0]->p[j][U[0]->NbColumns - 1]);
//     }
//     value_assign(Diag[0]->p[i][j], sum);
//   }
//   Matrix_Free(Uinv);
//   for (i = 0; i < U[0]->NbRows - 1; i++)
//     value_set_si(U[0]->p[i][U[0]->NbColumns - 1], 0);
//   for (i = 0; i < Diag[0]->NbRows - 1; i++) {
//     value_division(quo, Diag[0]->p[i][Diag[0]->NbColumns - 1],
//                    Diag[0]->p[i][i]);
//     value_modulus(rem, Diag[0]->p[i][Diag[0]->NbColumns - 1], Diag[0]->p[i][i]);

//     /* Apparently the % operator is strange when sign are different */
//     if (value_neg_p(rem)) {
//       value_addto(rem, rem, Diag[0]->p[i][i]);
//       value_decrement(quo, quo);
//     };
//     value_assign(Diag[0]->p[i][Diag[0]->NbColumns - 1], rem);
//     value_assign(V[0]->p[i][V[0]->NbColumns - 1], quo);
//   }
//   value_clear(sum);
//   value_clear(quo);
//   value_clear(rem);
//   return;
// } /* AffineSmith */

/*
 * Given a lattice 'A' and a boolean variable 'Forward', homogenize the lattice
 * if 'Forward' is True, otherwise if 'Forward' is False, de-homogenize the
 * lattice 'A'.
 * Algorithm:
 *            (1) If Forward == True
 *                Put the last row first.
 *                Put the last columns first.
 *            (2) Else
 *                Put the first row last.
 *                Put the first column last.
 *            (3) Return the result.
 */
Lattice *Homogenise(Lattice *A, Bool Forward) {

  Lattice *result;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered HOMOGENISE \n");
  fclose(fp);
#endif

  result = Matrix_Copy(A);
  if (Forward == True) {
    PutColumnFirst(result, A->NbColumns - 1);
    PutRowFirst(result, result->NbRows - 1);
  } else {
    PutColumnLast(result, 0);
    PutRowLast(result, 0);
  }
  return result;
} /* Homogenise */

/*
 * Given two lattices 'A' and 'B', verify if lattice 'A' is included in 'B' or
 * not. If 'A' is included in 'B' the 'A' intersection 'B', will be 'A'. So,
 * compute 'A' intersection 'B' and check if it is the same as 'A'.
 */
Bool LatticeIncludes(Lattice *A, Lattice *B) {

  Lattice *temp, *HA;
  Bool flag = False;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered LATTICE INCLUDES \n");
  fclose(fp);
#endif

  AffineHermite(A, &HA, NULL);
  temp = LatticeIntersection(B, HA);
  if(temp) {
    if(sameLattice(temp, HA))
      flag = True;
    Matrix_Free(temp);
  }

  Matrix_Free(HA);
  return flag;
} /* LatticeIncludes */

/*
 * Given two lattices 'A' and 'B', verify if 'A' and 'B' are the same lattice.
 * Algorithm:
 *           The Affine Hermite form of two full dimensional matrices are
 * unique. So, take the Affine Hermite form of both 'A' and 'B' and compare the
 * matrices. If they are equal, the function returns True, else it returns
 * False.
 */
Bool sameLattice(Lattice *A, Lattice *B) {

  Lattice *HA, *HB;
  int i, j;
  Bool result = True;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered SAME LATTICE \n");
  fclose(fp);
#endif

  if(A->NbRows != B->NbRows || A->NbColumns != B->NbColumns)
    return (False);

  AffineHermite(A, &HA, NULL);
  AffineHermite(B, &HB, NULL);

  for (i = 0; i < A->NbRows; i++)
    for (j = 0; j < A->NbColumns; j++)
      if (value_ne(HA->p[i][j], HB->p[i][j])) {
        result = False;
        break;
      }

  Matrix_Free(HA);
  Matrix_Free(HB);

  return result;
} /* sameLattice */

/*
 * Given a matrix 'A' and an integer 'dimension', do the following:
 * If dimension < A->dimension), output a (dimension * dimension) submatrix of
 * A. Otherwise the output matrix is [A 0][0 ID]. The order if the identity
 * matrix is (dimension - A->dimension). The input matrix is not necessarily
 * a Polylib matrix but the output is a polylib matrix.
 */
Lattice *ChangeLatticeDimension(Lattice *A, int dimension) {

  int i, j;
  Lattice *Result;

  Result = Matrix_Alloc(dimension, dimension);
  if (dimension <= A->NbRows) {
    for (i = 0; i < dimension; i++)
      for (j = 0; j < dimension; j++)
        value_assign(Result->p[i][j], A->p[i][j]);
    return Result;
  }
  for (i = 0; i < A->NbRows; i++)
    for (j = 0; j < A->NbRows; j++)
      value_assign(Result->p[i][j], A->p[i][j]);

  for (i = A->NbRows; i < dimension; i++)
    for (j = 0; j < dimension; j++) {
      value_set_si(Result->p[i][j], 0);
      value_set_si(Result->p[j][i], 0);
    }
  for (i = A->NbRows; i < dimension; i++)
    value_set_si(Result->p[i][i], 1);
  return Result;
} /* ChangeLatticeDimension */

/*
 * Given an affine lattice 'A', return a matrix of the linear part of the
 * lattice.
 */
Lattice *ExtractLinearPart(Lattice *A) {

  Lattice *Result;
  int i, j;
  Result = (Lattice *)Matrix_Alloc(A->NbRows - 1, A->NbColumns - 1);
  for (i = 0; i < A->NbRows - 1; i++)
    for (j = 0; j < A->NbColumns - 1; j++)
      value_assign(Result->p[i][j], A->p[i][j]);
  return Result;
} /* ExtractLinearPart */

// static Matrix *MakeDioEqforInter(Matrix *A, Matrix *B);

// /*
//  * Given two lattices 'A' and 'B', return the intersection of the two lattcies.
//  * The dimension of 'A' and 'B' should be the same.
//  * Algorithm:
//  *           (1) Verify if the lattcies 'A' and 'B' have the same affine part.
//  *               If they have same affine part, then only their Linear parts
//  *               need to be intersected. If they don't have the same affine
//  *               part then the affine part has to be taken into consideration.
//  *               For this, homogenise the lattices to get their Hermite Forms
//  *               and then find their intersection.
//  *
//  *           (2) Step(2) involves, solving the Diophantine Equations in order
//  *               to extract the intersection of the Lattices. The Diophantine
//  *               equations are formed taking into consideration whether the
//  *               affine part has to be included or not.
//  *
//  *           (3) Solve the Diophantine equations.
//  *
//  *           (4) Extract the necessary information from the result.
//  *
//  *           (5) If the lattices have different affine parts and they were
//  *               homogenised, the result is dehomogenised.
//  */
// Lattice *OldLatticeIntersection(Lattice *X, Lattice *Y) {

//   int i, j, exist;
//   Lattice *result = NULL, *U = NULL;
//   Lattice *A = NULL, *B = NULL, *H = NULL;
//   Matrix *fordio;
//   Vector *X1 = NULL;

// #ifdef DOMDEBUG
//   FILE *fp;
//   fp = fopen("_debug", "a");
//   fprintf(fp, "\nEntered LATTICEINTERSECTION \n");
//   fclose(fp);
// #endif

//   if (X->NbRows != X->NbColumns) {
//     fprintf(stderr, "\nIn LatticeIntersection : The Input Matrix X is a not a "
//                     "well defined Lattice\n");
//     return EmptyLattice(X->NbRows);
//   }

//   if (Y->NbRows != Y->NbColumns) {
//     fprintf(stderr, "\nIn LatticeIntersection : The Input Matrix Y is a not a "
//                     "well defined Lattice\n");
//     return EmptyLattice(X->NbRows);
//   }

//   if (Y->NbRows != X->NbRows) {
//     fprintf(stderr, "\nIn LatticeIntersection : the input lattices X and Y are "
//                     "of incompatible dimensions\n");
//     return EmptyLattice(X->NbRows);
//   }

//   if (isNormalLattice(X))
//     A = (Lattice *)Matrix_Copy(X);
//   else {
//     AffineHermite(X, &H, &U);
//     A = (Lattice *)Matrix_Copy(H);
//     Matrix_Free((Matrix *)H);
//     Matrix_Free((Matrix *)U);
//   }

//   if (isNormalLattice(Y))
//     B = (Lattice *)Matrix_Copy(Y);
//   else {
//     AffineHermite(Y, &H, &U);
//     B = (Lattice *)Matrix_Copy(H);
//     Matrix_Free((Matrix *)H);
//     Matrix_Free((Matrix *)U);
//   }

//   if ((isEmptyLattice(A)) || (isEmptyLattice(B))) {
//     result = EmptyLattice(X->NbRows);
//     Matrix_Free((Matrix *)A);
//     Matrix_Free((Matrix *)B);
//     return result;
//   }
//   fordio = MakeDioEqforInter(A, B);
//   Matrix_Free(A);
//   Matrix_Free(B);
//   exist = SolveDiophantine(fordio, (Matrix **)&U, &X1);
//   if (exist < 0) { /* Intersection is NULL */
//     result = (EmptyLattice(X->NbRows));
//     return result;
//   }

//   result = (Lattice *)Matrix_Alloc(X->NbRows, X->NbColumns);
//   for (i = 0; i < result->NbRows - 1; i++)
//     for (j = 0; j < result->NbColumns - 1; j++)
//       value_assign(result->p[i][j], U->p[i][j]);

//   for (i = 0; i < result->NbRows - 1; i++)
//     value_assign(result->p[i][result->NbColumns - 1], X1->p[i]);
//   for (i = 0; i < result->NbColumns - 1; i++)
//     value_set_si(result->p[result->NbRows - 1][i], 0);
//   value_set_si(result->p[result->NbRows - 1][result->NbColumns - 1], 1);

//   Matrix_Free((Matrix *)U);
//   Vector_Free(X1);
//   Matrix_Free(fordio);

//   AffineHermite(result, &H, &U);
//   Matrix_Free((Matrix *)result);
//   result = (Lattice *)Matrix_Copy(H);

//   Matrix_Free((Matrix *)H);
//   Matrix_Free((Matrix *)U);

//   /* Check whether the Lattice is NULL or not */

//   if (isEmptyLattice(result)) {
//     Matrix_Free((Matrix *)result);
//     return (EmptyLattice(X->NbRows));
//   }
//   return result;
// } /* LatticeIntersection */

// static Matrix *MakeDioEqforInter(Lattice *A, Lattice *B) {

//   Matrix *Dio;
//   int i, j;

// #ifdef DOMDEBUG
//   FILE *fp;
//   fp = fopen("_debug", "a");
//   fprintf(fp, "\nEntered MAKEDIOEQFORINTER \n");
//   fclose(fp);
// #endif

//   Dio = Matrix_Alloc(2 * (A->NbRows - 1) + 1, 3 * (A->NbColumns - 1) + 1);

//   for (i = 0; i < Dio->NbRows; i++)
//     for (j = 0; j < Dio->NbColumns; j++)
//       value_set_si(Dio->p[i][j], 0);

//   for (i = 0; i < A->NbRows - 1; i++) {
//     value_set_si(Dio->p[i][i], 1);
//     value_set_si(Dio->p[i + A->NbRows - 1][i], 1);
//   }
//   for (i = 0; i < A->NbRows - 1; i++)
//     for (j = 0; j < A->NbRows - 1; j++) {
//       value_oppose(Dio->p[i][j + A->NbRows - 1], A->p[i][j]);
//       value_oppose(Dio->p[i + (A->NbRows - 1)][j + 2 * (A->NbRows - 1)],
//                    B->p[i][j]);
//     }

//   /* Adding the affine part */

//   for (i = 0; i < A->NbColumns - 1; i++) {
//     value_oppose(Dio->p[i][Dio->NbColumns - 1], A->p[i][A->NbColumns - 1]);
//     value_oppose(Dio->p[i + A->NbRows - 1][Dio->NbColumns - 1],
//                  B->p[i][A->NbColumns - 1]);
//   }
//   value_set_si(Dio->p[Dio->NbRows - 1][Dio->NbColumns - 1], 1);
//   return Dio;
// } /* MakeDioEqforInter */

static void AddLattice(LatticeUnion *, Matrix *, Matrix *, Value, int);
LatticeUnion *SplitLattice(Matrix *, Matrix *, Matrix *);

/*
 * The function is transforming a lattice X in a union of lattices based on a
 starting lattice Y.
 * Note1: If the intersection of X and Y lattices is empty the result is identic
 with the first argument (X) because no operation can be made. *Note2: The
 function is available only for simple Lattices and not for a union of Lattices.
 *
 * Theorem : Given Two Lattices L1 and L2, (L2 subset of L1) there exists a
 *           Basis B = {b1, b2,..bn} of L1 and integers {a1, a2...,an} such
 *           that a1 divides a2, a2 divides a3 and so on and {a1b1, a2b2 ,...,
 *           .., anbn} is a Basis of L2. So given this theorem we can express
 *           the Lattice L1 in terms of Union of Lattices Involving L2, such
 *           that Lattice L1 = B1 = Union of (B2 + i1b1 + i2b2 + .. inbn) such
 *           that 0 <= i1 < a1; 0 <= i2 < a2; .......   0 <= in < an. We also
 *           know that A/B = A/(A Intersection B) and that (A Intersection B)
 *           is a subset of A. So, Making use of these two facts, we find the
 *           A/B. We Split The Lattice A in terms of Lattice (A Int B). From
 *           this Union of Lattices Delete the Lattice (A Int B).
 *
 *       Step 1:  Find Intersection = LatticeIntersection (A, B).
 *       Step 2:  Extract the Linear Parts of the Lattices A and Intersection.
 *                (while dealing with Basis we only deal with the Linear Parts)
 *       Step 3:  Let M1 = linear basis of A and M2 = linear basis of B.
 *                Let B1 and B2 be the Basis of A and B respectively,
 *                corresponding to the above Theorem.
 *                Then we Have B1 = M1 * U1 {a unimodular Matrix }
 *                and B2 = M2 * U2. M1 and M2 we know, they are the linear
 *                parts we obtained in Step 2. Our Task is now to find U1 and
 *                U2.
 *                We know that B1  * Delta = B2.
 *                i.e. M1 * U1 * Delta = M2 * U2
 *                or U1*Delta*U2Inverse = M1Inverse * M2.
 *                and Delta is the Diagonal Matrix which satisfies the
 *                above properties (in the Theorem).
 *                So Delta is nothing but the Smith Normal Form of
 *                M1Inverse * M2.
 *                So, first we have to find M1Inverse.
 *
 *                This Step, involves finding the Inverse of the Matrix M1.
 *                We find the Inverse using the Polylib function
 *                Matrix_Inverse. There is a catch here, the result of this
 *                function is an integral matrix, not necessarily the exact
 *                Inverse (since M1 need not be Unimodular), but a multiple
 *                of the actual inverse. The number by which we have to divide
 *                the matrix, is not obtained here as the input matrix is not
 *                a Polylib matrix { We input only the Linear part }. Later I
 *                give a way for finding that number.
 *
 *                M1Inverse = Matrix_Inverse ( M1 );
 *
 *      Step 4 :  MtProduct = Matrix_Product (M1Inverse, M2);
 *      Step 5 :  SmithNormalFrom (MtProduct, Delta, U, V);
 *                U1 = U and U2Inverse = V.
 *      Step 6 :  Find U2 = Matrix_Inverse  (U2inverse). Here there is no prob
 *                as U1 and its inverse are unimodular.
 *
 *      Step 7 :  Compute B1 = M1 * U1;
 *      Step 8 :  Compute B2 = M2 * U2;
 *      Step 9 :  Earlier when we computed M1Inverse, we knew that it was not
 *                the exact inverse but a multiple of it. Now we find the
 *                number, such that ( M1Inverse / number ) would give us the
 *                exact inverse of M1.
 *                We know that B1 * Delta = B2.
 *                Let k = B2[0][0] / B1[0][0].
 *                Let number = Delta[0][0]/k;
 *                This 'number' is the number we want.
 *                We Divide the matrix Delta by this number, to get the actual
 *                Delta such that B1 * Delta = B2.
 *     Step 10 :  Call Split Lattice (B1, B2, Delta ).
 *                This function returns the Union of Lattices in such a way
 *                that B2 is at the Head of this List.
 *
 *If the intersection between X and Y is empty then the result is NULL.
 */

LatticeUnion *Lattice2LatticeUnion(Lattice *X, Lattice *Y) {
  Lattice *B1 = NULL, *B2 = NULL, *newB1 = NULL, *newB2 = NULL,
          *Intersection = NULL;
  Matrix *U = NULL, *M1 = NULL, *M2 = NULL, *M1Inverse = NULL,
         *MtProduct = NULL;
  Matrix *Vinv, *V, *temp, *DiagMatrix;

  LatticeUnion *Head = NULL;
  int i;
  Value k;

  Intersection = LatticeIntersection(X, Y);
  #ifdef LATDIF_DEBUG
    fprintf(stderr,"Lattice intersection = ");
    Matrix_Print(stderr, P_VALUE_FMT, Intersection);
  #endif
  if (isEmptyLattice(Intersection)) {
    #ifdef LATDIF_DEBUG
      fprintf(stderr, "\nIn Lattice2LatticeUnion : the input lattices X and Y do "
                        "not have any common part\n");
    #endif
    Matrix_Free(Intersection);
    return NULL;
  }

  value_init(k);
  M1 = (Matrix *)ExtractLinearPart(X);
  M2 = (Matrix *)ExtractLinearPart(Intersection);
  #ifdef LATDIF_DEBUG
    fprintf(stderr, "M1 = ");
    Matrix_Print(stderr, P_VALUE_FMT, M1);
    fprintf(stderr, "M2 = ");
    Matrix_Print(stderr, P_VALUE_FMT, M2);
  #endif

  int maxdim = (M1->NbColumns > M1->NbRows)? M1->NbColumns : M1->NbRows;
  // compute left(!) inverse of M1 using a dirty patch
  M1Inverse = Matrix_Alloc(maxdim, maxdim);
  // temp = Matrix_Copy(M1);

  // complete temp with Id to make it square
  temp = Matrix_Alloc(maxdim, maxdim);
  for (i = 0; i < M1->NbRows; i++) {
    int j;
    for (j = 0; j < M1->NbColumns; j++)
      value_assign(temp->p[i][j], M1->p[i][j]);
    for( ; j < maxdim; j++)
      value_set_si(temp->p[i][j], (i==j));
  }
  for( ; i < temp->NbColumns ; i++ )
    for (int j = 0; j < M1->NbColumns; j++)
      value_set_si(temp->p[i][j], (i==j));

  Matrix_Inverse(temp, M1Inverse);
  Matrix_Free(temp);

  // get the right result: the rightmost columns of M1Inverse are just ignored :)
  M1Inverse->NbColumns = M1->NbRows;
  M1Inverse->NbRows = M1->NbColumns;
  temp = Matrix_Copy(M1Inverse);
  Matrix_Free(M1Inverse);
  M1Inverse = temp;

  #ifdef LATDIF_DEBUG
    fprintf(stderr, "M1Inverse = ");
    Matrix_Print(stderr, P_VALUE_FMT, M1Inverse);
  #endif
  
  MtProduct = Matrix_Alloc(M1Inverse->NbRows, M2->NbColumns);
  Matrix_Product(M1Inverse, M2, MtProduct);
  left_hermite(MtProduct, &DiagMatrix, &U, &Vinv);
  V = Matrix_Alloc(Vinv->NbRows, Vinv->NbColumns);
  Matrix_Inverse(Vinv, V);
  Matrix_Free(Vinv);
  B1 = Matrix_Alloc(M1->NbRows, U->NbColumns);
  B2 = Matrix_Alloc(M2->NbRows, V->NbColumns);
  Matrix_Product(M1, U, B1);
  Matrix_Product(M2, V, B2);
  Matrix_Free(M1);
  Matrix_Free(M2);
  value_division(k, B2->p[0][0], B1->p[0][0]);
  value_division(k, DiagMatrix->p[0][0], k);
  for (i = 0; i < DiagMatrix->NbRows; i++)
    value_division(DiagMatrix->p[i][i], DiagMatrix->p[i][i], k);

  #ifdef LATDIF_DEBUG
    fprintf(stderr, "B1 = ");
    Matrix_Print(stderr, P_VALUE_FMT, B1);
    fprintf(stderr, "B2 = ");
    Matrix_Print(stderr, P_VALUE_FMT, B2);
  #endif

  // previous method, supposed that B1 and B2 are square matrices:
  // newB1 = ChangeLatticeDimension(B1, B1->NbRows + 1);
  // Matrix_Free(B1);
  // newB2 = ChangeLatticeDimension(B2, B2->NbRows + 1);
  // Matrix_Free(B2);

  // New method: add a row and a column to B1 and B2:
  newB1 = Matrix_Alloc(B1->NbRows+1, B1->NbColumns+1);
  newB2 = Matrix_Alloc(B2->NbRows+1, B2->NbColumns+1);
  // copy the core of the matrix:
  for(i=0 ; i<B1->NbRows; i++) {
    for(int j=0 ; j<B1->NbColumns; j++) {
      value_assign(newB1->p[i][j], B1->p[i][j]);
      value_assign(newB2->p[i][j], B2->p[i][j]);
    }
  }
  // complete the last row with zeroes:
  for(int j=0 ; j<B1->NbColumns; j++) {
    value_set_si(newB1->p[B1->NbRows][j], 0);
    value_set_si(newB2->p[B1->NbRows][j], 0);
  }
  // complete last column (empty for B1, Intersection for B2):
  for (i = 0; i < newB2->NbRows; i++)
  {
    value_set_si(newB1->p[i][newB2->NbColumns - 1], (i==newB2->NbColumns - 1));
    value_assign(newB2->p[i][newB2->NbColumns - 1],
                 Intersection->p[i][Intersection->NbColumns - 1]);
  }

  #ifdef LATDIF_DEBUG
    fprintf(stderr, "newB2 = ");
    Matrix_Print(stderr, P_VALUE_FMT, newB2);
    fprintf(stderr, "Diag = ");
    Matrix_Print(stderr, P_VALUE_FMT, DiagMatrix);
  #endif

  Head = SplitLattice(newB1, newB2, DiagMatrix);
  Matrix_Free(newB1);
  Matrix_Free(MtProduct);
  Matrix_Free(M1Inverse);
  Matrix_Free(B1);
  Matrix_Free(B2);
  Matrix_Free(DiagMatrix);
  Matrix_Free(U);
  Matrix_Free(V);
  Matrix_Free(Intersection);
  value_clear(k);
  return Head;
}

/*
 * Return the Union of lattices that constitute the difference between
 * two single lattices: A - B.
 * The dimensions of A and B should be the same. 
 * Main algorithm: compute the intersection of A and B and take it out of A
 * If the difference is empty return NULL.
 * Allocates a LatticeUnion
 */
LatticeUnion *LatticeDifference(Lattice *A, Lattice *B) {

  Matrix *H, *X, *Y;
  int *pivots_columns;
  Lattice *Inter, *rest;
  LatticeUnion *Result;

  #ifdef DOMDEBUG
    FILE *fp;
    fp = fopen("_debug", "a");
    fprintf(fp, "\nEntered LATTICEDIFFERENCE \n");
    fclose(fp);
  #endif

  // Checking inputs:
  if (A->NbRows != B->NbRows) {
    errormsg1("LatticeDifference", "dimincomp",
      "input lattices A and B have incompatible dimensions (rows)");
    return NULL;
  }
  if (A->NbColumns != B->NbColumns) {
    errormsg1("LatticeDifference", "dimincomp",
      "input lattices A and B have incompatible dimensions (columns)");
    return NULL;
  }
  // normalize and create a copy A->X
  if (! isNormalLattice(A)) {
    AffineHermite(A, &H, NULL);
    X = H;
  }
  else {
    X = Matrix_Copy(A);
  }
  if (isEmptyLattice(X)) {
    Matrix_Free(X);
    return(NULL);
  }
  // no need to normalize B, it will just be used to compute
  // the (normalized) intersection with X
  // if (! isNormalLattice(B)) {
  //   AffineHermite(B, &H, NULL);
  //   Y = H;
  // }
  // else {
  //   Y = Matrix_Copy(B);
  // }
  #ifdef LATDIF_DEBUG
    fprintf(stderr, "--- Entering LatDiff ---\n"
        "A (normalized) = ");
    Matrix_Print(stderr, P_VALUE_FMT, X);
    fprintf(stderr, "B = ");
    Matrix_Print(stderr, P_VALUE_FMT, B);
  #endif

  // calculate the intersection between X and B
  Inter = LatticeIntersection(X, B);

  Result = NULL;
  #ifdef LATDIF_DEBUG
    fprintf(stderr, "Inter = ");
    Matrix_Print(stderr, P_VALUE_FMT, Inter);
  #endif
  if(!Inter){
    // if empty intersection return a copy of A (normalized)
    Result = LatticeUnion_Alloc();
    Result->M = X;
    return (Result);
  }

  // Prepare for main loop

  // rest is the rest of the lattice to be treated (intersection on first line(s), X on last line(s))
  rest = Matrix_Copy(X); // keep X safe, we still need it
  // get the positions of the pivots of (X and) Inter
  pivots_columns = malloc(sizeof(int) * A->NbRows);
  get_pivots_columns(Inter, pivots_columns);


  // -------------- MAIN LOOP --------------------

  // add each matrix with the line variant to the Result
  // and keep the intersection line variant in the rest list
  for (int line = 0; line < Inter->NbRows-1; line++) {
    // TODO: only consider the real pivots here, not lines below a previously treated pivot.

    #ifdef LATDIF_DEBUG
      fprintf(stderr, "+++ Enter main loop (%d)\n", line);
      fprintf(stderr, "+++ rest =\n");
      Matrix_Print(stderr, P_VALUE_FMT, rest);
    #endif

    // TODO: Inter can be removed below, just replace the line of rest by the one of the intersection before the function call and use this one.

    Result = generate_lattice_union_line(line, pivots_columns, X, Inter, rest, Result);
    #ifdef LATDIF_DEBUG
      fprintf(stderr, "+++ Intermediate result =\n");
      PrintLatticeUnion(stderr, P_VALUE_FMT, Result);
    #endif
  }

  // ------------ END MAIN LOOP --------------------

  // cleanup
  free(pivots_columns);
  Matrix_Free(rest);
  Matrix_Free(X);

  #ifdef LATDIF_DEBUG
    if(!Result)
      fprintf(stderr, "Empty Result\n");
    fprintf(stderr, "--- Exit LatDiff ---\n\n");
  #endif

  // No need to simplify since it is already computed in minimal form :)
  // if ((Result != NULL)) {
  //   Result = LatticeSimplify(Result);
  //   #ifdef LATDIF_DEBUG
  //     fprintf(stderr, "Simplified result = ");
  //     PrintLatticeUnion(stderr, P_VALUE_FMT, Result);
  //   #endif
  // }
  return Result;
} /* LatticeDifference */

/*
 * Given a Lattice 'B1' and a Lattice 'B2' and a Diagonal Matrix 'C' such that
 * 'B2' is a subset of 'B1' and C[0][0] divides C[1][1], C[1][1] divides C[2]
 * [2] and so on, output the list of matrices whose union is B1. The function
 * expresses the Lattice B1 in terms of B2 Unions of B1 = Union of {B2 + i0b0 +
 * i1b1 + .... + inbn} where 0 <= i0 < C[0][0]; 0 <= i1 < C[1][1] and so on and
 * {b0 ... bn} are the columns of Lattice B1. The list is so formed that the
 * Lattice B2 is the Head of the list.
 */
LatticeUnion *SplitLattice(Lattice *B1, Lattice *B2, Matrix *C) {

  int i;

  LatticeUnion *Head = NULL;
  Head = malloc(sizeof(LatticeUnion));
  Head->M = B2;
  Head->next = NULL;
  for (i = 0; i < C->NbRows; i++)
    AddLattice(Head, B1, B2, C->p[i][i], i);
  return Head;
} /* SplitLattice */

/*
 * Given lattices 'B1' and 'B2', an integer 'NumofTimes', a column number
 * 'Colnumber' and a pointer to a list of lattices, the function does the
 * following :-
 * For every lattice in the list, it adds a set of lattices such that the
 * affine part of the new lattices is greater than the original lattice by 0 to
 * NumofTimes-1 * {the (ColumnNumber)-th column of B1}.
 * Note :
 * Three pointers are defined to point at various points of the list. They are:
 * Head   -> It always points to the head of the list.
 * tail   -> It always points to the last element in the list.
 * marker -> It points to the element, which is the last element of the Input
 *           list.
 */
static void AddLattice(LatticeUnion *Head, Matrix *B1, Matrix *B2,
                       Value NumofTimes, int Colnumber) {

  LatticeUnion *temp, *tail, *marker;
  int j;
  Value i, maxval;

  // printf("Apply " P_VALUE_FMT "to column %d\n", NumofTimes, Colnumber);
  value_init(i);
  value_init(maxval);
  tail = Head;
  while (tail->next != NULL)
    tail = tail->next;
  marker = tail;

  for (temp = Head; temp != NULL; temp = temp->next) {
    for (value_set_si(i, 1); value_lt(i, NumofTimes); value_increment(i, i)) {
      Lattice *tempMatrix, *H;

      tempMatrix = Matrix_Copy(temp->M);
      for (j = 0; j < B2->NbRows; j++) {
        value_addmul(tempMatrix->p[j][B2->NbColumns - 1], i,
                     B1->p[j][Colnumber]);
        if(value_notzero_p(B1->p[j][Colnumber])) {
          value_multiply(maxval, NumofTimes, B1->p[j][Colnumber]);
          // if(value_neg_p(tempMatrix->p[j][B2->NbColumns - 1])) {

          // }
          // else {
          // value_print(stderr, P_VALUE_FMT, maxval);
          value_modulus(tempMatrix->p[j][B2->NbColumns - 1],
                          tempMatrix->p[j][B2->NbColumns - 1], maxval);
          // }
        }
      }
      tail->next = malloc(sizeof(LatticeUnion));
      AffineHermite(tempMatrix, &H, NULL);
      Matrix_Free(tempMatrix);
      tail->next->M = H;
      tail->next->next = NULL;
      tail = tail->next;
    }
    if (temp == marker)
      break;
  }
  value_clear(i);
  value_clear(maxval);
  return;
} /* AddLattice */

// /*
//  * Given a polyhedron 'A', store the Hermite basis 'B' and return the true
//  * dimension of the polyhedron 'A'.
//  * Algorithm :
//  *
//  *             1) First we find all the vertices of the Polyhedron A.
//  *                Now suppose the vertices are [v1, v2...vn], then
//  *                a particular set of vectors governing the space of A are
//  *                given by [v1-v2, v1-v3, ... v1-vn] (let us say V).
//  *                So we initially calculate these vectors.
//  *             2) Then there are the rays and lines which contribute to the
//  *                space in which A is going to lie.
//  *                So we append to the rays and lines. So now we get a matrix
//  *                {These are the rows} [ V ] [l1] [l2]...[lk]
//  *                where l1 to lk are either rays or lines of the Polyhedron A.
//  *             3) The above matrix is the set of vectors which determine
//  *                the space in which A is going to lie.
//  *                Using this matrix we find a Basis which is such that
//  *                the first 'm' columns of it determine the space of A.
//  *             4) But we also have to ensure that in the last 'n-m'
//  *                coordinates the Polyhedron is '0', this is done by
//  *                taking the image by B(inv) of A and finding the remaining
//  *                equalities, and composing it with the matrix B, so as
//  *                to get a new matrix which is the actual Hermite Basis of
//  *                the Polyhedron.
//  */
// int FindHermiteBasisofDomain(Polyhedron *A, Matrix **B) {

//   int i, j;
//   Matrix *temp, *temp1, *tempinv, *Newmat;
//   Matrix *vert, *rays, *result;
//   Polyhedron *Image;
//   int rank, equcount;
//   int noofvertices = 0, noofrays = 0;
//   int vercount, raycount;
//   Value lcm, fact;

// #ifdef DOMDEBUG
//   FILE *fp;
//   fp = fopen("_debug", "a");
//   fprintf(fp, "\nEntered FINDHERMITEBASISOFDOMAIN \n");
//   fclose(fp);
// #endif

//   POL_ENSURE_FACETS(A);
//   POL_ENSURE_VERTICES(A);

//   /* Checking is empty */
//   if (emptyQ(A)) {
//     B[0] = Identity(A->Dimension + 1);
//     return (-1);
//   }

//   value_init(lcm);
//   value_init(fact);
//   value_set_si(lcm, 1);

//   /* Finding the Vertices */
//   for (i = 0; i < A->NbRays; i++)
//     if ((value_notzero_p(A->Ray[i][0])) &&
//         value_notzero_p(A->Ray[i][A->Dimension + 1]))
//       noofvertices++;
//     else
//       noofrays++;

//   vert = Matrix_Alloc(noofvertices, A->Dimension + 1);
//   rays = Matrix_Alloc(noofrays, A->Dimension);
//   vercount = 0;
//   raycount = 0;

//   for (i = 0; i < A->NbRays; i++) {
//     if ((value_notzero_p(A->Ray[i][0])) &&
//         value_notzero_p(A->Ray[i][A->Dimension + 1])) {
//       for (j = 1; j < A->Dimension + 2; j++)
//         value_assign(vert->p[vercount][j - 1], A->Ray[i][j]);
//       value_lcm(lcm, lcm, A->Ray[i][j - 1]);
//       vercount++;
//     } else {
//       for (j = 1; j < A->Dimension + 1; j++)
//         value_assign(rays->p[raycount][j - 1], A->Ray[i][j]);
//       raycount++;
//     }
//   }

//   /* Multiplying the rows by the lcm */
//   for (i = 0; i < vert->NbRows; i++) {
//     value_division(fact, lcm, vert->p[i][vert->NbColumns - 1]);
//     for (j = 0; j < vert->NbColumns - 1; j++)
//       value_multiply(vert->p[i][j], vert->p[i][j], fact);
//   }

//   /* Drop the Last Columns */
//   temp = RemoveColumn(vert, vert->NbColumns - 1);
//   Matrix_Free(vert);

//   /* Getting the Vectors */
//   vert = Matrix_Alloc(temp->NbRows - 1, temp->NbColumns);
//   for (i = 1; i < temp->NbRows; i++)
//     for (j = 0; j < temp->NbColumns; j++)
//       value_subtract(vert->p[i - 1][j], temp->p[0][j], temp->p[i][j]);

//   Matrix_Free(temp);

//   /* Add the Rays and Lines */
//   /* Combined Matrix */
//   result = Matrix_Alloc(vert->NbRows + rays->NbRows, vert->NbColumns);
//   for (i = 0; i < vert->NbRows; i++)
//     for (j = 0; j < result->NbColumns; j++)
//       value_assign(result->p[i][j], vert->p[i][j]);

//   for (; i < result->NbRows; i++)
//     for (j = 0; j < result->NbColumns; j++)
//       value_assign(result->p[i][j], rays->p[i - vert->NbRows][j]);

//   Matrix_Free(vert);
//   Matrix_Free(rays);

//   rank = findHermiteBasis(result, &temp);
//   temp1 = ChangeLatticeDimension(temp, temp->NbRows + 1);

//   Matrix_Free(result);
//   Matrix_Free(temp);

//   /* Adding the Affine Part to take care of the Equalities */
//   temp = Matrix_Copy(temp1);
//   tempinv = Matrix_Alloc(temp->NbRows, temp->NbColumns);
//   Matrix_Inverse(temp, tempinv);
//   Matrix_Free(temp);
//   Image = DomainImage(A, tempinv, MAXNOOFRAYS);
//   Matrix_Free(tempinv);
//   Newmat = Matrix_Alloc(temp1->NbRows, temp1->NbColumns);
//   for (i = 0; i < rank; i++)
//     for (j = 0; j < Newmat->NbColumns; j++)
//       value_set_si(Newmat->p[i][j], 0);
//   for (i = 0; i < rank; i++)
//     value_set_si(Newmat->p[i][i], 1);
//   equcount = 0;
//   for (i = 0; i < Image->NbConstraints; i++)
//     if (value_zero_p(Image->Constraint[i][0])) {
//       for (j = 1; j < Image->Dimension + 2; j++)
//         value_assign(Newmat->p[rank + equcount][j - 1],
//                      Image->Constraint[i][j]);
//       ++equcount;
//     }
//   Domain_Free(Image);
//   for (i = 0; i < Newmat->NbColumns - 1; i++)
//     value_set_si(Newmat->p[Newmat->NbRows - 1][i], 0);
//   value_set_si(Newmat->p[Newmat->NbRows - 1][Newmat->NbColumns - 1], 1);
//   temp = Matrix_Alloc(Newmat->NbRows, Newmat->NbColumns);
//   Matrix_Inverse(Newmat, temp);
//   Matrix_Free(Newmat);
//   B[0] = Matrix_Alloc(temp1->NbRows, temp->NbColumns);

//   Matrix_Product(temp1, temp, B[0]);
//   Matrix_Free(temp1);
//   Matrix_Free(temp);
//   value_clear(lcm);
//   value_clear(fact);
//   return rank;
// } /* FindHermiteBasisofDomain */

// /*
//  * Return the image of a lattice 'A' by the invertible, affine, rational
//  * function 'M'.
//  */
// Lattice *LatticeImage(Lattice *A, Matrix *M) {

//   Lattice *Img, *temp, *Minv;

// #ifdef DOMDEBUG
//   FILE *fp;
//   fp = fopen("_debug", "a");
//   fprintf(fp, "\nEntered LATTICEIMAGE \n");
//   fclose(fp);
// #endif

//   if ((A->NbRows != M->NbRows) || (M->NbRows != M->NbColumns))
//     return (EmptyLattice(A->NbRows));

//   if (value_one_p(M->p[M->NbRows - 1][M->NbColumns - 1])) {
//     Img = Matrix_Alloc(M->NbRows, A->NbColumns);
//     Matrix_Product(M, A, Img);
//     return Img;
//   }
//   temp = Matrix_Copy(M);
//   Minv = Matrix_Alloc(temp->NbColumns, temp->NbRows);
//   Matrix_Inverse(temp, Minv);
//   Matrix_Free(temp);

//   Img = LatticePreimage(A, Minv);
//   Matrix_Free(Minv);
//   return Img;
// } /* LatticeImage */

// /*
//  * Return the preimage of a lattice 'L' by an affine, rational function 'G'.
//  * Algorithm:
//  *           (1) Prepare Diophantine equation :
//  *               [Gl -Ll][x y] = [Ga -La]{"l-linear, a-affine"}
//  *           (2) Solve the Diophantine equations.
//  *           (3) If there is solution to the Diophantine eq., extract the
//  *               general solution and the particular solution of x and that
//  *               forms the preimage of 'L' by 'G'.
//  */
// Lattice *LatticePreimage(Lattice *L, Matrix *G) {

//   Matrix *Dio, *U;
//   Lattice *Result;
//   Vector *X;
//   int i, j;
//   int rank;
//   Value divisor, tmp;

// #ifdef DOMDEBUG
//   FILE *fp;
//   fp = fopen("_debug", "a");
//   fprintf(fp, "\nEntered LATTICEPREIMAGE \n");
//   fclose(fp);
// #endif

//   /* Check for the validity of the function */
//   if (G->NbRows != L->NbRows) {
//     fprintf(stderr, "\nIn LatticePreimage: Incompatible types of Lattice and "
//                     "the function\n");
//     return (EmptyLattice(G->NbColumns));
//   }

//   value_init(divisor);
//   value_init(tmp);

//   /* Making Diophantine Equations [g -L] */
//   value_assign(divisor, G->p[G->NbRows - 1][G->NbColumns - 1]);
//   Dio = Matrix_Alloc(G->NbRows, G->NbColumns + L->NbColumns - 1);
//   for (i = 0; i < G->NbRows - 1; i++)
//     for (j = 0; j < G->NbColumns - 1; j++)
//       value_assign(Dio->p[i][j], G->p[i][j]);

//   for (i = 0; i < G->NbRows - 1; i++)
//     for (j = 0; j < L->NbColumns - 1; j++) {
//       value_multiply(tmp, divisor, L->p[i][j]);
//       value_oppose(Dio->p[i][j + G->NbColumns - 1], tmp);
//     }

//   for (i = 0; i < Dio->NbRows - 1; i++) {
//     value_multiply(tmp, divisor, L->p[i][L->NbColumns - 1]);
//     value_subtract(tmp, G->p[i][G->NbColumns - 1], tmp);
//     value_assign(Dio->p[i][Dio->NbColumns - 1], tmp);
//   }
//   for (i = 0; i < Dio->NbColumns - 1; i++)
//     value_set_si(Dio->p[Dio->NbRows - 1][i], 0);

//   value_set_si(Dio->p[Dio->NbRows - 1][Dio->NbColumns - 1], 1);
//   rank = SolveDiophantine(Dio, &U, &X);

//   if (rank == -1)
//     Result = EmptyLattice(G->NbColumns);
//   else {
//     Result = Matrix_Alloc(G->NbColumns, G->NbColumns);
//     for (i = 0; i < Result->NbRows - 1; i++)
//       for (j = 0; j < Result->NbColumns - 1; j++)
//         value_assign(Result->p[i][j], U->p[i][j]);

//     for (i = 0; i < Result->NbRows - 1; i++)
//       value_assign(Result->p[i][Result->NbColumns - 1], X->p[i]);
//     Matrix_Free(U);
//     Vector_Free(X);
//     for (i = 0; i < Result->NbColumns - 1; i++)
//       value_set_si(Result->p[Result->NbRows - 1][i], 0);
//     value_set_si(Result->p[i][i], 1);
//   }
//   Matrix_Free(Dio);
//   value_clear(divisor);
//   value_clear(tmp);
//   return Result;
// } /* LatticePreimage */

/*
 * Return True if the matrix 'm' is a valid lattice, otherwise return False.
 * Note: A valid lattice has the last row as [0 0 0 ... 1].
 */
Bool IsLattice(Matrix *m) {

  int i;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ISLATTICE \n");
  fclose(fp);
#endif

  /* Is it necessary to check if the lattice
     is fulldimensional or not here only? */

  if (m->NbRows != m->NbColumns)
    return False;

  for (i = 0; i < m->NbColumns - 1; i++)
    if (value_notzero_p(m->p[m->NbRows - 1][i]))
      return False;
  if (value_notone_p(m->p[i][i]))
    return False;
  return True;
} /* IsLattice */

// /*
//  *  Check whether the matrix 'm' is full row-rank or not.
//  */
// Bool isfulldim(Matrix *m) {

//   Matrix *h, *u;
//   int i;

//   /*
//      res = Hermite (m, &h, &u);
//      if (res != m->NbRows)
//      return False ;
//   */

//   Hermite(m, &h, &u);
//   for (i = 0; i < h->NbRows; i++)
//     if (value_zero_p(h->p[i][i])) {
//       Matrix_Free(h);
//       Matrix_Free(u);
//       return False;
//     }
//   Matrix_Free(h);
//   Matrix_Free(u);
//   return True;
// } /* isfulldim */

/*
 * This function takes as input a lattice list in which the lattices have the
 * same linear part, and almost the same affinepart, i.e. if A and B are two
 * of the lattices in the above lattice list and [a1, .. , an] and [b1 .. bn]
 * are the affineparts of A and B respectively, then for 0 < i < n ai = bi and
 * 'an' may not be equal to 'bn'. These are not the affine parts in the n-th
 * dimension, but the lattices have been tranformed such that the value of the
 * elment in the dimension on which we are simplifying is in the last row and
 * also the lattices are in a sorted order.
 *              This function also takes as input the dimension along which we
 * are simplifying and takes the diagonal element of the lattice along that
 * dimension and tries to find out the factors of that element and sees if the
 * list of lattices can be simplified using these factors. The output of this
 * function is the list of lattices in the simplified form and a flag to indic-
 * ate whether any form of simplification was actually done or not.
 */
static Bool Simplify(LatticeUnion **InputList, LatticeUnion **ResultList,
                     int dim) {

  int i;
  LatticeUnion *prev, *temp;
  factor allfac;
  Bool retval = False;
  int width;
  Value cnt, aux, k, fac, tmp, foobar; //, num

  if ((*InputList == NULL) || (InputList[0]->next == NULL))
    return False;

  value_init(aux);
  value_init(cnt);
  value_init(k);
  value_init(fac);
  // value_init(num);
  value_init(tmp);
  value_init(foobar);

  width = InputList[0]->M->NbRows - 1;
  allfac = allfactors(VALUE_TO_INT(InputList[0]->M->p[dim][dim]));
  value_set_si(cnt, 0);
  for (temp = InputList[0]; temp != NULL; temp = temp->next)
    value_increment(cnt, cnt);
  for (i = 0; i < allfac.count; i++) {
    value_set_si(foobar, allfac.fac[i]);
    value_division(aux, InputList[0]->M->p[dim][dim], foobar);
    if (value_ge(cnt, aux))
      break;
  }
  if (i == allfac.count) {
    value_clear(cnt);
    value_clear(aux);
    value_clear(k);
    value_clear(fac);
    // value_clear(num);
    value_clear(tmp);
    value_clear(foobar);
    free(allfac.fac);
    return False;
  }
  for (; i < allfac.count; i++) {
    Bool Present = False;
    value_set_si(k, 0);

    if (*InputList == NULL) {
      value_clear(cnt);
      value_clear(aux);
      value_clear(k);
      value_clear(fac);
      // value_clear(num);
      value_clear(tmp);
      value_clear(foobar);
      free(allfac.fac);
      return retval;
    }
    value_set_si(foobar, allfac.fac[i]);
    // value_division(num, InputList[0]->M->p[dim][dim], foobar);
    while (value_lt(k, foobar)) {
      Present = False;
      value_assign(fac, k);
      for (temp = *InputList; temp != NULL; temp = temp->next) {
        if (value_eq(temp->M->p[temp->M->NbRows - 1][temp->M->NbColumns - 1],
                     fac)) {
          value_set_si(foobar, allfac.fac[i]);
          value_addto(fac, fac, foobar);
          if (value_ge(fac, (*InputList)->M->p[dim][dim])) {
            Present = True;
            break;
          }
        }
        if (value_gt(temp->M->p[temp->M->NbRows - 1][temp->M->NbColumns - 1],
                     fac))
          break;
      }
      if (Present == True) {
        retval = True;
        if (*ResultList == NULL)
          *ResultList = temp = (LatticeUnion *)malloc(sizeof(LatticeUnion));
        else {
          for (temp = *ResultList; temp->next != NULL; temp = temp->next)
            ;
          temp->next = (LatticeUnion *)malloc(sizeof(LatticeUnion));
          temp = temp->next;
        }
        temp->M = Matrix_Copy(InputList[0]->M);
        temp->next = NULL;
        value_set_si(foobar, allfac.fac[i]);
        value_assign(temp->M->p[dim][dim], foobar);
        value_assign(temp->M->p[dim][width], k);
        value_set_si(temp->M->p[width][width], 1);

        /* Deleting the Lattices from the curlist */
        value_assign(tmp, k);
        prev = NULL;
        temp = InputList[0];
        while (temp != NULL) {
          if (value_eq(temp->M->p[width][width], tmp)) {
            if (temp == InputList[0]) {
              prev = temp;
              temp = InputList[0] = temp->next;
              Matrix_Free(prev->M);
              free(prev);
            } else {
              prev->next = temp->next;
              Matrix_Free(temp->M);
              free(temp);
              temp = prev->next;
            }
            value_set_si(foobar, allfac.fac[i]);
            value_addto(tmp, tmp, foobar);
          } else {
            prev = temp;
            temp = temp->next;
          }
        }
      }
      value_increment(k, k);
    }
  }
  free(allfac.fac);
  value_clear(cnt);
  value_clear(aux);
  value_clear(k);
  value_clear(fac);
  // value_clear(num);
  value_clear(tmp);
  value_clear(foobar);
  return retval;
} /* Simplify */

/*
 * This function is used in the qsort function in sorting the lattices. Given
 * two lattices 'A' and 'B', both in HNF, where A = [ [a11 0], [a21, a22, 0] .
 * .... [an1, .., ann] ] and B = [ [b11 0], [b21, b22, 0] ..[bn1, .., bnn] ],
 * then A < B, if there exists a pair <i,j> such that [aij < bij] and for every
 * other pair <i1, j1>, 0<=i1<i, 0<=j1<j [ai1j1 = bi1j1].
 */
static int LinearPartCompare(const void *A, const void *B) {

  Lattice **L1, **L2;
  int i, j;

  L1 = (Lattice **)A;
  L2 = (Lattice **)B;

  for (i = 0; i < L1[0]->NbRows - 1; i++)
    for (j = 0; j <= i; j++) {
      if (value_gt(L1[0]->p[i][j], L2[0]->p[i][j]))
        return 1;
      if (value_lt(L1[0]->p[i][j], L2[0]->p[i][j]))
        return -1;
    }
  return 0;
} /* LinearPartCompare */

/*
 * This function takes as input a List of Lattices and sorts them on the basis
 * of their Linear parts. It sorts in place, as a result of which the input
 * list is modified to the sorted order.
 */
static void LinearPartSort(LatticeUnion *Head) {

  int cnt;
  Lattice **Latlist;
  LatticeUnion *temp;

  cnt = 0;
  for (temp = Head; temp != NULL; temp = temp->next)
    cnt++;

  Latlist = (Lattice **)malloc(sizeof(Lattice *) * cnt);

  cnt = 0;
  for (temp = Head; temp != NULL; temp = temp->next)
    Latlist[cnt++] = temp->M;

  qsort(Latlist, cnt, sizeof(Lattice *), LinearPartCompare);

  cnt = 0;
  for (temp = Head; temp != NULL; temp = temp->next)
    temp->M = Latlist[cnt++];

  free(Latlist);
  return;
} /* LinearPartSort */

/*
 * This function is used in 'AfiinePartSort' in sorting the lattices with the
 * same linear part. GIven two lattices 'A' and 'B' with affineparts [a1 .. an]
 * and [b1 ... bn], then A < B if for some 0 < i <= n, ai < bi and for 0 < i1 <
 * i, ai1 = bi1.
 */
static int AffinePartCompare(const void *A, const void *B) {

  int i;
  Lattice **L1, **L2;

  L1 = (Lattice **)A;
  L2 = (Lattice **)B;

  for (i = 0; i < L1[0]->NbRows; i++) {
    if (value_gt(L1[0]->p[i][L1[0]->NbColumns - 1],
                 L2[0]->p[i][L2[0]->NbColumns - 1]))
      return 1;

    if (value_lt(L1[0]->p[i][L1[0]->NbColumns - 1],
                 L2[0]->p[i][L2[0]->NbColumns - 1]))
      return -1;
  }
  return 0;
} /* AffinePartCompare */

/*
 * This function takes a list of lattices with the same linear part and sorts
 * them on the basis of their affine part. The sorting is done in place.
 */
static void AffinePartSort(LatticeUnion *List) {

  int cnt;
  Lattice **LatList;
  LatticeUnion *tmp;

  cnt = 0;
  for (tmp = List; tmp != NULL; tmp = tmp->next)
    cnt++;

  LatList = malloc(sizeof(Lattice *) * cnt);

  cnt = 0;
  for (tmp = List; tmp != NULL; tmp = tmp->next)
    LatList[cnt++] = tmp->M;

  qsort(LatList, cnt, sizeof(Lattice *), AffinePartCompare);

  cnt = 0;
  for (tmp = List; tmp != NULL; tmp = tmp->next)
    tmp->M = LatList[cnt++];
  free(LatList);
  return;
} /* AffinePartSort */

static Bool AlmostSameAffinePart(LatticeUnion *A, LatticeUnion *B) {

  int i;

  if ((A == NULL) || (B == NULL))
    return False;

  for (i = 0; i < A->M->NbRows - 1; i++)
    if (value_ne(A->M->p[i][A->M->NbColumns - 1],
                 B->M->p[i][A->M->NbColumns - 1]))
      return False;
  return True;
} /* AlmostSameAffinePart */

/*
 * This function takes a list of lattices having the same linear part and tries
 * to simplify these lattices. This may not be the only way of simplifying the
 * lattices. The function returns a list of partially simplified lattices and
 * also a flag to tell whether any simplification was performed at all.
 */
static Bool AffinePartSimplify(LatticeUnion *curlist, LatticeUnion **newlist) {

  int i;
  Value aux;
  LatticeUnion *temp, *curr, *next;
  LatticeUnion *nextlist;
  Bool change = False, chng;

  if (curlist == NULL)
    return False;

  if (curlist->next == NULL) {
    curlist->next = newlist[0];
    newlist[0] = curlist;
    return False;
  }

  value_init(aux);
  for (i = 0; i < curlist->M->NbRows - 1; i++) {

    /* Interchanging the elements of the Affine part for easy computation
       of the sort (using qsort) */

    for (temp = curlist; temp != NULL; temp = temp->next) {
      value_assign(aux,
                   temp->M->p[temp->M->NbRows - 1][temp->M->NbColumns - 1]);
      value_assign(temp->M->p[temp->M->NbRows - 1][temp->M->NbColumns - 1],
                   temp->M->p[i][temp->M->NbColumns - 1]);
      value_assign(temp->M->p[i][temp->M->NbColumns - 1], aux);
    }
    AffinePartSort(curlist);
    nextlist = NULL;
    curr = curlist;
    while (curr != NULL) {
      next = curr->next;
      if (!AlmostSameAffinePart(curr, next)) {
        curr->next = NULL;
        chng = Simplify(&curlist, newlist, i);
        if (nextlist == NULL)
          nextlist = curlist;
        else {
          LatticeUnion *tmp;
          for (tmp = nextlist; tmp->next; tmp = tmp->next)
            ;
          tmp->next = curlist;
        }
        change = (Bool)(change | chng);
        curlist = next;
      }
      curr = next;
    }
    curlist = nextlist;

    /* Interchanging the elements of the Affine part for easy computation
       of the sort (using qsort) */

    for (temp = curlist; temp != NULL; temp = temp->next) {
      value_assign(aux,
                   temp->M->p[temp->M->NbRows - 1][temp->M->NbColumns - 1]);
      value_assign(temp->M->p[temp->M->NbRows - 1][temp->M->NbColumns - 1],
                   temp->M->p[i][temp->M->NbColumns - 1]);
      value_assign(temp->M->p[i][temp->M->NbColumns - 1], aux);
    }
    if (curlist == NULL)
      break;
  }
  if (*newlist == NULL)
    *newlist = nextlist;
  else {
    for (curr = *newlist; curr->next != NULL; curr = curr->next)
      ;
    curr->next = nextlist;
  }
  value_clear(aux);
  return change;
} /* AffinePartSimplify */

static Bool SameLinearPart(LatticeUnion *A, LatticeUnion *B) {

  int i, j;
  if ((A == NULL) || (B == NULL))
    return False;
  for (i = 0; i < A->M->NbRows - 1; i++)
    for (j = 0; j <= i; j++)
      if (value_ne(A->M->p[i][j], B->M->p[i][j]))
        return False;

  return True;
} /* SameLinearPart */

/*
 * Given a union of lattices, return a simplified list of lattices.
 */
LatticeUnion *LatticeSimplify(LatticeUnion *latlist) {
  // fprintf(stderr, "Entering LatticeSimplify with:");
  // PrintLatticeUnion(stderr, P_VALUE_FMT, latlist);
  LatticeUnion *curlist, *nextlist;
  LatticeUnion *curr, *next;
  Bool change = True;

  curlist = latlist;
  while (change == True) {
    change = False;
    LinearPartSort(curlist);
    curr = curlist;
    nextlist = NULL;
    while (curr != NULL) {
      next = curr->next;
      if (!SameLinearPart(curr, next)) {
        curr->next = NULL;
        change |= AffinePartSimplify(curlist, &nextlist);
        curlist = next;
      }
      curr = next;
    }
    curlist = nextlist;
  }
  return curlist;
} /* LatticeSimplify */

int intcompare(const void *a, const void *b) {

  int *i, *j;

  i = (int *)a;
  j = (int *)b;
  if (*i > *j)
    return 1;
  if (*i < *j)
    return -1;
  return 0;
} /* intcompare */

// compute the prime factors of n, including n itself if it is prime.
static factor prime_factors(int n)
{
  int tabsize = 10;
  int div = 2;
  int rest = n;
  factor result;

  // result init
  result.count = 0;
  result.fac = malloc(tabsize * sizeof(int));

  while(div*div <= rest)
  {
    if(rest % div == 0)
    {
      // double the size of the array if necessary
      if(result.count == tabsize)
      {
        tabsize *= 2;
        result.fac = realloc(result.fac, tabsize * sizeof(int));
      }
      // add div to the result
      result.fac[result.count++] = div;
      rest /= div;
    }
    else
      div += (1 + (1&div)); // 2, 3, 5, 7, 9, ...
  }
  if(rest != 1)
  {
    // add rest
    if(result.count == tabsize)
      {
        tabsize += 1;
        result.fac = realloc(result.fac, tabsize * sizeof(int));
      }
      result.fac[result.count++] = rest;
  }
  return result;
}

/*
 * Compute the prime factors of Value n, including n itself if it is prime.
 * reuses or allocates a Vector of Values
 * returns the number of values put into the result
 */
static int value_prime_factors(Value n, Vector **result) {

  if(!*result) {
    // allocate if NULL
    *result = Vector_Alloc(10);
  }

  int tabsize = 0; // current size
  Value div, rest, two, tmp;
  value_init(div);
  value_init(rest);
  value_init(two);
  value_init(tmp); // temp variable for tests

  value_set_si(div, 2);
  value_assign(rest, n);
  value_set_si(two, 2);

  while(True) {
    value_multiply(tmp, div, div);
    if(value_gt(tmp, rest)) // if(div * div > rest)
      break;                // exit while

    value_modulus(tmp, rest, div);
    if(value_zero_p(tmp)) { // if(rest % div == 0)

      if((*result)->Size == tabsize) { // increase vector size if needed
        *result = Vector_Realloc((*result), (*result)->Size*2);
      }
      // add div to result
      value_assign((*result)->p[tabsize], div); 
      tabsize++;

      value_division(rest, rest, div); // rest /= div;
    }
    else{
      // div = 2, 3, 5, 7, 9, 11, ...
      if(value_eq(div, two))        //  if(div == 2)
        value_set_si(div, 3);       //    div = 3;
      else
        value_addto(div, div, two); //    div += 2;
    }
  }

  // if something's left add it to the result:
  if(value_notone_p(rest)) {
    if((*result)->Size == tabsize){
      (*result) = Vector_Realloc((*result), (*result)->Size+1);
    }
    value_assign((*result)->p[tabsize], rest);
    tabsize++;
  }
  value_clear(rest);
  value_clear(div);
  value_clear(two);
  value_clear(tmp);

  return (tabsize);
}

static factor allfactors(int n)
{
  // compute the prime decomposition (including n if it is prime)
  factor primes = prime_factors(n);
  factor result;

  int size = (1<<(primes.count));   // max size of result = 2^(prime.count)
  result.count = 1;
  result.fac = malloc(size * sizeof(int));
  result.fac[0] = 1;

  for(int mask=1; mask < size; mask++) {
    // multiply the prime factors for the bits that are 1 in the mask :)
    int val = 1;
    for(int bits=0, rest=mask; rest ; bits++, rest>>=1) {
      if((rest & 1))
        val *= primes.fac[bits];
    }
    // add val to res.fac only if it is not already there
    // (this will suppress duplicates, in case a prime factor is there several times)
    if(val > result.fac[result.count-1]) {
      result.fac[result.count++] = val;
    }
  }
  // always remove the last one, that is n.
  result.count--;

  free(primes.fac);
  return result;
}



// static factor allfactors(int num) {

//   int i, j, tmp;
//   int noofelmts = 1;
//   int *list, *newlist;
//   int count;
//   factor result;
//   printf("%d\n",num);

//   list = (int *)malloc(sizeof(int));
//   list[0] = 1;

//   tmp = num;
//   for (i = 2; i <= polylib_sqrt(tmp); i++) {
//     if ((tmp % i) == 0) {
//       newlist = (int *)malloc(sizeof(int) * 2 * (noofelmts + 1));
//       for (j = 0; j < noofelmts; j++)
//         newlist[j] = list[j];
//       newlist[j] = i;
//       for (j = 0; j < noofelmts; j++)
//         newlist[j + noofelmts + 1] = i * list[j];
//       free(list);
//       list = newlist;
//       noofelmts = 2 * noofelmts + 1;
//       tmp = tmp / i;
//       i = 1;
//     }
//   }

//   if ((tmp != 0) && (tmp != num)) {
//     newlist = (int *)malloc(sizeof(int) * 2 * (noofelmts + 1));
//     for (j = 0; j < noofelmts; j++)
//       newlist[j] = list[j];
//     newlist[j] = tmp;
//     for (j = 0; j < noofelmts; j++)
//       newlist[j + noofelmts + 1] = tmp * list[j];
//     free(list);
//     list = newlist;
//     noofelmts = 2 * noofelmts + 1;
//   }
//   qsort(list, noofelmts, sizeof(int), intcompare);
//   count = 1;
//   for (i = 1; i < noofelmts; i++)
//     if (list[i] != list[i - 1])
//       list[count++] = list[i];
//   if (list[count - 1] == num)
//     count--;

//   result.fac = (int *)malloc(sizeof(int) * count);
//   result.count = count;
//   for (i = 0; i < count; i++)
//     result.fac[i] = list[i];
//   free(list);
//   return result;
// } /* allfactors */


/*
 * Takes into parameters two lattices A and B of the form:
 *  A =   A' | a      B =   B' | b
 *      0..0 | 1          0..0 | 1
 * 
 *  Copies them in a matrix (Tmp), used to calculate the left hermite,
 *  the Lattice of size A->Nbrows x min(A->NbCols, B->NbCols) is the
 *  intersection of A and B.
 * 
 * If the result is empty return NULL.
 * 
 */
Lattice* LatticeIntersection(Lattice* A, Lattice* B) {
  Lattice *Tmp, *H, *Res;
  if(A->NbRows != B->NbRows){
    errormsg1("LatticeIntersection", "dimincomp", "incompatible dimensions!");
    return NULL;
  }
  #ifdef LATINTER_DEBUG
  fprintf(stderr,"---Entering LatInter---\nMatrix A = ");
  Matrix_Print(stderr, P_VALUE_FMT, A);
  fprintf(stderr,"Matrix B = ");
  Matrix_Print(stderr, P_VALUE_FMT, B);
  #endif

  // Tmp will be in the form:
  // 
  //   1     0...0 |   1      0...0
  //   a      A'   |   b       B'
  // -------------------------------
  //   1     0...0 |    0 ..    0
  //   a      A'   |    0 ..    0
  
  Tmp = Matrix_Alloc(A->NbRows*2, A->NbColumns+B->NbColumns);

  if (!Tmp) {
    errormsg1("LatticeIntersection", "outofmem", "Not enough memory space!");
    return NULL;
  }
  
  //copying A in Tmp:
  
  // initalizing the top-left 1
  value_assign(Tmp->p[0][0], A->p[A->NbRows-1][A->NbColumns-1]);
  value_assign(Tmp->p[A->NbRows][0], A->p[A->NbRows-1][A->NbColumns-1]);
  
  //copy of the constant vector a
  for(int i=1; i<A->NbRows; i++) {
    value_assign(Tmp->p[i][0], A->p[i-1][A->NbColumns-1]);
    value_assign(Tmp->p[i+A->NbRows][0], A->p[i-1][A->NbColumns-1]);
  }
  //copy of the matrix kernel A'
  for(int i = 1 ; i < A->NbRows; i++) {
    for(int j = 1; j < A->NbColumns; j++){
      value_assign(Tmp->p[i][j], A->p[i-1][j-1]);
      value_assign(Tmp->p[i+A->NbRows][j], A->p[i-1][j-1]);
    }
  }

  // copying B into tmp:
  value_assign(Tmp->p[0][A->NbColumns], B->p[B->NbRows-1][B->NbColumns-1]);
  
  //the constant b (last col of lattice)
  for (int i = 1; i < B->NbRows; i++) {
    value_assign(Tmp->p[i][A->NbColumns], B->p[i-1][B->NbColumns-1]);
  }
  
  for (int i = 1; i < B->NbRows; i++){
    for (int j = 1; j < B->NbColumns; j++) {
      value_assign(Tmp->p[i][j+A->NbColumns], B->p[i-1][j-1]);
    }
  }
  
  #ifdef LATINTER_DEBUG
    fprintf(stderr,"H init = ");
    Matrix_Print(stderr,P_VALUE_FMT, Tmp);
  #endif
  

  // left_hermite of the TMP
  // H is the matrix that contains the solution. it is of the form:
  // 
  // H =   D  |   0          D is a square matrix
  //     -----------------
  //       X  |  1 0.0
  //          |  r  R
  // 
  // with  R    r
  //      0..0  1   being our result
  // if the number above r is not 1 then the intersection is not integer
  // (no solution to the intersection)

  left_hermite(Tmp, &H, NULL, NULL);


  #ifdef LATINTER_DEBUG
    fprintf(stderr,"\nH = ");
    Matrix_Print(stderr,P_VALUE_FMT,H);
  #endif
  Matrix_Free(Tmp);

  // get the result. if the value on top-left of R is not 1 then we have an empty solution.

  // what is the number of columns of zeros on the first NbRows Rows of H?
  // the matrix has A->NbColumns + B-> NbColumns columns.
  int nbcol = 0;
  for(int col_num = H->NbColumns-1 ; col_num >= 0; col_num--) {
    int i;
    for(i = 0; i < A->NbRows; i++) {
      if(value_notzero_p(H->p[i][col_num]))
        break;
    }
    if(i != A->NbRows) {
      // there is a non-zero value
      break;
    }
    nbcol++;
  }

  if(value_notone_p(H->p[A->NbRows][H->NbColumns-nbcol])) {
    #ifdef LATINTER_DEBUG
      fprintf(stderr,"\nEmpty intersection\n");
    #endif
    Matrix_Free(H);
    return NULL;
  }

  Res = Matrix_Alloc(A->NbRows, nbcol);
  for (int i = 0; i < Res->NbRows; i++) {
    for (int j = 0; j < Res->NbColumns; j++) {
        value_assign(Res->p[i][j], H->p[i + H->NbRows - Res->NbRows][j + H->NbColumns - Res->NbColumns]);
    }
  }
  Matrix_Free(H);
  // put Res in the proper affine form (didn't want to write a loop, had a pre-written function.)
  Matrix_Move_Homogeneous_Dim_Last(Res);

  #ifdef LATINTER_DEBUG
    fprintf(stderr, "\nLatticeIntersection result = ");
    Matrix_Print(stderr, P_VALUE_FMT, Res);
    fprintf(stderr,"---Exiting LatInter---\n\n");
  #endif

  return Res;
}

/*
 * Moves the constant part (last line and last row) as first line and row
 * of the matrix.
 * This is useful to compute the HNF and keeping the affine part as top-left
 * non-nul result. The same function can be called again to get the result
 * of affine HNF.
 *  A =  A'  | c     ->    z | 0..0
 *      0..0 | z           c |  A'
 */
void Matrix_Move_Homogeneous_Dim_First(Matrix *A) {
  if(A->NbRows == 0 || A->NbColumns == 0)
    return;

  Value tmp;
  value_init(tmp);
  // puts the last column first:
  for(int i=0; i<A->NbRows; i++) {
    // on row i
    value_assign(tmp, A->p[i][A->NbColumns-1]); // tmp = last col value
    for(int j = A->NbColumns-1; j > 0; j--) {
      value_assign(A->p[i][j], A->p[i][j-1]);  // [j] <- [j-1]
    }
    value_assign(A->p[i][0], tmp);  // [0]<- tmp
  }
  // then puts the last row first:
  for(int j = 0; j < A->NbColumns; j++) {
    value_assign(tmp, A->p[A->NbRows-1][j]); // tmp = last row value
    for(int i = A->NbRows-1; i > 0; i--) {
      value_assign(A->p[i][j], A->p[i-1][j]); // [i] <- [i-1]
    }
    value_assign(A->p[0][j], tmp); // [0]<- tmp
  }
  value_clear(tmp);
} /* Matrix_Move_Homogeneous_Dim_First */


/*
 * Moves the constant part on a homogenous matrix (first line and first row) as last line and last row
 * of the matrix.
 * This is useful to compute the HNF and keeping the affine part as top-left
 * non-nul result. The same function can be called again to get the result
 * of affine HNF.
 *  A =  A'  | c     ->    z | 0..0
 *      0..0 | z           c |  A'
 */
void Matrix_Move_Homogeneous_Dim_Last(Matrix *A) {
  if(A->NbRows == 0 || A->NbColumns == 0)
    return;

  Value tmp;
  value_init(tmp);
  // puts the first col in the end
  for (int i = 0; i < A->NbRows; i++) {
    value_assign(tmp,A->p[i][0]); // tmp = first col value
    for (int j = 0; j < A->NbColumns-1; j++) {
      value_assign(A->p[i][j],A->p[i][j+1]); // [i] <- [i+1]
    }
    value_assign(A->p[i][A->NbColumns-1],tmp); //[last] <- tmp
  }
  for (int j = 0; j < A->NbColumns; j++) {
    value_assign(tmp,A->p[0][j]); // tmp first row value
    for (int i = 0; i < A->NbRows-1; i++) {
      value_assign(A->p[i][j],A->p[i+1][j]);
    }
    value_assign(A->p[A->NbRows-1][j],tmp);
  }
  value_clear(tmp);
} /* Matrix_Move_Homogeneous_Dim_Last */

// not used anymore:
// Vector* get_pivots(Matrix* A){
//   Vector* pivot = Vector_Alloc(A->NbRows);
//   int j=0;
//   for(int i=0; i < A->NbRows ; i++) {
//     if(value_zero_p(A->p[i][j])){
//       value_set_si(pivot->p[i],1);
//     }else {
//       value_assign(pivot->p[i],A->p[i][j]);
//       j++;
//     }
//   }
//   return pivot;
// }

/*
 * get the column numbers of the pivots.
 * since the matrix is not necessarily square, retreive the right pivot
 * column number for each line number.
 * Fills the array columns (needs to be allocated by caller)
 */
static void get_pivots_columns(Matrix* A, int *columns)
{
  int col = 0;
  for(int i=0; i < A->NbRows ; i++) {
    if(i>0 && value_zero_p(A->p[i][col])){
      // zero in this column: take the previous one
      columns[i] = col-1;
    }
    else {
      // there's a non zero value on this column, take it and increase col.
      columns[i] = col;
      col++;
    }
  }
}

/*
 * Generate all variants of the line 'line_nb' in lattice matrix A
 * - add all variants not intersecting Intersection to the result
 * - replace the corresponding line of the rest L by the intersection
 *
 * The intersection line is used as a basis reference lattice, and all
 * variants of the corresponding line in A are generated:
 *    if Intersection contains line "*..* p 0..0 c"
 *    and the corresponding pivot of A is pA
 *    then generate new lines :  *..* p 0..0 c+pA; *..* p 0..0 c+2pA;
 *    *..* p 0..0 c+3pA; ...  (optim.: can be decomposed in prime factors)
 *
 *    Adjust the constants of the lines below, that depend on p,
 *
 * Add all newly generated lattices (with the new lines) to Result.
 */
static LatticeUnion *generate_lattice_union_line(int line_nb, int *pivots_columns,
    Lattice *A, Lattice *Intersection, Lattice *rest, LatticeUnion *Result)
{
  // Value cst; // loop index is a Value = iteration * pivot
  // Value iteration; // this is the iteration number
  Value step, multiply, modulo, ratio, tmp;
  Vector *prime_factors = NULL; // Vector of Values, reuse memory several times (from previous step).
  int num_factors;
  int pivot_col = pivots_columns[line_nb];

  value_init(step);
  value_init(multiply);
  value_init(modulo);
  value_init(ratio);
  value_init(tmp); // for computing tests on values


  // simplify directly when building

  // replace the current line in the rest with the one from the intersection
  // will consider this line as basis for generating all variants
  for(int i = 0; i < rest->NbColumns; i++) {
    value_assign(rest->p[line_nb][i], Intersection->p[line_nb][i]);
  }
  // no need to update the constant of lines below the pivot here

  // get the ratio between A and rest, to be used as multiplier for every generated new line
  value_division(ratio, rest->p[line_nb][pivot_col], A->p[line_nb][pivot_col]); // diag inter / diag A

  #ifdef LATDIF_DEBUG
    fprintf(stderr, "Considering line %d. Rest pivot = ", line_nb);
    value_print(stderr, P_VALUE_FMT, rest->p[line_nb][pivot_col]);
    fprintf(stderr, " Ratio = ");
    value_print(stderr, P_VALUE_FMT, ratio);
    fprintf(stderr, "\n");
  #endif

  // consider the decomposition in prime factors of the "pivot" = ratio (inters. pivot / A pivot):
  // if the "pivot" is 15, will take out the right p%3==0/1/2 and p%5==0/1/2/3/4
  // only one case p%15=c (the intersection) will not enter these (combination of) cases :)
  // can be empty, if p=1 then size=0 and the whole loop is skipped.
  // if a prime factor appears multiple times, multiply it: (2,2,2) -> (2,4,8)
  num_factors = value_prime_factors(ratio, &prime_factors);  // prime factors of pivot (ratio).

  for(int p = 0; p < num_factors; p++) {
    // consider the prime factor: prime_factors->p[p]

    // if the previous prime factor is the same:
    // example with 2:
    // get the cases p%4==0/1 and 2/3 after we took out i%2=0 or 1 in the previous step
    // next step get p%8==0+hit or 4+hit (hit is a potential intersection hit constant)

    // example with: 3*3*3
    // step 0: m=3, it=1 (init=0):
    // 3i + 0/1/2 -> if 3i+1 is the intersection, take out 3i+0 and 3i+2 (add to result).
    // Next step will not consider 3 * this (so 9i+0/3/6 and 9i+2/5/8), and just scan this:
    // step 1: m=9, it=3, (init=1):
    // 9i + 1/4/7 -> if 9i+7 is the intersection, take out 9i+1 and 9i+4
    // step 2: m=27, it=9, (init=7):
    // 27i+ 7/16/25

    // general case:
    // - if same primer factor as previously:
    //    * iterator step = previous multiplier
    //    * multiply = prime factor * previous multiplier
    //   else (new multiplier):
    //    * iteration step = 1 (* A initial pivot)
    //    * multiply = prime factor (* A initial pivot)
    //      (A initial pivot always divides the pivot of the intersection)
    // - init loop value = intersection constant % iterator step

    // check if same prime factor as before
    if(p>0 && value_eq(prime_factors->p[p-1], prime_factors->p[p])) {
      value_assign(step, multiply);
      value_multiply(multiply, multiply, prime_factors->p[p]);
                  // multiply = multiply*prime_factors->p[p];
    }
    else {
      value_assign(step, A->p[line_nb][pivot_col]);
      value_multiply(multiply, prime_factors->p[p], A->p[line_nb][pivot_col]);
    }
    #ifdef LATDIF_DEBUG
      fprintf(stderr, "multiply = ");
      value_print(stderr, P_VALUE_FMT, multiply);
      fprintf(stderr, ", step = ");
      value_print(stderr, P_VALUE_FMT, step);
      fprintf(stderr, "\n");
    #endif

    // Iterate on each possible 'modulo':
    // from a possible intersection value, to 'multiply', with 'step'
    for(value_modulus(modulo, Intersection->p[line_nb][Intersection->NbColumns-1], step);
        value_lt(modulo, multiply);
        value_addto(modulo, modulo, step)) {
      // consider line: multiply * x + modulo
      #ifdef LATDIF_DEBUG
        fprintf(stderr, "  -> considering line: ");
        value_print(stderr, P_VALUE_FMT, multiply);
        fprintf(stderr, " * x + ");
        value_print(stderr, P_VALUE_FMT, modulo);
      #endif

      value_modulus(tmp, Intersection->p[line_nb][Intersection->NbColumns-1], multiply);
      if(value_eq(tmp, modulo)) {
        // no need to do anything there, this modulo hits the intersection
        // and will be considered at next loop iteration or in the rest :)
        #ifdef LATDIF_DEBUG
          fprintf(stderr, " -> part of the intersection, ignoring\n");
        #endif
      }
      else {
        // this line does not hit the intersection, add it to the result.
        #ifdef LATDIF_DEBUG
          fprintf(stderr, " -> add it to result\n");
        #endif

        LatticeUnion *newResult = LatticeUnion_Alloc();
        Matrix *newLat = Matrix_Copy(rest); // get a copy of rest

        // link newResult to Result,
        newResult->M = NULL;
        newResult->next = Result;
        Result = newResult;

        // then update current line
        // new pivot:
        value_assign(newLat->p[line_nb][pivot_col], multiply);
        // new constant:
        value_assign(newLat->p[line_nb][newLat->NbColumns-1], modulo);

        // // adjust the coefficients below the changed line:
        // // take them from the intersection!
        // for(int cc = 0; cc < pivot_col; cc++) {
        //   for(int ll = line_nb+1; ll < A->NbRows; ll++) {
        //     value_assign(newLat->p[ll][cc], Intersection->p[ll][cc]);
        //   }
        // }
        // adjust the rows below because they change depending on the changed pivot&constant above:
        // if a coefficient below the pivot is not zero, multiply the coefficient by ratio
        // recompute the constant accordingly (adding (modulo/step))   ->  * ratio ??????
        for(int ll = line_nb+1; ll < A->NbRows; ll++) {
          if(value_notzero_p(newLat->p[ll][pivot_col])) {

            // new coefficient
            // value_multiply(newLat->p[ll][pivot_col], newLat->p[ll][pivot_col], ratio);
            // new coefficient: set it to the coef of the intersection
            value_assign(newLat->p[ll][pivot_col], Intersection->p[ll][pivot_col]);

            // new constant
            value_division(tmp, modulo, step); // iteration number 0/1/2/3/... one of them is the intersection
            value_multiply(tmp, tmp, A->p[ll][pivot_col]); // multiplied by old value below pivot
            value_addto(newLat->p[ll][newLat->NbColumns-1],
                        newLat->p[ll][newLat->NbColumns-1], tmp);
          }
        }
        // tranforms the new lattice back to HNF and store it into Result
        AffineHermite(newLat, &Result->M, NULL);
        Matrix_Free(newLat);
      }
    }
  }

  // adjust the column below the pivot: take it from the intersection
  for(int ll = line_nb+1; ll < A->NbRows; ll++) {
    value_assign(rest->p[ll][pivot_col], Intersection->p[ll][pivot_col]);
  }

  
  // // OLD VERSION:
  // // Add each possible alternate line to the Result
  // for(value_assign(cst, pivotA), value_set_si(iteration, 1);
  //     value_lt(cst, diagInter->p[line_nb]);
  //     value_addto(cst, cst, pivotA), value_add_int(iteration, iteration, 1))
  // {
  //   LatticeUnion *newResult = LatticeUnion_Alloc();
  //   Matrix *newLat = Matrix_Copy(rest); // from rest
  //   // store and link the new lattice in first position of the Result list
  //   newResult->M = newLat;
  //   newResult->next = *Result;
  //   *Result = newResult;

  //   // now update this line
  //   value_addto(newLat->p[line_nb][newLat->NbColumns-1], 
  //               newLat->p[line_nb][newLat->NbColumns-1], cst);

  //   value_modulus(newLat->p[line_nb][newLat->NbColumns-1], 
  //                 newLat->p[line_nb][newLat->NbColumns-1], diagInter->p[line_nb]);

  //   // ajust the constants below this row because they change depending on the changed line above ^^
  //   // if the corresponding value below the pivot (= multiplier) is not zero
  //   for(int ll = line_nb+1; ll < A->NbRows; ll++) {
  //     // scan the lines below searching for non-zero values

  //     // TODO: bug here, line_nb is not necessarily the column of the pivot (matrix can be non square)
  //     if(value_notzero_p(A->p[ll][line_nb])) {
  //       // add this value to the constant of line ll:
  //       // the multiplier  A->p[ll][line_nb] * the iteration number
  //       value_addmul(newLat->p[ll][newLat->NbColumns-1], iteration, A->p[ll][line_nb]);
  //       value_modulus(newLat->p[ll][newLat->NbColumns-1], 
  //                     newLat->p[ll][newLat->NbColumns-1], diagInter->p[ll]);
  //     }
  //   }
  // }

  // cleanup
  if(prime_factors)
    Vector_Free(prime_factors);
  value_clear(tmp);
  value_clear(ratio);
  value_clear(modulo);
  value_clear(multiply);
  value_clear(step);
  return (Result);
}

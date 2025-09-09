#include <polylib/polylib.h>
#include <stdlib.h>

// #define LATINTER_DEBUG 1
// #define LATDIF_DEBUG 1


static void get_pivots_columns(Matrix* A, int *columns);
static int value_prime_factors(Value n, Vector **result);
static LatticeUnion *generate_lattice_union_line(
  int line_nb, int *pivots_columns, Matrix *A, Matrix *Intersection, Matrix *L,
  LatticeUnion *Result);


  /*
 * Print the contents of a list of Lattices 'Head'
 */
void PrintLatticeUnion(FILE *fp, char *format, LatticeUnion *Head) {

  LatticeUnion *temp;

  for (temp = Head; temp != NULL; temp = temp->next)
    Matrix_Print(fp, format, temp->M);
  return;
} /* PrintLatticeUnion */

/*
 * Free the memory allocated for a list of lattices 'Head'
 */
void LatticeUnion_Free(LatticeUnion *Head) {
  while (Head) {
    LatticeUnion *temp;
    temp = Head;
    Head = temp->next;
    Matrix_Free(temp->M);
    free(temp);
  }
  return;
} /* LatticeUnion_Free */


/*
 * Allocate a head for a list of Lattices
 */
LatticeUnion *LatticeUnion_Alloc(void) {

  LatticeUnion *temp;

  temp = malloc(sizeof(LatticeUnion));
  temp->M = NULL;
  temp->next = NULL;
  return temp;
} /* LatticeUnion_Alloc */


/*
 * Return True if lattice 'A' is empty, otherwise return False.
 */
Bool isEmptyLattice(Matrix *A) {
  return(A == NULL || A->NbColumns == 0);
} /* isEmptyLattice */


/*
 * Move the constant part (last line and last row) as first line and row
 * of the matrix.
 * This is useful to compute a HNF and keep the affine part (as top-left
 * non-nul result).
 *  A =  A'  | c     ->    z | 0..0
 *      0..0 | z           c |  A'
 * Usage: modifies A in place.
 */
void Matrix_Move_Homogeneous_Dim_First(Matrix *A)
{
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
 * Moves the constant part of a transformed lattice (first line and first row)
 * as last line and row.
 * This is useful to convert back the affine part at bottom-right after a HNF.
 *  A =  A'  | c     <-    z | 0..0
 *      0..0 | z           c |  A'
 * Usage: modifies A in place.
 */
void Matrix_Move_Homogeneous_Dim_Last(Matrix *A)
{
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


/*
 * Check if the matrix 'A' is in affine Hermite normal form.
 */
Bool isNormalLattice(Matrix *A)
{
  // Check the affine-homogeneous lattice (last line/column is first):
  // Let A =  L   l
  //         0..0 1
  //
  // A' =  0..0 1
  //         L  l
  //
  // check if A' verifies:
  // - first element of column (= pivot) greater than zero
  // - all elements left of pivot lower than pivot.
  // 
  // Example:
  //    1 0 0 0
  //    * 0 0 0
  //    < + 0 0
  //    < < + 0
  //    * * * 0
  // all < of a line are lower than the + (pivot) of this line
  // * is anything.
  // a column of zero is valid on the right of L.

  int previous_nnl = -1;

  if(value_notone_p(A->p[A->NbRows-1][A->NbColumns-1])) {
    // the bottom right value should be 1.
    return (False);
  }

  for (int j = 0; j < A->NbColumns - 1; j++) {
    // consider column j of L
    int nnl; // line number of the pivot (non-null line)

    // find the line number of the first non-null element
    for(nnl = 0; nnl<A->NbRows; nnl++) {
      if(value_notzero_p(A->p[nnl][j]))
        break;
    }
    if(nnl == A->NbRows) {
      // this is a column of zeros, valid, just continue.
      previous_nnl = nnl;
      continue;
    }
    if(nnl <= previous_nnl) {
      // there is a non-zero value higher than the previous column
      return(False);
    }
    previous_nnl = nnl;

    // The pivot is: A->p[nnl][j]
    // check that the pivot is positive
    if(value_neg_p(A->p[nnl][j])) {
      return(False);
    }

    // check that the values left of pivot are lower than pivot
    for(int i = 0; i < j; i++) {
      if (value_ge(A->p[nnl][i], A->p[nnl][j])) {
        return(False);
      }
    }
    // check that the constant in l is lower than pivot
    if (value_ge(A->p[nnl][A->NbColumns-1], A->p[nnl][j])) {
      return(False);
    }
  }
  return(True);
} /* isNormalLattice */

/*
 * Return the affine Hermite normal form of the affine lattice 'A'. The
 * affine Hermite form if a lattice is stored in 'H' and the unimodular
 * matrix corresponding to A = H U is stored in the matrix 'U' (if not NULL)
 *
 * Algorithm:
 *     1) move the homogeneous dimensions first (on top-left) in A
 *     2) compute the left_hermite of the matrix A' = H' U'
 *     3) move back the homogeneous dimensions (bottom-right) in A, H and U.
 * Works also on non square matrices (lattices having less row than columns)
 *
 * U can be NULL (will be ignored)
 * *H and *U can be NULL (will be allocated by left_hermite)
 */
void AffineHermite(Matrix *A, Matrix **H, Matrix **U)
{
  // for left hermite to include the constant, move it on top-left:
  Matrix_Move_Homogeneous_Dim_First(A);
  left_hermite(A, H, U, NULL);
  Matrix_Move_Homogeneous_Dim_Last(*H);
  if(U)
    Matrix_Move_Homogeneous_Dim_Last(*U);
  Matrix_Move_Homogeneous_Dim_Last(A); // restore A as it was

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
// void AffineSmith(Matrix *A, Matrix **U, Matrix **V, Matrix **Diag) {

//   Matrix *temp;
//   Matrix *Uinv;
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
 * Given two lattices 'A' and 'B', verify if lattice 'A' is included in 'B' or
 * not. If 'A' is included in 'B' the 'A' intersection 'B', will be 'A'. So,
 * compute 'A' intersection 'B' and check if it is the same as 'A'.
 */
Bool LatticeIncludes(Matrix *A, Matrix *B) {

  Matrix *temp, *HA;
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
Bool sameLattice(Matrix *A, Matrix *B) {

  Matrix *HA, *HB;
  int i, j;
  Bool result = True;

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
 * Return the Union of lattices that constitute the difference between
 * two single lattices: A - B.
 * 
 * if A = NULL it's the universe.
 * The dimensions of A and B should be the same. 
 * Main algorithm: compute the intersection of A and B and take it out of A
 * If the difference is empty return NULL.
 * Allocates a LatticeUnion.
 */
LatticeUnion *LatticeDifference(Matrix *A, Matrix *B) {

  Matrix *H, *X;
  int *pivots_columns;
  Matrix *Inter, *rest;
  LatticeUnion *Result = NULL;

  if(B->NbRows == 1) {
    return(NULL);
  }
  // Checking inputs:
  if(!A) {
    Value gcd;
    value_init(gcd);
    A = Matrix_Copy(B);
    // Divide each column of A by its gcd to get the minimal lattice
    // spreading the same space than B
    for(int j = 0; j < A->NbColumns-1; j++) {
      int i;
      // initial gcd value
      for(i = 0; i < A->NbRows; i++) {
        if(value_notzero_p(A->p[i][j])) {
          value_assign(gcd, A->p[i][j]);
          break;
        }
      }
      // complete gcd computation
      for(i = i+1; i < A->NbRows; i++) {
        value_gcd(gcd, gcd, A->p[i][j]);
      }

      if(value_notone_p(gcd)) {
        // divide the column by its gcd:
        for(i = 0; i < A->NbRows; i++) {
          value_division(A->p[i][j], A->p[i][j], gcd);
        }
      }
    }
    // AffineHermite(A, &X, NULL); // not necessary
    // Matrix_Free(A);
    X = A;
    value_clear(gcd);
  }
  else {
    if (A->NbRows != B->NbRows) {
      errormsg1("LatticeDifference", "dimincomp",
        "input lattices A and B have incompatible dimensions (rows)");
      return NULL;
    }
    // NbColumn can be different if there is a column of zero in one of them

    if(! isNormalLattice(A)) {
      AffineHermite(A, &H, NULL);
      X = H;
    }
    else {
      X = Matrix_Copy(A);
    }
  }
  if(isEmptyLattice(X)) {
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
        "A (normalized) = X = ");
    Matrix_Print(stderr, P_VALUE_FMT, X);
    fprintf(stderr, "B = ");
    Matrix_Print(stderr, P_VALUE_FMT, B);
  #endif

  // calculate the intersection between X and B
  Inter = LatticeIntersection(X, B);
  if(!Inter) {
    #ifdef LATDIF_DEBUG
      fprintf(stderr, "Empty intersection, returning A\n");
    #endif
    // if empty intersection return a copy of A (normalized)
    Result = LatticeUnion_Alloc();
    Result->M = X;
    return (Result);
  }
  #ifdef LATDIF_DEBUG
    fprintf(stderr, "Inter = ");
    Matrix_Print(stderr, P_VALUE_FMT, Inter);
  #endif

  // if Inter has only one column, there is a problem in the loop:
  // line 0 would have no pivot.
  if(Inter->NbColumns == 1) {
    Matrix_Free(Inter);
    Matrix_Free(X);
    return(NULL);
  }

  // Prepare for main loop:
  // rest will be the rest of the lattice X to be treated
  // (Intersection on first line(s)/column(s), X on bottom-right)
  rest = Matrix_Copy(X);
  // get the positions of the pivots of (X and) Inter
  pivots_columns = malloc(sizeof(int) * X->NbRows);
  get_pivots_columns(Inter, pivots_columns);

  // -------------- MAIN LOOP --------------------

  // add each matrix with the line variant to the Result
  for (int line = 0; line < Inter->NbRows-1; line++) {
    if(line > 0 && pivots_columns[line] == pivots_columns[line-1]) {
      // only consider the *real* pivots here,
      // ignore lines below a previously treated pivot.
      continue;
    }

    #ifdef LATDIF_DEBUG
      fprintf(stderr, "+++ Enter main loop (%d)\n", line);
      fprintf(stderr, "+++ rest =\n");
      Matrix_Print(stderr, P_VALUE_FMT, rest);
    #endif

    Result = generate_lattice_union_line(line, pivots_columns, X, Inter,
                rest, Result);
    #ifdef LATDIF_DEBUG
      fprintf(stderr, "+++ Intermediate result =\n");
      PrintLatticeUnion(stderr, P_VALUE_FMT, Result);
    #endif
  }

  // ------------ END MAIN LOOP --------------------
  #ifdef LATDIF_DEBUG
    if(!Result)
      fprintf(stderr, "Empty Result\n");
    fprintf(stderr, "--- Exit LatDiff ---\n\n");
  #endif

  // cleanup
  free(pivots_columns);
  Matrix_Free(Inter);
  Matrix_Free(rest);
  Matrix_Free(X);

  return Result;
} /* LatticeDifference */

// Tried to use this, but Urt has zero columns too :(
// use: left_hermite(Matrix *M, Matrix **Hp, Matrix **Qp, Matrix **Up)
//  |A  B| . |Ult Urt| =  M U = H = |D    0  |
//  |A  0|   |Ulb Urb|              |X  inter|  <- A Urt = inter
// (then do the preimage of P by Urt)

/*
 * Compute the intersection between two lattices.
 * If the result is empty return NULL.
 *
 * Algorithm:
 * Let:
 *  A =   A' | a      B =   B'  | b
 *      0..0 | 1           0..0 | 1
 * 
 * Build matrix Tmp as:
 *   1   0...0 |   1    0...0
 *   a    A'   |   b     B'
 * ------------+--------------
 *   1   0...0 |    0 .. 0
 *   a    A'   |    0 .. 0
 * 
 * Then computes H = left Hermite of Tmp
 * H is of the form:
 * H =   D  |     0            D is a square matrix
 *     -----+-----------
 *       X  |  1 0 ... 0
 *          |  r    R
 * 
 * with   R   | r
 *      0...0 | 1   being our result
 *
 * if the number above r is not 1 then the intersection is not integer
 * (there is no solution to the intersection)
 */
Matrix* LatticeIntersection(Matrix* A, Matrix* B)
{
  Matrix *Tmp, *H, *Res;
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
  
  Tmp = Matrix_Alloc(A->NbRows*2, A->NbColumns + B->NbColumns);

  // copying A in Tmp (twice):
  value_assign(Tmp->p[0][0], A->p[A->NbRows-1][A->NbColumns-1]);
  value_assign(Tmp->p[A->NbRows][0], A->p[A->NbRows-1][A->NbColumns-1]);
  // constant vector a
  for(int i = 1; i < A->NbRows; i++) {
    value_assign(Tmp->p[i][0], A->p[i-1][A->NbColumns-1]);
    value_assign(Tmp->p[i+A->NbRows][0], A->p[i-1][A->NbColumns-1]);
  }
  // matrix kernel A'
  for(int i = 1 ; i < A->NbRows; i++) {
    for(int j = 1; j < A->NbColumns; j++){
      value_assign(Tmp->p[i][j], A->p[i-1][j-1]);
      value_assign(Tmp->p[i+A->NbRows][j], A->p[i-1][j-1]);
    }
  }

  // copying B into Tmp:
  value_assign(Tmp->p[0][A->NbColumns], B->p[B->NbRows-1][B->NbColumns-1]);
  // constant vector b
  for (int i = 1; i < B->NbRows; i++) {
    value_assign(Tmp->p[i][A->NbColumns], B->p[i-1][B->NbColumns-1]);
  }
  // matrix kernel B'
  for (int i = 1; i < B->NbRows; i++){
    for (int j = 1; j < B->NbColumns; j++) {
      value_assign(Tmp->p[i][j+A->NbColumns], B->p[i-1][j-1]);
    }
  }
  #ifdef LATINTER_DEBUG
    fprintf(stderr,"H init = ");
    Matrix_Print(stderr, P_VALUE_FMT, Tmp);
  #endif

  // // TRIED:
  // Matrix *U = NULL;
  // left_hermite(Tmp, &H, NULL, &U);

  left_hermite(Tmp, &H, NULL, NULL);

  

  #ifdef LATINTER_DEBUG
    fprintf(stderr,"\nH = ");
    Matrix_Print(stderr,P_VALUE_FMT,H);
  #endif
  Matrix_Free(Tmp);

  // what is the number of columns of zeros on the first NbRows rows of H?
  // the matrix has A->NbColumns + B-> NbColumns columns.
  int nbcol = 0;
  for(int col_num = H->NbColumns-1 ; col_num >= 0; col_num--) {
    int i;
    for(i = 0; i < A->NbRows; i++) {
      if(value_notzero_p(H->p[i][col_num]))
        break;
    }
    if(i != A->NbRows) { // there is a non-zero value
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
        value_assign(Res->p[i][j],
          H->p[i + H->NbRows - Res->NbRows][j + H->NbColumns - Res->NbColumns]);
    }
  }
  Matrix_Free(H);
  Matrix_Move_Homogeneous_Dim_Last(Res);
  #ifdef LATINTER_DEBUG
    fprintf(stderr, "\nLatticeIntersection result = ");
    Matrix_Print(stderr, P_VALUE_FMT, Res);
    // // TRIED:
    // // Get Urt such that: A Urt = Res
    // Matrix *Urt = NULL;
    // Matrix_Print(stdout, P_VALUE_FMT, U);
    // fprintf(stderr,"A Urt = Inter.\n Urt = ");
    // Matrix_subMatrix(U, 0, U->NbColumns-Res->NbColumns, Res->NbRows, U->NbColumns, &Urt);
    // Matrix_Move_Homogeneous_Dim_Last(Urt);
    // Matrix_Print(stdout, P_VALUE_FMT, Urt);

    fprintf(stderr,"---Exiting LatInter---\n\n");
  #endif

  return(Res);
}


// Utilities for LatticeDifference:

/*
 * Get the column numbers of the pivots.
 * since the matrix is not necessarily square, retrieve the right pivot
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
} /* get_pivots_columns */


/*
 * Compute the prime factors of Value n, including n itself if it is prime.
 * reuses or allocates a Vector of Values
 * returns the number of values put into the result
 *
 * *result is a vector of Values, that can be larger than the return value
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
        *result = Vector_Realloc((*result), (*result)->Size * 2);
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


/*
 * Generate all variants of the pivot on line 'line_nb' in lattice matrix A
 * - add all variants not intersecting Intersection to Result
 * - replace the corresponding line of the rest by the intersection
 * - replace the corresponding column of the rest by the intersection
*
 * The intersection line is used as a basis reference lattice, and all
 * variants of the corresponding line in A are generated:
 *    if Intersection contains line "*..* p 0..0 c"
 *    and the corresponding pivot of A is pA
 *    then generate new lines :  *..* p 0..0 c+pA; *..* p 0..0 c+2pA;
 *    *..* p 0..0 c+3pA; ...
 *   (optimized: p is decomposed in prime factors)
 *
 * Adjust the lines below (that need to be updated since p changes):
 *   - coefficients below the pivot p
 *   - constants for the non-zero coefficients
 *
 * Add all newly generated lattices to Result, and return the new Result.
 */
static LatticeUnion *generate_lattice_union_line(int line_nb, int *pivots_columns,
            Matrix *A, Matrix *Intersection, Matrix *rest, LatticeUnion *Result)
{
  Value step, multiply, modulo, ratio, tmp;
  Vector *prime_factors = NULL; // Vector of Values, reuse memory several times
                                // (from previous step).
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

  // get the ratio between A and rest, to be used as multiplier for every
  // generated new line
  value_division(ratio, rest->p[line_nb][pivot_col], A->p[line_nb][pivot_col]);

  #ifdef LATDIF_DEBUG
    fprintf(stderr, "Considering line %d. Rest pivot = ", line_nb);
    value_print(stderr, P_VALUE_FMT, rest->p[line_nb][pivot_col]);
    fprintf(stderr, " Ratio = ");
    value_print(stderr, P_VALUE_FMT, ratio);
    fprintf(stderr, "\n");
  #endif

  // consider the decomposition in prime factors of the "pivot" = ratio
  // (inters. pivot / A pivot):
  // if the "pivot" is 15, will take out the right p%3==0/1/2 and p%5==0/1/2/3/4
  // only one case p%15=c (the intersection) will not enter these (combination
  // of) cases :)
  // can be empty, if p=1 then size=0 and the whole loop is skipped.
  // if a prime factor appears multiple times, multiply-accumulate:
  // (2,2,2) -> (2,4,8)
  num_factors = value_prime_factors(ratio, &prime_factors);
                                          // prime factors of pivot ratio.

  // scan the prime factors: prime_factors->p[p].
  for(int p = 0; p < num_factors; p++) {
    // if the previous prime factor is the same:
    // example with: 3*3*3
    // - step 0: m=3, it=1 (init=0):
    // 3i + 0/1/2 -> if 3i+1 is the intersection, take out 3i+0 and 3i+2
    // (add to result).
    // Next step will not consider 3 * these (so 9i+0/3/6 and 9i+2/5/8).
    // - step 1: m=9, it=3, (init=1):
    // 9i + 1/4/7 -> if 9i+7 is the intersection, take out 9i+1 and 9i+4
    // - step 2: m=27, it=9, (init=7):
    // 27i+ 7/16/25

    // general case:
    // - if same prime factor as previously:
    //    * iterator step = previous multiplier
    //    * multiply = prime factor * previous multiplier
    //   else (new multiplier):
    //    * iteration step = 1 (* initial pivot of A)
    //    * multiply = prime factor (* initial pivot of A)
    //      (initial pivot of A always divides the pivot of the intersection)
    if(p>0 && value_eq(prime_factors->p[p-1], prime_factors->p[p])) {
      value_assign(step, multiply);
                  // step = multiply (from previous iteration)
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
    // from a possible intersection value, to 'multiply', with step 'step'
    // -> init loop value = intersection constant % iterator step
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

      value_modulus(tmp, Intersection->p[line_nb][Intersection->NbColumns-1],
                    multiply);
      if(value_eq(tmp, modulo)) {
        // no need to do anything there, this modulo hits the intersection
        // and will be considered in the rest :)
        #ifdef LATDIF_DEBUG
          fprintf(stderr, " -> part of the intersection, ignoring\n");
        #endif
      }
      else {
        LatticeUnion *newResult;
        // this line does not hit the intersection, add it to the result.
        #ifdef LATDIF_DEBUG
          fprintf(stderr, " -> add it to result\n");
        #endif

        Matrix *newLat = Matrix_Copy(rest); // get a copy of rest
        // and update current line. New pivot:
        value_assign(newLat->p[line_nb][pivot_col], multiply);
        // New constant:
        value_assign(newLat->p[line_nb][newLat->NbColumns-1], modulo);

        // adjust the rows below: they change depending on the changed pivot
        // and constant: if a coefficient below the pivot is not zero, set
        // it to the intersection coef., and recompute the constant
        // accordingly (adding (modulo/step))
        for(int ll = line_nb+1; ll < A->NbRows; ll++) {
          if(value_notzero_p(newLat->p[ll][pivot_col])) {
            // new coefficient: set it to the one of the intersection
            value_assign(newLat->p[ll][pivot_col],
                         Intersection->p[ll][pivot_col]);
            // adjust constant:
            value_division(tmp, modulo, step); // iteration number 0/1/2/...
            value_addmul(newLat->p[ll][newLat->NbColumns-1],
                         tmp, A->p[ll][pivot_col]); // multiplied by pivot A
          }
        }

        // link newResult to Result
        newResult = LatticeUnion_Alloc();
        newResult->M = NULL; // set by Hermite
        newResult->next = Result;
        Result = newResult;
        // transforms the new lattice to HNF and store it into Result
        AffineHermite(newLat, &Result->M, NULL);
        Matrix_Free(newLat);
      }
    }
  }

  // adjust the column below the pivot in rest (from the intersection)
  for(int ll = line_nb+1; ll < A->NbRows; ll++) {
    value_assign(rest->p[ll][pivot_col], Intersection->p[ll][pivot_col]);
  }

  // cleanup
  if(prime_factors)
    Vector_Free(prime_factors);
  value_clear(tmp);
  value_clear(ratio);
  value_clear(modulo);
  value_clear(multiply);
  value_clear(step);

  return (Result);
} /* generate_lattice_union_line */

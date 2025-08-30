#include <polylib/polylib.h>
#include <stdlib.h>

// DEFINITIONS USED BELOW:
//
// - a single LBL (sLBL function prefix) is a single lattice associated to
//   a polyhedral domain - can be a union of polyhedra
//
// - an LBL is a chained list of single LBLs, with possibly multiple
//   different lattices using the ->next structure element.
//
// They use the same type (LBL *), but next is not used in the first case.
//
// See the header of Zpolyhedron.h for more information


// debug this file:
// #define DEBUG
#ifdef DEBUG
  #define LAT_TEST 1
  #define CANONICAL_DEBUG 1
  #define INTERSECTION_DEBUG 1
  #define DIFFERENCE_DEBUG 1
#endif
  #define CANONICAL_DEBUG 1

static LBL *sLBL_Intersection(LBL *, LBL *);
static LBL *sLBL_Copy(LBL *A);
static void sLBL_Free(LBL *L);
static LBL *sLBL_Difference(LBL *, LBL *);
static LBL *sLBL_Image(LBL *, Matrix *);
static LBL *sLBL_Preimage(LBL *, Matrix *);
static void sLBL_Print(FILE *fp, const char *format, LBL *A);
static void sLBL_Canonical(LBL* A);
static LBL *FindLatticePred(Matrix *L, LBL *A);
static LBL *LBL_sLBL_Difference(LBL* A, LBL* B);
static int count_zeroCols (Matrix* M);

// typedef struct forsimplify {
//   Polyhedron *Pol;
//   LatticeUnion *LatUni;
//   struct forsimplify *next;
// } ForSimplify;


/*
 * Returns True if 'A' is empty, otherwise returns False.
 * A can be a non-simplified list of empty LBLs
 */
Bool isEmptyLBL(LBL *A)
{
  if (A == NULL)
    return True;
  if (emptyQ(A->P)) {
    // check the emptiness of next
    return(isEmptyLBL(A->next));
  }
  return False;
} /* isEmptyLBL */


/*
 * Given a matrix 'Lat' and a domain 'Domain', allocate space and return
 * the LBL corresponding to the image of the integer points of 'Poly'
 * by the affine function 'Lat', in canonical form (HNF, no equalities in
 * Domain).
 */
LBL *LBLAlloc(Matrix *Lat, Polyhedron *Domain)
{
  LBL *A;

  POL_ENSURE_FACETS(Domain);
  POL_ENSURE_VERTICES(Domain);

  if (Lat->NbColumns != Domain->Dimension + 1) {
    errormsg1("LBLAlloc", "dimincomp",
      "the Lattice and the Polyhedron are not compatible to form a LBL");
    return NULL;
  }

  A = malloc(sizeof(LBL));
  if (!A) {
    errormsg1("LBLAlloc", "outofmem", "Out of Memory");
    return NULL;
  }
  A->next = NULL;
  A->P = Domain_Copy(Domain);
  A->Lat = Matrix_Copy(Lat);

  CanonicalLBL(A);
  return A;
} /* LBLAlloc */


/*
 * Free the memory used by the single LBL 'L'
 * Internal function, users should use LBLFree.
 */
static void sLBL_Free(LBL *L)
{
  if (L == NULL)
    return;
  if(L->Lat)
    Matrix_Free(L->Lat);
  if(L->P)
    Domain_Free(L->P);
  free(L);
  return;
} /* sLBL_Free */


/*
 * Free the memory used by the LBL 'L'
 */
void LBLFree(LBL *L)
{
  if (L == NULL)
    return;
  LBLFree(L->next);
  sLBL_Free(L);
} /* LBLFree */


/*
 * Return a copy of the single LBL 'A'.
 * Internal function, users should use LBLCopy.
 */
static LBL *sLBL_Copy(LBL *A)
{
  return (LBLAlloc(A->Lat, A->P));
} /* sLBL_Copy */


/*
 * Return a copy of the LBL 'L'
 */
LBL *LBLCopy(LBL *L)
{
  LBL *copy;
  copy = sLBL_Copy(L);

  if (L->next != NULL)
    copy->next = LBLCopy(L->next);
  return copy;
} /* LBLCopy */


/*
 * Concatenate the LBLs 'A' and 'B',
 * and return a pointer to the new LBL.
 * Consumes the memory of A and of B (no need to free) to build
 * the result. Internal function only, do not use to build unions.
 */
static LBL *LBL_concatenate(LBL *A, LBL *B)
{
  LBL *tmp;

  if (isEmptyLBL(A)) {
    LBLFree(A);
    return (B);
  }
  if (isEmptyLBL(B)) {
    LBLFree(B);
    return (A);
  }

  for(tmp = A; tmp->next; tmp = tmp->next)
    ;
  tmp->next = B;
  
  return (A);
} /* LBL_concatenate */


/*
 * Return the empty Z-polyhedron of dimension 'dimension'
 * Lat = (0 ... 0 1)^T
 * P = empty polyhedron of dimension 0
 */
LBL *EmptyLBL(int dimension)
{
  LBL *A;

  A = malloc(sizeof(LBL));
  if(!A) {
    errormsg1("EmptyLBL", "outofmem", "Out of Memory");
  }

  A->Lat = Matrix_Alloc(dimension+1, 1);
  for(int j = 0 ; j < dimension; j++) {
      value_set_si(A->Lat->p[0][j], 0);
  }
  value_set_si(A->Lat->p[0][dimension], 1);

  A->P = Empty_Polyhedron(0);
  A->next = NULL;

  return (A);
} /* EmptyLBL */


/*
 * Given LBLs A and B, return True if A is included in B,
 * otherwise return False.
 */
Bool LBLIncludes(LBL *A, LBL *B)
{
  Bool ret = False;
  LBL *diff;

  // TODO: can we do better on ZDomains?

  diff = LBLDifference(A, B);
  if(isEmptyLBL(diff)) {
    ret = True;
  }
  LBLFree(diff);

  return ret;
} /* LBLIncludes */


/*
 * Print the contents of a single LBL 'A'
 */
static void sLBL_Print(FILE *fp, const char *format, LBL *A)
{
  if (A == NULL)
    return;
  fprintf(fp, "LBL: Dimension %d \n", A->Lat->NbRows - 1);

  if(emptyQ(A->P)) {
    fprintf(fp, "\n<empty>>\n");
  }
  else {
    fprintf(fp, "\nLATTICE: \n");
    Matrix_Print(fp, format, A->Lat);
    Polyhedron_Print(fp, format, A->P);
  }
  return;
} /* sLBL_Print */


/*
 * Print the contents of an LBL 'A'
 */
void LBLPrint(FILE *fp, const char *format, LBL *A)
{
  for( ; A; A = A->next) {
    sLBL_Print(fp, format, A);
    if(A->next)
      fprintf(fp, "\nUNION ");
  }
} /* LBLPrint */


/*
 * Return the LBL union of the LBLs 'A' and 'B'. The dimensions of
 * 'A' and 'B' must be equal.
 * All the LBLs of the resulting union are in Canonical form.
 */
LBL *LBLUnion(LBL *A, LBL *B)
{
  LBL *Result = NULL, *tmp;

  // copy A and B, concatenate, Canonicalize, and return :)

  Result = LBL_concatenate(LBLCopy(A), LBLCopy(B));

  CanonicalLBL(Result);

  return Result;
} /* LBLUnion */

/*
 * Return the intersection of the LBLs 'A' and 'B'.
 * The dimensions of 'A' and 'B' must be equal.
 */
LBL *LBLIntersection(LBL *A, LBL *B)
{
  LBL *Result = NULL, *tempA = NULL, *tempB = NULL;

  if (A->Lat->NbRows != B->Lat->NbRows) {
    errormsg1("LBLIntersection", "dimincomp", "incompatible dimensions between domains");
    return (NULL);
  }

  for (tempA = A; tempA; tempA = tempA->next) {
    for (tempB = B; tempB; tempB = tempB->next) {
      LBL *Inter;
      Inter = sLBL_Intersection(tempA, tempB);
      if(Inter) {
        Result = LBL_concatenate(Inter, Result);
      }
    }
  }

  if (!Result)
    return (EmptyLBL(A->Lat->NbRows - 1));

  CanonicalLBL(Result);

  return (Result);
} /* LBLIntersection */


/*
 * Return the difference of the LBLs 'A' - 'B' in canonical form.
 * The dimensions of 'A' and 'B' must be equal. Note that the
 * difference of two single LBLs can be a union of LBLs
 */
LBL *LBLDifference(LBL *A, LBL *B)
{ 
  LBL *res;

  if (A->Lat->NbRows != B->Lat->NbRows) {
    errormsg1("LBLDifference", "dimincomp", "incompatible dimensions between domains");
    return (NULL);
  }
  
  res = LBLCopy(A);
  // remove all single LBLs composing B from a copy of A:
  for (LBL *tempB = B; tempB; tempB = tempB->next) {
    LBL *tmp;
    tmp = LBL_sLBL_Difference(res, tempB);
    LBLFree(res);
    res = tmp;
  }

  if (!res)
    return (EmptyLBL(A->Lat->NbRows - 1));

  CanonicalLBL(res);

  return (res);
} /* LBLDifference */

/*
 * Return the image of the LBL 'A' under the affine transformation function
 * 'Func'. The number of columns of the function must be equal to the number
 * of rows in the matrix representing the lattice of 'A'.
 * Note:: Image((Z1 U Z2), F) = Image(Z1,F) U Image(Z2 U F).
 */
LBL *LBLImage(LBL *A, Matrix *Func)
{
  LBL *Result = NULL;

  for (LBL *temp = A; temp; temp = temp->next) {
    LBL *Im;
    Im = sLBL_Image(temp, Func);
    Result = LBL_concatenate(Im, Result);
  }
  if (Result == NULL)
    return EmptyLBL(A->Lat->NbRows - 1);

  CanonicalLBL(Result);

  return Result;
} /* LBLImage */

/*
 * Return the preimage of the Z-domain 'A' under the invertible, affine, ratio-
 * nal transformation 'Func'. The number of rows of the matrix representing
 * the function 'Func' must be equal to the number of rows of the matrix repr-
 * senting the lattice of 'A'.
 */
LBL *LBLPreimage(LBL *A, Matrix *Func) {

  LBL *Result = NULL;

  for (LBL *temp = A; temp; temp = temp->next) {
    LBL *B;
    B = sLBL_Preimage(temp, Func);
    Result = LBL_concatenate(B, Result);
  }

  if (Result == NULL)
    return (EmptyLBL(Func->NbColumns - 1));

  CanonicalLBL(Result);
  return Result;
} /* LBLPreimage */

/*
 * Return the LBL intersection of the single-LBLs 'A' and 'B'.
 * The result is always a single LBL, NULL if empty.
 * 
 * LInter is the intersection of the two lattices of A and B.
 * If LInter is empty, we return NULL.
 * Otherwise, we calculate PInter = intersection of the rational hulls of
 * A and B. We calculate P = Preimage of PInter by LInter and finally we
 * build the result LBL (Linter, P), in canonical form.
 *
 * USAGE: A and B's first Lattice considered only (no chained list),
 *        but can contain a polyhedral domain.
 */
static LBL *sLBL_Intersection(LBL *A, LBL *B) {

  LBL *Result = NULL;
  Matrix *LInter;
  Polyhedron *PInter, *ImageA, *ImageB, *PreImage;

  LInter = LatticeIntersection(A->Lat, B->Lat);

  if (isEmptyLattice(LInter)) {
    Matrix_Free(LInter);
    return (NULL);
  }

  ImageA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
  ImageB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
  PInter = DomainIntersection(ImageA, ImageB, MAXNOOFRAYS);
  // Although PInter can be an over approximation (a rational convex hull of
  // the resulting LBL), its preimage by LInter will constrain it to be the
  // right LBL, including potential "holes".
  // TODO: proof ^^
  // (??) does this work when there are columns of zeros in Linter?
  // will equalities be removed from the result?
  // if it does not, need to build explicitly.

  if (emptyQ(PInter))
    Result = NULL;
  else {
    PreImage = DomainPreimage(PInter, LInter, MAXNOOFRAYS);
    Result = LBLAlloc(LInter, PreImage);
    Domain_Free(PreImage);
  }

  Matrix_Free(LInter);
  Domain_Free(PInter);
  Domain_Free(ImageB);
  Domain_Free(ImageA);

  return (Result);
} /* sLBL_Intersection */


/*
 * Return the difference A - B
 * between a (union of) LBL(s) 'A' and a single LBL 'B'.
 * Algo: remove B from each part of A, and build a list of the result.
 *
 * USAGE: only the first lattice of B is considered
 *       (even if next is not NULL).
 * Creates a new allocated LBL, not necessarily in canonical form
 */
static LBL *LBL_sLBL_Difference(LBL* A, LBL* B)
{
  LBL *Result = NULL;

  for(LBL *z = A; z; z = z->next) {
    LBL *diff;

    diff = sLBL_Difference(z, B);
    // simple concatenate of diff and result (not canonical)
    Result = LBL_concatenate(diff, Result);
  }

  // Result contains every piece of the solution,
  // but it is not necessarily in canonical form (will be done be callee)
  return Result;
} /* LBL_sLBL_Difference */


/*
 * Return the difference of two single LBLs A and B.
 * Inspired from the method Gautam describes in his thesis,
 * modified to handle LBLs.
 * A and B are single LBLs, but the return value can be a union of LBLs!
 * Creates a new allocated LBL union
 *
 * USAGE: only the first lattice of A and B is considered (no union),
 *        but A and B can contain a coordinate polyhedral domain (in ->P).
 * Internal function, users should use LBLDifference.
 */

static LBL *sLBL_Difference(LBL* A, LBL* B)
{
  LBL *Result = NULL, *Final_Result; // U. of LBLs
  LBL *Ainter, *Binter; // single LBL
  LatticeUnion *LatDiff;
  Polyhedron *imA, *imB, *preimA, *ImDiff, *ImInter; // polyhedral domains

  if (A->Lat->NbRows != B->Lat->NbRows) {
    errormsg1("sLBL_Difference", "dimincomp", "incompatible dimensions");
    return(NULL);
  }

  // treat the simple case where the LBLs do not intersect
  Binter = LBLIntersection(A, B); // reused below
  if(isEmptyLBL(Binter)) {
    // if B does not intersect A, return A.
    #ifdef DIFFERENCE_DEBUG
      fprintf(stderr, "Binter=(A inter B) is empty, so B does not intersect A, we return A\n");
    #endif
    LBLFree(Binter);
    return(LBLCopy(A));
  }

  // Separate the computation in 3 phases:
  // 0. compute the difference of the image polyhedra P_A \ P_B (=temp) and
  //    add it to the solution LBL (using lattice L_A).
  //    This is an over-approximation of A (but not B)
  // 1. compute the rest where the intersection of P_A and P_B have same
  //    dimensions
  // 2. intersect the result with A to get rid of the over-approximations

  // [STEP 0 (includes Gautam's Step 2)]
  imA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
  imB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
  ImDiff = DomainDifference(imA, imB, MAXNOOFRAYS);
  #ifdef DIFFERENCE_DEBUG
    fprintf(stderr, "ImDiff (hull of A that does not cover B) = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, ImDiff);
  #endif

  // Add (A->Lat, A - hull(B)) to the result:
  if (!emptyQ(ImDiff)) {
    Polyhedron *RedPolyDiff;
    RedPolyDiff = DomainPreimage(ImDiff, A->Lat, MAXNOOFRAYS);
    // NOTICE: this is an over-approximation of A
    Result = LBLAlloc(A->Lat, RedPolyDiff);
    #ifdef DIFFERENCE_DEBUG
      fprintf(stderr, "Adding this to the temporary result: ");
      LBLPrint(stderr, P_VALUE_FMT, Result);
    #endif
    Domain_Free(RedPolyDiff);
  }

  // compute the images intersection of A and B
  ImInter = DomainIntersection(imA, imB, MAXNOOFRAYS);
  #ifdef DIFFERENCE_DEBUG
    fprintf(stderr, "ImInter (hull of A inter B) = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, ImInter);
  #endif
  
  // TODO: can be simplified, Ainter not really needed!

  // compute the part of A that intersects the hull of B in the image space
  preimA = DomainPreimage(ImInter, A->Lat, MAXNOOFRAYS);
  Ainter = LBLAlloc(A->Lat, preimA);
  // NOTICE: this Ainter can be a over-approximation of A

  Domain_Free(preimA);
  Domain_Free(ImDiff);
  Domain_Free(imA);
  Domain_Free(imB);

  // now Ainter and Binter have same lattices and polyhedra dimensions
  #ifdef DIFFERENCE_DEBUG
    fprintf(stderr, "-- [STEP1] now we compute the intersection on same lattice dimensions\n");
    fprintf(stderr, "Ainter = ");
    LBLPrint(stderr, P_VALUE_FMT, Ainter);
    fprintf(stderr, "and Binter = ");
    LBLPrint(stderr, P_VALUE_FMT, Binter);
  #endif

  // LatDiff (union of lattices) is the difference : (A->Lat) - (B->Lat) of same dimensions
  LatDiff = LatticeDifference(Ainter->Lat, Binter->Lat); 
  #ifdef DIFFERENCE_DEBUG
    if(!LatDiff)
      fprintf(stderr, "Empty Lattice difference\n");
  #endif

  // [STEP 1 of Gautam]:
  // Add all Z-polyhedra applying the (list of) lattice difference on ImInter
  for(LatticeUnion *tmp = LatDiff; tmp; tmp = tmp->next) {
    LBL *Ztmp;
    #ifdef DIFFERENCE_DEBUG
      fprintf(stderr, "Considering Lat diff: ");
      Matrix_Print(stderr, P_VALUE_FMT, tmp->M);
    #endif
    Ztmp = malloc(sizeof(*Ztmp));
    Ztmp->next = Result;
    Ztmp->Lat = tmp->M;
    Ztmp->P = DomainPreimage(ImInter, tmp->M, MAXNOOFRAYS);
    // NOTICE: this can be an over-approximation of the exact LBL

    Result = Ztmp;
  }
  // free LatticeUnion remaining memory (M has been reused as a lattice of Result)
  while(LatDiff) {
    LatticeUnion *next = LatDiff->next;
    free(LatDiff);
    LatDiff = next;
  }

  Domain_Free(ImInter);
  LBLFree(Ainter);
  LBLFree(Binter);

  if(!Result) {
    #ifdef DIFFERENCE_DEBUG
      fprintf(stderr, "-- result = (NULL)\n");
    #endif
    return(NULL);
  }

  #ifdef DIFFERENCE_DEBUG
    fprintf(stderr, "-- temporary over-approximation of result = ");
    LBLPrint(stderr, P_VALUE_FMT, Result);
  #endif
  // intersect the result with A to get the exact LBL.
  Final_Result = LBLIntersection(Result, A);
  LBLFree(Result);

  return(Final_Result);

} /* sLBL_Difference */


/*
 * Return the image of the single LBL 'A' under the affine function 'Func'
 * 
 * Algorithm:
 * - Multiply Lat by Func,
 * - Canonicalize the result (done by LBLAlloc)
 */
static LBL *sLBL_Image(LBL *A, Matrix *Func)
{
  Matrix *newL;
  LBL *result;

  if ((Func->NbColumns != A->Lat->NbRows)) {
    errormsg1("sLBL_Image", "dimincomp", "Incompatible dimensions");
    return NULL;
  }

  newL = Matrix_Alloc(Func->NbRows, A->Lat->NbColumns);
  Matrix_Product(Func, A->Lat, newL);
  result = LBLAlloc(newL, A->P);

  Matrix_Free(newL);
  return(result);
} /* sLBL_Image */

/*
 * Return the preimage of the single LBL 'Z' under an affine
 * transformation function 'G'. The number of rows of matrix 'G' must
 * be equal to the number of rows of the matrix representing the
 * lattice of Z.
 * Algorithm:
 * - build the LBL { z' | Lz = Gz', z \in A->P },
 * - remove z by normalizing the result
 */
static LBL *sLBL_Preimage(LBL *Z, Matrix *G)
{
  LBL *Result;
  Polyhedron *P, *newP;
  Matrix *Con;

  if(G->NbRows != Z->Lat->NbRows) {
    // G z' = L z
    errormsg1("sLBLPreimage", "dimincomp", "incompatible dimensions");
    return(NULL);
  }

  // d is the dimension of Z.
  // d' is the number of columns of G = the dimension of the result
  // build the Z-polyhedron = { z' | with P in dimension d + d'
  // such that L z = G z' }
  // then eliminate z by simplifying the result

  // the lattice is spreading z'
  // homogeneous d and d', homogeneous sum is d+d'-1
  int d = Z->Lat->NbRows;
  int dp = G->NbColumns;

  //          z'    z    cst
  // newL =   Id |  0   | 0 |
  //           0 |  0   | 1 |
  Matrix *newL = Matrix_Alloc(dp, d+dp-1);
  for(int i = 0; i < newL->NbRows; i++) {
    for(int j = 0; j < newL->NbColumns; j++) {
      if(j == i && i != newL->NbRows-1) { // left diagonal
        value_set_si(newL->p[i][j], 1);
      }
      else {
        value_set_si(newL->p[i][j], 0);
      }
    }
  }
  value_set_si(newL->p[newL->NbRows-1][newL->NbColumns-1], 1);

  // add the extra dimension on P
  newP = align_context(Z->P, d+dp-2, MAXNOOFRAYS);

  // build the extra constraint to be added to newP: G z' = L z
  // con =    0 |     |     |
  //          . |  G  | -L  | (g-l)
  //          0 |     |     |
  Con = Matrix_Alloc(d-1, d+dp-1+1);
  // copy G
  for(int i = 0; i < Con->NbRows; i++) {
    value_set_si(Con->p[i][0], 0); // equality
    for(int j = 0; j < dp-1; j++) {
      value_assign(Con->p[i][j+1], G->p[i][j]);
    }
    // constant
    value_assign(Con->p[i][Con->NbColumns-1], G->p[i][G->NbColumns-1]);
  }
  // copy NEGATIVE L
  for(int i = 0; i < Con->NbRows; i++) {
    for(int j = 0; j < d-1; j++) {
      value_oppose(Con->p[i][j+dp-1+1], Z->Lat->p[i][j]);
    }
    // substract constant l from g
    value_substract(Con->p[i][Con->NbColumns-1], Con->p[i][Con->NbColumns-1], Z->Lat->p[i][Z->Lat->NbColumns-1]);
  }

  P = DomainAddConstraints(newP, Con, MAXNOOFRAYS);
  Matrix_Free(Con);
  Domain_Free(newP);

  Result = LBLAlloc(newL, P);
  Domain_Free(P);
  Matrix_Free(newL);

  return(Result);
} /* sLBLPreimage */


// /*
//  * Given a Z-polyhedron 'A' in which the Lattice is not integral, return the
//  * Z-polyhedron which contains all the integral points in the input lattice.
//  */
// LBL *IntegraliseLattice(LBL *A) {

//   LBL *Result;
//   Matrix *M = NULL, *Id;
//   Polyhedron *Im = NULL, *Preim = NULL;

// #ifdef DOMDEBUG
//   FILE *fp;
//   fp = fopen("_debug", "a");
//   fprintf(fp, "\nEntered INTEGRALISELATTICE\n");
//   fclose(fp);
// #endif

//   Im = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
//   Id = Identity(A->Lat->NbRows);
//   M = LatticeImage(Id, A->Lat);
//   if (isEmptyLattice(M))
//     Result = EmptyLBL(A->Lat->NbRows - 1);
//   else {
//     Preim = DomainPreimage(Im, M, MAXNOOFRAYS);
//     Result = LBLAlloc(M, Preim);
//   }
//   Matrix_Free(M);
//   Domain_Free(Im);
//   Domain_Free(Preim);
//   return Result;
// } /* IntegraliseLattice */

// /*
//  * Return the simplified representation of the Z-domain 'ZDom'. It attempts to
//  * convexize unions of polyhedra when they correspond to the same lattices and
//  * to simplify union of lattices when they correspond to the same polyhdera.
//  */
// LBL *ZDomainSimplify(LBL *ZDom) {

//   LBL *Ztmp, *Result;
//   ForSimplify *Head, *Prev, *Curr;
//   LBL *ZDomHead, *Emp;

//   if (ZDom == NULL) {
//     fprintf(stderr, "\nError in ZDomainSimplify - ZDomHead = NULL\n");
//     return NULL;
//   }
//   if (ZDom->next == NULL)
//     return (LBL_Copy(ZDom));
//   Emp = EmptyLBL(ZDom->Lat->NbRows - 1);
//   ZDomHead = LBLUnion(ZDom, Emp);
//   LBL_Free(Emp);
//   Head = NULL;
//   Ztmp = ZDomHead;
//   do {
//     Polyhedron *Img;
//     Img = DomainImage(Ztmp->P, Ztmp->Lat, MAXNOOFRAYS);
//     for (Curr = Head; Curr != NULL; Curr = Curr->next) {
//       Polyhedron *Diff1;
//       Bool flag = False;

//       Diff1 = DomainDifference(Img, Curr->Pol, MAXNOOFRAYS);
//       if (emptyQ(Diff1)) {
//         Polyhedron *Diff2;

//         Diff2 = DomainDifference(Curr->Pol, Img, MAXNOOFRAYS);
//         if (emptyQ(Diff2))
//           flag = True;
//         Domain_Free(Diff2);
//       }
//       Domain_Free(Diff1);
//       if (flag == True) {
//         LatticeUnion *temp;

//         temp = malloc(sizeof(LatticeUnion));
//         temp->M = Matrix_Copy(Ztmp->Lat);
//         temp->next = Curr->LatUni;
//         Curr->LatUni = temp;
//         break;
//       }
//     }
//     if (Curr == NULL) {
//       Curr = malloc(sizeof(ForSimplify));
//       Curr->Pol = Domain_Copy(Img);
//       Curr->LatUni = malloc(sizeof(LatticeUnion));
//       Curr->LatUni->M = Matrix_Copy(Ztmp->Lat);
//       Curr->LatUni->next = NULL;
//       Curr->next = Head;
//       Head = Curr;
//     }
//     Domain_Free(Img);
//     Ztmp = Ztmp->next;
//   } while (Ztmp != NULL);

//   for (Curr = Head; Curr != NULL; Curr = Curr->next)
//     Curr->LatUni = LatticeSimplify(Curr->LatUni);
//   Result = NULL;
//   for (Curr = Head; Curr != NULL; Curr = Curr->next) {
//     LatticeUnion *L;
//     for (L = Curr->LatUni; L != NULL; L = L->next) {
//       Polyhedron *Preim;
//       LBL *Zpol;

//       Preim = DomainPreimage(Curr->Pol, L->M, MAXNOOFRAYS);
//       Zpol = LBLAlloc(L->M, Preim);
//       Zpol->next = Result;
//       Result = Zpol;
//       Domain_Free(Preim);
//     }
//   }
//   Curr = Head;
//   while (Curr != NULL) {
//     Prev = Curr;
//     Curr = Curr->next;
//     LatticeUnion_Free(Prev->LatUni);
//     Domain_Free(Prev->Pol);
//     free(Prev);
//   }
//   return Result;
// } /* ZDomainSimplify */

// /*
//  * 
//  *
// */

// LBL *SplitLBL(LBL *ZPol, Matrix *B) {

//   Matrix *H, *U1, *X, *Y;
//   LBL *zpnew, *Result;
//   LatticeUnion *Head = NULL, *tempHead = NULL;

// #ifdef DOMDEBUG
//   FILE *fp;
//   fp = fopen("_debug", "a");
//   fprintf(fp, "\nEntered SplitLBL \n");
//   fclose(fp);
// #endif

//   if (B->NbRows != B->NbColumns) {
//     fprintf(
//         stderr,
//         "\n SplitLBL : The Input Matrix B is not a proper Lattice \n");
//     return NULL;
//   }

//   if (ZPol->Lat->NbRows != B->NbRows) {
//     fprintf(stderr,
//             "\nSplitLBL : The Lattice in LBL and B have ");
//     fprintf(stderr, "incompatible dimensions \n");
//     return NULL;
//   }

//   if (isNormalLattice(ZPol->Lat) != True) {
//     AffineHermite(ZPol->Lat, &H, &U1);
//     X = Matrix_Copy(H);
//     Matrix_Free(U1);
//     Matrix_Free(H);
//   } else
//     X = Matrix_Copy(ZPol->Lat);

//   if (isNormalLattice(B) != True) {
//     AffineHermite(B, &H, &U1);
//     Y = Matrix_Copy(H);
//     Matrix_Free(H);
//     Matrix_Free(U1);
//   } else
//     Y = Matrix_Copy(B);
//   if (isEmptyLattice(X)) {
//     return NULL;
//   }

//   Head = Lattice2LatticeUnion(X, Y);

//   /* If the spliting operation can't be done the result is the original
//    * Zplyhedron. */

//   if (Head == NULL) {
//     Matrix_Free(X);
//     Matrix_Free(Y);
//     return LBL_Copy(ZPol);
//   }

//   Result = NULL;

//   while (Head) {
//     tempHead = Head;
//     Head = Head->next;
//     zpnew = LBLAlloc(tempHead->M, ZPol->P);
//     Result = ZDconcatenate(zpnew, Result);
//     tempHead->next = NULL;
//     Matrix_Free(tempHead->M);
//     free(tempHead);
//   }
//   Matrix_Free(X);
//   Matrix_Free(Y);
//   return Result;
// } /* SplitLBL */

/*
 * get the matrix of equalities from a polyhedron
 * (without the first columns of 0's)
 */
static Matrix *get_equalities(Polyhedron *P)
{
  // Eq is the matrix of equations of P (including the constant)
  Matrix* Eq = Matrix_Alloc(P->NbEq, P->Dimension+1);
  // get equalities (first rows of P->Constraint)
  for(int i=0; i<P->NbEq; i++) {
    for(int j=0; j<Eq->NbColumns; j++) {
      value_assign(Eq->p[i][j], P->Constraint[i][j+1]);
    }
  }
  return (Eq);
} /* get_equalities */

/*
 * compare a matrix of equalities to the one of a polyhedron P
 */
static Bool same_equalities(Matrix *Eq, Polyhedron *P)
{
  if(P->NbEq != Eq->NbRows)
    return (False);

  for(int i=0; i<P->NbEq && i<Eq->NbRows; i++) {
    for(int j=0; j<Eq->NbColumns; j++) {
      if(value_ne(Eq->p[i][j], P->Constraint[i][j+1]))
        return (False);
    }
  }
  return (True);
} /* same_equalities */

/*
 * Try to remove the equalities from A->P in the single LBL A.
 * In place. A->P is a domain.
 */
static void sLBL_Remove_Equalities(LBL *A, Matrix *Equalities)
{
  if (A->P->Dimension > 0 && A->P->NbEq != 0) {
    Matrix *ker=NULL, *H = NULL, *NewL;
    Matrix *eq_hermite = NULL, *eq_U = NULL;

    // remove equalities in domain P and change Lat to spread the original space
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "P has equalities\n");
      fprintf(stderr, "Equality matrix (including constants): ");
      Matrix_Print(stderr, P_VALUE_FMT, Equalities);
    #endif

    // compute the kernel of Equalities matrix.
    left_hermite(Equalities, &eq_hermite, NULL, &eq_U);
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "hermite eq = ");
      Matrix_Print(stderr, P_VALUE_FMT, eq_hermite);
      fprintf(stderr, "hermite U = ");
      Matrix_Print(stderr, P_VALUE_FMT, eq_U);
      #endif
    // the kernel of Equalities is the last n-NbEq columns of eq_U
    Matrix_subMatrix(eq_U, 0, A->P->NbEq, eq_U->NbRows, eq_U->NbColumns, &ker);
    Matrix_Free(eq_hermite);
    Matrix_Free(eq_U);

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "ker of eq: ");
      Matrix_Print(stderr, P_VALUE_FMT, ker);
    #endif
    
    AffineHermite(ker, &H, NULL);
    Matrix_Free(ker);
    // We know that: Eq . H = Eq . Ker(Eq) = 0

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "Matrix H: ");
      Matrix_Print(stderr, P_VALUE_FMT, H);
      fprintf(stderr, "Lattice of A: ");
      Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
    #endif

    // TODO: complete the matrix such that it is full row
    //       with (0..0 1 0..0)^T columns on the right ???
    // is that right? is that enough?

    
    // TODO: CHECK THAT THIS IS CORRECT!
    // if the bottom right value of H is not one, this means that
    // the transformation matrix is not integer but rational.
    // just make it integer to eliminate rational points.
    // (see ZImPre3 for an example where this is necessary)
    value_set_si(H->p[H->NbRows-1][H->NbColumns-1], 1);


    // #ifdef CANONICAL_DEBUG
    // {
    //   Value det;
    //   value_init(det);
    //   value_set_si(det, 1);
    //   // check if the invariant factor is 1
    //   for (int col = 0; col < H->NbColumns - 1; col++) {
    //     int lin;
    //     for(lin = 0; lin<H->NbRows-1; lin++) {
    //       if(value_notzero_p(H->p[lin][col])) {
    //           value_multiply(det, det, H->p[lin][col]);
    //         break;
    //       }
    //     }
    //   }
    //   fprintf(stderr, "invariant factor = ");
    //   value_print(stderr, P_VALUE_FMT, det);
    //   fprintf(stderr, "\n");
    //   value_clear(det);
    // }
    // #endif


    // NewL = L . H
    NewL = Matrix_Alloc(A->Lat->NbRows, H->NbColumns);
    Matrix_Product(A->Lat, H, NewL);
    // NewP = H^{-1} . P
    Polyhedron* NewP = DomainPreimage(A->P, H, MAXNOOFRAYS);   // TODO: wrong since H is not unimodular!

    // update A
    Domain_Free(A->P);
    Matrix_Free(A->Lat);

    A->P = NewP;
    A->Lat = NewL;
    Matrix_Free(H);
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "New Lat: ");
      Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
      fprintf(stderr, "New P: ");
      Polyhedron_Print(stderr, P_VALUE_FMT, A->P);
    #endif
  }
  else { // P contains no equalities
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "P has no equalities.\n");
    #endif
  }
} /* sLBL_Remove_Equalities */

/*
 * Set the affine function A->Lat to normal form in single LBL 'A'.
 * In place. A->P is a domain.
 */
static void sLBL_Lat_Normalize(LBL *A)
{
  // check if A->Lat is in Hermite form
  if(!isNormalLattice(A->Lat)) {
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "A is not HNF\n");
    #endif
    Matrix* U = NULL;
    Matrix* H = NULL;

    // for left hermite to include the constant:
    Matrix_Move_Homogeneous_Dim_First(A->Lat);
    // to compute HNF of the lattice (constant part left-top)
    // We will use U = Q^{-1}, such that LU = H.
    left_hermite(A->Lat, &H, NULL, &U);
    // Move the constant back to right-bottom
    Matrix_Move_Homogeneous_Dim_Last(H);
    Matrix_Move_Homogeneous_Dim_Last(U);

    // set the new Lat matrix as H
    Matrix_Free(A->Lat);
    A->Lat = H;

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "New Lat (HNF): ");
      Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
    #endif

    // Now update of A->P using the preimage by U (unimodular)
    Polyhedron *NewP = DomainPreimage(A->P, U, MAXNOOFRAYS);
    Domain_Free(A->P);
    A->P = NewP;
    Matrix_Free(U);

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "New P: ");
      Polyhedron_Print(stderr, P_VALUE_FMT, A->P);
    #endif
  } // A->Lat in canonical form
  else {
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "A is HNF.\n");
    #endif
  }
} /* sLBL_Lat_Normalize */


/*
 * A is empty, check if A has no successor in the list and is canonical.
 * Can relink A->next to current A if needed, and returns True.
 * returns False if A->next did not change.
 * In place.
 */
static Bool sLBL_Verify_Empty(LBL *A)
{
  // A is empty
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Empty LBL\n");
  #endif

  // if there is something linked to an empty LBL, need to replace the
  // current LBL with the next LBL: replace A with next and free A->next
  if(A->next) {
    LBL *nextA = A->next;

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "... But the next one is not empty, relinking\n");
    #endif
    Domain_Free(A->P);
    Matrix_Free(A->Lat);
    A->P = nextA->P;
    A->Lat = nextA->Lat;
    A->next = nextA->next;
    free(nextA);

    // A changed
    return(True);
  }

  // A is empty and alone.
  // Verify that it is canonical and return.
  if(A->P->Dimension > 0) {
    int dimension = A->Lat->NbRows;
    Domain_Free(A->P);
    Matrix_Free(A->Lat);
    A->P = Empty_Polyhedron(0);
    A->Lat = Matrix_Alloc(dimension, 1);
    for(int j=0 ; j < dimension-1; j++) {
      value_set_si(A->Lat->p[0][j], 0);
    }
    value_set_si(A->Lat->p[0][dimension-1], 1);
  }

  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Returning, A is empty: ");
    LBLPrint(stderr, P_VALUE_FMT, A);
  #endif

  return(False);
}


/*
 * Remove the columns of zeros from A->Lat.
 * In place. A->P is a domain.
 */
static void sLBL_Lat_Remove_Zeros(LBL *A)
{
  // TODO: this is the naive version,
  // need to consider integer existential variable elimination and dark shadow!

  Polyhedron *NewP;
  int nbZeros = count_zeroCols(A->Lat);
  if(nbZeros) {
    Matrix* Transformation = Matrix_Alloc(A->Lat->NbColumns-nbZeros, A->Lat->NbColumns);
    for (int  i = 0; i < Transformation->NbRows; i++) {
      for (int j = 0; j < Transformation->NbColumns; j++) {
        if(i==j && i!=Transformation->NbRows-1) {
          value_set_si(Transformation->p[i][j], 1);
        }
        else {
          value_set_si(Transformation->p[i][j], 0);
        }
      }
    }
    value_set_si(Transformation->p[Transformation->NbRows-1][Transformation->NbColumns-1], 1);

    NewP = DomainImage(A->P, Transformation, MAXNOOFRAYS);
    Domain_Free(A->P);
    A->P = NewP;
    Matrix_Free(Transformation);
    Matrix* NewL = Matrix_Alloc(A->Lat->NbRows,A->Lat->NbColumns-nbZeros);
    for (int  i = 0; i < NewL->NbRows; i++) {
      for (int j = 0; j < NewL->NbColumns; j++) {
        if(j < NewL->NbColumns-1) {
          value_assign(NewL->p[i][j],A->Lat->p[i][j]);
        }else{
          value_assign(NewL->p[i][j],A->Lat->p[i][A->Lat->NbColumns-1]);
        }
      }
    }
    Matrix_Free(A->Lat);
    A->Lat = NewL;
    #ifdef CANONICAL_DEBUG
      LBLPrint(stderr, P_VALUE_FMT, A);
    #endif
  }
} /* sLBL_Lat_Remove_Zeros */


/*
 * The function takes the head of LBL 'A'
 * --- single lattice function, domain (list of polyhedra) ---
 * and modifies it in place to be in canonical form (A->Lat in HNF and no
 * equalities in A->P)
 * IN PLACE: modifies A itself
 * 
 * WARNING: this function modifies the head single LBL of A to build a union:
 * it may add or remove the LBLs stored in A->next
 */
static void sLBL_Canonical(LBL* A)
{
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Entering sLBL_Canonical\n");
    fprintf(stderr, "--------- Input Lat: ");
    Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
    fprintf(stderr, "--------- Input P: ");
    Polyhedron_Print(stderr, P_VALUE_FMT, A->P);
  #endif
  
  if (A->P->Dimension+1 != A->Lat->NbColumns) {
    errormsg1("sLBL_Canonical", "dimincomp", "incompatible dimensions");
    return;
  }

  if(emptyQ(A->P)) {
    if(sLBL_Verify_Empty(A)) {
      // head was empty and has been replaced (with next).
      // need to canonicalize the (new) current LBL itself
      // (since the caller will not rescan it!):
      sLBL_Canonical(A);
    }
  }

  // change P such that all polyhedra in this domain have the same set of
  // equalities, that is, the equalities of the first one.
  // all the other ones are added to a new LBL, linked to LBL A (in A->next)
  // A->next will be treated at next step by the caller
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Checking for equalites in P\n");
  #endif
  LBL *new = NULL;
  Matrix * Equalities = get_equalities(A->P); // get eq from the first one
  Polyhedron *nextpp, *prevpp = A->P; // keep a ref to the previous to relink

  for(Polyhedron *pp = A->P->next; pp; prevpp = pp, pp = nextpp) {
    // check that the equalities of pp->Constraints are the same as the ones
    // of matrix Equalities.
    if(!same_equalities(Equalities, pp)) {
      // if not, get pp out.
      if(!new) {
        new = malloc(sizeof(*new));
        if (!new) {
          errormsg1("sLBL_Canonical", "outofmem", "Out of Memory");
          return;
        }
        new->P = NULL;
        new->Lat = Matrix_Copy(A->Lat);
      }
      // remove pp from the list A->P, and get the right next iteration
      nextpp = prevpp->next = pp->next;
      // add pp to new->P
      pp->next = new->P;
      new->P = AddPolyToDomain(pp, new->P);
    }
    else {
      nextpp = pp->next; // next polyhedron of the domain
    }
  }
  if(new) {
    // include new in the LBL list A
    new->next = A->next;
    A->next = new;
  }
  // Now all polyhedra of domain A->P have the same equalities

  // We can try to remove some equalities from A->P
  sLBL_Remove_Equalities(A, Equalities);
  Matrix_Free(Equalities);

  // check that the result is not empty after removing equalities
  // (example ZAlloc1b.in needs that)
  if(emptyQ(A->P)) {
    if(sLBL_Verify_Empty(A)) {
      // head was empty and has been replaced (with next).
      // need to canonicalize the (new) current LBL itself
      // (since the caller will not rescan it!):
      sLBL_Canonical(A);
    }
  }

  // Normalize the affine lattice A->Lat
  sLBL_Lat_Normalize(A);

  // Remove the columns of zeros from A->Lat
  // sLBL_Lat_Remove_Zeros(A);

  return;
} /* sLBL_Canonical */


/*
 * The function takes an LBL 'A' and transforms it into its canonical form:
 *  - all lattices in HNF, and
 *  - no equalities in all polyhedral domains.
 * Performs the operation IN PLACE (modifies A)
 * 
 * USAGE NOTICE:
 * If A is a real LBL and not a Z-polyhedron, there will remain zero-columns
 * on the right of the lattice matrice(s).
 * To transform a union of LBLs into a union of Z-polyhedra, you need to
 * call LBL2ZDomain().
 */
void CanonicalLBL(LBL *A) {

  // here, just transform every LBL of the list individually
  // careful, this may add a new LBL to the list A itself
  // (after tmp) but they will be scanned by this loop :)
  for(LBL *tmp = A; tmp; tmp = tmp->next) {
    sLBL_Canonical(tmp);
  }

  // check if a Lat is present twice in A, and if it is, union this
  // polyhedron to the existing one and remove the second reference
  for(LBL *tmp = A; tmp; tmp = tmp->next) {
    LBL *ZZ;
    if((ZZ = FindLatticePred(tmp->Lat, tmp))) {
      LBL *remove;
      Polyhedron *nextpp;
      // add all polyhedra of the domain ZZ->next->P to tmp
      // consumes ZZ->next->P.
      Polyhedron *pp = ZZ->next->P;
      while(pp) {
        nextpp = pp->next;
        pp->next = NULL;
        // this consumes pp, so need to get next before
        tmp->P = AddPolyToDomain(pp, tmp->P);
        pp = nextpp;
      }
      // remove ZZ->next by changing the linked list
      remove = ZZ->next;
      ZZ->next = ZZ->next->next;
      Matrix_Free(remove->Lat);
      // remove->P has been reused
      free(remove);
    }
  }
} /* CanonicalLBL */


/*
 * Find if a given lattice is present in a LBL.
 * Returns the address of the *previous* LBL
 * (such that ZZ->next->Lat == L).
 * NULL if not found.
 */
static LBL *FindLatticePred(Matrix *L, LBL *A) {
  LBL* tmp;

  for(tmp = A; tmp->next; tmp=tmp->next) {
    if(sameLattice(L, tmp->next->Lat)) {
      return (tmp);
    }
  }
  return (NULL);
} /* FindLatticePred */


/*
 * count the number of columns of zeros on the right of the linear part
 * of a lattice function
 */
static int count_zeroCols (Matrix* M)
{
  int nb=0;
  for (int j = M->NbColumns-2 ; j >= 0; j--) {
    Bool isZero=True;
    for (int i = 0; i < M->NbRows; i++) {
      if (value_notzero_p(M->p[i][j])){
        isZero=False;
        break;
      }
    }
    if(!isZero)
      break;
    nb++;
  }
  return nb;
}

/*
 * Build a union of Z-domains from a union of LBLs.
 * 
 * The resulting Z-domain union lattices do not contain any zero columns.
 * 
 * The result is proved to be finite, but can pretty easily explode in
 * complexity, for example with a very pointy initial LBL
 * -- take something based on (2^16 x - (2^16-1) y) for example.
 */
LBL *LBL2ZDomain(LBL *A)
{
  // TODO: eliminate existential variables :)

  return (NULL);
}

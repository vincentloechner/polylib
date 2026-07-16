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


// debug all functions:
// #define DEBUG
#ifdef DEBUG
#define CANONICAL_DEBUG 1
#define INTERSECTION_DEBUG 1
#define DIFFERENCE_DEBUG 1
#define LBLDIFF_DEBUG 1
#define COMP_DEBUG 1
#define HOLES_DEBUG 1
#define SIMPLIFY_DEBUG 1
#define SIMPLIFY2_DEBUG 1
#define IMAGE_DEBUG 1
#endif

static LBL *LBLConcatenate(LBL *A, LBL *B);
static LBL *sLBLIntersection(LBL *, LBL *);
static LBL *sLBLCopy(LBL *A);
static void sLBLFree(LBL *L);
static LBL *sLBLComplement(LBL *A);
static Bool LBL_simple_inclusion_check(LBL *A, LBL *B);
// static LBL *sLBL_Difference(LBL *, LBL *);
static LBL *sLBLImage(LBL *, Matrix *);
static LBL *sLBLPreimage(LBL *, Matrix *);
static void sLBLCanonical(LBL *A);
// static LBL *LBL_sLBL_Difference(LBL *A, LBL *B);
static Polyhedron *sLBLCompute_holes(LBL *A, Polyhedron **pExact);
static Bool polyhedron_int_solution(Polyhedron *scan, Value *val, int position);
static Polyhedron *Domain_Remove_Integer_Empty(Polyhedron *D);
static Polyhedron *domain_project(Polyhedron *P, int eliminate);
static Polyhedron *domain_insert_dim(Polyhedron *D, int move);
static void LBL_Remove_Empty(LBL *A);
static void sLBLMake_lattice_equal_to(LBL *A, Matrix *ref);


/*
 * Returns True if 'A' is empty, otherwise returns False.
 * A can be a non-simplified list of empty LBLs
 */
Bool isEmptyLBL(LBL *A)
{
  if(A == NULL)
    return True;
  if(A->P == NULL)
    return True;
  if(emptyQ(A->P)) {
    // this one is empty.
    // return the emptiness of next
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

  if(Domain && (Lat->NbColumns != Domain->Dimension + 1)) {
    errormsg1("LBLAlloc", "dimincomp",
      "the Lattice and the Polyhedron are not compatible to form a LBL");
    return NULL;
  }
  if(Domain && emptyQ(Domain)) {
    Domain = NULL;
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
static void sLBLFree(LBL *L)
{
  if (L == NULL)
    return;
  if(L->Lat)
    Matrix_Free(L->Lat);
  if(L->P)
    Domain_Free(L->P);
  free(L);
  return;
} /* sLBLFree */


/*
 * Free the memory used by the LBL 'L'
 */
void LBLFree(LBL *L)
{
  if (L == NULL)
    return;
  LBLFree(L->next);
  sLBLFree(L);
} /* LBLFree */


/*
 * Return a copy of the single LBL 'A'.
 * Internal function, users should use LBLCopy.
 */
static LBL *sLBLCopy(LBL *A)
{
  return (LBLAlloc(A->Lat, A->P));
} /* sLBLCopy */


/*
 * Return a copy of the LBL 'L'
 */
LBL *LBLCopy(LBL *L)
{
  LBL *copy;
  copy = sLBLCopy(L);

  if (L->next != NULL)
    copy->next = LBLCopy(L->next);
  return copy;
} /* LBLCopy */


/*
 * Concatenate the LBLs 'A' and 'B',
 * and return a pointer to the new LBL.
 * Consumes the memory of A and of B (no need to free) to build
 * the result. Internal function only, do not use to build unions.
 * No simplification, just join them.
 */
static LBL *LBLConcatenate(LBL *A, LBL *B)
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

  // go to the end of A and link B there:
  for(tmp = A; tmp->next; tmp = tmp->next)
    ;
  tmp->next = B;
  
  return (A);
} /* LBLConcatenate */


/*
 * Return the empty LBL of dimension 'dimension'
 * Lat = (0 ... 0 1)^T
 * P = empty polyhedron of dimension 0
 */
LBL *EmptyLBL(int dimension)
{
  LBL *A;

  A = malloc(sizeof(LBL));
  if(!A) {
    errormsg1("EmptyLBL", "outofmem", "Out of Memory");
    return(NULL);
  }

  A->Lat = Matrix_Alloc(dimension+1, 1);
  for(int j = 0 ; j < dimension; j++) {
      value_set_si(A->Lat->p[0][j], 0);
  }
  value_set_si(A->Lat->p[0][dimension], 1);

  // A->P = Empty_Polyhedron(0);
  A->P = NULL;
  A->next = NULL;

  return (A);
} /* EmptyLBL */


/*
 * Return the universe Z-polyhedron of dimension 'dimension'
 * Lat = Id
 * P = Universe_Polyhedron()
 */
LBL *UniverseLBL(int dimension)
{
  LBL *A;

  A = malloc(sizeof(LBL));
  if(!A) {
    errormsg1("EmptyLBL", "outofmem", "Out of Memory");
    return(NULL);
  }
  A->Lat = NULL;
  Matrix_identity(dimension + 1, &(A->Lat));
  A->P = Universe_Polyhedron(dimension);
  A->next = NULL;

  return (A);
} /* UniverseLBL */


/*
 * Simple inclusion test, return True if all sLBLs of A are also part of B.
 *
 * Check if all sLBLs of A are present in B (same Lat matrix, domain covered)
 */
static Bool LBL_simple_inclusion_check(LBL *A, LBL *B)
{
  // check if all pieces of A are also present in B
  for( ; A; A = A->next) {
    LBL *tmpB;
    Polyhedron *ddiff;

    // For now this test is very rudimentary!
    // 
    // TODO: the lattice search could search for a lattice included in the
    // other, not necessarily equal! But then we have to transform the
    // coordinate polyhedra to be in the same space. Other parts of A could
    // also be covered by other parts of B.

    for(tmpB = B; tmpB; tmpB = tmpB->next) {
      if(isEqualLattice(A->Lat, tmpB->Lat))
        break; // found
    }
    if(! tmpB)
      return(False);  // did not find A->Lat in B
  

    // check if A->P is included in tmpB->P
    if(! tmpB->P)
      return(False);  // B is empty

    // if tmpB->P is a single polyhedron, a simple inclusion check of all
    // parts of A->P is enough:
    if(tmpB->P->next == NULL) {
      for(Polyhedron *Apart = A->P; Apart; Apart = Apart->next) {
        if(! PolyhedronIncludes(tmpB->P, Apart))
        {
          // tmpB->P does not include (cover) Apart
          return(False);
        }
      }
    }
    else {
      // if tmpB->P is a domain, need to explicitly compute the difference and
      // check its emptiness:
      // (A->P - tmpB->P) should be empty.
      ddiff = DomainDifference(A->P, tmpB->P, MAXNOOFRAYS);
      if(! emptyQ(ddiff)) {
        Domain_Free(ddiff);
        return(False);  // not included
      }
      Domain_Free(ddiff);
    }
    // success for this part of A.
  } // next part of A

  // every part of A was found in B
  return(True);
}


/*
 * Given LBLs A and B, return True if A is included in B,
 * otherwise return False.
 */
Bool LBLIncluded(LBL *A, LBL *B)
{
  Bool ret = False;
  LBL *diff;

  if(LBL_simple_inclusion_check(A, B)) {
    return(True);
  }

  // Could we do better on ZDomains?
  // the answer is no: the complicated part of the difference computation is
  // not executed anyway when the input LBLs are ZDomains.

  diff = LBLDifference(A, B);

  // check if diff is integer empty
  // search an integer solution, stop as soon as found
  for(LBL *tmp = diff; tmp; tmp = tmp->next) {
    if((tmp->P = Domain_Remove_Integer_Empty(tmp->P))) {
      break;
    }
  }

  if(isEmptyLBL(diff)) {
    ret = True;
  }
  LBLFree(diff);

  return ret;
} /* LBLIncluded */


/*
 * Return True if the point pt is included in LBL A,
 * otherwise return False.
 * 
 * pt should be an integer array of Values of dimension (A->Lat->NbRows - 1)
 */
Bool LBLContainsPoint(LBL *A, Value *pt)
{

  // scan each simple LBL of A:
  for( ; A ; A = A->next) {
    // build the LBL { AL z |  AL z = pt, z \in AP } by adding the equalities
    // {AL z = pt} to AP, normalize the LBL,
    // and check if the resulting LBL has a solution

    int dimLBL = A->Lat->NbRows - 1;
    int hdimP = A->P->Dimension + 1;  // homogeneous dimension of A->P
                                      // == A->Lat->NbColumns
    LBL *inter;
    Matrix *Eq;
    Polyhedron *newP;
    Bool empty;

    // Eq = [ 0 |  Al  |  Ac-pt ]
    // with A = [ Al | Ac ]  (linear part/constant part)
    Eq = Matrix_Alloc(dimLBL, hdimP + 1);
    for(int d = 0; d < dimLBL; d++) {
      value_set_si(Eq->p[d][0], 0);
      Vector_Copy(A->Lat->p[d], &Eq->p[d][1], hdimP);
      value_substract(Eq->p[d][hdimP], Eq->p[d][hdimP], pt[d]);
    }

    // build new LBL inter
    newP = DomainAddConstraints(A->P, Eq, MAXNOOFRAYS);
    inter = LBLAlloc(A->Lat, newP);

    // simplify and check emptiness
    LBLSimplifyEmpty(inter);
    empty = isEmptyLBL(inter);

    // free memory
    LBLFree(inter);
    Domain_Free(newP);
    Matrix_Free(Eq);

    // early exit if found
    if(!empty) {
      return True;
    }
  }  // continue to next LBL of A

  return False;
} /* LBLContainsPoint */


/*
 * Print the contents of a single LBL 'A'
 */
static void sLBLPrint(FILE *fp, const char *format, LBL *A)
{
  fprintf(fp, "LBL: Dimension %d\n", A->Lat->NbRows - 1);
  if(!A->P ) {
    fprintf(fp, "\n<empty>\n");
  }
  else {
    fprintf(fp, "\nLATTICE:\n");
    Matrix_Print(fp, format, A->Lat);
    Polyhedron_Print(fp, format, A->P);
  }
} /* LBLPrint */


/*
 * Print the contents of an LBL 'A'
 */
void LBLPrint(FILE *fp, const char *format, LBL *A)
{
  if(!A) {
    fprintf(fp, "\n<empty>\n");
    return;
  }

  sLBLPrint(fp, format, A);
  for( A = A->next; A; A = A->next) {
    fprintf(fp, "\nUNION ");
    sLBLPrint(fp, format, A);
  }
} /* LBLPrint */


/*
 * Return the LBL union of the LBLs 'A' and 'B'. The dimensions of
 * 'A' and 'B' must be equal.
 * All the LBLs of the resulting union are in Canonical form.
 */
LBL *LBLUnion(LBL *A, LBL *B)
{
  LBL *Result = NULL;

  // copy A and B, concatenate, Canonicalize, and return :)

  Result = LBLConcatenate(LBLCopy(A), LBLCopy(B));

  CanonicalLBL(Result);

  return Result;
} /* LBLUnion */

/*
 * Return the intersection of the LBLs 'A' and 'B'.
 * The dimensions of 'A' and 'B' must be equal.
 * 
 * Algorithm:
 * intersect each piece of A with each piece of B and union all results
 */
LBL *LBLIntersection(LBL *A, LBL *B)
{
  LBL *Result = NULL, *tempA, *tempB;

  if (A->Lat->NbRows != B->Lat->NbRows) {
    errormsg1("LBLIntersection", "dimincomp",
      "incompatible dimensions between domains");
    return (NULL);
  }

  for (tempA = A; tempA; tempA = tempA->next) {
    for (tempB = B; tempB; tempB = tempB->next) {
      LBL *Inter;
      Inter = sLBLIntersection(tempA, tempB);
      if(Inter) {
        Result = LBLConcatenate(Inter, Result);
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
 * difference of two single LBLs can be a union of LBLs.
 * 
 * Algorithm:
 * Let I = A inter B
 * - check if I is empty -> result = A
 * - check if A included in I -> result = empty
 * - general case:
 *   successively remove each single LBL of I from A, return the rest.
 * 
 */
LBL *LBLDifference(LBL *A, LBL *B)
{ 
  LBL *inter, *res = NULL;
  
  #ifdef LBLDIFF_DEBUG
  fprintf(stderr, "-----Entering LBLDiff-----\n");
  fprintf(stderr, "---- A = ");
  LBLPrint(stderr, P_VALUE_FMT, A);
  fprintf(stderr, "---- B = ");
  LBLPrint(stderr, P_VALUE_FMT, B);
  #endif

  if(!A) {
    return(NULL);
  }
  if(!B) {
    return(LBLCopy(A));
  }
  if (A->Lat->NbRows != B->Lat->NbRows) {
    errormsg1("LBLDifference", "dimincomp",
        "incompatible dimensions between domains");
    return (NULL);
  }

  // compute the intersection of A and B
  inter = LBLIntersection(A, B);
  if(isEmptyLBL(inter)) {
    LBLFree(inter);
    return(LBLCopy(A));
  }

  // simple check if A is included in inter, then inter == A and result = empty.
  if(LBL_simple_inclusion_check(A, inter)) {
    LBLFree(inter);
    return(EmptyLBL(A->Lat->NbRows - 1));
  }

  // compute res = A - inter
  // initialize result: A
  res = A;
  // remove all single LBLs composing inter from res:
  for (LBL *tmpi = inter; tmpi; tmpi = tmpi->next) {
    LBL *diff, *comp;

    // compute res - tmpi:
    comp = sLBLComplement(tmpi);
    diff = LBLIntersection(res, comp);
    LBLFree(comp);

    if(res != A)
      LBLFree(res);        // free previous res
    res = diff;            // new res = diff

    // early exit if empty
    if(isEmptyLBL(res))
      break;
  }

  LBLFree(inter);

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
    Im = sLBLImage(temp, Func);
    Result = LBLConcatenate(Im, Result);
  }
  if (Result == NULL)
    return EmptyLBL(A->Lat->NbRows - 1);

  CanonicalLBL(Result);

  return Result;
} /* LBLImage */

/*
 * Return the preimage of the LBL 'A' under the affine transformation 'Func'. 
 * The number of rows of the matrix representing the function 'Func' must be
 *  equal to the number of rows of the matrix representing the lattice of 'A'.
 */
LBL *LBLPreimage(LBL *A, Matrix *Func) {

  LBL *Result = NULL;

  for (LBL *temp = A; temp; temp = temp->next) {
    LBL *B;
    B = sLBLPreimage(temp, Func);
    Result = LBLConcatenate(B, Result);
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
 * USAGE: A and B's first Lattice considered only (no chained list),
 *        but can contain a polyhedral domain in ->P.
 * 
 * Algorithm:
 * - IF the input LBLs are Z-Polyhedra, we can simply compute:
 *    LInter is the intersection of the two lattices AL and BL.
 *    Compute PI = the intersection of the images of AP by AL, and BP by BL
 *    Build the result Z-polyhedron (Linter, preimage of PI by Linter)
 * - ELSE:
 *    build explicit equalities between points of A and B and simplify.
 *      build the LBL { AL z |  BL z' = AL z, z \in AP, z' \in BP },
 *      and remove z' by normalizing the result
 */
static LBL *sLBLIntersection(LBL *A, LBL *B)
{
  LBL *Result = NULL;
  Matrix *LInter;
  Polyhedron *PInter, *ImageA, *ImageB, *PreImage;

  #ifdef INTERSECTION_DEBUG
    fprintf(stderr, "-- Entering sLBLIntersection\nA = ");
    sLBLPrint(stderr, P_VALUE_FMT, A);
    fprintf(stderr, "B = ");
    sLBLPrint(stderr, P_VALUE_FMT, B);
  #endif
  if(isEmptyLBL(A) || isEmptyLBL(B))
  {
    return(EmptyLBL(A->Lat->NbRows - 1));
  }
  LInter = LatticeIntersection(A->Lat, B->Lat);
  if(isEmptyLattice(LInter)) {
    Matrix_Free(LInter);
    #ifdef INTERSECTION_DEBUG
    fprintf(stderr, "Empty Lattice intersection, result = <empty>\n");
    fprintf(stderr, "-- exit sLBLIntersection\n");
    #endif
    return (NULL);
  }
  #ifdef INTERSECTION_DEBUG
  // fprintf(stderr, "Lattice intersection = LInter = ");
  // Matrix_Print(stderr, P_VALUE_FMT, LInter);
  #endif

  if(LatCountZeroCols(A->Lat) == 0 && LatCountZeroCols(B->Lat) == 0 &&
    LatCountZeroCols(LInter) == 0)
  {
    // This works only IF there are no columns of zeros in the LBLs:
    // they are Z-polyhedra

    ImageA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
    ImageB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
    PInter = DomainIntersection(ImageA, ImageB, MAXNOOFRAYS);
    Domain_Free(ImageB);
    Domain_Free(ImageA);
    #ifdef INTERSECTION_DEBUG
    // fprintf(stderr, "imageA inter imageB = PInter = ");
    // Polyhedron_Print(stderr, P_VALUE_FMT, PInter);
    #endif
    if (emptyQ(PInter))
      Result = NULL;
    else {
      PreImage = DomainPreimage(PInter, LInter, MAXNOOFRAYS);
      Result = LBLAlloc(LInter, PreImage);
      Domain_Free(PreImage);
    }

    Matrix_Free(LInter);
    Domain_Free(PInter);
    #ifdef INTERSECTION_DEBUG
    fprintf(stderr, "Z-polyhedra, simplified intersection = ");
    if(!Result)
      fprintf(stderr, "empty\n");
    else
      LBLPrint(stderr, P_VALUE_FMT, Result);
    fprintf(stderr, "-- exit sLBLIntersection\n");
    #endif

    return (Result);
  }
  else {
    // there's an LBL, need to build the exact LBL manually:
    // { AL z |  BL z' = AL z, z \in AP, z' \in BP }
    // { (0 AL) z' |  (-BL AL Al-Bl)  z' = 0,  z' \in BP
    //          z  |                  z        z  \in AP  }
    int extra_max_rows = 0, extra_B_row = A->Lat->NbRows - 1;
    Matrix *newL = NULL, *extra;
    Polyhedron *newP = NULL, *AP_aligned;
    Matrix_Free(LInter);

    //  (size) BLCols ALCols 1
    //          z'      z   cst
    // newL =   0   |  AL  | Al |
    //          0   |  0   | 1  |   <- this row is in AL already
    newL = Matrix_Alloc(A->Lat->NbRows,
      A->Lat->NbColumns + B->Lat->NbColumns - 1);
    for(int i = 0; i < newL->NbRows; i++) {
      Vector_Set(newL->p[i], 0, B->Lat->NbColumns-1);
      Vector_Copy(A->Lat->p[i], &newL->p[i][B->Lat->NbColumns-1],
        A->Lat->NbColumns);
    }

    // newP:

    // scan the polyhedra of domains A->P and B->P, and build their constraints
    // intersections.
    // Build the constraints:
    //  0/1    z'      z      cst
    //  ap0    0      AP       Ap     # from A->P    -> AP_aligned
    //   0    -BL     AL      Al-Bl   # equalities   \ extra
    //  bp0    BP      0       Bp     # from B->P    /

    // start with a domain of the right dimension (expand dimension of A)
    AP_aligned = align_context(A->P, A->P->Dimension + B->P->Dimension,
      MAXNOOFRAYS);
    #ifdef INTERSECTION_DEBUG
    fprintf(stderr, "AP_aligned = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, AP_aligned);
    #endif
    // extra will be the matrix containing the extra constraints (including
    // equalities: AL z = BL z', initialized once before scanning the polyhedra)
    // count max nbrows of extra
    for(Polyhedron *BP = B->P; BP; BP = BP->next) {
      if(A->Lat->NbRows + BP->NbConstraints > extra_max_rows) {
        extra_max_rows = A->Lat->NbRows + BP->NbConstraints;
      }
    }
    extra = Matrix_Alloc(extra_max_rows, AP_aligned->Dimension + 2);
    Vector_Set(extra->p[0], 0, extra->NbRows * extra->NbColumns); // init 0
    // initialize |0 -BL  AL  Al-Bl| in extra
    for(int row = 0; row < extra_B_row; row++) {
      // equality (0 is set already)
      Vector_Oppose(&B->Lat->p[row][0], &extra->p[row][1],
        B->Lat->NbColumns - 1); // -BL
      Vector_Copy  (&A->Lat->p[row][0], &extra->p[row][B->Lat->NbColumns],
        A->Lat->NbColumns - 1); // AL
      value_substract(extra->p[row][extra->NbColumns - 1], // constant
        A->Lat->p[row][A->Lat->NbColumns - 1],
        B->Lat->p[row][B->Lat->NbColumns - 1]);
    }
    // scan the intersections of each BP with each AP and build union newP
    for(Polyhedron *BP = B->P; BP; BP = BP->next) {
      Polyhedron *P;
      // complement extra with the Constraints of BP
      for(int con = 0; con < BP->NbConstraints; con++) {
        // constraint BP
        Vector_Copy(BP->Constraint[con], extra->p[extra_B_row + con],
          BP->Dimension + 1);
        // + constant Bp
        value_assign(extra->p[extra_B_row + con][extra->NbColumns-1],
          BP->Constraint[con][BP->Dimension+1]);
      }
      extra->NbRows = extra_B_row + BP->NbConstraints;
      #ifdef INTERSECTION_DEBUG
      fprintf(stderr, "extra = ");
      Matrix_Print(stderr, P_VALUE_FMT, extra);
      #endif

      for(Polyhedron *AP = AP_aligned; AP; AP = AP->next) {
        P = AddConstraints(extra->p[0], extra->NbRows, AP, MAXNOOFRAYS);
        #ifdef INTERSECTION_DEBUG
        fprintf(stderr, "Adding P = ");
        Polyhedron_Print(stderr, P_VALUE_FMT, P);
        #endif
        if(emptyQ(P)) {
          Polyhedron_Free(P);
        }
        else {
          newP = AddPolyToDomain(P, newP); // consumes P and newP
        }
      }
    }
    Matrix_Free(extra);
    Domain_Free(AP_aligned);

    // Use newL and newP to build result
    if(newP) {
      Result = LBLAlloc(newL, newP);
    }
    Matrix_Free(newL);
    Domain_Free(newP);
    #ifdef INTERSECTION_DEBUG
    fprintf(stderr, "Manual built intersection = ");
    LBLPrint(stderr, P_VALUE_FMT, Result);
    fprintf(stderr, "-- exit sLBLIntersection\n");
    #endif

    return(Result);
  }

} /* sLBLIntersection */


// /*
//  * Return the difference A - B
//  * between a union of LBLs 'A' and a single LBL 'B'.
//  * 
//  * Algo: remove B from each part of A, and build a list of the resulting LBLs
//  *
//  * USAGE: only the first lattice of B is considered
//  *       (even if next is not NULL).
//  * Creates a new allocated LBL, not necessarily in canonical form
//  */
// static LBL *LBL_sLBL_Difference(LBL *A, LBL *B)
// {
//   LBL *Result = NULL;

//   for(; A; A = A->next) {
//     LBL *diff;

//     #ifdef LBLDIFF_DEBUG
//     fprintf(stderr, "  LBL_sLBL_diff: Removing B from (a part of) A. A = ");
//     sLBLPrint(stderr, P_VALUE_FMT, A);
//     fprintf(stderr, "  LBL_sLBL_diff: Removing B from (a part of) A. B = ");
//     sLBLPrint(stderr, P_VALUE_FMT, B);
//     #endif

//     diff = sLBL_Difference(A, B);  // A - B
//     #ifdef LBLDIFF_DEBUG
//     fprintf(stderr, "  LBL_sLBL_diff: diff = ");
//     LBLPrint(stderr, P_VALUE_FMT, diff);
//     #endif

//     // simple concatenate of diff and result (not canonical)
//     Result = LBLConcatenate(diff, Result);
//   }

//   // Result contains every piece of the solution,
//   // but it is not necessarily in canonical form (will be done be callee)
//   return Result;
// } /* LBL_sLBL_Difference */


// /*
//  * Compute the complement of sLBL A: all points z such that z is not in A.
//  *
//  * Algorithm:
//  * Let L = A->Lat, P = A->P.
//  * complement(A) = Universe() - A = union of:
//  *   1- LBL (Z^d, complement hull(A)), with hull(A) = image by L of P
//  *   2- LBL ((Z^d - L), hull(A)) ---- or ((Z^d - L), universe())
//  *   3- holes of A
//  *      if L has no zero columns -> empty
//  *      = L z' such that there exist no z in A->P such that L z' = L z
//  *      -> need exact shadow
//  */
// static LBL *sLBLComplement2(LBL *A)
// {
//   // testing a new version, just add dimensions and 0 columns in L, and
//   // compute the complement of the coordinate polyhedron
//   // and add the holes of A at the end.
//   Matrix *id = NULL;
//   Polyhedron *univ, *comp, *holes;

//   // compute holes of A:
//   holes = sLBLCompute_holes(A, NULL);

//   Matrix_identity(A->Lat->NbRows, &id);
//   A = sLBLCopy(A);
//   sLBLMake_lattice_equal_to(A, id);
//   Matrix_Free(id);
//   fprintf(stderr, "A lattice equal to Id = ");
//   sLBLPrint(stderr, P_VALUE_FMT, A);

//   // // remove the zero dimensions of A->P:
//   // nz = LatCountZeroCols(A->Lat);
//   // Polyhedron *newP = A->P;
//   // for(int dim = 0; dim < nz; dim++) {

//   // }
//   // fprintf(stderr, "A expanded with lines in 0-dims = ");
//   // sLBLPrint(stderr, P_VALUE_FMT, A);

//   // Compute Universe - A->P
//   univ = Universe_Polyhedron(A->P->Dimension);
//   comp = DomainDifference(univ, A->P, MAXNOOFRAYS);
//   Domain_Free(univ);

//   fprintf(stderr, "LINKING: holes(A) = ");
//   Polyhedron_Print(stderr, P_VALUE_FMT, holes);
//   fprintf(stderr, "         with comp = ");
//   Polyhedron_Print(stderr, P_VALUE_FMT, comp);

//   // just link holes at the end of comp (they are separated)
//   Polyhedron *compEnd = comp;
//   while(compEnd->next)
//     compEnd = compEnd->next;
//   compEnd->next = holes;
//   Domain_Free(A->P);
//   A->P = comp;

//   CanonicalLBL(A);
//   fprintf(stderr, "comp(A) = ");
//   sLBLPrint(stderr, P_VALUE_FMT, A);

//   return(A);
// }
static LBL *sLBLComplement(LBL *A)
{
  LBL *Result = NULL;
  Polyhedron *Univ, *hullA, *comp_hullA;
  LatticeUnion *LatDiff;
  int nbzeros;
  #ifdef COMP_DEBUG
  fprintf(stderr, "\n-- Entering sLBLComplement. A = ");
  sLBLPrint(stderr, P_VALUE_FMT, A);
  #endif
  
  // STEP 1: complement of hull(A)
  hullA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
  Univ = Universe_Polyhedron(A->Lat->NbRows - 1);
  comp_hullA = DomainDifference(Univ, hullA, MAXNOOFRAYS);
  #ifdef COMP_DEBUG
  fprintf(stderr, "STEP 1: Adding hull complement polyhedron: ");
  Polyhedron_Print(stderr, P_VALUE_FMT, comp_hullA);
  #endif
  if (!emptyQ(comp_hullA)) {
    Matrix *Id;
    Id = Identity_Matrix(comp_hullA->Dimension + 1);
    Result = LBLAlloc(Id, comp_hullA);
    Matrix_Free(Id);
  }
  Domain_Free(comp_hullA);
  Domain_Free(Univ);

  // STEP 2: lattice differences (not L) on hull(A) /or Universe
  if(hullA && !emptyQ(hullA))
  {
    LatDiff = LatticeDifference(NULL, A->Lat);
    #ifdef COMP_DEBUG
    fprintf(stderr, "\nSTEP 2: LatDiff =\n");
    PrintLatticeUnion(stderr, P_VALUE_FMT, LatDiff);
    #endif
    // Add all Z-polyhedra to Result: the list of lattices on Universe
    for(LatticeUnion *lat = LatDiff; lat; lat = lat->next) {
      LBL *Ztmp;
      #ifdef COMP_DEBUG
      fprintf(stderr, "Considering Lat diff: ");
      Matrix_Print(stderr, P_VALUE_FMT, lat->M);
      #endif
      Ztmp = malloc(sizeof(LBL));
      Ztmp->Lat = lat->M;
      // can use universe, no need to restrict on the more complicated hullA:
      // Ztmp->P = DomainPreimage(hullA, lat->M, MAXNOOFRAYS);
      Ztmp->P = Universe_Polyhedron(lat->M->NbColumns - 1);
      // remove obvious simplification?
      // -> not necessary since preimage by integer function.
      // Ztmp->P = DomainConstraintSimplify(Ztmp->P, MAXNOOFRAYS);
      #ifdef COMP_DEBUG
      fprintf(stderr, "Adding: ");
      sLBLPrint(stderr, P_VALUE_FMT, Ztmp);
      #endif
      Ztmp->next = Result;
      Result = Ztmp;
    }
    // free LatticeUnion remaining memory (M has been reused as a lattice)
    while(LatDiff) {
      LatticeUnion *next = LatDiff->next;
      free(LatDiff);
      LatDiff = next;
    }
  }
  Domain_Free(hullA);

  // STEP 3: holes (if there are zero columns in Lat)
  if((nbzeros = LatCountZeroCols(A->Lat))) {
    // there are potential holes
    Matrix *newL;
    Polyhedron *holes;
    holes = sLBLCompute_holes(A, NULL);
    #ifdef COMP_DEBUG
    fprintf(stderr, "\nSTEP 3 adding holes = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, holes);
    #endif

    if(holes && !emptyQ(holes))
    {
      newL = RemoveNColumns(A->Lat, A->Lat->NbColumns-1-nbzeros, nbzeros);
    
      Result = LBLConcatenate(LBLAlloc(newL, holes), Result);
      Matrix_Free(newL);
    }
    Domain_Free(holes);
  }

  // CanonicalLBL(Result);

  // Don't need to simplify (remove integer-empty polyhedra)
  // LBLSimplifyEmpty(Result);
  #ifdef COMP_DEBUG
  fprintf(stderr, "\n-- sLBLComplement final result (normalized) = ");
  LBLPrint(stderr, P_VALUE_FMT, Result);
  #endif

  return (Result);
} /* sLBLComplement */


/*
 * complement of LBL A.
 *
 * Intersection of the complements of the single LBLs of A
 */
LBL *LBLComplement(LBL *A)
{
  LBL *Result;
  Result = sLBLComplement(A);
  for(LBL *tmp = A->next; tmp && Result; tmp = tmp->next) {
    LBL *comp, *inter;
    comp = sLBLComplement(tmp);
    inter = LBLIntersection(Result, comp);
    LBLFree(Result);
    LBLFree(comp);
    Result = inter;
  }

  CanonicalLBL(Result);

  return(Result);
} /* LBLComplement */


// /*
//  * Return the difference of two single LBLs A - B.
//  * A and B are single LBLs, but the return value can be a union of LBLs!
//  * Creates a new allocated LBL union
//  *
//  * USAGE: only the first lattice of A and B is considered (no union),
//  *        but A and B can contain several coordinate polyhedra (in ->P).
//  * 
//  * Algorithm:
//  * -> New version: compute A inter complement(B).
//  * -> Former version inspired from the method Gautam describes in his thesis,
//  * modified to handle LBLs.
//  */
// static LBL *sLBL_Difference(LBL* A, LBL* B)
// {
//   LBL *Result, *Bcomp; // union of LBLs
//   LBL *Binter; // intersection = single LBL.

//   #ifdef DIFFERENCE_DEBUG
//   fprintf(stderr, "-- Entering sLBL_Difference. A = ");
//   sLBLPrint(stderr, P_VALUE_FMT, A);
//   #endif
//   if (A->Lat->NbRows != B->Lat->NbRows) {
//     errormsg1("sLBL_Difference", "dimincomp", "incompatible dimensions");
//     return(NULL);
//   }

//   // treat the simple case where the LBLs do not intersect
//   Binter = sLBLIntersection(A, B); // reused below
//   #ifdef DIFFERENCE_DEBUG
//   fprintf(stderr, "Binter = ");
//   LBLPrint(stderr, P_VALUE_FMT, Binter);
//   #endif
//   if(isEmptyLBL(Binter)) {
//     // if B does not intersect A, return A.
//     #ifdef DIFFERENCE_DEBUG
//     fprintf(stderr,
//       "Binter=(A inter B) is empty, so B does not intersect A, we return A\n");
//     #endif
//     LBLFree(Binter);
//     return(LBLCopy(A));
//   }

//   // // Separate the computation in 3 phases:
//   // // 0. compute the difference of the image polyhedra P_A \ P_B (=ImDiff) and
//   // //    add it to the solution LBL (with lattice L_A).
//   // //    This can be an over-approximation of A if A->Lat has zero columns
//   // //    (but not of B)
//   // // 1. compute the rest where the intersection of P_A and P_B have same
//   // //    dimensions (required for lattice difference)
//   // // 2. intersect the result with A to get rid of the over-approximations

//   // LBL *Result = NULL, *Final_Result; // U. of LBLs
//   // LBL *Ainter, *Binter; // single LBL
//   // LatticeUnion *LatDiff;
//   // Polyhedron *imA, *imB, *preimA, *ImDiff, *ImInter; // polyhedral domains

//   // // [STEP 0 (includes Gautam's Step 2)]
//   // imA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
//   // imB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
//   // ImDiff = DomainDifference(imA, imB, MAXNOOFRAYS);
//   // #ifdef DIFFERENCE_DEBUG
//   //   fprintf(stderr, "ImDiff (hull of A that does not cover B) = ");
//   //   Polyhedron_Print(stderr, P_VALUE_FMT, ImDiff);
//   // #endif

//   // // Add (A->Lat, A->P - hull(B)) to the result:
//   // if (!emptyQ(ImDiff)) {
//   //   Polyhedron *RedPolyDiff;
//   //   RedPolyDiff = DomainPreimage(ImDiff, A->Lat, MAXNOOFRAYS);
//   //   // NOTICE: this can be an over-approximation of A
//   //   Result = LBLAlloc(A->Lat, RedPolyDiff);
//   //   #ifdef DIFFERENCE_DEBUG
//   //     fprintf(stderr, "Adding this to the temporary result: ");
//   //     LBLPrint(stderr, P_VALUE_FMT, Result);
//   //   #endif
//   //   Domain_Free(RedPolyDiff);
//   // }

//   // // compute the images intersection of A and B
//   // ImInter = DomainIntersection(imA, imB, MAXNOOFRAYS);
//   // #ifdef DIFFERENCE_DEBUG
//   //   fprintf(stderr, "ImInter (hull of A inter B) = ");
//   //   Polyhedron_Print(stderr, P_VALUE_FMT, ImInter);
//   // #endif
  
//   // // can be simplified, Ainter not really needed!

//   // // compute the part of A that intersects the hull of B in the image space
//   // preimA = DomainPreimage(ImInter, A->Lat, MAXNOOFRAYS);
//   // Ainter = LBLAlloc(A->Lat, preimA);
//   // // NOTICE: this Ainter can be a over-approximation of A

//   // Domain_Free(preimA);
//   // Domain_Free(ImDiff);
//   // Domain_Free(imA);
//   // Domain_Free(imB);

//   // // now Ainter and Binter have same lattices and polyhedra dimensions
//   // #ifdef DIFFERENCE_DEBUG
//   //   fprintf(stderr,
//   //     "-- [STEP1] now we compute the intersection on same lattice dimensions\n");
//   //   fprintf(stderr, "Ainter = ");
//   //   LBLPrint(stderr, P_VALUE_FMT, Ainter);
//   //   fprintf(stderr, "and Binter = ");
//   //   LBLPrint(stderr, P_VALUE_FMT, Binter);
//   // #endif

//   // // LatDiff (union of lattices) is the difference : (A->Lat) - (B->Lat) of
//   // // same dimensions
//   // LatDiff = LatticeDifference(Ainter->Lat, Binter->Lat); 
//   // #ifdef DIFFERENCE_DEBUG
//   //   if(!LatDiff)
//   //     fprintf(stderr, "Empty Lattice difference\n");
//   // #endif

//   // // [STEP 1 of Gautam]:
//   // // Add all Z-polyhedra applying the (list of) lattice difference on ImInter
//   // for(LatticeUnion *tmp = LatDiff; tmp; tmp = tmp->next) {
//   //   LBL *Ztmp;
//   //   #ifdef DIFFERENCE_DEBUG
//   //     fprintf(stderr, "Considering Lat diff: ");
//   //     Matrix_Print(stderr, P_VALUE_FMT, tmp->M);
//   //   #endif
//   //   Ztmp = malloc(sizeof(*Ztmp));
//   //   Ztmp->next = Result;
//   //   Ztmp->Lat = tmp->M;
//   //   Ztmp->P = DomainPreimage(ImInter, tmp->M, MAXNOOFRAYS);
//   //   // NOTICE: this can be an over-approximation of A (but not of B)

//   //   Result = Ztmp;
//   // }
//   // // free LatticeUnion remaining memory (M has been reused as a lattice of
//   // // Result)
//   // while(LatDiff) {
//   //   LatticeUnion *next = LatDiff->next;
//   //   free(LatDiff);
//   //   LatDiff = next;
//   // }

//   // // could also consider the intersection of lattices, where some points of
//   // // lattice B->Lat could have no integer antecedent in B->P and should
//   // // be kept in the result A - B:
//   // // Add the holes of B (that can be included in A but not in B).


//   // Domain_Free(ImInter);
//   // LBLFree(Ainter);
//   // LBLFree(Binter);

//   // if(!Result) {
//   //   #ifdef DIFFERENCE_DEBUG
//   //     fprintf(stderr, "-- result = (NULL)\n");
//   //   #endif
//   //   return(NULL);
//   // }

//   // #ifdef DIFFERENCE_DEBUG
//   //   fprintf(stderr, "-- temporary over-approximation of result = ");
//   //   LBLPrint(stderr, P_VALUE_FMT, Result);
//   // #endif
//   // // intersect the result with A to get the exact LBL in case there was an
//   // // over-approximation of A before.
//   // Final_Result = LBLIntersection(Result, A);
//   // LBLFree(Result);
//   // return(Final_Result);

//   // Which one to use to compute the complement, Binter or B?
//   // which one is simpler? B is larger... but Binter is part of A
//   // Binter is probably better to prepare for the intersection
//   Bcomp = sLBLComplement(Binter);
//   #ifdef DIFFERENCE_DEBUG
//   fprintf(stderr, "Difference = intersection (between A and) Bcomp = ");
//   LBLPrint(stderr, P_VALUE_FMT, Bcomp);
//   #endif
//   LBL *nextA = A->next; A->next = NULL; // unlink A from its next
//   Result = LBLIntersection(Bcomp, A);
//   A->next = nextA;

//   LBLFree(Binter);
//   LBLFree(Bcomp);

//   #ifdef DIFFERENCE_DEBUG
//   fprintf(stderr, "Difference = ");
//   LBLPrint(stderr, P_VALUE_FMT, Result);
//   #endif

//   return(Result);
// } /* sLBL_Difference */


/*
 * Return the image of the single LBL 'A' under the affine function 'Func'
 * 
 * Algorithm:
 * - Multiply Lat by Func,
 * - Canonicalize the result (done by LBLAlloc)
 */
static LBL *sLBLImage(LBL *A, Matrix *Func)
{
  Matrix *newL;
  LBL *result;

  if ((Func->NbColumns != A->Lat->NbRows)) {
    errormsg1("sLBLImage", "dimincomp", "Incompatible dimensions");
    return NULL;
  }

  // if there is an eliminated column of zeros in the result:
  // LBLAlloc will ensure that eliminated points in the projection
  // are integer.

  newL = Matrix_Alloc(Func->NbRows, A->Lat->NbColumns);
  Matrix_Product(Func, A->Lat, newL);
  result = LBLAlloc(newL, A->P);

  Matrix_Free(newL);
  return(result);
} /* sLBLImage */

/*
 * Return the preimage of the single LBL 'Z' under an affine
 * transformation function 'G'. The number of rows of matrix 'G' must
 * be equal to the number of rows of the lattice of Z.
 * Algorithm:
 * - if G is invertible, compute the LBL {G^{-1} L, D},
 * - else, build the LBL { z' | G z' = L z, z \in A->P, z' free},
 *         and remove z by normalizing the result
 */
static LBL *sLBLPreimage(LBL *Z, Matrix *G)
{
  LBL *Result;
  Polyhedron *P, *newP;
  Matrix *Con;

  #ifdef IMAGE_DEBUG
  fprintf(stderr, "===== entering sLBLPreimage =====\n Z = ");
  sLBLPrint(stderr, P_VALUE_FMT, Z);
  fprintf(stderr, "G = ");
  Matrix_Print(stderr, P_VALUE_FMT, G);
  #endif

  if(G->NbRows != Z->Lat->NbRows) {
    // G z' = L z
    errormsg1("sLBLPreimage", "dimincomp", "incompatible dimensions");
    return(NULL);
  }

  // first try if G is invertible
  if(G->NbColumns == G->NbRows) {
    const int dim = G->NbColumns;
    Matrix *tmp, *inv;
    tmp = Matrix_Copy(G);
    inv = Matrix_Alloc(dim, dim);
    if(Matrix_Inverse(tmp, inv) && value_one_p(inv->p[dim-1][dim-1])) {
      // reuse tmp to compute the product Inv Lat
      #ifdef IMAGE_DEBUG
      fprintf(stderr, "G is invertible\nInverse =");
      Matrix_Print(stderr, P_VALUE_FMT, inv);
      #endif
      Matrix_Product(inv, Z->Lat, tmp);
      Result = LBLAlloc(tmp, Z->P);
      Matrix_Free(tmp);
      Matrix_Free(inv);
      return(Result);
    }
    Matrix_Free(tmp);
    Matrix_Free(inv);
  }

  #ifdef IMAGE_DEBUG
  fprintf(stderr, "G is not invertible\n");
  #endif
  // need to build the LBL { z' | G z' = L z, z \in A->P, z' free},

  // d is the dimension of A->P (nb columns of Lat)
  // d' is the number of columns of G = the dimension of the result
  // build the Z-polyhedron = { z' | with P in dimension d + d'
  // such that L z = G z' }
  // then eliminate z by simplifying the result

  // the lattice is spreading z'
  // homogeneous d and d', homogeneous sum is d+d'-1
  int d = Z->Lat->NbColumns;
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
  #ifdef IMAGE_DEBUG
  fprintf(stderr, "newL = ");
  Matrix_Print(stderr, P_VALUE_FMT, newL);
  #endif

  // add the extra dimension on P (first dimensions!)
  newP = align_context(Z->P, d+dp-2, MAXNOOFRAYS);
  #ifdef IMAGE_DEBUG
  fprintf(stderr, "newP = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, newP);
  #endif


  // build the extra constraint to be added to newP: G z' = L z
  // con =    0 |     |     |
  //          . |  G  | -L  | (g-l)
  //          0 |     |     |
  Con = Matrix_Alloc(G->NbRows - 1, d+dp-1+1);
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
    value_substract(Con->p[i][Con->NbColumns-1],
              Con->p[i][Con->NbColumns-1], Z->Lat->p[i][Z->Lat->NbColumns-1]);
  }
  #ifdef IMAGE_DEBUG
  fprintf(stderr, "Con = ");
  Matrix_Print(stderr, P_VALUE_FMT, Con);
  #endif

  P = DomainAddConstraints(newP, Con, MAXNOOFRAYS);
  Matrix_Free(Con);
  Domain_Free(newP);

  Result = LBLAlloc(newL, P);
  Domain_Free(P);
  Matrix_Free(newL);

  return(Result);
} /* sLBLPreimage */


// typedef struct forsimplify {
//   Polyhedron *Pol;
//   LatticeUnion *LatUni;
//   struct forsimplify *next;
// } ForSimplify;
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
//     return (LBLCopy(ZDom));
//   Emp = EmptyLBL(ZDom->Lat->NbRows - 1);
//   ZDomHead = LBLUnion(ZDom, Emp);
//   LBLFree(Emp);
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
//     return LBLCopy(ZPol);
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



/***************************************************************************/
/*  Utility functions for LBL simplification and normalization             */
/***************************************************************************/

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
 * Change A->P such that all polyhedra in this domain have the same set of
 * equalities, that is, the equalities of the first polyhedron of this domain.
 * All the other ones are added to a new LBL, linked to LBL A (in A->next).
 */
static Matrix *sLBLHomogenize_equalities(LBL *A)
{
  #ifdef CANONICAL_DEBUG
  fprintf(stderr, "Checking for equalites in P\n");
  #endif
  LBL *new = NULL;
  Matrix *Equalities = get_equalities(A->P); // get eq from the first one
  Polyhedron *nextpp, *prevpp, *pp;

  prevpp = A->P;
  pp = A->P->next;
  while(pp)
  {
    // check if the equalities of pp->Constraints are the same as the ones
    // of matrix Equalities.
    nextpp = pp->next;
    // nextpp = next polyhedron of the domain (pp can be relinked below:)
    if(!same_equalities(Equalities, pp)) {
      // if not, get pp out.
      if(!new) {
        new = malloc(sizeof(LBL));
        if (!new) {
          errormsg1("sLBLHomogenize_equalities", "outofmem", "Out of Memory");
          return(NULL);
        }
        new->P = NULL;
        new->Lat = Matrix_Copy(A->Lat);
      }
      // remove pp from the list A->P, and get the right next iteration
      prevpp->next = pp->next;
      // add pp to new->P
      pp->next = new->P;
      new->P = pp;
      // prevpp does not change
    }
    else {
      // move on, keep a pointer to the previous one to relink easily
      prevpp = pp;
    }
    pp = nextpp;
  }

  if(new) {
    // include new in the LBL list A
    new->next = A->next;
    A->next = new;
    #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Unified equalites in A->P!\n - First A->P = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, A->P);
    fprintf(stderr, " - A Next ->P = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, A->next->P);
    #endif
  }
  // Now all polyhedra of domain A->P have the same equalities
  return(Equalities);
}


/*
 * Simplify the equalities from A->P, in the single LBL A.
 * In place. A->P is a domain (all polyhedra have the same set of equalities).
 * 
 * Modifies A->Lat and A->P.
 */
static void sLBLSimplify_equalities(LBL *A, Matrix *Equalities)
{
  Matrix *H = NULL, *NewL;
  Matrix *eq_hermite = NULL, *ker = NULL;
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "P has equalities\n");
    fprintf(stderr, "Equality matrix (including constants): ");
    Matrix_Print(stderr, P_VALUE_FMT, Equalities);
  #endif

  // compute the kernel of Equalities matrix using Hermite:
  left_hermite(Equalities, &eq_hermite, NULL, &ker);
  Matrix_Free(eq_hermite);

  // the kernel of Equalities is the last (NbEq) columns of ker
  //            -----NbEq------
  // ker = *..* k_0..k_{NbEq-1} c
  // move them left, matrix ker becomes:
  //       -----NbEq------
  // ker = k_0..k_{NbEq-1} c *..*
  for(int i = 0; i < ker->NbRows; i++) {
    for(int j = 0; j < ker->NbColumns - A->P->NbEq; j++) {
      value_assign(ker->p[i][j], ker->p[i][j + A->P->NbEq]);
    }
  }
  // set the right number of columns
  ker->NbColumns = ker->NbColumns - A->P->NbEq;
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "ker of Eq: ");
    Matrix_Print(stderr, P_VALUE_FMT, ker);
  #endif

  // Compute H = affine HNF of ker:
  AffineHermite(ker, &H, NULL);
  Matrix_Free(ker);
  // We know that: Eq . H = Eq . Ker(Eq) = 0
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Matrix H: ");
    Matrix_Print(stderr, P_VALUE_FMT, H);
  #endif

  // The result is just empty (because it is rational) when the bottom-right
  // value of H is not one.
  // ----> But this never happens since H is HNF!
  // if(value_notone_p(H->p[H->NbRows-1][H->NbColumns-1])) {
  //   Domain_Free(A->P);
  //   A->P = NULL;      // empty
  // }
  // else
  {
    // NewL = L . H
    NewL = Matrix_Alloc(A->Lat->NbRows, H->NbColumns);
    Matrix_Product(A->Lat, H, NewL);

    // NewP = H^{-1} . P
    Polyhedron* NewP = DomainPreimage(A->P, H, MAXNOOFRAYS);
    // H is not necessarily unimodular, but it has multiplied Lat: this
    // newP could have less points than the original one, which is correct
    // since some non-integer solutions to the equalities have been removed.

    // update A
    Domain_Free(A->P);
    Matrix_Free(A->Lat);
    A->P = NewP;
    A->Lat = NewL;
  }

  Matrix_Free(H);
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "New Lat: ");
    Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
    fprintf(stderr, "New P: ");
    Polyhedron_Print(stderr, P_VALUE_FMT, A->P);
  #endif
} /* sLBLSimplify_equalities */


/*
 * compute the inside of polyhedron P that can be projected (along dim) to
 * get the dark shadow.
 * 
 * consider P a single polyhedron, even if P->next is set.
 * return NULL if dark == input P
 * 
 * Args:
 * - P: polyhedron (ignore ->next)
 * - dim: 0 <= dim < P->Dimension
 */
static Polyhedron *polyhedron_dark_source(Polyhedron *P, int dim)
{
  // check if that dimension is constrained in P
  // - if it is not constrained, ignore.
  // - if it is positive constrained,
  //   scan all positive constraints on i_dim of the form:
  //     {... + alpha . i_dim + ... + c >= 0}, with alpha > 0
  //   and add this constraint to P:
  //     {... + alpha . i_dim + ... + c + alpha-1 >= 0}
  // - if negative constrained do the opposite
  //   (take the smallest number of constraints between them)

  int pos_constrained = 0, // number of alpha's > 1
      neg_constrained = 0, // number of alpha's < -1
      eq_constrained = 0;  // equality (both pos and neg)
  Polyhedron *result;

  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Entering dark_source. dimension: %d\n", dim);
    fprintf(stderr, "Polyhedron: ");
    Polyhedron_Print(stderr, P_VALUE_FMT, P);
  #endif
  // count constraints abs(alpha) > 1 on this variable (dim)
  for(int c = 0; c < P->NbConstraints; c++) {
    if(value_zero_p(P->Constraint[c][dim+1])) {
      // alpha = 0
      continue;
    }
    if(value_zero_p(P->Constraint[c][0])) {
      // should not happen, equalities have been removed before
      eq_constrained++;
    }
    else if(value_pos_p(P->Constraint[c][dim+1]) &&
      value_notone_p(P->Constraint[c][dim+1])) {
      // alpha > 1 (strict)
      pos_constrained++;
    }
    else if(value_neg_p(P->Constraint[c][dim+1]) &&
      value_notmone_p(P->Constraint[c][dim+1])) {
      // alpha < -1 (strict)
      neg_constrained++;
    }
  }

  if(eq_constrained > 0 || pos_constrained == 0 || neg_constrained == 0) {
    // at least one side of the polyhedron is integer or open to infinite,
    // so the dark source is just equal to P
    return(NULL);
  }
  else {
    // pos_constrained > 0 and neg_constrained > 0.
    Matrix *extra;
    int nb_extra = 0;
    if(pos_constrained < neg_constrained) {
      // consider it pos_constrained (less extra to add)
      extra = Matrix_Alloc(pos_constrained, P->Dimension + 2);
      // add extra constraints on positive ones
      for(int c = 0; c < P->NbConstraints; c++) {
        if(value_zero_p(P->Constraint[c][dim+1]) ||
          value_zero_p(P->Constraint[c][0])) {
          continue;
        }
        if(value_pos_p(P->Constraint[c][dim+1]) &&
          value_notone_p(P->Constraint[c][dim+1])) { // alpha > 1 (strict)
          // from constraint:
          // {... + alpha . i_dim + ... + c >= 0}
          // add constraint:
          // {... + alpha . i_dim + ... + c - (alpha-1) >= 0}
          Vector_Copy(P->Constraint[c], extra->p[nb_extra], P->Dimension + 2);
          // constant update: - alpha + 1
          value_substract(extra->p[nb_extra][extra->NbColumns-1],
                          extra->p[nb_extra][extra->NbColumns-1],
                          extra->p[nb_extra][dim+1]);
          value_add_int(extra->p[nb_extra][extra->NbColumns-1],
                        extra->p[nb_extra][extra->NbColumns-1], 1);
          nb_extra++;
        }
      }
    }
    else { // neg_constrained
      extra = Matrix_Alloc(neg_constrained, P->Dimension + 2);
      for(int c = 0; c < P->NbConstraints; c++) {
        if(value_zero_p(P->Constraint[c][dim+1]) ||
          value_zero_p(P->Constraint[c][0])) {
          continue;
        }
        if(value_neg_p(P->Constraint[c][dim+1]) &&
          value_notmone_p(P->Constraint[c][dim+1])) { // alpha < -1 (strict)
          // from constraint:
          // {... + alpha . i_dim + ... + c >= 0}
          // add constraint:
          // {... + alpha . i_dim + ... + c - (-alpha-1) >= 0}
          Vector_Copy(P->Constraint[c], extra->p[nb_extra], P->Dimension + 2);
          // constant update: + alpha + 1
          value_addto(extra->p[nb_extra][extra->NbColumns-1],
                      extra->p[nb_extra][extra->NbColumns-1],
                      extra->p[nb_extra][dim+1]);
          value_add_int(extra->p[nb_extra][extra->NbColumns-1],
                        extra->p[nb_extra][extra->NbColumns-1], 1);
          nb_extra++;
        }
      }
    }
    // add the extra constraints
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "Adding constraints: ");
      Matrix_Print(stderr, P_VALUE_FMT, extra);
    #endif
    result = AddConstraints(extra->p[0], nb_extra, P, MAXNOOFRAYS);
    Matrix_Free(extra);
  }

  return (result);
} /* polyhedron_dark_source */


/*
 * compute the projection of domain P along dimension eliminate.
 * 
 * P is a polyhedral domain
 */
static Polyhedron *domain_project(Polyhedron *P, int eliminate)
{
  Polyhedron *Pext;
  Matrix *ray;
  // # 1 ..........  elim+1 elim+2... dim cst
  // --- elim+1 ---   xx    ---  dim-elim ---
  // 
  // new homogeneous dimension of cons and ray matrices: dim + 1.
  const int rest = P->Dimension - eliminate;

  if(!P || emptyQ(P)) {
    return(NULL);
  }
  #ifdef CANONICAL_DEBUG
  fprintf(stderr, "Entering exact shadow -dimension %d- of:", eliminate);
  Polyhedron_Print(stderr, P_VALUE_FMT, P);
  #endif

  // add ray: a line along dimension dim
  ray = Matrix_Alloc(1, P->Dimension + 2);
  Vector_Set(ray->p[0], 0, P->Dimension + 2);
  value_set_si(ray->p[0][eliminate+1], 1);

  Pext = DomainAddRays(P, ray, MAXNOOFRAYS);
  Matrix_Free(ray);


  // just remove the dimension from Pext
  for(Polyhedron *tmp = Pext; tmp; tmp = tmp->next) {
    // eliminate dimension in constraint matrix
    // (keep memory alignment of the Constraints vector of values)
    for(int c = 0; c < tmp->NbConstraints; c++) {
      if(c != 0)
        Vector_Copy(tmp->Constraint[c]+0, tmp->Constraint[c]-c, eliminate + 1);
      Vector_Copy(tmp->Constraint[c]+eliminate+2, tmp->Constraint[c]-c+eliminate+1, rest);
      tmp->Constraint[c] -= c;
    }
    // eliminate dimension in ray matrix
    // (keep memory alignment of the Ray vector of values)
    for(int r = 0; r < tmp->NbRays; r++) {
      if(r != 0)
        Vector_Copy(tmp->Ray[r]+0, tmp->Ray[r]-r, eliminate + 1);
      Vector_Copy(tmp->Ray[r]+eliminate+2, tmp->Ray[r]-r+eliminate+1, rest);
      tmp->Ray[r] -= r;
    }
    tmp->Dimension--;
  
    // remove the null ray that is somewhere...
    for(int r = 0; r < tmp->NbRays; r++) {
      int i;
      for(i = 0; i <= tmp->Dimension; i++) {
        if(value_notzero_p(tmp->Ray[r][i]))
          break;
      }
      if(i == tmp->Dimension + 1) {
        // r is the null ray, erase it with the end of the ray matrix.
        Vector_Copy(tmp->Ray[r+1], tmp->Ray[r], (tmp->NbRays - r - 1) * (tmp->Dimension + 2));
        tmp->NbBid--;
        tmp->NbRays--;
      }
    }
  }
  #ifdef CANONICAL_DEBUG
  fprintf(stderr, "projected result P = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, Pext);
  #endif

  return(Pext);
} /* domain_project */

// // previous version (much slower on large polyhedra!)
// static Polyhedron *domain_project1(Polyhedron *P, int dim)
// {
//   Matrix *T; // transformation: Id without the dim column
//   Polyhedron *image;
//   #ifdef CANONICAL_DEBUG
//     fprintf(stderr, "Entering exact shadow -dimension %d- of:", dim);
//     Polyhedron_Print(stderr, P_VALUE_FMT, P);
//   #endif

//   T = Matrix_Alloc(P->Dimension, P->Dimension + 1);
//   Vector_Set(T->p_Init, 0, T->p_Init_size);
//   for(int i = 0; i < P->Dimension; i++) {
//     if(i >= dim) {
//       value_set_si(T->p[i][i+1], 1);
//     }
//     else {
//       value_set_si(T->p[i][i], 1);
//     }
//   }

//   image = DomainImage(P, T, MAXNOOFRAYS);
//   #ifdef CANONICAL_DEBUG
//   fprintf(stderr, "projected result P = ");
//   Polyhedron_Print(stderr, P_VALUE_FMT, image);
//   #endif
//   Matrix_Free(T);

//   return (image);
// } /* domain_project */


/*
 * compute the dark shadow of domain P projected along dimension dim.
 * 
 * P is a polyhedral domain
 */
static Polyhedron *domain_dark_shadow(Polyhedron *P, int dim)
{
  Polyhedron *dark = NULL;
  for(Polyhedron *pp = P; pp; pp = pp->next) {
    Polyhedron *pp_inside, *pp_shadow;

    pp_inside = polyhedron_dark_source(pp, dim);
    if(pp_inside) {
      pp_shadow = domain_project(pp_inside, dim);
    }
    else {
      // pp_inside == pp
      Polyhedron *ppnext = pp->next;
      pp->next = NULL;
      pp_shadow = domain_project(pp, dim);
      pp->next = ppnext;
    }
    dark = AddPolyToDomain(pp_shadow, dark);
    if(pp_inside) {
      Polyhedron_Free(pp_inside);
    }
  }
  return(dark);
} /* domain_dark_shadow */


/*
 * Generate a polyhedron of just one point
 */
Polyhedron *GenPoly(int dim, Value *val)
{
  Matrix *rays;
  Polyhedron *res;

  rays = Matrix_Alloc(1, dim + 2);
  Vector_Copy(val, rays->p[0], dim + 1);
  value_set_si(rays->p[0][0], 1);        // vertex
  value_set_si(rays->p[0][dim + 1], 1);  // denominator = 1
  res = Rays2Polyhedron(rays, MAXNOOFRAYS);
  #ifdef HOLES_DEBUG
  fprintf(stderr, "*** generating polyhedron for vertex ***: ");
  Matrix_Print(stderr, P_VALUE_FMT, rays);
  #endif

  Matrix_Free(rays);
  return(res);
}


/*
 * Scan a single (hopefully small) pointy polyhedron R, and check if the
 * scan hits an integer point.
 * R is the intersection of:
 * - the shadow (rest) to be scanned whole,
 * - AP, the whole domain to check for an integer solution
 * Returns the union of polyhedra in the rest that are not holes.
 * 
 * Uses the allocated array of values val of dimension R->dimension+2
 * position = index position of current loop index
 * (from 1 to R->Dimension)
 * 
 * Recursive,
 * - if scanning the shadow (position <= dimrest), scan all possible values
 *   and recursive call to scan, accumulating all results together
 * - if the shadow is scanned already (position > dimrest),
 *   recursive scan R, but stop as soon as an integer solution is found
 * 
 * val must be set to 0's where not used, val[Dimension+1] must be set to 1.
 */
Polyhedron *Scan_RestAP(Polyhedron *R, Value *val, int position, int dimrest)
{
  Polyhedron *Result = NULL;
  Value LB, UB;

  //           -----------------R------------------
  //           -------rest-------  ----------------
  // val = 0   *    ...         *  *      ...     *  1
  // idx:  0   1    ...   dimrest       ...   R->Dimension
  //           |-> position...
  
  if(!R) {
    // end here, it is a hit!
    // generate polyhedron of the point = value val
    return(GenPoly(dimrest, val));
  }

  #ifdef HOLES_DEBUG
  fprintf(stderr, "Enter Scan_Rest, position = %2d.", position);
  fprintf(stderr, "val = (");
  for(int i=0; i <= R->Dimension + 1; i++) {
    value_print(stderr, P_VALUE_FMT, val[i]);
    if(i == dimrest || i == 0 || i == R->Dimension)
      fprintf(stderr, "::");
  }
  fprintf(stderr, ")\n");
  #endif
  value_init(LB);
  value_init(UB);

  if(position <= dimrest) {
    // scan Rest
    if(lower_upper_bounds(position, R, val, &LB, &UB) != 0) {
      errormsg1("Scan_RestAP", "infinite", "trying to scan infinity!");
      return(NULL);
    }
    // #ifdef HOLES_DEBUG
    // fprintf(stderr, "scanR - position %d - looping from ", position);
    // value_print(stderr, P_VALUE_FMT, LB);
    // fprintf(stderr, "to ");
    // value_print(stderr, P_VALUE_FMT "\n", UB);
    // #endif

    // loop LB -> UB
    for(; value_le(LB, UB); value_increment(LB, LB)) {
      value_assign(val[position], LB);

      Result = AddPolyToDomain(Scan_RestAP(R->next, val, position + 1,
        dimrest), Result);
      // and continue with all other values
    }
    // reset value[position] for next scans
    value_set_si(val[position], 0);
  }
  else {
    // scan AP
    if(lower_upper_bounds(position, R, val, &LB, &UB) != 0) {
      // infinite AP, has an int solution.
      return(GenPoly(dimrest, val));
    }
    // #ifdef HOLES_DEBUG
    // fprintf(stderr, "scanAP - position %d - looping from ", position);
    // value_print(stderr, P_VALUE_FMT, LB);
    // fprintf(stderr, "to ");
    // value_print(stderr, P_VALUE_FMT "\n", UB);
    // #endif

    for(; value_le(LB, UB); value_increment(LB, LB)) {
      value_assign(val[position], LB);

      if((Result = Scan_RestAP(R->next, val, position + 1, dimrest))) {
        // it's a hit: stop here and add this point to result :)
        value_clear(UB);
        value_clear(LB);
        value_set_si(val[position], 0);
        // early exit
        return(Result);
      }
    }
    // reset value for next scans
    value_set_si(val[position], 0);
  }

  value_clear(UB);
  value_clear(LB);
  return(Result);
} /* Scan_RestAP */


/*
 * Bound the single polyhedron P by a box containing an integer point
 * 
 * returns a newly allocated polyhedron, or P if there are no lines/rays
 */
static Polyhedron *bound_polyhedron(Polyhedron *P)
{
  int num_lr; // number of lines+rays
  const int dim = P->Dimension;
  Value ONE;
  Matrix *new_rays;
  Polyhedron *res;

  // count number of lines and rays
  num_lr = 0;
  for(int r = 0; r < P->NbRays; r++) {
    if(value_zero_p(P->Ray[r][0]) || value_zero_p(P->Ray[r][dim + 1])) {
      num_lr++;
    }
  }
  if(num_lr == 0) {
    return(P);
  }

  value_init(ONE);
  value_set_si(ONE, 1);

  
  // build the matrix of vertices of the bounded P:
  // original vertices + each line/ray added to each of them.
  new_rays = Matrix_Alloc((P->NbRays - num_lr) * (num_lr + 1), dim + 2);
  new_rays->NbRows = 0;

  // scan the vertices (check each of them if it is a vertex or line/ray):
  for(int v = 0; v < P->NbRays; v++) {
    if(value_zero_p(P->Ray[v][0]) || value_zero_p(P->Ray[v][dim + 1])) {
      // it's a line/ray
      continue;
    }
    // P->Ray[v] is a vertex.
    // copy the original vertex v:
    Vector_Copy(P->Ray[v], new_rays->p[new_rays->NbRows], dim + 2);
    new_rays->NbRows++;
    // add each line/ray to v (taking v's divisor into account)
    for(int l = 0; l < P->NbRays; l++) {
      if(value_zero_p(P->Ray[l][0]) || value_zero_p(P->Ray[l][dim + 1])) {
        // it's a line/ray
        Vector_Combine(P->Ray[v], P->Ray[l], new_rays->p[new_rays->NbRows], ONE,
          P->Ray[v][dim + 1], dim + 2);
          // value_assign(new_rays->p[new_rays->NbRows][dim + 1], P->Ray[v][dim + 1]);
        value_set_si(new_rays->p[new_rays->NbRows][0], 1);
        new_rays->NbRows++;
      }
    }
  }
  #ifdef SIMPLIFY_DEBUG
  fprintf(stderr, "new_rays = ");
  Matrix_Print(stderr, P_VALUE_FMT, new_rays);
  #endif

  res = Rays2Polyhedron(new_rays, MAXNOOFRAYS);
  Matrix_Free(new_rays);
  value_clear(ONE);
  return(res);
}

/*
 * Compute the coordinate polyhedron containing the holes of the single LBL A.
 *
 * Usage: set *pExact to the exact shadow if the pointer is not NULL
 *        (can be used by caller)
 *
 * Algo:
 * - compute the domain rest = (exact shadow - dark shadow)
 * - scan all its integer points and verify for each point:
 *      if there is an integer point in the intersection with the coordinate
 *      polyhedron, add it to the polyhedral domain not_a_hole
 * - return (rest - not_a_hole)
 */
static Polyhedron *sLBLCompute_holes(LBL *A, Polyhedron **pExact)
{
  int nbzeros;
  Polyhedron *exact = A->P, *dark = A->P; // initialize with P then project
  Polyhedron *rest, *rest_AP; // exact shadow - dark shadow (polyhedral domain)
  Polyhedron *tmp, *U0, *not_a_hole = NULL, *holes;
  Vector *v;


  #ifdef HOLES_DEBUG
  fprintf(stderr, "\n-- Entering compute holes. A = ");
  sLBLPrint(stderr, P_VALUE_FMT, A);
  #endif
  nbzeros = LatCountZeroCols(A->Lat);
  if(nbzeros == 0) {
    return(NULL);
  }
  for(int z = 0; z < nbzeros; z++) {
    Polyhedron *d, *e; // shadow polyhedra after eliminating the column
    int col = A->Lat->NbColumns - 2 - z;

    d = domain_dark_shadow(dark, col);
    e = domain_project(exact, col);
    if(dark != A->P) {    // (do not free the original A->P first ones)
      Domain_Free(dark);  // no longer need the previous calculated ones
      Domain_Free(exact);
    }
    dark = d;
    exact = e;
  }
  #ifdef HOLES_DEBUG
  fprintf(stderr, "\nDark = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, dark);
  fprintf(stderr, "Exact = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, exact);
  #endif

  // rest is the polyhedral domain (exact - dark) in origin-nbzero col space
  rest = DomainDifference(exact, dark, MAXNOOFRAYS);
  Domain_Free(dark);
  if(pExact) {
    *pExact = exact;  // keep a copy of exact shadow (caller can reuse it)
  }
  else {
    Domain_Free(exact);
  }

  // make rest disjoint to scan each point once
  tmp = Disjoint_Domain(rest, 0, MAXNOOFRAYS);
  Domain_Free(rest);
  rest = tmp;

  // simplify obvious non integer cases
  rest = DomainConstraintSimplify(rest, MAXNOOFRAYS);

  // exit if rest is empty
  if(!rest || emptyQ(rest)) {
    Domain_Free(rest);
    return(NULL);
  }
  #ifdef HOLES_DEBUG
  fprintf(stderr, "Rest = (Exact - Dark) = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, rest);
  #endif

  // // PREPARE SCAN:
  // vector of fixed values
  v = Vector_Alloc(A->P->Dimension + 2);
  // universe (dim 0)
  U0 = Universe_Polyhedron(0);


  // scan each piece of rest_AP
  for(Polyhedron *R = rest; R; R = R->next) {
    Polyhedron *nextR, *inter, *bounded_rest;
    Polyhedron *new_not_a_hole = NULL;
    // nullify next R to ensure we work on a single rest
    nextR = R->next; // save and
    R->next = NULL;  // unlink next

    // Check if R is bounded. If it is not, just make it a bounded box
    // (add the lines/rays to the vertices to get new vertices)
    // the holes are necessarily regular in a not bounded piece of R.
    // keep memory of the lines/rays in R to add them back to the solution
    bounded_rest = bound_polyhedron(R);
  
    // rest_AP = bounded_rest dimension expanded to A->P
    rest_AP = bounded_rest;
    for(int d = R->Dimension + 1; d <= A->P->Dimension; d++) {
      tmp = domain_insert_dim(rest_AP, d);
      if(rest_AP != bounded_rest)
        Domain_Free(rest_AP);
      rest_AP = tmp;
    }
    // intersect with A->P
    inter = DomainIntersection(rest_AP, A->P, MAXNOOFRAYS);
    inter = DomainConstraintSimplify(inter, MAXNOOFRAYS);
    if(rest_AP != bounded_rest)
      Domain_Free(rest_AP);

    // scan the pieces of inter = (R inter A->P), compute holes
    while(inter) {
      Polyhedron *scanR, *nextI;
      nextI = inter->next;
      inter->next = NULL;

      // prepare to scan the points of rest
      scanR = Polyhedron_Scan(inter, U0, MAXNOOFRAYS);
  
      // init vector to (0...0 1)
      Vector_Set(v->p, v->Size-1, 0);
      value_set_si(v->p[v->Size-1], 1);
  
      // scan
      #ifdef HOLES_DEBUG
      fprintf(stderr, "------- Calling Scan_Rest -------\n");
      fprintf(stderr, " scanR = ");
      Polyhedron_Print(stderr, P_VALUE_FMT, scanR);
      #endif
      // scan: search points that are not holes
      // (*original* R dimension passed to the function)
      tmp = Scan_RestAP(scanR, v->p, 1, R->Dimension);
      // update new_not_a_hole
      while(tmp) {
        Polyhedron *next = tmp->next;
        tmp->next = NULL;
        new_not_a_hole = AddPolyToDomain(tmp, new_not_a_hole);
        tmp = next;
      }
  
      Domain_Free(scanR);
      Polyhedron_Free(inter);
      inter = nextI;
    } // while(inter)

    #ifdef HOLES_DEBUG
    fprintf(stderr, "got not holes = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, new_not_a_hole);
    fprintf(stderr, "------- End Scan_Rest -------\n");
    #endif

    // if R was unbounded
    if(bounded_rest != R) {
      Matrix *lines_rays;
      Domain_Free(bounded_rest);

      if(new_not_a_hole) {
        // add lines/rays of rest to new_not_a_hole
        lines_rays = Matrix_Alloc(R->NbRays, R->Dimension + 2);
        lines_rays->NbRows = 0;
        for(int r = 0; r < R->NbRays; r++) {
          if(value_zero_p(R->Ray[r][0])
            || value_zero_p(R->Ray[r][R->Dimension + 1])) {
              Vector_Copy(R->Ray[0], lines_rays->p[lines_rays->NbRows], (R->Dimension + 2));
              lines_rays->NbRows++;
            }
        }
        tmp = DomainAddRays(new_not_a_hole, lines_rays, MAXNOOFRAYS);
        Domain_Free(new_not_a_hole);
        new_not_a_hole = tmp;
        Matrix_Free(lines_rays);
      }
      #ifdef HOLES_DEBUG
      fprintf(stderr, "not holes was unbounded, new_not_a_hole = ");
      Polyhedron_Print(stderr, P_VALUE_FMT, new_not_a_hole);
      #endif
    }
    // not_a_hole = link new_not_a_hole and not_a_hole
    if(new_not_a_hole) {
      Polyhedron *end = new_not_a_hole;
      while(end->next)
        end = end->next;
      end->next = not_a_hole;
      not_a_hole = new_not_a_hole;
    }

    R->next = nextR;
  } // for(R)

  Vector_Free(v);
  Domain_Free(U0);
  #ifdef HOLES_DEBUG
  fprintf(stderr, "------- end loop on rest -------\n");
  fprintf(stderr, "not holes in (exact-dark) = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, not_a_hole);
  #endif

  // build final domain: (rest - not_a_hole)
  holes = DomainDifference(rest, not_a_hole, MAXNOOFRAYS);
  holes = DomainConstraintSimplify(holes, MAXNOOFRAYS);
  holes = Domain_Remove_Integer_Empty(holes);
  Domain_Free(rest);
  Domain_Free(not_a_hole);
  if(!holes || emptyQ(holes)) {
    #ifdef HOLES_DEBUG
    fprintf(stderr, "sLBLCompute_holes returning: <NULL>\n");
    #endif
    Domain_Free(holes);
    return(NULL); // no holes
  }

  #ifdef HOLES_DEBUG
  fprintf(stderr, "------- end sLBLCompute_holes -------\n");
  fprintf(stderr, "returning holes = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, holes);
  #endif
  return(holes);
} /* sLBLCompute_holes */


/*
 * Try to eliminate the zero columns of lattice A->Lat through
 * projection.
 * 
 * Successive (along each zero-dimension)
 *   elimination if exact shadow == dark shadow
 * Restart again from last column after a successful column elimination.
 */
static void sLBL_Simplify_Zero_Dimensions(LBL *A)
{
  Bool modified = False; // True if something was projected
  #ifdef CANONICAL_DEBUG
  fprintf(stderr, "Entering sLBL_Simplify_Zero_Dimensions\n");
  #endif

  // scan the zero columns on the right of A->Lat
  for (int col = A->Lat->NbColumns-2 ; col >= 0; col--) {
    int i;
    for (i = 0; i < A->Lat->NbRows; i++) {
      if (value_notzero_p(A->Lat->p[i][col])) {
        break;
      }
    }
    if(i == A->Lat->NbRows) {
      // col is empty
      // Compute the dark shadow and the exact shadow.
      // If dark shadow is in exact shadow: can project out the dimension
      // else: keep

      // What if some polyhedra of the union can be simplified and others
      // cannot? should we separate them or just stay at the domain level?
      // -> staying at the domain level keeps the number of lattices low,
      //    so it's probably best, even if they have zero-columns...
      Polyhedron *dark = domain_dark_shadow(A->P, col);
      // compute exact projection of P and check if dark covers exact:
      Polyhedron *exact = domain_project(A->P, col);
      Polyhedron *diff;
      diff = DomainDifference(exact, dark, MAXNOOFRAYS); // diff = exact - dark

      if(! emptyQ(diff)) {
        // try to remove obvious integer-empty solutions.
        diff = DomainConstraintSimplify(diff, MAXNOOFRAYS);
      }
      // could check if diff has no integer solution... but this is complex.
      // Reserved for LBLSimplify()

      if(emptyQ(diff)) {
        // if exact - dark = 0, project out the column :)
        Matrix *newL;
        #ifdef CANONICAL_DEBUG
        fprintf(stderr, "Exact == Dark. Removing column %d of Lat\n", col);
        #endif

        Domain_Free(A->P);
        A->P = exact;
        // remove column from A->Lat
        newL = RemoveColumn(A->Lat, col);
        Matrix_Free(A->Lat);
        A->Lat = newL;
        if(col != A->Lat->NbColumns-2) {
          // one of the "inner" columns was eliminated, check at the end if
          // one of the outer one can be eliminated now: call this same
          // function at the end to try again.
          modified = True;
        }
      }
      else {
        #ifdef CANONICAL_DEBUG
        fprintf(stderr, "Exact != Dark. Keeping column %d of Lat\n", col);
        #endif
        Domain_Free(exact);
      }
      Domain_Free(diff);
      Domain_Free(dark);
    }
    else {
      // non empty column, everything on the left is also non empty, exit loop.
      break;
    }
  }
  #ifdef CANONICAL_DEBUG
  fprintf(stderr, "Exiting sLBL_Simplify_Zero_Dimensions\n");
  LBLPrint(stderr, P_VALUE_FMT, A);
  #endif
  if(modified) {
    // an inner column was eliminated, we can try again to eliminate
    // another one
    sLBL_Simplify_Zero_Dimensions(A);
  }
} /* sLBL_Simplify_Zero_Dimensions */


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

    // to compute HNF of the lattice (constant part moved left/top)
    Matrix_Move_Homogeneous_Dim_First(A->Lat);

    // We will use U of left Hermite, such that LU = H.
    left_hermite(A->Lat, &H, NULL, &U);

    // Move the constant back to right/bottom in H and U.
    Matrix_Move_Homogeneous_Dim_Last(H);
    Matrix_Move_Homogeneous_Dim_Last(U);

    // set the new Lat matrix as H
    Matrix_Free(A->Lat);
    A->Lat = H;

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "New Lat (HNF): ");
      Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
      fprintf(stderr, "U = ");
      Matrix_Print(stderr, P_VALUE_FMT, U);
    #endif

    // Now update of A->P using the preimage by U (is unimodular)
    Polyhedron *NewP = DomainPreimage(A->P, U, MAXNOOFRAYS);
    Domain_Free(A->P);
    A->P = NewP;
    Matrix_Free(U);

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "New P: ");
      Polyhedron_Print(stderr, P_VALUE_FMT, A->P);
    #endif
    // A->Lat is now in canonical form
  }
  else {
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "A is HNF.\n");
    #endif
  }  
} /* sLBL_Lat_Normalize */


/*
 * Remove the empty sLBL's from a list of LBLs.
 * In place.
 */
static void LBL_Remove_Empty(LBL *A)
{
  LBL *current, *previous;
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Entering Remove_Empty\n");
  #endif

  if(!A)
    return;

  // Scan from A->next, and relink previous to next if empty found
  previous = A;
  current = A->next;
  while(current) {
    if(!current->P || emptyQ(current->P)) {
      // remove current
      #ifdef CANONICAL_DEBUG
        fprintf(stderr, "Found empty sLBL, relinking previous to next\n");
      #endif
      previous->next = current->next; // relink previous
      Domain_Free(current->P);        // free current
      Matrix_Free(current->Lat);
      free(current);
      current = previous->next;       // new current (previous does not change)
    }
    else {
      // advance to next
      previous = current;
      current = current->next;
    }
  }

  // At the end, if there is something linked to an empty LBL, need to replace
  // the content of head A with the one of A->next (and free A->next)
  if(A->next && (!A->P || emptyQ(A->P))) {
    LBL *nextA;
    #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Found empty sLBL at head, replacing head with next\n");
    #endif
    nextA = A->next;
    Domain_Free(A->P);
    Matrix_Free(A->Lat);
    *A = *nextA;
    free(nextA);
  }

  // If A is empty, make sure it is canonical
  if(A->P && emptyQ(A->P)) {
    // not canonical if A->P is not NULL
    Domain_Free(A->P);
    A->P = NULL;
    // if empty, we do not care about the number of columns of A
    // (was:)
    // int dimension = A->Lat->NbRows;
    // if(A->Lat->NbColumns != 1) {
    //   Matrix_Free(A->Lat);
    //   A->Lat = Matrix_Alloc(dimension, 1);
    //   for(int j=0 ; j < dimension-1; j++) {
    //     value_set_si(A->Lat->p[0][j], 0);
    //   }
    //   value_set_si(A->Lat->p[0][dimension-1], 1);
    // }
  }
  #ifdef CANONICAL_DEBUG
  fprintf(stderr, "Exit Remove_Empty, A = ");
  LBLPrint(stderr, P_VALUE_FMT, A);
  #endif
} /* LBL_Remove_Empty */


/*
 * check for an integer solution of a polyhedron scan, return True if found
 * else, return False.
 * scan should be a bounded polyhedron scan.
 * 
 * early exit: return True as soon as an integer point is found.
 * The found integer solution is in the array val.
 */
static Bool polyhedron_int_solution(Polyhedron *scan, Value *val, int position)
{
  Value LB, UB;

  #ifdef SIMPLIFY_DEBUG2
  fprintf(stderr, "  Entering polyhedron_int_solution: position = %d, val = [",
    position);
  for(int p=1 ; p<=position; p++)
    value_print(stderr, P_VALUE_FMT, val[p]);
  fprintf(stderr, "]\n");
  #endif
  if(! scan) {
    // found an integer point, exit
    return(True);
  }

  value_init(LB);
  value_init(UB);
  if(lower_upper_bounds(position, scan, val, &LB, &UB)) {
    #ifdef SIMPLIFY_DEBUG2
    fprintf(stderr, "  polyhedron_int_solution: lower_upper_bounds not zero\n");
    #endif
    // infinite, should never happen here (should not call the scan)
    fprintf(stderr, "polyhedron_int_solution: infinite lower/upper bound\n");
    value_clear(UB);
    value_clear(LB);
    return(True); // True if infinity
  }
  else {
    #ifdef SIMPLIFY_DEBUG2
    fprintf(stderr, "  polyhedron_int_solution: LB/UB = ");
    value_print(stderr, P_VALUE_FMT, LB);
    value_print(stderr, P_VALUE_FMT, UB);
    fprintf(stderr, "\n");
    #endif
    // LB increment till UB:
    for(; value_le(LB, UB); value_increment(LB, LB)) {
      value_assign(val[position], LB);
      if(polyhedron_int_solution(scan->next, val, position + 1)) {
        // early exit as soon as an integer point is found
        value_clear(UB);
        value_clear(LB);
        return(True);
      }
    }
    // reset val[position] to 0 after the scan, to prepare for next scans
    value_set_si(val[position], 0);
  }

  // cleanup
  value_clear(UB);
  value_clear(LB);
  return(False);
}

/*
 * Remove polyhedra that have no integer solutions from a domain.
 * Return a polyhedral domain, reusing the memory of D (do not free)
 * 
 * Note: DomainConstraintSimplify() has been called already before entering
 * this function.
 * 
 * Careful with unbounded polyhedra in this domain, we compute a bounding box
 * for them
 */
static Polyhedron *Domain_Remove_Integer_Empty(Polyhedron *D)
{
  Polyhedron *result = NULL, *universe = NULL;
  Vector *vec = NULL;

  #ifdef SIMPLIFY_DEBUG
  fprintf(stderr, "--- Entering Domain_Remove_Integer_Empty\n");
  #endif

  // this should always be done already:
  // D = DomainConstraintSimplify(D, MAXNOOFRAYS);

  // scan each polyhedron, if no integer solution eliminate it
  // (do not copy to result, free memory)
  while(D) {
    Polyhedron *scan = NULL, *next;
    int ray;
    Bool int_solution_found = False;
    int unbounded = 0;

    next = D->next;
    D->next = NULL;

    if(emptyQ(D)) { // D is rational-empty, eliminate.
      Polyhedron_Free(D);
      D = next;
      continue;
    }
    #ifdef SIMPLIFY_DEBUG
    fprintf(stderr, "Checking integer emptiness of: ");
    Polyhedron_Print(stderr, P_VALUE_FMT, D);
    #endif

    for(ray = 0; ray < D->NbRays; ray++) {
      // check if there's an integer vertex
      if(value_notzero_p(D->Ray[ray][0]) &&
         value_one_p(D->Ray[ray][D->Dimension + 1])) {
          // One of the vertices is integer then it's not empty, algo done
          int_solution_found = True;
        break;
      }
      // check if there's a line or ray:
      if(value_zero_p(D->Ray[ray][0]) ||
         value_zero_p(D->Ray[ray][D->Dimension + 1])) {
        unbounded = 1;
      }
    }

    // Here, check if dark shadow to 0-dim space is non empty.
    Polyhedron *dark = D;
    // compute dark shadow (source, no need to project) along every dimensions:
    for(int d = 0; d < D->Dimension; d++) {
      Polyhedron *temp = polyhedron_dark_source(dark, d);
      if(temp) {
        if(dark != D)
          Polyhedron_Free(dark);
        dark = temp;
      }
    }
    // is the dark shadow empty ?
    if(! emptyQ(dark)) {
      #ifdef SIMPLIFY_DEBUG
      fprintf(stderr, "The dark shadow to 0-dim is not empty! dark = ");
      Polyhedron_Print(stderr, P_VALUE_FMT, dark);
      fprintf(stderr, "Integer solution found = True");
      #endif
      int_solution_found = True;
    }
    if(dark != D) {
      Polyhedron_Free(dark);
    }


    if(!int_solution_found) {
      // if there is no obvious integer solution, try a scan.

      // but first, bound the polyhedron if it was not.
      Polyhedron *Dbounded = (unbounded)?bound_polyhedron(D):D;
      #ifdef SIMPLIFY_DEBUG
      if(unbounded) {
        fprintf(stderr, "bounded D= ");
        Polyhedron_Print(stderr, P_VALUE_FMT, Dbounded);
      }
      #endif

      // allocate memory if not done yet
      if(!vec)
      {
        vec = Vector_Alloc(Dbounded->Dimension + 2);
        Vector_Set(vec->p, 0, Dbounded->Dimension + 1);
        value_set_si(vec->p[Dbounded->Dimension + 1], 1);
        universe = Universe_Polyhedron(0);
      }

      scan = Polyhedron_Scan(Dbounded, universe, MAXNOOFRAYS);
      #ifdef SIMPLIFY_DEBUG
      fprintf(stderr, " scan = ");
      Polyhedron_Print(stderr, P_VALUE_FMT, scan);
      #endif
      int_solution_found = polyhedron_int_solution(scan, vec->p, 1);

      if(unbounded)
        Polyhedron_Free(Dbounded);
      Domain_Free(scan);

      // If found, then vec->p contains an integer solution of D,
      // use it to simplify D
      // -> at least we can set an integer bound on its first dimension
      // -> or can we force the vertex to be in D?
      //    splitting the polyhedron in parts? -> can be very complex in
      //    higher dimensions, don't!
      if(int_solution_found) {
        Polyhedron *tmp;
        // at least we can cut the first dimension of D to the lower bound
        // reuse vec to build the constraint on first dimension to be added:
        value_set_si(vec->p[0], 1);                     // inequality
        value_oppose(vec->p[vec->Size - 1], vec->p[1]); // constant = -vec[1]
        value_set_si(vec->p[1], 1);                     // vec[1] = 1
        Vector_Set(vec->p + 2, 0, vec->Size - 3);       // 0's everywhere else
        tmp = AddConstraints(vec->p, 1, D, MAXNOOFRAYS);
        Polyhedron_Free(D);
        D = tmp;
        value_set_si(vec->p[vec->Size - 1], 1); // restore vec constant
      }
    }

    if(int_solution_found) {
      // link polyhedron D to result
      D->next = result;
      result = D;
      #ifdef SIMPLIFY_DEBUG
      fprintf(stderr, "-> not empty\n");
      #endif
    }
    else {
      // free the polyhedron (next was saved above)
      Polyhedron_Free(D);
      #ifdef SIMPLIFY_DEBUG
      fprintf(stderr, "-> empty\n");
      #endif
    }

    D = next;
  } // while(D)

  if(vec) {
    // free memory used for scan, if allocated (vec and universe together)
    Polyhedron_Free(universe);
    Vector_Free(vec);
  }

  #ifdef SIMPLIFY_DEBUG
  fprintf(stderr, "--- Exit Domain_Remove_Integer_Empty with result = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, result);
  #endif
  return(result);
}

/*
 * Return a new domain that is D expanded such that the dimension 'move' is
 * moved to a new dimension (right) and the existing one is left empty
 * (add a line along this dimension)
 * 1 <= move <= D->Dimension + 1.
 * 
 * In D we have constraints/rays:
 * flag d_1 d_2 ... move-1 move move+1 ... d_dim constant
 * the function will build:
 * flag d_1 d_2 ... move-1  0   move+1 ... d_dim move constant
 * 
 * if move == D->Dimension + 1, this function will just just add a dimension.
*/
static Polyhedron *domain_insert_dim(Polyhedron *D, int move)
{
  Polyhedron *R = NULL; // result
  #ifdef SIMPLIFY2_DEBUG
  fprintf(stderr, "Entering domain_insert_dim. move = %d\n", move);
  // Polyhedron_Print(stderr, P_VALUE_FMT, D);
  #endif
  for(Polyhedron *p = D; p; p = p->next) {
    Polyhedron *new;
    new = Polyhedron_Alloc(p->Dimension + 1, p->NbConstraints, p->NbRays + 1);
    new->NbBid = p->NbBid + 1;
    new->NbEq = p->NbEq;
    new->flags = p->flags;
    // copy+expand Constraint matrix:
    for(int c = 0; c < p->NbConstraints; c++) {
      // linear part
      Vector_Copy(p->Constraint[c], new->Constraint[c], p->Dimension + 1);
      // constant
      value_assign(new->Constraint[c][new->Dimension + 1],
        p->Constraint[c][p->Dimension + 1]);
      // new dim and move
      value_assign(new->Constraint[c][new->Dimension], p->Constraint[c][move]);
      value_set_si(new->Constraint[c][move], 0);
      }
    // add new line (0 0 ... 0 1 0 ... 0) in new->Ray[0]
    Vector_Set(new->Ray[0], 0, new->Dimension + 2);
    value_set_si(new->Ray[0][move], 1);
    // copy+expand Ray matrix (in new->Ray[1+])
    for(int r = 0; r < p->NbRays; r++) {
      // linear part
      Vector_Copy(p->Ray[r], new->Ray[r+1], p->Dimension + 1);
      // constant
      value_assign(new->Ray[r+1][new->Dimension + 1],
        p->Ray[r][p->Dimension + 1]);
      // new dim and move
      value_assign(new->Ray[r+1][new->Dimension], p->Ray[r][move]);
      value_set_si(new->Ray[r+1][move], 0);
    }
    // link new to result
    new->next = R;
    R = new;
  }
  #ifdef SIMPLIFY2_DEBUG
  // fprintf(stderr, "Exiting domain_insert_dim. Result = ");
  // Polyhedron_Print(stderr, P_VALUE_FMT, R);
  #endif

  return(R);
} /* domain_insert_dim */


/*
 * Fonction to expand the lattice of LBL A such that it is equal to the ref
 * lattice (ignore zero columns). The ref lattice includes the lattice
 * A->Lat.
 * This function will expand the domain D to the required dimension and add
 * equalities, such that the same LBL points are spread (does the opposite
 * of remove_equalities)
 * In place: modifies A.
 */
static void sLBLMake_lattice_equal_to(LBL *A, Matrix *ref)
{
  // A->Lat is included in ref, so we know that
  // the pivots of A are multiple of the pivots of ref.

  // some examples:
  /*
  merging A = LBL: Dimension 2

  LATTICE:
  3 3
    6    0    4  -> 6i+0j+4,  transf. into:
                    2i'+0j+0 and \exists i' such that 2i' = 6i+4
    0    1    0
    0    0    1
  [included in] tmp = LBL: Dimension 2

  LATTICE:
  3 3
    2    0    0
    0    1    0
    0    0    1
  */

  // a more complex example:
  /*
  LATTICE:
  4 4
    2    0    0    0  -> 2i + 0j + 0, transf. into:
                         i' + 0j + 0 and \exists i' with i'=2i
    0    0    0    0
    1    3    0    1  -> i + 3j + 1, transf. into:
                             1j'+ 0 -> \exists j' with j'=i+3j+1
    0    0    0    1

  [included in] tmp = LBL: Dimension 3

  LATTICE:
  4 4
    1    0    0    0

    0    0    0    0
    0    1    0    0
    0    0    0    1
  */

  // and a difficult one:
  /*
  LATTICE:
  4 4
    2    0    0    0  -> same i = i'
    5   15    0   10  -> 5i+15j+10 transf. into:
                             5j'    \exists j', 5j'=5i+15j+10
                                    and a domain transfo
    13   15   66   54 -> 13i + 15j + 66k + 54 transf. into:
                           52i + 27j' + 66k' + 0 -> same k = k'

                           (27j' = 27i + 3*27j + 54 so the line becomes
                            79i + 81j + 66k + 54, do modulo 66 =
                            13i + 15j + 66k + 54 that is correct.
                            )
    0    0    0    1

  [included in] tmp = LBL: Dimension 3
  LATTICE:
  4 4
    2    0    0    0
    0    5    0    0
    52   27   66    0
    0    0    0    1

  */

  // dimension reduction:
  /*
  LATTICE:
  3 2
    0    3  -> add column (0 0 0)^T and add new dimension with equality {3i=3}
    1    0                              (same procedure as above!)
    0    1
  [included in] tmp = LBL: Dimension 2
  LATTICE:
  3 3
    3    0    0
    0    1    0
    0    0    1
  */
  // A included in ref
  // the pivot of A->Lat is a multiple of the pivot of ref (if exists)

  // search the pivot in ref, it's not necessarily the same columns than
  // A->Lat since nbcolumns(A->Lat) can be different than nbcolumns(ref)

  int nb_added_dims = 0;
  Matrix *L = A->Lat;
  #ifdef SIMPLIFY2_DEBUG
  fprintf(stderr, "Entering sLBLMake_lattice_equal_to\n");
  #endif

  for(int col = 0; col < ref->NbColumns - 1; col++) {
    int row;
    Polyhedron *expand, *new;

    // search pivot of column col
    for(row = 0; row < ref->NbRows ; row++) {
      if(value_notzero_p(ref->p[row][col])) {
        break;
      }
    }
    if(row == ref->NbRows) {
      // no more pivots (zero column), exit the loop.
      break;
    }

    // no corresponding pivot in L!
    if(col == L->NbColumns - 1 || value_zero_p(L->p[row][col])) {
      // need to add a column in L
      nb_added_dims--;

      Matrix *newL;
      newL = Matrix_Alloc(L->NbRows, L->NbColumns + 1);
      for(int r = 0; r < L->NbRows; r++) {
        // column number col is just a zero column
        Vector_Copy(L->p[r], newL->p[r], col); // 0 to col-1
        value_set_si(newL->p[r][col], 0);      // col
        Vector_Copy(L->p[r] + col, newL->p[r] + col + 1, L->NbColumns - col);
      }
      Matrix_Free(A->Lat);
      A->Lat = newL;
      L = newL;
    }
    // now pivot is in position [row][col] of both matrices.
    // if same row, just ignore, goto next
    else {
      int c;
      for(c = 0; c <= col ; c++) {
        if(value_ne(L->p[row][c], ref->p[row][c])) {
          break;
        }
      }
      if(c == col + 1) {
        // all values are equal
        // if constant is equal too, continue to next pivot
        if(value_eq(L->p[row][L->NbColumns-1], ref->p[row][ref->NbColumns-1]))
        {
          continue;
        }
      }
    }

    nb_added_dims++;
    // update A->P (if not NULL)
    // expand one dimension of the domain (the pivot one)
    // and add equality L[row] = ref[row]
    // the original dimension should be moved to the right (eliminated) column,
    // such that the new lattice is just ref
    if(A->P)
    {
      Matrix *new_equality = Matrix_Alloc(1, A->P->Dimension + 2 + 1);

      expand = domain_insert_dim(A->P, col + 1);
  
      // build equality L[row] = ref[row],
      // (the dimension col+1 of L moved to the new dimension)
      Vector_Set(new_equality->p[0], 0, new_equality->NbColumns);
      // (1) copy linear part of L:
      Vector_Copy(L->p[row], new_equality->p[0] + 1, L->NbColumns - 1);
      // constant
      value_assign(new_equality->p[0][new_equality->NbColumns-1],
        L->p[row][L->NbColumns-1]);
      // move the pivot to the new dimension
      value_assign(new_equality->p[0][new_equality->NbColumns-2],
        new_equality->p[0][col+1]);
      // and replace it with 0.
      value_set_si(new_equality->p[0][col+1], 0);
      // (2) substract linear part of ref[row] from the equality:
      for(int c = 0; c < ref->NbColumns-1 && c < new_equality->NbColumns-1; c++) {
        value_substract(new_equality->p[0][c+1], new_equality->p[0][c+1],
          ref->p[row][c]);
      }
      // substract constant:
      value_substract(new_equality->p[0][new_equality->NbColumns-1],
        new_equality->p[0][new_equality->NbColumns-1],
        ref->p[row][ref->NbColumns-1]);
      #ifdef SIMPLIFY2_DEBUG
      fprintf(stderr, "Adding equality: ");
      Matrix_Print(stderr, P_VALUE_FMT, new_equality);
      #endif
  
      // build the final domain and update A->P
      new = DomainAddConstraints(expand, new_equality, MAXNOOFRAYS);
      // Domain_Free(expand);
      Domain_Free(A->P);
      A->P = new;
      Matrix_Free(new_equality);
    }
  }

  // update A->Lat:
  // replace Lat with ref
  // and add zero columns for the extra dimensions
  Matrix *Lat;
  // Lat = Matrix_Alloc(A->Lat->NbRows, A->P->Dimension + 1);
  // new nbcolumns = nbcolumns of A->Lat + number of added dimensions
  Lat = Matrix_Alloc(A->Lat->NbRows, A->Lat->NbColumns + nb_added_dims);
  for(int r = 0; r < Lat->NbRows; r++) {
    if(Lat->NbColumns > ref->NbColumns) {
      // linear part:
      Vector_Copy(ref->p[r], Lat->p[r], ref->NbColumns - 1);
      // add zero columns:
      Vector_Set(Lat->p[r] + ref->NbColumns - 1, 0,
        Lat->NbColumns - ref->NbColumns);
    }
    else {
      // in case ref has more zero columns than the new Lat:
      Vector_Copy(ref->p[r], Lat->p[r], Lat->NbColumns - 1);
    }
    // constant:
    value_assign(Lat->p[r][Lat->NbColumns - 1], ref->p[r][ref->NbColumns - 1]);
  }
  Matrix_Free(A->Lat);
  A->Lat = Lat;
  #ifdef SIMPLIFY2_DEBUG
  fprintf(stderr, "New Lat: ");
  Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
  #endif

  // just normalize lattice (+update domain) without removing equalities
  sLBL_Lat_Normalize(A);

  #ifdef SIMPLIFY2_DEBUG
  fprintf(stderr, "Exit sLBLMake_lattice_equal_to. A = ");
  sLBLPrint(stderr, P_VALUE_FMT, A);
  #endif
} /* sLBLMake_lattice_equal_to */


/***************************************************************************/
/*       CanonicalLBL, LBL2ZDomain, and LBLSimplify                        */
/***************************************************************************/

/*
 * Modify the single LBL 'A' to be in canonical form:
 * A->Lat in HNF and no equalities in A->P.
 * Also tries to remove the columns of zeros from A->Lat if possible:
 * do the projection along those dimensions and eliminate only if
 * dark shadow = exact shadow
 *
 * USAGE: in place, modifies A itself.
 * 
 * This function can leave an arbitrary number of empty LBLs in the list
 * (will be removed later on).
 * Can modify A->next to insert new sub-LBLs.
 */
static void sLBLCanonical(LBL *A)
{
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Entering sLBLCanonical\n");
    fprintf(stderr, "--------- Input Lat: ");
    Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
    fprintf(stderr, "--------- Input P: ");
    Polyhedron_Print(stderr, P_VALUE_FMT, A->P);
  #endif
  
  if (A->P) {
    if( A->P->Dimension+1 != A->Lat->NbColumns) {
      errormsg1("sLBLCanonical", "dimincomp", "incompatible dimensions");
      return;
    }
  }

  // ************************
  // STEP 1: normalize A->Lat
  // ************************

  // Normalize the affine lattice of A (update A->Lat and A->P)
  sLBL_Lat_Normalize(A);

  // simplify non-integer constraints such that they intersect at least one
  // integer point (to avoid some infinite integer-empty polyhedra and
  // obviously empty polyhedra)
  A->P = DomainConstraintSimplify(A->P, MAXNOOFRAYS);

  // check emptiness
  if(!A->P || emptyQ(A->P)) {
    Domain_Free(A->P);
    A->P = NULL;
    // set the lattice to a single empty column
    int dimension = A->Lat->NbRows - 1;
    if(A->Lat->NbColumns > 1) {
      Matrix_Free(A->Lat);
      A->Lat = Matrix_Alloc(dimension + 1, 1);
      for(int j = 0 ; j < dimension; j++) {
        value_set_si(A->Lat->p[0][j], 0);
      }
      value_set_si(A->Lat->p[0][dimension], 1);
    }
    return;
  }

  // ***********************************
  // STEP 2: remove equalities from A->P
  // ***********************************

  // homogenize the equalities of A->P: ensure that all polyhedra of the
  // domain verify the same set of equalities
  Matrix *Equalities = sLBLHomogenize_equalities(A);

  if (A->P->Dimension > 0 && A->P->NbEq != 0) {
    // Remove the equalities from A->P
    sLBLSimplify_equalities(A, Equalities);

    // some equalities were eliminated. Do we need to start again from scratch?
    // Lat and P changed, but Lat is kept in HNF, so no.
    // Matrix_Free(Equalities);
    // sLBLCanonical(A);
    // return;
  }
  Matrix_Free(Equalities);
  if(!A->P) {
    // empty, done.
    return;
  }

  // ****************************************
  // STEP 3: eliminate zero columns of A->Lat
  // ****************************************

  // Remove the columns of zeros from A->Lat if possible:
  // do the projection along the zero-dimensions,
  // eliminate only if dark shadow \in exact shadow
  // (do not convert into ZDomains/check for integer-emptiness of polyhedra)
  sLBL_Simplify_Zero_Dimensions(A);

  return;
} /* sLBLCanonical */


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
void CanonicalLBL(LBL *A)
{
  #ifdef CANONICAL_DEBUG
  fprintf(stderr, "Entering CanonicalLBL.\nA =");
  LBLPrint(stderr, P_VALUE_FMT, A);
  #endif

  // transform every LBL of the list individually
  for(LBL *tmp = A; tmp; tmp = tmp->next) {
    sLBLCanonical(tmp);
  }

  // remove empty LBLs from the list, keep at least something in A
  LBL_Remove_Empty(A);

  #ifdef CANONICAL_DEBUG
  fprintf(stderr, "sLBLCanonical and RemoveEmpty done.\nA =");
  LBLPrint(stderr, P_VALUE_FMT, A);
  #endif
  
  // check if a lattice is present twice in A, and if it is, union the second
  // polyhedral domain with the first one and remove the second sLBL
  for(; A; A = A->next)
  {
    LBL *previous = A, *current = A->next;
    while(current) {
      if(isEqualLattice(A->Lat, current->Lat)) {
        // move current->P to A->P, remove current, and relink previous to next
        Polyhedron *pp = current->P;
        while(pp) {
          Polyhedron *nextpp = pp->next;
          pp->next = NULL;
          A->P = AddPolyToDomain(pp, A->P);
          pp = nextpp;
        }
        Matrix_Free(current->Lat);
        previous->next = current->next;
        free(current);
        current = previous->next;
        // and previous does not change
      }
      else {
        previous = current;
        current = current->next;
      }
    }
  }
} /* CanonicalLBL */


/*
 * Transform a single LBL into a list of Z-domains
 *
 * Remove zero columns from the lattice, build a union of Z-polyhedra
 * 
 */
static LBL *sLBL2ZDomain(LBL *A)
{
  int nbzeros;
  LBL *Result;

  if((nbzeros = LatCountZeroCols(A->Lat))) {
    // there are potential holes
    Matrix *newL;
    Polyhedron *holes, *not_holes, *exact;
    holes = sLBLCompute_holes(A, &exact);

    not_holes = DomainDifference(exact, holes, MAXNOOFRAYS);
    #ifdef HOLES_DEBUG
    fprintf(stderr, "-- in LBL2ZDomain - Not holes = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, not_holes);
    #endif
    Domain_Free(holes);
    Domain_Free(exact);

    // build result LBL
    newL = RemoveNColumns(A->Lat, A->Lat->NbColumns-1-nbzeros, nbzeros);
    Result = LBLAlloc(newL, not_holes);
    Matrix_Free(newL);
    Domain_Free(not_holes);
  }
  else {
    Result = sLBLCopy(A);
  }
  return(Result);
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
  LBL *Result = sLBL2ZDomain(A);
  for(LBL *Z = A->next; Z; Z = Z->next)
  {
    LBL *tmp;
    tmp = sLBL2ZDomain(Z);
    Result = LBLConcatenate(tmp, Result);
  }
  CanonicalLBL(Result);
  LBLSimplifyEmpty(Result);
  return(Result);
}


/*
 * Remove integer empty single LBLs from an LBL, in place.
 *
 * Algorithm:
 * Check all polyhedra for integer-emptiness, remove if a polyhedron does not
 *    contain any integer point
 */
void LBLSimplifyEmpty(LBL *A)
{
  // just for testing, call simplify from simplifyEmpty:
  #ifdef SIMPLIFY2_DEBUG
  LBLSimplify(A);
  #endif

  for(LBL *tmp = A; tmp; tmp = tmp->next) {
    tmp->P = Domain_Remove_Integer_Empty(tmp->P);
  }
  LBL_Remove_Empty(A); // remove emptied LBLs from the list

  CanonicalLBL(A);
} /* LBLSimplifyEmpty */

/*
 * Simplify an LBL, in place.
 *
 * Algorithm:
 * Merge all the lattices that can be merged. Ideas... :
 *    a- merge same lattices, ignore zero columns (extend dimension of P) -> easy.
 *    b- check lattice inclusion and merge?
 *       (involves a-) but how to merge?
 *         only if poly included too?
 *         or separate the part that is not included?
 *         -> or at least their hulls should intersect?
 *       if merge, transform the lattices to be equal:
 *       add equalities to P in order to generate the original points
 *         of L x | x \in P
 *    c- sort lattices by linear part...
 *            and try to merge same ones..., when domains overlap?
 *             or split domains such that a part of it does overlap?
 *       example: can (2i+0) be merged with (2i+1)? -> yes if P is same (but this is unlikely)
 * let A = im((2i+0),P), and B = im((2i+1),P') -> if A = B then (i+0), preim((i+0), A) -> Z-pol only

  another idea: make all lattices = (Id 0), adding existential variables in the domains to generate the right lbls
  then merge everything together
  -> can be pretty complex..., but obviously covers all cases!

  lattice inclusion test covers this case too, if the lattice Id is in the LBL list
  and the hulls intersect.

 * Then:
 * - fuse/simplify all adjacent polyhedral domains (complement of simplify of complement)
 * - CanonicalLBL to remove equalities
 */
void LBLSimplify(LBL *A)
{
  // LBL A should be canonical (all functions return a canonical LBL)
  if(isEmptyLBL(A))
    return;

  // TESTING, add lattice Id as first lattice in LBL A (with an empty domain).
  // this ensures that all lattices are merged together!

  // LBL *first;
  // first = malloc(sizeof(LBL));
  // // copy A in first
  // *first = *A;
  // // set A to (Z^d, empty)
  // A->next = first;
  // A->Lat = Identity_Matrix(A->Lat->NbRows);
  // A->P = NULL;

  // check if a lattice of A is included in another, merge them if yes.
  // warning, need to modify A while scanning it: check inclusion both ways
  // and always merge with the first one (suppress the last one)
  for(LBL *tmp = A; tmp; tmp = tmp->next) {
    LBL *current = tmp->next, *prec = tmp;
    while(current) {
      int flag = 0; // if flag, current included in tmp

      if(isSameLatticeSpace(current->Lat, tmp->Lat)) {
        // the two lattices do not have the same zero columns but are equal
        flag = 2;
      }
      else if(LatticeIncluded(current->Lat, tmp->Lat)) {
        // current is included in tmp
        flag = 1;
      }
      else if(LatticeIncluded(tmp->Lat, current->Lat)) {
        Matrix *L;
        Polyhedron *P;
        // tmp is included in current
        flag = 1;
        // exchange tmp and current, so that current is included in tmp.
        L = tmp->Lat;            P = tmp->P;
        tmp->Lat = current->Lat; tmp->P = current->P;
        current->Lat = L;        current->P = P;
      }
      if(flag) {
        #ifdef SIMPLIFY2_DEBUG
        fprintf(stderr, "merging current = ");
        sLBLPrint(stderr, P_VALUE_FMT, current);
        fprintf(stderr, "[included in] tmp = ");
        sLBLPrint(stderr, P_VALUE_FMT, tmp);
        #endif
        if(flag == 1) {
          // make the lattice of current equal to the one of tmp
          // apart from zero columns
          // (current is included in tmp)

          sLBLMake_lattice_equal_to(current, tmp->Lat);
        }

        // after this, merge the zero columns such that the two lattices are
        // perfectly equal
        // can require to swap dimensions to align same existential variables,
        // or to increase dimension to take in all of them (lazy) <- do that

        int current_nbzero, tmp_nbzero;
        Polyhedron *newP;
        Matrix *newL;
        // expand domains current and tmp:
        current_nbzero = LatCountZeroCols(current->Lat);
        tmp_nbzero = LatCountZeroCols(tmp->Lat);
        if(current->P) {
          for(int i = 0; i < tmp_nbzero; i++) {
            // add tmp_nbzero dimensions to current->P
            newP = domain_insert_dim(current->P, current->P->Dimension + 1);
            Domain_Free(current->P);
            current->P = newP;
          }
        }
        if(tmp->P) {
          for(int i = 0; i < current_nbzero; i++) {
            // add current_nbzero dimensions to tmp->P
            assert(tmp->P);
            newP = domain_insert_dim(tmp->P, tmp->P->Dimension + 1);
            Domain_Free(tmp->P);
            tmp->P = newP;
          }  
        }

        // merge current in tmp:
        // - update tmp->P
        newP = DomainUnion(tmp->P, current->P, MAXNOOFRAYS);
        Domain_Free(tmp->P);
        tmp->P = newP;

        // - update tmp->Lat
        newL = Matrix_Alloc(tmp->Lat->NbRows,
          tmp->Lat->NbColumns + current_nbzero);
        for(int r = 0; r < newL->NbRows; r++) {
          // linear part
          Vector_Copy(tmp->Lat->p[r], newL->p[r], tmp->Lat->NbColumns - 1);
          // zeros
          Vector_Set(newL->p[r] + tmp->Lat->NbColumns - 1, 0, current_nbzero);
          // constant
          value_assign(newL->p[r][newL->NbColumns - 1],
            tmp->Lat->p[r][tmp->Lat->NbColumns - 1]);
        }
        Matrix_Free(tmp->Lat);
        tmp->Lat = newL;
        #ifdef SIMPLIFY2_DEBUG
        fprintf(stderr, "[new merged] tmp = ");
        sLBLPrint(stderr, P_VALUE_FMT, tmp);
        #endif

        // remove current, update prec->next and current
        Domain_Free(current->P);
        Matrix_Free(current->Lat);
        prec->next = current->next;
        free(current);
        current = prec->next;
      }
      else {
        // simple advance prec and current
        prec = current;
        current = current->next;
      }
    }
  }
  // now the the polyhedra have their equalities back and are merged together,
  // what to do with them??


  // THIS IS EASY BUT DOES NOT HELP :(
  for(LBL *tmp = A; tmp; tmp = tmp->next) {
    tmp->P = DomainConstraintSimplify(tmp->P, MAXNOOFRAYS);
  }


  // // THIS IS TOO COMPLEX! -> ZDiff11 takes ages...
  // Polyhedron *universe;
  // universe = Universe_Polyhedron(A->P->Dimension);

  // for(LBL *tmp = A; tmp; tmp = tmp->next) {
  //   Polyhedron *newP;

  //   // compute complement -> complement
  //   // first complement
  //   newP = DomainDifference(universe, tmp->P, MAXNOOFRAYS);
  //   Domain_Free(tmp->P);
  //   tmp->P = newP;
  
  //   // disjoint of complement + simplify
  //   // TOO COMPLEX:
  //   // tmp->P = Disjoint_Domain(tmp->P, 0, MAXNOOFRAYS);
  //   // tmp->P = DomainConstraintSimplify(tmp->P, MAXNOOFRAYS);
  
  //   // and complement + simplify:
  //   newP = DomainDifference(universe, tmp->P, MAXNOOFRAYS);
  //   Domain_Free(tmp->P);
  //   tmp->P = newP;
  //   tmp->P = DomainConstraintSimplify(tmp->P, MAXNOOFRAYS);
  // }


  // // THIS IS NOT WHAT WE WANT!
  // // but could be useful to count?
  // for(LBL *tmp = A; tmp; tmp = tmp->next) {
  //   Polyhedron *newP;

  //   // disjoint + simplify
  //   newP = Disjoint_Domain(tmp->P, 0, MAXNOOFRAYS);
  //   Domain_Free(tmp->P);
  //   tmp->P = DomainConstraintSimplify(newP, MAXNOOFRAYS);
  // }

  // TODO: simplify the polyhedral domains stored in A to get a minimal
  // and disjoint(?) form


  #ifdef SIMPLIFY2_DEBUG
  fprintf(stderr, "\n Exit LBLSimplify. A (before simplify) = ");
  LBLPrint(stderr, P_VALUE_FMT, A);
  #endif
  CanonicalLBL(A);

  // avoid recursive loop, just for testing:
  #ifndef SIMPLIFY2_DEBUG
  LBLSimplifyEmpty(A);
  #endif
} /* LBLSimplify */


/***************************************************************************/
/*                      LBLDisjointUnion                                   */
/***************************************************************************/

/*
 * Compute the disjoint union of LBL A, by performing succesive intersections
 * and differences.
 */
LBL *LBLDisjointUnion(LBL *A)
{
  if(!A)
    return(NULL);

  // TODO: infinite loop for ZDisj11b.in

  // TODO: split lattices first, then polyhedra in them

  // need a ZDomain version?
  // A = LBL2ZDomain(A);

  // get a copy of the first domain of A
  LBL *res = sLBLCopy(A);

  // then compute the inter/diff of each sLBL stored in A
  for(LBL *tmpA = A->next; tmpA; tmpA = tmpA->next) {
    LBL *next, *inter, *diffA, *diffB, *newres;

    // unlink next
    next = tmpA->next;
    tmpA->next = NULL;

    // compute intersections and differences
    inter = LBLIntersection(res, tmpA);
    diffA = LBLDifference(tmpA, res);
    diffB = LBLDifference(res, tmpA);

    newres = LBLConcatenate(LBLConcatenate(inter, diffA), diffB);
    CanonicalLBL(newres);
    LBLFree(res);
    res = newres;

    // relink next
    tmpA->next = next;
  }

  return(res);
}

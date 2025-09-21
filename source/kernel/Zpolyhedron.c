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
#define COMP_DEBUG 1
#define HOLES_DEBUG 1
#endif

static LBL *sLBL_Intersection(LBL *, LBL *);
static LBL *sLBL_Copy(LBL *A);
static void sLBL_Free(LBL *L);
static LBL *sLBL_Difference(LBL *, LBL *);
static LBL *sLBL_Image(LBL *, Matrix *);
static LBL *sLBL_Preimage(LBL *, Matrix *);
static void sLBL_Canonical(LBL* A);
static LBL *FindLatticePred(Matrix *L, LBL *A);
static LBL *LBL_sLBL_Difference(LBL* A, LBL* B);
static int count_zeroCols (Matrix* M);
static Polyhedron *sLBL_compute_holes(LBL *A, Polyhedron **pExact);

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
    return(NULL);
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

  // TODO: could check lattice union inclusion first, if not included return fail.

  diff = LBLDifference(A, B);
  if(isEmptyLBL(diff)) {
    ret = True;
  }
  LBLFree(diff);

  return ret;
} /* LBLIncludes */


/*
 * Print the contents of an LBL 'A'
 */
void LBLPrint(FILE *fp, const char *format, LBL *A)
{
  for( ; A; A = A->next) {
    fprintf(fp, "LBL: Dimension %d\n", A->Lat->NbRows - 1);
    if(emptyQ(A->P)) {
      fprintf(fp, "\n<empty>>\n");
    }
    else {
      fprintf(fp, "\nLATTICE:\n");
      Matrix_Print(fp, format, A->Lat);
      Polyhedron_Print(fp, format, A->P);
    }
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
  LBL *Result = NULL;

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
    errormsg1("LBLIntersection", "dimincomp",
      "incompatible dimensions between domains");
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
    errormsg1("LBLDifference", "dimincomp",
        "incompatible dimensions between domains");
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


void sLBL_Print(FILE *out, char *fmt, LBL *A)
{
  LBL *next = A->next;
  A->next=NULL;
  LBLPrint(out, fmt, A);
  A->next = next;
}


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
static LBL *sLBL_Intersection(LBL *A, LBL *B) {

  LBL *Result = NULL;
  Matrix *LInter;
  Polyhedron *PInter, *ImageA, *ImageB, *PreImage;

  #ifdef INTERSECTION_DEBUG
    fprintf(stderr, "-- Entering sLBL_Intersection\nA = ");
    sLBL_Print(stderr, P_VALUE_FMT, A);
    fprintf(stderr, "B = ");
    sLBL_Print(stderr, P_VALUE_FMT, B);
  #endif
  LInter = LatticeIntersection(A->Lat, B->Lat);
  if (isEmptyLattice(LInter)) {
    Matrix_Free(LInter);
    #ifdef INTERSECTION_DEBUG
    fprintf(stderr, "Empty Lattice intersection, result = <empty>\n");
    fprintf(stderr, "-- exit sLBL_Intersection\n");
    #endif
    return (NULL);
  }
  #ifdef INTERSECTION_DEBUG
  fprintf(stderr, "Lattice intersection = LInter = ");
  Matrix_Print(stderr, P_VALUE_FMT, LInter);
  #endif

  if(count_zeroCols(A->Lat) == 0 && count_zeroCols(B->Lat) == 0 &&
    count_zeroCols(LInter) == 0)
  {
    // This works only IF there are no columns of zeros in the LBLs:
    // they are Z-polyhedra

    ImageA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
    ImageB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
    PInter = DomainIntersection(ImageA, ImageB, MAXNOOFRAYS);
    Domain_Free(ImageB);
    Domain_Free(ImageA);
    #ifdef INTERSECTION_DEBUG
    fprintf(stderr, "imageA inter imageB = PInter = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, PInter);
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
    LBLPrint(stderr, P_VALUE_FMT, Result);
    fprintf(stderr, "-- exit sLBL_Intersection\n");
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
      if(A->Lat->NbRows + B->P->NbConstraints > extra_max_rows) {
        extra_max_rows = A->Lat->NbRows + B->P->NbConstraints;
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
      extra->NbRows = extra_B_row + B->P->NbConstraints;
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
    Result = LBLAlloc(newL, newP);
    Matrix_Free(newL);
    Domain_Free(newP);
    #ifdef INTERSECTION_DEBUG
    fprintf(stderr, "Manual built intersection = ");
    sLBL_Print(stderr, P_VALUE_FMT, Result);
    fprintf(stderr, "-- exit sLBL_Intersection\n");
    #endif

    return(Result);
  }

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
 * Compute the complement of LBL A: all points z such that z is not in A.
 *
 * Algorithm:
 * Let L = A->Lat, P = A->P.
 * complement(A) = Universe() - A = union of:
 *   1- LBL (Z^d, complement hull(A)), with hull(A) = image by L of P
 *   2- LBL ((Z^d - L), hull(A)) ---- or ((Z^d - L), universe())
 *   3- holes of A
 *      if L has no zero columns -> empty
 *      = L z' such that there exist no z in A->P such that L z' = L z
 *      -> need exact shadow
 */
static LBL *sLBLComplement(LBL *A)
{
  LBL *Result = NULL;
  Polyhedron *Univ, *hullA, *comp_hullA;
  LatticeUnion *LatDiff;
  int nbzeros;
  #ifdef COMP_DEBUG
  fprintf(stderr, "\n-- Entering sLBLComplement. A = ");
  LBLPrint(stderr, P_VALUE_FMT, A);
  #endif
  
  // STEP 1: complement of hull(A)
  hullA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
  Univ = Universe_Polyhedron(hullA->Dimension);
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

  // STEP 2: lattice differences (not L) on hull(A)
  LatDiff = LatticeDifference(NULL, A->Lat);
  #ifdef COMP_DEBUG
  fprintf(stderr, "\nLatDiff: ");
  PrintLatticeUnion(stderr, P_VALUE_FMT, LatDiff);
  #endif
  // Add all Z-polyhedra to Result, applying the list of lattices on hullA
  for(LatticeUnion *lat = LatDiff; lat; lat = lat->next) {
    LBL *Ztmp;
    #ifdef COMP_DEBUG
    fprintf(stderr, "Considering Lat diff: ");
    Matrix_Print(stderr, P_VALUE_FMT, lat->M);
    #endif
    Ztmp = malloc(sizeof(LBL));
    Ztmp->Lat = lat->M;
    Ztmp->P = DomainPreimage(hullA, lat->M, MAXNOOFRAYS);
    // remove obvious simplification?
    // -> not necessary since preimage by integer function.
    // Ztmp->P = DomainConstraintSimplify(Ztmp->P, MAXNOOFRAYS);
    #ifdef COMP_DEBUG
    Ztmp->next = NULL;
    fprintf(stderr, "Adding: ");
    LBLPrint(stderr, P_VALUE_FMT, Ztmp);
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
  Domain_Free(hullA);

  // STEP 3: holes
  if((nbzeros = count_zeroCols(A->Lat))) {
    // there are potential holes
    Matrix *newL;
    Polyhedron *holes = sLBL_compute_holes(A, NULL);
    #ifdef COMP_DEBUG
    fprintf(stderr, "STEP 3 adding holes = ");
    PolyhedronPrint(stderr, P_VALUE_FMT, holes);
    #endif
    newL = RemoveNColumns(A->Lat, A->Lat->NbColumns-1-nbzeros, nbzeros);
    
    Result = LBL_concatenate(LBLAlloc(newL, holes), Result);
    Matrix_Free(newL);
    Polyhedron_Free(holes);
  }

  CanonicalLBL(Result);
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
  for(LBL *tmp = A->next; tmp; tmp = tmp->next) {
    LBL *comp, *inter;
    comp = sLBLComplement(tmp);
    inter = LBLIntersection(Result, comp);
    LBLFree(Result);
    Result = inter;
  }

  return(Result);
} /* LBLComplement */


/*
 * Return the difference of two single LBLs A and B.
 * A and B are single LBLs, but the return value can be a union of LBLs!
 * Creates a new allocated LBL union
 *
 * USAGE: only the first lattice of A and B is considered (no union),
 *        but A and B can contain a coordinate polyhedral domain (in ->P).
 * Internal function, users should use LBLDifference.
 * 
 * Algorithm:
 * -> New version: compute A inter complement(B).
 * -> Former version inspired from the method Gautam describes in his thesis,
 * modified to handle LBLs.
 */
static LBL *sLBL_Difference(LBL* A, LBL* B)
{
  LBL *Result, *Binter, *Bcomp;

  #ifdef DIFFERENCE_DEBUG
  fprintf(stderr, "-- Entering sLBL_Difference. A = ");
  LBLPrint(stderr, P_VALUE_FMT, A);
  #endif
  if (A->Lat->NbRows != B->Lat->NbRows) {
    errormsg1("sLBL_Difference", "dimincomp", "incompatible dimensions");
    return(NULL);
  }

  // treat the simple case where the LBLs do not intersect
  Binter = sLBL_Intersection(A, B); // reused below
  #ifdef DIFFERENCE_DEBUG
  fprintf(stderr, "Binter = ");
  LBLPrint(stderr, P_VALUE_FMT, Binter);
  #endif
  if(isEmptyLBL(Binter)) {
    // if B does not intersect A, return A.
    #ifdef DIFFERENCE_DEBUG
    fprintf(stderr,
      "Binter=(A inter B) is empty, so B does not intersect A, we return A\n");
    #endif
    LBLFree(Binter);
    return(LBLCopy(A));
  }

  // // Separate the computation in 3 phases:
  // // 0. compute the difference of the image polyhedra P_A \ P_B (=ImDiff) and
  // //    add it to the solution LBL (with lattice L_A).
  // //    This can be an over-approximation of A if A->Lat has zero columns
  // //    (but not of B)
  // // 1. compute the rest where the intersection of P_A and P_B have same
  // //    dimensions (required for lattice difference)
  // // 2. intersect the result with A to get rid of the over-approximations

  // LBL *Result = NULL, *Final_Result; // U. of LBLs
  // LBL *Ainter, *Binter; // single LBL
  // LatticeUnion *LatDiff;
  // Polyhedron *imA, *imB, *preimA, *ImDiff, *ImInter; // polyhedral domains

  // // [STEP 0 (includes Gautam's Step 2)]
  // imA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
  // imB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
  // ImDiff = DomainDifference(imA, imB, MAXNOOFRAYS);
  // #ifdef DIFFERENCE_DEBUG
  //   fprintf(stderr, "ImDiff (hull of A that does not cover B) = ");
  //   Polyhedron_Print(stderr, P_VALUE_FMT, ImDiff);
  // #endif

  // // Add (A->Lat, A->P - hull(B)) to the result:
  // if (!emptyQ(ImDiff)) {
  //   Polyhedron *RedPolyDiff;
  //   RedPolyDiff = DomainPreimage(ImDiff, A->Lat, MAXNOOFRAYS);
  //   // NOTICE: this can be an over-approximation of A
  //   Result = LBLAlloc(A->Lat, RedPolyDiff);
  //   #ifdef DIFFERENCE_DEBUG
  //     fprintf(stderr, "Adding this to the temporary result: ");
  //     LBLPrint(stderr, P_VALUE_FMT, Result);
  //   #endif
  //   Domain_Free(RedPolyDiff);
  // }

  // // compute the images intersection of A and B
  // ImInter = DomainIntersection(imA, imB, MAXNOOFRAYS);
  // #ifdef DIFFERENCE_DEBUG
  //   fprintf(stderr, "ImInter (hull of A inter B) = ");
  //   Polyhedron_Print(stderr, P_VALUE_FMT, ImInter);
  // #endif
  
  // // (TODO) can be simplified, Ainter not really needed!

  // // compute the part of A that intersects the hull of B in the image space
  // preimA = DomainPreimage(ImInter, A->Lat, MAXNOOFRAYS);
  // Ainter = LBLAlloc(A->Lat, preimA);
  // // NOTICE: this Ainter can be a over-approximation of A

  // Domain_Free(preimA);
  // Domain_Free(ImDiff);
  // Domain_Free(imA);
  // Domain_Free(imB);

  // // now Ainter and Binter have same lattices and polyhedra dimensions
  // #ifdef DIFFERENCE_DEBUG
  //   fprintf(stderr,
  //     "-- [STEP1] now we compute the intersection on same lattice dimensions\n");
  //   fprintf(stderr, "Ainter = ");
  //   LBLPrint(stderr, P_VALUE_FMT, Ainter);
  //   fprintf(stderr, "and Binter = ");
  //   LBLPrint(stderr, P_VALUE_FMT, Binter);
  // #endif

  // // LatDiff (union of lattices) is the difference : (A->Lat) - (B->Lat) of
  // // same dimensions
  // LatDiff = LatticeDifference(Ainter->Lat, Binter->Lat); 
  // #ifdef DIFFERENCE_DEBUG
  //   if(!LatDiff)
  //     fprintf(stderr, "Empty Lattice difference\n");
  // #endif

  // // [STEP 1 of Gautam]:
  // // Add all Z-polyhedra applying the (list of) lattice difference on ImInter
  // for(LatticeUnion *tmp = LatDiff; tmp; tmp = tmp->next) {
  //   LBL *Ztmp;
  //   #ifdef DIFFERENCE_DEBUG
  //     fprintf(stderr, "Considering Lat diff: ");
  //     Matrix_Print(stderr, P_VALUE_FMT, tmp->M);
  //   #endif
  //   Ztmp = malloc(sizeof(*Ztmp));
  //   Ztmp->next = Result;
  //   Ztmp->Lat = tmp->M;
  //   Ztmp->P = DomainPreimage(ImInter, tmp->M, MAXNOOFRAYS);
  //   // NOTICE: this can be an over-approximation of A (but not of B)

  //   Result = Ztmp;
  // }
  // // free LatticeUnion remaining memory (M has been reused as a lattice of
  // // Result)
  // while(LatDiff) {
  //   LatticeUnion *next = LatDiff->next;
  //   free(LatDiff);
  //   LatDiff = next;
  // }

  // // (TODO) also consider the intersection of lattices, where some points of
  // // lattice B->Lat could have no integer antecedent in B->P and should
  // // be kept in the result A - B:
  // // Add the holes of B (that can be included in A but not in B).


  // Domain_Free(ImInter);
  // LBLFree(Ainter);
  // LBLFree(Binter);

  // if(!Result) {
  //   #ifdef DIFFERENCE_DEBUG
  //     fprintf(stderr, "-- result = (NULL)\n");
  //   #endif
  //   return(NULL);
  // }

  // #ifdef DIFFERENCE_DEBUG
  //   fprintf(stderr, "-- temporary over-approximation of result = ");
  //   LBLPrint(stderr, P_VALUE_FMT, Result);
  // #endif
  // // intersect the result with A to get the exact LBL in case there was an
  // // over-approximation of A before.
  // Final_Result = LBLIntersection(Result, A);
  // LBLFree(Result);
  // return(Final_Result);

  // TODO: which one to use, Binter or B?
  // which one is simpler? B is larger... Binter is part of A
  Bcomp = LBLComplement(Binter);
  #ifdef DIFFERENCE_DEBUG
  fprintf(stderr, "Difference = intersection between Bcomp = ");
  LBLPrint(stderr, P_VALUE_FMT, Bcomp);
  fprintf(stderr, "[Difference = intersection between Bcomp] and A = ");
  LBLPrint(stderr, P_VALUE_FMT, A);
  #endif
  Result = LBLIntersection(Bcomp, A);

  LBLFree(Binter);
  LBLFree(Bcomp);

  return(Result);
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

  // if there is an eliminated column of zeros in the result:
  // LBLAlloc will ensure that eliminated points in the projection
  // are integer.

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
 * - build the LBL { z' | L z = G z', z \in A->P, z' free},
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

  // add the extra dimension on P (first dimensions!)
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
    value_substract(Con->p[i][Con->NbColumns-1],
              Con->p[i][Con->NbColumns-1], Z->Lat->p[i][Z->Lat->NbColumns-1]);
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
 * Change A->P such that all polyhedra in this domain have the same set of
 * equalities, that is, the equalities of the first polyhedron of this domain.
 * All the other ones are added to a new LBL, linked to LBL A (in A->next).
 */
static Matrix *sLBL_Homogenize_Equalities(LBL *A)
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
          errormsg1("sLBL_Canonical", "outofmem", "Out of Memory");
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
 * Returns True if A is modified
 */
static Bool sLBL_Simplify_Equalities(LBL *A, Matrix *Equalities)
{
  if (A->P->Dimension > 0 && A->P->NbEq != 0) {
    Matrix *H = NULL, *NewL;
    Matrix *eq_hermite = NULL, *ker = NULL;

    // simplify Lat with the equalities of domain P
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

    // previous method was:
    // if the bottom right value of H is not one, this means that
    // the transformation matrix is not integer but rational.
    // just make it integer to eliminate rational points.
    // (see ZImPre3 for an example where this is necessary)
    // value_set_si(H->p[H->NbRows-1][H->NbColumns-1], 1);

    // The result is just empty (because it is rational) when the bottom-right
    // value of H is not one.
    if(value_notone_p(H->p[H->NbRows-1][H->NbColumns-1])) {
      Domain_Free(A->P);
      A->P = NULL;  // will be fixed by caller
    }
    else {
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
    return(True);
  }
  else { // P contains no equalities
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "P has no equalities.\n");
    #endif
    return(False);
  }
} /* sLBL_Simplify_Equalities */


/*
 * compute the inside of polyhedron P that can be projected (along dim) to
 * get the dark shadow.
 * 
 * consider P a single polyhedron, even if P->next is set.
 * return NULL if dark == input P
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
 * compute the projection of domain P along dimension dim.
 * 
 * P is a polyhedral domain
 */
static Polyhedron *domain_project(Polyhedron *P, int dim)
{
  Matrix *T; // transformation: Id without the dim column
  Polyhedron *image;
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Entering exact shadow -dimension %d- of:", dim);
    Polyhedron_Print(stderr, P_VALUE_FMT, P);
  #endif

  T = Matrix_Alloc(P->Dimension, P->Dimension + 1);
  Vector_Set(T->p_Init, 0, T->p_Init_size);
  for(int i = 0; i < P->Dimension; i++) {
    if(i >= dim) {
      value_set_si(T->p[i][i+1], 1);
    }
    else {
      value_set_si(T->p[i][i], 1);
    }
  }

  image = DomainImage(P, T, MAXNOOFRAYS);
  #ifdef CANONICAL_DEBUG
  fprintf(stderr, "projected result P = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, image);
  #endif
  Matrix_Free(T);

  return (image);
} /* domain_project */


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

  rays = Matrix_Alloc(1, dim+2);
  Vector_Copy(val, rays->p[0], dim+1);
  value_set_si(rays->p[0][dim+1], 1);
  res = Rays2Polyhedron(rays, MAXNOOFRAYS);
  #ifdef HOLES_DEBUG
  fprintf(stderr, "generating polyhedron for vertex : ");
  Matrix_Print(stderr, P_VALUE_FMT, rays);
  #endif

  Matrix_Free(rays);
  return(res);
}


/*
 * Scan a single small pointy polyhedron R, and check for each of these
 * values if the scan hits an integer point.
 * 'scan' contains the whole domain top scan (A->P), R is the first-dimension
 * part of the (exact-dark) shadow.
 * Returns the union of polyhedra that verify this condition (not a hole)
 * uses the allocated array of values val of dimension scan->dimension+1
 * position = index position of current loop index (starting at 1 up to
 * scan->Dimension)
 * 
 * recursive,
 * - if R is scanned already (end of val init position),
 *   generate a for loop on all possible values on scan, and check if there
 *   is an integer solution at the end, early exit if found.
 * - if position <= R->Dimension, scan all possible values in R and recursive
 *   call on scan->next
 * 
 * val[hdim] must be one.
 */
Polyhedron *Scan_Rest(Polyhedron *scan, Polyhedron *R, Value *val,
  int position, int dim_R, Polyhedron *Result)
{
  Value LB, UB;

  value_set_si(val[position], 0); // ensure no previous value is assigned there
  #ifdef HOLES_DEBUG
  fprintf(stderr, "Enter Scan_Rest, position = %d\n", position);
  fprintf(stderr, "val = (");
  for(int i=0; i<=((scan)?(dim_R+scan->Dimension):position); i++) {
    if(i == position || i == 1)
      fprintf(stderr, "::");
    value_print(stderr, P_VALUE_FMT, val[i]);
  }
  fprintf(stderr, ")\n");
  #endif
  if(!R && !scan) {
      // end here, it is a hit!
      // generate polyhedron of the point = value val
      return(GenPoly(dim_R, val));
  }
  value_init(LB);
  value_init(UB);
  if(!R) {
    // no more R to scan, scan the 'scan' polytope, setting the values of
    // lower/upper bounds to val up to position.
    // loop and recursive call to scan->next, exit as soon as we found a hit!

    if(lower_upper_bounds(position, scan, val, &LB, &UB) != 0) {
      // should never happen (?)
      // errormsg1("Scan_Rest", "infinitepoly",
      //   "trying to scan an infinite (A->P) domain");

      // it's a hit if infinity!
      return(GenPoly(dim_R, val));
    }
    #ifdef HOLES_DEBUG
    fprintf(stderr, "position %d - looping from ", position);
    value_print(stderr, P_VALUE_FMT, LB);
    fprintf(stderr, "to ");
    value_print(stderr, P_VALUE_FMT "\n", UB);
    #endif
    // loop LB -> UB
    for(; value_le(LB, UB); value_increment(LB, LB)) {
      Polyhedron *res;
      // use LB in val[position] and scan next dimensions
      value_assign(val[position], LB);
      // recursive call (no accumulation of results here: Result=NULL)
      if((res = Scan_Rest(scan->next, NULL, val, position+1, dim_R, NULL))) {
        // it's a hit!
        value_clear(UB);
        value_clear(LB);
        return(AddPolyToDomain(res, Result));
        // if Result != NULL it's the initial call and res will be added to
        // Result, else this will just return res.
      }
    }
    // it's a hole!
    // goto the end, free memory and return Result unchanged.
  }
  else {
    // scaning R, recursive call to next dimension,
    // and accumulate the results in the Result union.
    value_set_si(val[R->Dimension+1], 1); // for lower_upper_bound
    if(lower_upper_bounds(position, R, val, &LB, &UB)) {
      // problem: infinity somewhere!
      errormsg1("Scan_Rest", "infinitepoly",
        "trying to scan an infinite (exact - dark) domain");
      return(Result);
    }
    #ifdef HOLES_DEBUG
    fprintf(stderr, "position %d - looping from ", position);
    value_print(stderr, P_VALUE_FMT, LB);
    fprintf(stderr, "to ");
    value_print(stderr, P_VALUE_FMT "\n", UB);
    #endif
    // loop: LB from LB to UB
    for(; value_le(LB, UB); value_increment(LB, LB)) {
      value_assign(val[position], LB);
      // accumulate all results
      Result = Scan_Rest(scan->next, R->next, val, position+1, dim_R, Result);
    }
  }
  value_clear(UB);
  value_clear(LB);
  return(Result);
} /* Scan_Rest */


/*
 * Compute the coordinate polyhedron containing the holes of the single LBL A.
 *
 * Usage: set *pExact to the exact shadow if the pointer is not NULL
 *
 * Algo:
 * - compute the domain (exact shadow - dark shadow)
 * - scan all its integer points and verify for each point:
 *      if there is an integer point in the origin intersection with the LBL
 *      add it to the polyhedral domain not_a_hole
 * - return (exact shadow - dark shadow) - not_a_hole
 */
static Polyhedron *sLBL_compute_holes(LBL *A, Polyhedron **pExact)
{
  int nbzeros;
  Polyhedron *exact = A->P, *dark = A->P; // initialize with P then project
  Polyhedron *rest; // exact shadow - dark shadow (polyhedral domain)
  Polyhedron *tmp, *AP, *U0, *not_a_hole = NULL, *holes;
  Vector *v;

  #ifdef HOLES_DEBUG
  fprintf(stderr, "Entering compute holes. A = ");
  LBLPrint(stderr, P_VALUE_FMT, A);
  #endif
  nbzeros = count_zeroCols(A->Lat);
  for(int z = 0; z < nbzeros; z++) {
    Polyhedron *d, *e; // shadow polyhedra after eliminating the column
    int col = A->Lat->NbColumns - 2 - z;

    d = domain_dark_shadow(dark, col);
    e = domain_project(exact, col);
    if(dark != A->P) {
      Domain_Free(dark);  // no longer need the previous calculated ones
      Domain_Free(exact); // (keep the original A->P first ones)
    }
    dark = d;
    exact = e;
  }

  // rest is the polyhedral domain (exact - dark) in origin-nbzero col space
  rest = DomainDifference(exact, dark, MAXNOOFRAYS);
  Domain_Free(dark);
  if(pExact) {
    *pExact = exact;
  }
  else {
    Domain_Free(exact);
  }
  // simplify obvious non integer cases
  rest = DomainConstraintSimplify(rest, MAXNOOFRAYS);
  if(emptyQ(rest)) {
    Domain_Free(rest);
    return(NULL);
  }

  // PREPARE SCAN:
  // make rest a disjoint domain, to avoid scanning a point multiple times
  tmp = Disjoint_Domain(rest, 0, MAXNOOFRAYS);
  Domain_Free(rest);
  rest = tmp;
  #ifdef HOLES_DEBUG
  fprintf(stderr, "disjoint (exact - dark) = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, rest);
  #endif
  // disjoint domain A->P (same)
  AP = Disjoint_Domain(A->P, 0, MAXNOOFRAYS);
  v = Vector_Alloc(A->P->Dimension + 2);
  // and universe (dim 0)
  U0 = Universe_Polyhedron(0);

  // need to:
  // - scan the points that can be holes (domain_disjoint + polyhedron scan
  //   of each + lower_upper_bound to scan)
  // - add them to the polyhedron list of hits if they are not holes
  while(AP) {
    // AP is the polyhedron that needs to be scanned
    // rest is the context
    Polyhedron *scanAP, *nextAP;
    nextAP = AP->next; // save and
    AP->next = NULL;   // unlink next

    // polyhedron scan does not work on a domain (need to nullify next)
    scanAP = Polyhedron_Scan(AP, U0, MAXNOOFRAYS);
    #ifdef HOLES_DEBUG
    fprintf(stderr, "Scanning:");
    Polyhedron_Print(stderr, P_VALUE_FMT, scanAP);
    #endif

    for(Polyhedron *R=rest; R; R = R->next) {
      Polyhedron *scanR, *nextR;

      nextR = R->next; // save and
      R->next = NULL;  // unlink next

      // polyhedron scan does not work on a domain (need to nullify next)
      scanR = Polyhedron_Scan(R, U0, MAXNOOFRAYS);
      R->next = nextR; // relink next

      Vector_Set(v->p, v->Size-1, 0);
      value_set_si(v->p[v->Size-1], 1);

      #ifdef HOLES_DEBUG
      fprintf(stderr, "------- Calling Scan_Rest -------");
      fprintf(stderr, "R = ");
      Polyhedron_Print(stderr, P_VALUE_FMT, scanR);
      fprintf(stderr, "scan = ");
      Polyhedron_Print(stderr, P_VALUE_FMT, scanAP);
      #endif
      // scan and update not_a_hole (add points that are not holes)
      not_a_hole = Scan_Rest(scanAP, scanR, v->p, 1, scanR->Dimension,
        not_a_hole);

      Domain_Free(scanR);
    }

    Domain_Free(scanAP);
    Polyhedron_Free(AP);
    AP = nextAP; // continue with next disjoint part
  }
  Vector_Free(v);
  Domain_Free(U0);
  #ifdef HOLES_DEBUG
  fprintf(stderr, "not holes in (exact-dark) = ");
  Polyhedron_Print(stderr, P_VALUE_FMT, not_a_hole);
  #endif

  // build final domain: (rest - not_a_hole)
  holes = DomainDifference(rest, not_a_hole, MAXNOOFRAYS);
  Domain_Free(rest);
  Domain_Free(not_a_hole);
  if(emptyQ(holes)) {
    #ifdef HOLES_DEBUG
    fprintf(stderr, "sLBL_compute_holes returning: <NULL>");
    #endif
    Domain_Free(holes);
    return(NULL); // no holes
  }

  return(holes);
} /* sLBL_compute_holes */


/*
 * Try to eliminate the zero columns of lattice A->Lat through
 * projection.
 * 
 * Eliminate only if exact shadow == dark shadow along each dimension.
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
      // TODO: for now, stay at the domain level.
      Polyhedron *dark = domain_dark_shadow(A->P, col);
      // compute exact projection of P and check if dark covers exact:
      Polyhedron *exact = domain_project(A->P, col);
      Polyhedron *diff;
      diff = DomainDifference(exact, dark, MAXNOOFRAYS);

      if(! emptyQ(diff)) {
        // try to remove obvious integer-empty solutions.
        // is this useful in some case?
        diff = DomainConstraintSimplify(diff, MAXNOOFRAYS);
        // TODO: could check if diff has no integer solution...
      }

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
          // one of the outer one can be eliminated now... (call the function
          // again at the end)
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
    // another one (that has been scaned before).
    sLBL_Simplify_Zero_Dimensions(A);
  }
}

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
static Bool sLBL_Remove_Empty(LBL *A)
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
  if(!A->P || A->P->Dimension > 0) {
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


// previous method, not considering existential integer variables...
// /*
//  * Remove the columns of zeros from A->Lat.
//  * In place. A->P is a domain.
//  * 
//  * This is equivalent to removing an existential variable: need to verify that
//  * there is an integer solution in the removed dimension
//  */
// static void sLBL_Lat_Remove_Zeros(LBL *A)
// {
//   // Could use DomainConstraintSimplify() to eliminate obvious empty case

//   Polyhedron *NewP;
//   int nbZeros = count_zeroCols(A->Lat);
//   if(nbZeros) {
//     // check which dimensions can be eliminated
//     Bool *elim = malloc(sizeof(Bool)*nbZeros); // dimensions to eliminate
//     int nbelim = 0; // number of dimensions to eliminate

//     for(int dim = 0; dim < nbZeros; dim++) {
//       int position = dim + A->Lat->NbColumns - nbZeros;
//       elim[dim] = True;
//       // if there is an equality with a coefficient different than +/- 1,
//       // the dimension cannot be eliminated
//       for(int j = 0; j < A->P->NbEq; j++) {
//         if(value_notone_p(A->P->Constraint[j][position]) &&
//             value_notmone_p(A->P->Constraint[j][position]) ) {
//           elim[dim] = False;
//           break;
//         }
//       }
//     }

//     // Now transform the domain
//     Matrix *Transformation;
//     Transformation = Matrix_Alloc(A->Lat->NbColumns - nbelim,
//                                   A->Lat->NbColumns);
//     // Id on top-left
//     for (int  i = 0; i < Transformation->NbRows; i++) {
//       for (int j = 0; j < Transformation->NbColumns; j++) {
//         if(i==j && i!=Transformation->NbRows-1) {
//           value_set_si(Transformation->p[i][j], 1);
//         }
//         else {
//           value_set_si(Transformation->p[i][j], 0);
//         }
//       }
//     }
//     // 1 on bottom-right
//     value_set_si(
//       Transformation->p[Transformation->NbRows-1][Transformation->NbColumns-1],
//       1);

//     NewP = DomainImage(A->P, Transformation, MAXNOOFRAYS);
//     Domain_Free(A->P);
//     A->P = NewP;
//     Matrix_Free(Transformation);

//     // Take the first columns of Lat
//     Matrix* NewL = Matrix_Alloc(A->Lat->NbRows, A->Lat->NbColumns-nbZeros);
//     for (int  i = 0; i < NewL->NbRows; i++) {
//       for (int j = 0; j < NewL->NbColumns; j++) {
//         if(j < NewL->NbColumns-1) {
//           value_assign(NewL->p[i][j], A->Lat->p[i][j]);
//         }
//         else {
//           value_assign(NewL->p[i][j], A->Lat->p[i][A->Lat->NbColumns-1]);
//         }
//       }
//     }
//     Matrix_Free(A->Lat);
//     A->Lat = NewL;

//     #ifdef CANONICAL_DEBUG
//       LBLPrint(stderr, P_VALUE_FMT, A);
//     #endif
//   }
// } /* sLBL_Lat_Remove_Zeros */


/*
 * Modify the single LBL 'A' (next is ignored) to be in canonical form:
 * A->Lat in HNF and no equalities in A->P.
 * Also tries to remove the columns of zeros from A->Lat if possible:
 * do the projection along those dimensions and eliminate only if
 * dark shadow = exact shadow
 *
 * USAGE: in place, modifies A itself
 * 
 * IMPORTANT: this function modifies the head single LBL of A to build a union,
 * it may add or remove the LBLs stored in A->next
 */
static void sLBL_Canonical(LBL* A)
{
  int simplified;
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

  // ************************
  // STEP 1: normalize A->Lat
  // ************************

  // Normalize the affine lattice A->Lat (and update A->P)
  sLBL_Lat_Normalize(A);

  // simplify non-integer constraints such that they intersect at least one
  // integer point (to avoid infinite empty integer polyhedra)
  A->P = DomainConstraintSimplify(A->P, MAXNOOFRAYS);

  // check emptyness
  if(emptyQ(A->P)) {
    if(sLBL_Remove_Empty(A)) {
      // head was empty and has been replaced with A->next.
      // need to canonicalize the (new) current LBL itself
      // (so the caller does not need to rescan it!)
      sLBL_Canonical(A);
    }
  }
  // now we have a non-empty A->P and A->Lat is canonical

  // ***********************************
  // STEP 2: remove equalities from A->P
  // ***********************************

  // homogenize the equalities of A->P: ensure that all polyhedra of the
  // domain verify the same set of equalities
  Matrix *Equalities = sLBL_Homogenize_Equalities(A);

  // We can remove the equalities from A->P
  simplified = sLBL_Simplify_Equalities(A, Equalities);
  Matrix_Free(Equalities);

  // If some equalities were eliminated start again from scratch!
  // (Lat and P changed and could be further simplified)
  if(!A->P || simplified) {
    sLBL_Canonical(A);
  }

  // ****************************************
  // STEP 3: eliminate zero columns of A->Lat
  // ****************************************

  // Remove the columns of zeros from A->Lat if possible
  // do the projection along those dimensions,
  // eliminate only if dark shadow \in exact shadow
  sLBL_Simplify_Zero_Dimensions(A);

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
void CanonicalLBL(LBL *A)
{
  // here, just transform every LBL of the list individually
  // careful, this may add a new LBL to the list A itself
  // (after tmp) but they will be scanned by this loop :)
  for(LBL *tmp = A; tmp; tmp = tmp->next) {
    sLBL_Canonical(tmp);
  }
  
  // check if a lattice is present twice in A, and if it is, union the other
  // polyhedral domain with this one and remove the second reference
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
static int count_zeroCols(Matrix* M)
{
  int nb = 0;
  for (int j = M->NbColumns-2; j >= 0; j--) {
    Bool isZero = True;
    for(int i = 0; i < M->NbRows; i++) {
      if(value_notzero_p(M->p[i][j])) {
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
 * Transform a single LBL into a list of Z-domains
 *
 * Remove zero columns from the lattice, build a union of Z-polyhedra
 * 
 */
static LBL *sLBL2ZDomain(LBL *A)
{
  int nbzeros;
  LBL *Result;

  if((nbzeros = count_zeroCols(A->Lat))) {
    // there are potential holes
    Matrix *newL;
    Polyhedron *holes, *not_holes, *exact;
    holes = sLBL_compute_holes(A, &exact);

    not_holes = DomainDifference(exact, holes, MAXNOOFRAYS);
    Polyhedron_Free(holes);
    Polyhedron_Free(exact);

    // build result LBL
    newL = RemoveNColumns(A->Lat, A->Lat->NbColumns-1-nbzeros, nbzeros);
    Result = LBLAlloc(newL, not_holes);
    Matrix_Free(newL);
    Polyhedron_Free(not_holes);
  }
  else {
    Result = sLBL_Copy(A);
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
  LBL *Result = NULL;
  for(LBL *Z = A; Z; Z = Z->next)
  {
    LBL *tmp;
    tmp = sLBL2ZDomain(Z);
    Result = LBL_concatenate(tmp, Result);
  }
  CanonicalLBL(Result);
  return(Result);
}

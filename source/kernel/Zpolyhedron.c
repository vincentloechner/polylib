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
  #define CANONICAL_DEBUG 1
  #define INTERSECTION_DEBUG 1
  #define DIFFERENCE_DEBUG 1
  #define COMP_DEBUG 1
#endif
#define COMP_DEBUG 1

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
 *      build the LBL { AL z | AL z = BL z', z \in AP, z' \in BP },
 *      and remove z' by normalizing the result
 */
static LBL *sLBL_Intersection(LBL *A, LBL *B) {

  LBL *Result = NULL;
  Matrix *LInter;
  Polyhedron *PInter, *ImageA, *ImageB, *PreImage;

  #ifdef INTERSECTION_DEBUG
    fprintf(stderr, "entering sLBL_Intersection\nA = ");
    sLBL_Print(stderr, P_VALUE_FMT, A);
    fprintf(stderr, "B = ");
    sLBL_Print(stderr, P_VALUE_FMT, B);
  #endif
  LInter = LatticeIntersection(A->Lat, B->Lat);
  if (isEmptyLattice(LInter)) {
    Matrix_Free(LInter);
    return (NULL);
  }

  if(count_zeroCols(A->Lat) == 0 && count_zeroCols(B->Lat) == 0 &&
    count_zeroCols(LInter) == 0)
  {
    // This works only IF there are no columns of zeros in the LBLs:
    // they are Z-polyhedra

    ImageA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
    ImageB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
    PInter = DomainIntersection(ImageA, ImageB, MAXNOOFRAYS);
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
    #ifdef INTERSECTION_DEBUG
      fprintf(stderr, "Z-polyhedra, simplified intersection = ");
      sLBL_Print(stderr, P_VALUE_FMT, Result);
    #endif

    return (Result);
  }
  else {
    // build the LBL { AL z | AL z = BL z', z \in AP, z' \in BP }

    int extra_max_rows = 0, extra_B_row = A->Lat->NbRows - 1;
    Matrix *newL, *extra;
    Polyhedron *newP = NULL, *AP_aligned;
    Matrix_Free(LInter);

    //          z     z'   cst
    // newL =   AL |  0   | Al |
    //           0 |  0   | 1  |   <- this row is in AL already
    newL = Matrix_Alloc(A->Lat->NbRows, A->Lat->NbColumns + B->Lat->NbColumns - 1);
    for(int i = 0; i < newL->NbRows; i++) {
      for(int j = 0; j < newL->NbColumns; j++) {
        if(j < A->Lat->NbColumns - 1) { // left AL
          value_assign(newL->p[i][j], A->Lat->p[i][j]);
        }
        else if(j == newL->NbColumns - 1) {
          value_assign(newL->p[i][j], A->Lat->p[i][A->Lat->NbColumns - 1]);
        }
        else {
          value_set_si(newL->p[i][j], 0);
        }
      }
    }

    // newP:

    // scan the polyhedra of domains A->P and B->P, and build their constraints
    // intersections.
    // Build the constraints:
    //  0/1    z      z'      cst
    //  ap0    AP      0       Ap     # from A->P    -> AP_aligned
    //   0     AL     -BL      Al-Bl  # equalities   \ extra
    //  bp0    0      BP       Bp     # from B->P    /
    
    // start with a domain of the right dimension (expand dimension of A)
    AP_aligned = align_context(A->P, A->P->Dimension + B->P->Dimension,
      MAXNOOFRAYS);
    // extra will be the matrix containing the extra constraints (including
    // equalities: AL z = BL z', initialized once)
    for(Polyhedron *BP = B->P; BP; BP = BP->next) {
      if(A->Lat->NbRows + B->P->NbConstraints > extra_max_rows) {
        extra_max_rows = A->Lat->NbRows + B->P->NbConstraints;
      }
    }
    extra = Matrix_Alloc(extra_max_rows, AP_aligned->Dimension + 2);
    // initialize |0 AL  -BL  constant| in extra
    for(int row = 0; row < extra_B_row; row++) {
      value_set_si(extra->p[row][0], 0); // equality
      Vector_Copy  (&A->Lat->p[row][0], &extra->p[row][1],
        A->Lat->NbColumns - 1);
      Vector_Oppose(&B->Lat->p[row][0], &extra->p[row][A->Lat->NbColumns],
        B->Lat->NbColumns - 1);
      value_substract(extra->p[row][extra->NbColumns - 1],
        A->Lat->p[row][A->Lat->NbColumns - 1],
        B->Lat->p[row][B->Lat->NbColumns - 1]);
    }
    // scan the intersections and build union newP
    for(Polyhedron *AP = AP_aligned; AP; AP = AP->next) {
      for(Polyhedron *BP = B->P; BP; BP = BP->next) {
        Polyhedron *P;
        // complement extra with the Constraints of BP
        for(int con = 0; con < BP->NbConstraints; con++) {
          // 0/1 bp0
          value_assign(extra->p[extra_B_row + con][0], BP->Constraint[con][0]);
          // constraint BP + constant Bp
          Vector_Copy(&BP->Constraint[con][1],
            &extra->p[extra_B_row + con][A->Lat->NbColumns],
            BP->Dimension + 1);
        }
        extra->NbRows = extra_B_row + B->P->NbConstraints;
        #ifdef INTERSECTION_DEBUG
          fprintf(stderr, "extra = ");
          Matrix_Print(stderr, P_VALUE_FMT, extra);
        #endif
  
        P = AddConstraints(&extra->p[0][0], extra->NbRows,
          AP_aligned, MAXNOOFRAYS);
        newP = AddPolyToDomain(P, newP); // consumes P and newP
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


// TODO:
LBL *LBLComplement(LBL *A)
{
  // Let L = A->Lat, P = A->P.
  // complement(A) = Universe() - A = union of:
  //   - LBL (Z^d, not hull(A))
  //   - LBL ((Z^d - L), universe polyhedron)
  //   - holes of A -> how to compute these????
  //       = L z' such that there exist no z in A->P such that L z' = L z
  //       = {L . z' | L z' != L z, z \in P, AND z INTEGER }
  //       = {L . z' | L z' > L z, z \in P } union
  //         {L . z' | L z' < L z, z \in P }
  //             + integer z -> need something similar to exact shadow
  //     if L invertible -> empty.

  // we could also just build:
  //   - { Id z' | z' != L z, z in P }
  // does this work??
  // (?) no because missing condition: z INTEGER
  // but yes because LBLAlloc will ensure no rational points are eliminated,
  // and keep inequalities and zero columns in Lat if necessary.

  // Build this:
  //    { z' | L z > z', z in P } union
  //    { z' | z' < L z, z in P }

  // more exactly:
  //   (Id  0)  (z')  |  (-Id  L) . (z') !=  0,   (z') in P_extended
  //            (z )  |             (z )          (z )
  // z' has L->NbRows dimension

  LBL *Result;
  Polyhedron *PResult = NULL;
  Matrix *ZId = NULL;
  Matrix *extra; // extra constraints (allocated once, reused)
  int nb_col = A->Lat->NbColumns + A->Lat->NbRows - 1; // new lattice nb col
  int nb_rows = A->Lat->NbRows; // nb lattice nb rows
  #ifdef COMP_DEBUG
  fprintf(stderr, "Entering LBLComplement. A = ");
  LBLPrint(stderr, P_VALUE_FMT, A);
  #endif

  ZId = Matrix_Alloc(nb_rows, nb_col);
  for(int i = 0; i < nb_rows; i++) {
    for(int j = 0; j < nb_col; j++) {
      if((i == j && i != nb_rows-1) || (j == nb_col-1 && i == nb_rows-1))
        value_set_si(ZId->p[i][j], 1);
      else
        value_set_si(ZId->p[i][j], 0);
    }
  }
  #ifdef COMP_DEBUG
  fprintf(stderr, "ZId = ");
  Matrix_Print(stderr, P_VALUE_FMT, ZId);
  #endif
  extra = Matrix_Alloc(nb_rows, nb_col + 1); // PResult constraints dim.

  for(Polyhedron *AP = A->P; AP; AP=AP->next) {
    Polyhedron *Pextended;

    Pextended = align_context(AP, AP->Dimension + nb_rows - 1, MAXNOOFRAYS);
    #ifdef COMP_DEBUG
    fprintf(stderr, "Pextended = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, Pextended);
    #endif

    // build constraints z' > L z, pos is the number of equalities (first):
    //     z'(above pos) = L z (above pos)
    //     z'_pos <= L z_pos (-1 for strict)
    for(int pos = 0; pos < A->Lat->NbRows - 1; pos++) {
      Polyhedron *P;
      // constraints above pos have been initialized correctly already.

      // POSITIVE: (-Id  L) (z' z)^T > 0 (>= 1. --> -1 on line)
      // inequality
      value_set_si(extra->p[pos][0], 1);
      // L (linear part)
      Vector_Copy(&A->Lat->p[pos][0], &extra->p[pos][nb_rows],
        A->Lat->NbColumns-1);
      // -Id
      Vector_Set(&extra->p[pos][1], 0, nb_rows-1);
      value_set_si(extra->p[pos][pos+1], -1); // (-1) * z'_pos
      // constant - 1
      value_sub_int(extra->p[pos][extra->NbColumns-1],
        A->Lat->p[pos][A->Lat->NbColumns-1], 1);
      #ifdef COMP_DEBUG
      fprintf(stderr, "extra = ");
      Matrix_Print(stderr, P_VALUE_FMT, extra);
      #endif
      P = AddConstraints(extra->p[0], pos+1, Pextended, MAXNOOFRAYS);
      #ifdef COMP_DEBUG
      fprintf(stderr, "Adding P = ");
      Polyhedron_Print(stderr, P_VALUE_FMT, P);
      #endif
      PResult = AddPolyToDomain(P, PResult);

      // NEGATIVE: (Id  -L) (z' z)^T > 0 (>= 1. --> -1 on line)
      // inequality
      value_set_si(extra->p[pos][0], 1);
      // -L (linear part)
      Vector_Oppose(&A->Lat->p[pos][0], &extra->p[pos][nb_rows],
        A->Lat->NbColumns-1);
      // Id
      value_set_si(extra->p[pos][pos+1], 1); // (1) * z'_pos
      // -constant - 1
      value_oppose(extra->p[pos][extra->NbColumns-1],
        A->Lat->p[pos][A->Lat->NbColumns-1]);
      value_sub_int(extra->p[pos][extra->NbColumns-1],
        extra->p[pos][extra->NbColumns-1], 1);
      #ifdef COMP_DEBUG
      fprintf(stderr, "extra = ");
      Matrix_Print(stderr, P_VALUE_FMT, extra);
      #endif
      P = AddConstraints(extra->p[0], pos+1, Pextended, MAXNOOFRAYS);
      #ifdef COMP_DEBUG
      fprintf(stderr, "Adding P = ");
      Polyhedron_Print(stderr, P_VALUE_FMT, P);
      #endif
      PResult = AddPolyToDomain(P, PResult);

      // TRANSFORM THE INEQUALITY INTO AN EQUALITY (set constant right)
      value_set_si(extra->p[pos][0], 0);
      value_oppose(extra->p[pos][extra->NbColumns-1],
        A->Lat->p[pos][A->Lat->NbColumns-1]);
    }
  }
  #ifdef COMP_DEBUG
  fprintf(stderr, "PResult below. Trying to LBLAlloc next\n");
  Polyhedron_Print(stderr, P_VALUE_FMT, PResult);
  #endif
  Result = LBLAlloc(ZId, PResult);
  return (Result);
}

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
      fprintf(stderr,
      "Binter=(A inter B) is empty, so B does not intersect A, we return A\n");
    #endif
    LBLFree(Binter);
    return(LBLCopy(A));
  }

  // Separate the computation in 3 phases:
  // 0. compute the difference of the image polyhedra P_A \ P_B (=ImDiff) and
  //    add it to the solution LBL (with lattice L_A).
  //    This can be an over-approximation of A if A->Lat has zero columns
  //    (but not of B)
  // 1. compute the rest where the intersection of P_A and P_B have same
  //    dimensions (required for lattice difference)
  // 2. intersect the result with A to get rid of the over-approximations

  // Similar to computing (complement(B))  \inter  A

  // [STEP 0 (includes Gautam's Step 2)]
  imA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
  imB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
  ImDiff = DomainDifference(imA, imB, MAXNOOFRAYS);
  #ifdef DIFFERENCE_DEBUG
    fprintf(stderr, "ImDiff (hull of A that does not cover B) = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, ImDiff);
  #endif

  // Add (A->Lat, A->P - hull(B)) to the result:
  if (!emptyQ(ImDiff)) {
    Polyhedron *RedPolyDiff;
    RedPolyDiff = DomainPreimage(ImDiff, A->Lat, MAXNOOFRAYS);
    // NOTICE: this can be an over-approximation of A
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
    fprintf(stderr,
      "-- [STEP1] now we compute the intersection on same lattice dimensions\n");
    fprintf(stderr, "Ainter = ");
    LBLPrint(stderr, P_VALUE_FMT, Ainter);
    fprintf(stderr, "and Binter = ");
    LBLPrint(stderr, P_VALUE_FMT, Binter);
  #endif

  // LatDiff (union of lattices) is the difference : (A->Lat) - (B->Lat) of
  // same dimensions
  LatDiff = LatticeDifference(Ainter->Lat, Binter->Lat); 
  #ifdef DIFFERENCE_DEBUG
    if(!LatDiff)
      fprintf(stderr, "Empty Lattice difference\n");
  #endif

  // TODO: consider the intersection of lattices, where some points of lattice
  // B->Lat could have no integer antecedent in B->P (???) and should be kept
  // in the result A - B


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
  // free LatticeUnion remaining memory (M has been reused as a lattice of
  // Result)
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
  // intersect the result with A to get the exact LBL in case there was an
  // over-approximation before.
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

    // TODO: verify that.
    // The result is just empty when the bottom-right value of H is not one!
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
      Matrix_Print(stderr, P_VALUE_FMT, NewL);
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


// TODO: what about higher dimensions than dim? Take into account or not ?
// (probably yes)

/*
 * compute the inside of polyhedron P that can be projected (along dim) to
 * get the dark shadow.
 * 
 * consider P a single polyhedron, even if P->next is set.
 * returns NULL if dark shadow = exact shadow.
 */
static Polyhedron *polyhedron_dark_source(Polyhedron *P, int dim)
{
  // check if that dimension is constrained in P
  // - if it is not constrained, just add an equality {i_dim = 0}.
  // - if it is positive constrained,
  //   scan all positive constraints on i_dim of the form:
  //     {... + alpha . i_dim + ... + c >= 0}, with alpha > 0
  //   and add this constraint to P:
  //     {... + alpha . i_dim + ... + c + alpha-1 >= 0}
  // - if negative constrained do the opposite
  //   (take the smallest between them)

  int pos_constrained = 0, // number of alpha's > 0
      neg_constrained = 0, // number of alpha's < 0
      eq_constrained = 0;  // equality (both pos and neg), keep
  Polyhedron *result;

  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Entering dark_source. dimension: %d\n", dim);
    fprintf(stderr, "Polyhedron: ");
    Polyhedron_Print(stderr, P_VALUE_FMT, P);
  #endif
  // count constraints on this variable (dim)
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
    // one side is integer (or open to infinite), can be safely ignored :)
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
}


// TODO: check if all above dimensions should also be projected, or not ?!

/*
 * compute the projection of domain P along dimension dim.
 * 
 * P is a domain
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
  Matrix_Free(T);
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "projected result: ");
    Polyhedron_Print(stderr, P_VALUE_FMT, image);
  #endif

  return (image);
}

/*
 * Check that the zero columns of lattice A->Lat have one single integer
 * point in A->P as predecessor.
 * 
 * case 1. If A->Lat has empty columns on the right, and that they that do not
 * appear in the constraints of A->P, then A->P is under-constrained. We just
 * add an equality to P and it will get simplified automatically :)
 * case 2. If there are more than one integer point in A->P that map to the
 * same point by the A->Lat function, the redundant ones need to be eliminated
 * 
 */
static void sLBL_Simplify_Zero_Dimensions(LBL *A)
{
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
      // check if that dimension (corresponding to variable i_col) is
      // constrained in all polyhedra of domain A->P.
      // - if it is not constrained, just add an equality {i_col = 0}.
      // - if it is positive constrained, set the thickness to one:
      //   scan all positive constraints on i_col of the form:
      //     {... + alpha . i_col + ... >= 0}, with alpha > 0
      //   and add this constraint to P:
      //     {... + alpha . i_col + ... + alpha-1} <= 0 (oppose!)
      // - TODO: if it is not positive constrained but only negative
      //   constrained, what to do?
      // - what if there is a conflict between different members of the union?

      // easier method:
      // build the domain transformed to start at i_col = 0 (each constraint
      // of each polyhedron)
      // then add constraint i_col <= 1.
      // does not work, non unimodular transformation.

      // solution 3: compute the dark shadow and the exact shadow.
      // if dark shadow projection is in exact shadow: can project
      // else: keep

      // TODO: what if some polyhedra of the union can be projected and others
      // cannot? should we separate them or just stay at the domain level?
      Polyhedron *dark = NULL;
      for(Polyhedron *pp = A->P; pp; pp = pp->next) {
        Polyhedron *pp_inside, *pp_shadow;

        pp_inside = polyhedron_dark_source(pp, col);
        if(pp_inside) {
          pp_shadow = domain_project(pp_inside, col);
        }
        else {
          Polyhedron *ppnext = pp->next;
          pp->next = NULL;
          pp_shadow = domain_project(pp, col);
          pp->next = ppnext;
        }
        dark = AddPolyToDomain(pp_shadow, dark);
        if(pp_inside) {
          Polyhedron_Free(pp_inside);
        }
      }
      // compute exact projection of P and check if dark covers exact:
      Polyhedron *exact = domain_project(A->P, col);
      Polyhedron *diff;
      diff = DomainDifference(exact, dark, MAXNOOFRAYS);

      // is this useful? removes obvious integer empty solutions.
      diff = DomainConstraintSimplify(diff, MAXNOOFRAYS);
      // TODO: should check if diff has no integer solution... ???

      if(emptyQ(diff)) {
        // if exact - dark = 0, can project :)
        Matrix *newL;
        Domain_Free(A->P);
        A->P = exact;

// TODO: should create a new LBL here
        // remove column from A->Lat
        newL = RemoveColumn(A->Lat, col);
        Matrix_Free(A->Lat);
        A->Lat = newL;
      }
      else {
        Domain_Free(exact);
      }
      Domain_Free(diff);
      Domain_Free(dark);
    }
    else {
      // non empty column, everything on the left is also non empty, exit.
      break;
    }
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


/*
 * Remove the columns of zeros from A->Lat.
 * In place. A->P is a domain.
 * 
 * This is equivalent to removing an existential variable: need to verify that
 * there is an integer solution in the removed dimension
 */
static void sLBL_Lat_Remove_Zeros(LBL *A)
{
  // TODO: need to consider integer existential variable elimination and dark
  // shadow!

  // Could use DomainConstraintSimplify() to eliminate obvious empty case


  Polyhedron *NewP;
  int nbZeros = count_zeroCols(A->Lat);
  if(nbZeros) {
    // check which dimensions can be eliminated
    Bool *elim = malloc(sizeof(Bool)*nbZeros); // dimensions to eliminate
    int nbelim = 0; // number of dimensions to eliminate

    for(int dim = 0; dim < nbZeros; dim++) {
      int position = dim + A->Lat->NbColumns - nbZeros;
      elim[dim] = True;
      // if there is an equality with a coefficient different than +/- 1,
      // the dimension cannot be eliminated
      for(int j = 0; j < A->P->NbEq; j++) {
        if(value_notone_p(A->P->Constraint[j][position]) &&
            value_notmone_p(A->P->Constraint[j][position]) ) {
          elim[dim] = False;
          break;
        }
      }
    }


    // Now transform the domain
    Matrix *Transformation;
    Transformation = Matrix_Alloc(A->Lat->NbColumns - nbelim,
                                  A->Lat->NbColumns);
    // Id on top-left
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
    // 1 on bottom-right
    value_set_si(
      Transformation->p[Transformation->NbRows-1][Transformation->NbColumns-1],
      1);

    NewP = DomainImage(A->P, Transformation, MAXNOOFRAYS);
    Domain_Free(A->P);
    A->P = NewP;
    Matrix_Free(Transformation);

    // Take the first columns of Lat
    Matrix* NewL = Matrix_Alloc(A->Lat->NbRows, A->Lat->NbColumns-nbZeros);
    for (int  i = 0; i < NewL->NbRows; i++) {
      for (int j = 0; j < NewL->NbColumns; j++) {
        if(j < NewL->NbColumns-1) {
          value_assign(NewL->p[i][j], A->Lat->p[i][j]);
        }
        else {
          value_assign(NewL->p[i][j], A->Lat->p[i][A->Lat->NbColumns-1]);
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

  // Normalize the affine lattice A->Lat (can also update A->P)
  sLBL_Lat_Normalize(A);

  // simplify non-integer constraints such that they intersect at least one
  // integer point (to avoid infinite empty integer polyhedra for example)
  A->P = DomainConstraintSimplify(A->P, MAXNOOFRAYS);

  // check emptyness
  if(emptyQ(A->P)) {
    if(sLBL_Remove_Empty(A)) {
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

  // We can remove the equalities from A->P
  simplified = sLBL_Simplify_Equalities(A, Equalities);
  Matrix_Free(Equalities);

  // If some equalities were eliminated start again from scratch!
  // (Lat and P changed and could be further simplified)
  if(!A->P || simplified) {
    sLBL_Canonical(A);
  }

  // Remove the columns of zeros from A->Lat if possible
  // do the projection along those dimensions,
  // eliminate only if dark shadow = exact shadow
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

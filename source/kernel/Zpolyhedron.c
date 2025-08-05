#include <polylib/polylib.h>
#include <stdlib.h>

// DEFINITIONS USED IN ALL COMMENTS BELOW:
//
// - a ZPolyhedron is a single lattice associated to a polyhedral domain
//   (a polyhedral domain can be a union of polyhedra)
//
// - a ZDomain is a chained list of ZPolyhedra (with possibly multiple lattices;
//   using the ->next structure element).


// debug this file:
// #define DEBUG
#ifdef DEBUG
  #define LAT_TEST 1
  #define CANONICAL_DEBUG 1
  #define INTERSECTION_DEBUG 1
  #define DIFF_DEBUG 1
#endif

static ZPolyhedron *ZPolyhedronIntersection(ZPolyhedron *, ZPolyhedron *);
static ZPolyhedron *ZPolyhedron_Copy(ZPolyhedron *A);
static void ZPolyhedron_Free(ZPolyhedron *Zpol);
static ZPolyhedron *ZPolyhedronDifference(ZPolyhedron *, ZPolyhedron *);
static ZPolyhedron *ZPolyhedronImage(ZPolyhedron *, Matrix *);
static ZPolyhedron *ZPolyhedronPreimage(ZPolyhedron *, Matrix *);
static void ZPolyhedronPrint(FILE *fp, const char *format, ZPolyhedron *A);
static void Canonical_ZPolyhedron_Gautam(ZPolyhedron* A);
static ZPolyhedron *FindLatticePred(Lattice *L, ZPolyhedron *A);
static ZPolyhedron *ZD_ZP_Difference(ZPolyhedron* A, ZPolyhedron* B);
// static Bool ZPolyhedronIncludes(ZPolyhedron *A, ZPolyhedron *B);
static int count_zeroCols (Matrix* M);

typedef struct forsimplify {
  Polyhedron *Pol;
  LatticeUnion *LatUni;
  struct forsimplify *next;
} ForSimplify;

/*
 * Returns True if 'Zpol' is empty, otherwise returns False
 * ZPol can be a non-simplified list of empty ZDomain
 */
Bool isEmptyZDomain(ZPolyhedron *Zpol) {

  if (Zpol == NULL)
    return True;
  if (emptyQ(Zpol->P)) {
    // check the emptiness of next
    return(isEmptyZDomain(Zpol->next));
  }
  return False;
} /* isEmptyZDomain */

/*
 * Given a Lattice 'Lat' and a Domain 'Domain', allocate space, and return the
 * Z-polyhedron corresponding to the image of the integer points of 'Poly'
 * by the lattice 'Lat', in canonical form (HNF, no equalities)
 */
ZPolyhedron *ZPolyhedronAlloc(Lattice *Lat, Polyhedron *Domain) {

  ZPolyhedron *A;

  POL_ENSURE_FACETS(Domain);
  POL_ENSURE_VERTICES(Domain);

  if (Lat->NbColumns != Domain->Dimension + 1) {
    errormsg1("ZPolyhedronAlloc", "dimincomp",
      "the Lattice and the Polyhedron are not compatible to form a ZPolyhedron");
    return NULL;
  }

  A = malloc(sizeof(ZPolyhedron));
  if (!A) {
    errormsg1("ZPolyhedronAlloc", "outofmem", "Out of Memory");
    return NULL;
  }
  A->next = NULL;
  A->P = Domain_Copy(Domain);
  A->Lat = Matrix_Copy(Lat);

  Canonical_ZDomain(A);
  return A;
} /* ZPolyhedronAlloc */

/*
 * Free the memory used by the Z-domain 'Head'
 */
void ZDomain_Free(ZPolyhedron *Head) {

  if (Head == NULL)
    return;
  ZDomain_Free(Head->next);
  ZPolyhedron_Free(Head);
} /* ZDomain_Free */

/*
 * Free the memory used by the Z-polyhedron 'Zpol'
 */
static void ZPolyhedron_Free(ZPolyhedron *Zpol) {

  if (Zpol == NULL)
    return;
  if(Zpol->Lat)
    Matrix_Free(Zpol->Lat);
  if(Zpol->P)
    Domain_Free(Zpol->P);
  free(Zpol);
  return;
} /* ZPolyhderon_Free */

/*
 * Return a copy of the Z-domain 'Head'
 */
ZPolyhedron *ZDomain_Copy(ZPolyhedron *Head) {

  ZPolyhedron *Zpol;
  Zpol = ZPolyhedron_Copy(Head);

  if (Head->next != NULL)
    Zpol->next = ZDomain_Copy(Head->next);
  return Zpol;
} /* ZDomain_Copy */

/*
 * Return a copy of the Z-polyhedron 'A'
 */
static ZPolyhedron *ZPolyhedron_Copy(ZPolyhedron *A) {

  ZPolyhedron *Zpol;

  Zpol = ZPolyhedronAlloc(A->Lat, A->P);
  return Zpol;
} /* ZPolyhedron_Copy */

/*
 * Add the Z-Domain 'Zpol' as first list element to the Z-domain 'Result'
 * and return a pointer to the new Z-domain.
 * *Consumes the memory* of Result and of Zpol (no need to free) to build
 * the result.
 */
static ZPolyhedron *ZDconcatenate(ZPolyhedron *Zpol, ZPolyhedron *Result) {

  if (isEmptyZDomain(Zpol)) {
    ZDomain_Free(Zpol);
    return Result;
  }
  if (isEmptyZDomain(Result)) {
    ZDomain_Free(Result);
    return Zpol;
  }

  // go to end of Zpol and concatenate Result there
  ZPolyhedron *tmp = Zpol;
  while(tmp->next)
    tmp = tmp->next;
  tmp->next = Result;
  
  return Zpol;
} /* ZDconcatenate */


/*
 * Return the empty Z-polyhedron of dimension 'dimension'
 */
ZPolyhedron *EmptyZPolyhedron(int dimension) {

  ZPolyhedron *Zpol;
  Lattice *E;
  Polyhedron *P;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered EMPTYZPOLYHEDRON\n");
  fclose(fp);
#endif

  E = Matrix_Alloc(dimension+1,1);
  for(int j=0 ; j<dimension; j++) {
      value_set_si(E->p[0][j], 0);
  }
  value_set_si(E->p[0][dimension], 1);

  P = Empty_Polyhedron(0);

  Zpol = ZPolyhedronAlloc(E, P);
  Matrix_Free((Matrix *)E);
  Domain_Free(P);
  return Zpol;
} /* EmptyZPolyhedron */

/*
 * Given Z-domains A and B, return True if A is included in B,
 * otherwise return False.
 */
Bool ZDomainIncludes(ZPolyhedron *A, ZPolyhedron *B) {

  Bool ret = False;

  #ifdef DIFF_DEBUG
    fprintf(stderr, "entering Includes. Checking that A =");
    ZPolyhedronPrint(stderr, P_VALUE_FMT, A);
    fprintf(stderr, "is included in B = ");
    ZPolyhedronPrint(stderr, P_VALUE_FMT, B);
  #endif

  ZPolyhedron *diff;
  diff = ZDomainDifference(A, B);
  if(isEmptyZDomain(diff)) {
    ret = True;
  }
  #ifdef DIFF_DEBUG
    fprintf(stderr, "diff = ");
    ZDomainPrint(stderr, P_VALUE_FMT, diff);
  #endif
  ZDomain_Free(diff);
  return ret;
} /* ZDomainIncludes */

/*
 * Given Z-polyhedra 'A' and 'B', return True if 'A' is directly present in 'B',
 * otherwise return False
 */
// static Bool old_ZPolyhedronIncludes(ZPolyhedron *A, ZPolyhedron *B) {

//   Polyhedron *Diff = NULL;
//   Bool retval = False;
// #ifdef DOMDEBUG
//   FILE *fp;
//   fp = fopen("_debug", "a");
//   fprintf(fp, "\nEntered ZPOLYHEDRONINCLUDES\n");
//   fclose(fp);
// #endif
  
//   // fast first check: the lattices intersect
//   if (LatticeIncludes(A->Lat, B->Lat) == True) {
//     Polyhedron *ImageA, *ImageB;

//     // fast second check: the rational images intersect
//     ImageA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
//     ImageB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);

//     Diff = DomainDifference(ImageA, ImageB, MAXNOOFRAYS);
//     if (emptyQ(Diff))
//       retval = True;
//     else
//     {
//       // TODO: need to compute the exact difference here in case we are missing something

//     }
//     Domain_Free(ImageA);
//     Domain_Free(ImageB);
//     Domain_Free(Diff);
//   }
//   return retval;
// } /* ZPolyhedronIncludes */

/*
 * Print the contents of a Z-domain 'A'
 */
void ZDomainPrint(FILE *fp, const char *format, ZPolyhedron *A) {
#ifdef DOMDEBUG
  FILE *fp1;
  fp1 = fopen("_debug", "a");
  fprintf(fp1, "\nEntered ZDOMAINPRINT\n");
  fclose(fp1);
#endif

  for( ; A; A=A->next) {
    ZPolyhedronPrint(fp, format, A);
    if(A->next)
      fprintf(fp, "\nUNION ");
  }
} /* ZDomainPrint */

/*
 * Print the contents of a ZPolyhderon 'A'
 */
static void ZPolyhedronPrint(FILE *fp, const char *format, ZPolyhedron *A) {
  if (A == NULL)
    return;
  fprintf(fp, "ZPOLYHEDRON: Dimension %d \n", A->Lat->NbRows - 1);
  fprintf(fp, "\nLATTICE: \n");
  Matrix_Print(fp, format, (Matrix *)A->Lat);
  Polyhedron_Print(fp, format, A->P);
  return;
} /* ZPolyhedronPrint */

/*
 * Return the Z-domain union of the Z-domain 'A' and 'B'. The dimensions of the
 * Z-domains 'A' and 'B' must be equal. All the Z-polyhedra of the resulting
 * union are expected to be in Canonical forms.
 */
ZPolyhedron *ZDomainUnion(ZPolyhedron *A, ZPolyhedron *B) {

  ZPolyhedron *Result = NULL, *tempA , * tempB;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZDOMAINUNION\n");
  fclose(fp);
#endif

  // copy A and B, concatenate, Canonicalize, and return :)

  tempA = ZDomain_Copy(A); 
  tempB = ZDomain_Copy(B);

  Result = ZDconcatenate(tempA, tempB);
  Canonical_ZDomain(Result);
  return Result;
} /* ZDomainUnion */

/*
 * Return the Z-domain intersection of the Z-domains 'A' and 'B'.The dimensions
 * of domains 'A' and 'B' must be equal.
 */
ZPolyhedron *ZDomainIntersection(ZPolyhedron *A, ZPolyhedron *B) {

  ZPolyhedron *Result = NULL, *tempA = NULL, *tempB = NULL;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZDOMAININTERSECTION\n");
  fclose(fp);
#endif

  for (tempA = A; tempA != NULL; tempA = tempA->next)
    for (tempB = B; tempB != NULL; tempB = tempB->next) {
      ZPolyhedron *Zpol;
      Zpol = ZPolyhedronIntersection(tempA, tempB);
      Result = ZDconcatenate(Zpol, Result);
    }
  if (Result == NULL)
    return EmptyZPolyhedron(A->Lat->NbRows - 1);

  Canonical_ZDomain(Result);
  return (Result);
} /* ZDomainIntersection */

/*
 * Return the Z-domain difference of the domains (A - B) in canonical form.
 * The dimensions of the Z-domains A and B must be equal. Note that the
 * difference of two Z-polyhedra is a Union of Z-polyhedra.
 * Algorithm: 
 * 
 */
ZPolyhedron *ZDomainDifference(ZPolyhedron *A, ZPolyhedron *B) { 

  ZPolyhedron *tempB;
  ZPolyhedron *res;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZDOMAINDIFFERENCE\n");
  fclose(fp);
#endif

  if (A->Lat->NbRows != B->Lat->NbRows) {
    errormsg1("ZDomainDifference", "dimincomp", "incompatible dimensions between domains");
    return (NULL);
  }
  
  res = ZPolyhedron_Copy(A);
  // remove all Z-polyhedra of B from Z-domain A:
  for (tempB = B; tempB; tempB = tempB->next) {
    ZPolyhedron *tmp;
    tmp = ZD_ZP_Difference(res, tempB);
    ZDomain_Free(res);
    res = tmp;
  }

  if (res == NULL)
    return (EmptyZPolyhedron(A->Lat->NbRows - 1));
  Canonical_ZDomain(res);

  return (res);
} /* ZDomainDifference */

/*
 * Return the image of the Z-domain 'A' under the invertible, affine, rational
 * transformation function 'Func'. The matrix representing the function 'Func'
 * must be non-singular and the number of rows of the function must be equal
 * to the number of rows in the matrix representing the lattice of 'A'.
 * Note:: Image((Z1 U Z2),F) = Image(Z1,F) U Image(Z2 U F).
 */
ZPolyhedron *ZDomainImage(ZPolyhedron *A, Matrix *Func) {

  ZPolyhedron *Result = NULL;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZDOMAINIMAGE\n");
  fclose(fp);
#endif

  for (ZPolyhedron *temp = A; temp; temp = temp->next) {
    ZPolyhedron *Zpol;
    Zpol = ZPolyhedronImage(temp, Func);
    Result = ZDconcatenate(Zpol, Result);
  }
  if (Result == NULL)
    return EmptyZPolyhedron(A->Lat->NbRows - 1);
  Canonical_ZDomain(Result);
  return Result;
} /* ZDomainImage */

/*
 * Return the preimage of the Z-domain 'A' under the invertible, affine, ratio-
 * nal transformation 'Func'. The number of rows of the matrix representing
 * the function 'Func' must be equal to the number of rows of the matrix repr-
 * senting the lattice of 'A'.
 */
ZPolyhedron *ZDomainPreimage(ZPolyhedron *A, Matrix *Func) {

  ZPolyhedron *Result = NULL;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZDOMAINPREIMAGE\n");
  fclose(fp);
#endif

  for (ZPolyhedron *temp = A; temp; temp = temp->next) {
    ZPolyhedron *Zpol;
    Zpol = ZPolyhedronPreimage(temp, Func);
    Result = ZDconcatenate(Zpol, Result);
  }

  if (Result == NULL)
    return (EmptyZPolyhedron(Func->NbColumns - 1));

  Canonical_ZDomain(Result);
  return Result;
} /* ZDomainPreimage */

/*
 * Return the Z-polyhedron intersection of the Z-polyhedra 'A' and 'B'.
 * We are based on the intersection of the two lattices of the polyhedra, named LInter.
 * If LInter is empty, we return the empty Zpolyhedron.
 * Otherwise, we calculate the intersection of the polyhedral images of A and B (PInter).
 * We calculate the Preimage of PInter by LInter and finally we allocate the result,
 * a Zpolyhedron allocated in canonical form.
 *
 *  /!\ USAGE: A and B contain a single Lattice, but can contain a polyhedral domain.
 */
static ZPolyhedron *ZPolyhedronIntersection(ZPolyhedron *A, ZPolyhedron *B) {

  ZPolyhedron *Result = NULL;
  Lattice *LInter;
  Polyhedron *PInter, *ImageA, *ImageB, *PreImage;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZPOLYHEDRONINTERSECTION\n");
  fclose(fp);
#endif

  LInter = LatticeIntersection(A->Lat, B->Lat);

  if (isEmptyLattice(LInter)) {
    Matrix_Free(LInter);
    return (EmptyZPolyhedron(A->Lat->NbRows - 1));
  }

  // TODO: this can be an over-approximation! -> need to handle LBLs!
  
  ImageA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
  ImageB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
  PInter = DomainIntersection(ImageA, ImageB, MAXNOOFRAYS);

  if (emptyQ(PInter))
    Result = EmptyZPolyhedron(LInter->NbRows - 1);
  else {
    PreImage = DomainPreimage(PInter, LInter, MAXNOOFRAYS);
    Result = ZPolyhedronAlloc(LInter, PreImage);
    Domain_Free(PreImage);
  }

  Matrix_Free(LInter);
  Domain_Free(PInter);
  Domain_Free(ImageB);
  Domain_Free(ImageA);

  return Result;
} /* ZPolyhedronIntersection */


/*
 * Return the difference A - B, between a ZDomain A and a ZPolyhedron B.
 * A can contain a list of lattices, B has a single lattice.
 * Algo: remove B from each part of A, and build a list of the result.
 *
 * /!\ USAGE: only the first lattice of B is considered (no list/Zdomains),
 * Creates a new allocated ZDomain, not necessarily in canonical form
 */
static ZPolyhedron *ZD_ZP_Difference(ZPolyhedron* A, ZPolyhedron* B) {
  ZPolyhedron *Result = NULL;

  for(ZPolyhedron *zp = A; zp; zp = zp->next) {
    ZPolyhedron *diff;

    diff = ZPolyhedronDifference(zp, B);
    // concatenate diff and result (not canonical)
    Result = ZDconcatenate(diff, Result);
  }

  // Result contains every piece of the solution,
  // but it is not necessarily in canonical form (will be handled be callee)
  return Result;
}

/*
 * Return the difference of two Z-polyhedra A and B.
 * inspired from the method Gautam describes in his thesis,
 * modified to handle LBLs.
 * A and B are Zpolyhedra, but the return value is a ZDomain!
 * Creates a new allocated ZDomain
 *
 * /!\ USAGE: only the first lattice of A and B is considered (no list/Zdomains),
 *            but A and B can contain a coordinate polyhedral domain (in ->Pol).
 */

static ZPolyhedron *ZPolyhedronDifference(ZPolyhedron* A, ZPolyhedron* B) {
  ZPolyhedron *Result = NULL, *Final_Result; // ZDomains
  ZPolyhedron *Ainter, *Binter; // ZPolyhedra
  LatticeUnion *LatDiff;
  Polyhedron *imA, *imB, *preimA, *ImDiff, *ImInter; // polyhedral domains

  if (A->Lat->NbRows != B->Lat->NbRows) {
    errormsg1("ZPolyhedronDifference", "dimincomp", "incompatible dimensions");
    return(NULL);
  }
  #ifdef DIFF_DEBUG
    fprintf(stderr, "-------Entering ZPDifference-------\nA = ");
    ZPolyhedronPrint(stderr, P_VALUE_FMT, A);
    fprintf(stderr, "and B = ");
    ZPolyhedronPrint(stderr, P_VALUE_FMT, B);
  #endif

  // treat the simple case where the Zpolyhedra do not intersect
  // the exact intersection Binter will be reused later
  Binter = ZPolyhedronIntersection(A, B);
  if(isEmptyZDomain(Binter)) {
    // if B does not intersect A, return A.
    #ifdef DIFF_DEBUG
      fprintf(stderr, "Binter=(A inter B) is empty, so B does not intersect A, we return A\n");
    #endif
    ZPolyhedron_Free(Binter);
    return(ZPolyhedron_Copy(A));
  }

  // Separate the computation in 3 phases:
  // 0. compute the difference of the image polyhedra P_A \ P_B (=temp) and
  //    add it to the solution Zpolyhedron (using lattice L_A).
  //    This is an over-approximation of A (but not B)
  // 1. compute the rest where the intersection of P_A and P_B have same
  //    dimensions
  // 2. intersect the result with A to get rid of the over-approximations

  // [STEP 0 (includes Gautam's Step 2)]
  imA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
  imB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
  ImDiff = DomainDifference(imA, imB, MAXNOOFRAYS);
  #ifdef DIFF_DEBUG
    fprintf(stderr, "ImDiff (hull of A that does not cover B) = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, ImDiff);
  #endif

  // Add (A->Lat, A - hull(B)) to the result:
  if (!emptyQ(ImDiff)) {
    Polyhedron *RedPolyDiff;
    RedPolyDiff = DomainPreimage(ImDiff, A->Lat, MAXNOOFRAYS);
    // NOTICE: this is an over-approximation of A
    Result = ZPolyhedronAlloc(A->Lat, RedPolyDiff);
    #ifdef DIFF_DEBUG
      fprintf(stderr, "Adding this to the temporary result: ");
      ZDomainPrint(stderr, P_VALUE_FMT, Result);
    #endif
    Domain_Free(RedPolyDiff);
  }

  // compute the images intersection of A and B
  ImInter = DomainIntersection(imA, imB, MAXNOOFRAYS);
  #ifdef DIFF_DEBUG
    fprintf(stderr, "ImInter (hull of A inter B) = ");
    Polyhedron_Print(stderr, P_VALUE_FMT, ImInter);
  #endif
  
  // compute the part of A that intersects the hull of B in the image space
  preimA = DomainPreimage(ImInter, A->Lat, MAXNOOFRAYS);
  Ainter = ZPolyhedronAlloc(A->Lat, preimA);
  // NOTICE: this Ainter can be a over-approximation of A

  Domain_Free(preimA);
  Domain_Free(ImDiff);
  Domain_Free(imA);
  Domain_Free(imB);

  // now Ainter and Binter have same lattices and polyhedra dimensions
  #ifdef DIFF_DEBUG
    fprintf(stderr, "-- [STEP1] now we compute the intersection on same lattice dimensions\n");
    fprintf(stderr, "Ainter = ");
    ZPolyhedronPrint(stderr, P_VALUE_FMT, Ainter);
    fprintf(stderr, "and Binter = ");
    ZPolyhedronPrint(stderr, P_VALUE_FMT, Binter);
  #endif

  // LatDiff (union of lattices) is the difference : (A->Lat) - (B->Lat) of same dimensions
  LatDiff = LatticeDifference(Ainter->Lat, Binter->Lat); 
  #ifdef DIFF_DEBUG
    if(!LatDiff)
      fprintf(stderr, "Empty Lattice difference\n");
  #endif

  // [STEP 1 of Gautam]:
  // Add all Z-polyhedra applying the (list of) lattice difference on ImInter
  for(LatticeUnion *tmp = LatDiff; tmp; tmp = tmp->next) {
    ZPolyhedron *Ztmp;
    #ifdef DIFF_DEBUG
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
  ZPolyhedron_Free(Ainter);
  ZPolyhedron_Free(Binter);

  if(Result == NULL) {
    #ifdef DIFF_DEBUG
      fprintf(stderr, "-- result = (NULL)\n");
    #endif
    return(EmptyZPolyhedron(A->Lat->NbRows - 1));
  }

  #ifdef DIFF_DEBUG
    fprintf(stderr, "-- temporary over-approximation of result = ");
    ZPolyhedronPrint(stderr, P_VALUE_FMT, Result);
  #endif
  // intersect the result with A to get the exact LBL.
  Final_Result = ZDomainIntersection(Result, A);
  ZDomain_Free(Result);

  return(Final_Result);

} /* ZPolyhedronDifference */



/*
 * Return the image of the Z-polyhedron 'ZPol' under the invertible, affine,
 * rational transformation function 'Func'. 
 * 
 * Algorithm:
 * - Multiply Lat by Func,
 * - Canonicalize the result (done in ZPAlloc)
 */
static ZPolyhedron *ZPolyhedronImage(ZPolyhedron *ZPol, Matrix *Func) {

  Matrix *newL;
  ZPolyhedron *result;

  if ((Func->NbColumns != ZPol->Lat->NbRows)) {
    errormsg1("ZPolyhedronImage", "dimincomp", "Incompatible dimensions");
    return NULL;
  }

  newL = Matrix_Alloc(Func->NbRows, ZPol->Lat->NbColumns);
  Matrix_Product(Func, ZPol->Lat, newL);
  result = ZPolyhedronAlloc(newL, ZPol->P);

  Matrix_Free(newL);
  return(result);
} /* ZPolyhedronImage */

/*
 * Return the preimage of the Z-polyhedron 'Zpol' under an affine
 * transformation function 'G'. The number of rows of matrix 'G' must
 * be equal to the number of rows of the matrix representing the
 * lattice of Zpol
 * Algorithm:
 * - build the Z-polyhedron { z' | Lz = Gz', z \in Zpol.P },
 * - remove z by normalizing the result (remove equalities)
 */
static ZPolyhedron *ZPolyhedronPreimage(ZPolyhedron *Z, Matrix *G) {

  ZPolyhedron *Result;
  Polyhedron *P, *newP;
  Matrix *Con;

  if(G->NbRows != Z->Lat->NbRows) {
    // G z' = L z
    errormsg1("ZPolyhedronPreimage", "dimincomp", "incompatible dimensions");
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

  Result = ZPolyhedronAlloc(newL, P);
  Domain_Free(P);
  Matrix_Free(newL);

  return(Result);
} /* ZPolyhedronPreimage */

// /*
//  * Return the Z-polyhderon 'Zpol' in canonical form: 'Result' (for the Z-poly-
//  * hedron in canonical form) and Basis 'Basis' (for the basis with respect to
//  * which 'Result' is in canonical form.
//  */
// void CanonicalForm(ZPolyhedron *Zpol, ZPolyhedron **Result, Matrix **Basis) {

//   Matrix *B1 = NULL, *B2 = NULL, *T1, *B2inv;
//   int i, l1, l2;
//   Value tmp;
//   Polyhedron *Image, *ImageP;
//   Matrix *H, *U, *temp, *Hprime, *Uprime, *T2;

// #ifdef DOMDEBUG
//   FILE *fp;
//   fp = fopen("_debug", "a");
//   fprintf(fp, "\nEntered CANONICALFORM\n");
//   fclose(fp);
// #endif

//   if (isEmptyZDomain(Zpol)) {
//     Basis[0] = Identity(Zpol->Lat->NbRows);
//     Result[0] = ZDomain_Copy(Zpol);
//     return;
//   }
//   value_init(tmp);
//   l1 = FindHermiteBasisofDomain(Zpol->P, &B1);
//   Image = DomainImage(Zpol->P, Zpol->Lat, MAXNOOFRAYS);
//   l2 = FindHermiteBasisofDomain(Image, &B2);

//   if (l1 != l2)
//     fprintf(stderr, "In CNF : Something wrong with the Input Zpolyhedra \n");

//   B2inv = Matrix_Alloc(B2->NbRows, B2->NbColumns);
//   temp = Matrix_Copy(B2);
//   Matrix_Inverse(temp, B2inv);
//   Matrix_Free(temp);

//   temp = Matrix_Alloc(B2inv->NbRows, Zpol->Lat->NbColumns);
//   T1 = Matrix_Alloc(temp->NbRows, B1->NbColumns);
//   Matrix_Product(B2inv, Zpol->Lat, temp);
//   Matrix_Product(temp, B1, T1);
//   Matrix_Free(temp);

//   T2 = ChangeLatticeDimension(T1, l1);
//   temp = ChangeLatticeDimension(T2, T2->NbRows + 1);

//   /* Adding the affine part */
//   for (i = 0; i < l1; i++)
//     value_assign(temp->p[i][temp->NbColumns - 1], T1->p[i][T1->NbColumns - 1]);

//   AffineHermite(temp, &H, &U);
//   Hprime = ChangeLatticeDimension(H, Zpol->Lat->NbRows);

//   /* Exchanging the Affine part */
//   for (i = 0; i < l1; i++) {
//     value_assign(tmp, Hprime->p[i][Hprime->NbColumns - 1]);
//     value_assign(Hprime->p[i][Hprime->NbColumns - 1],
//                  Hprime->p[i][H->NbColumns - 1]);
//     value_assign(Hprime->p[i][H->NbColumns - 1], tmp);
//   }
//   Uprime = ChangeLatticeDimension(U, Zpol->Lat->NbRows);

//   /* Exchanging the Affine part */
//   for (i = 0; i < l1; i++) {
//     value_assign(tmp, Uprime->p[i][Uprime->NbColumns - 1]);
//     value_assign(Uprime->p[i][Uprime->NbColumns - 1],
//                  Uprime->p[i][U->NbColumns - 1]);
//     value_assign(Uprime->p[i][U->NbColumns - 1], tmp);
//   }
//   Polyhedron_Free(Image);
//   Matrix_Free(B2inv);
//   B2inv = Matrix_Alloc(B1->NbRows, B1->NbColumns);
//   Matrix_Inverse(B1, B2inv);
//   ImageP = DomainImage(Zpol->P, B2inv, MAXNOOFRAYS);
//   Matrix_Free(B2inv);
//   Image = DomainImage(ImageP, Uprime, MAXNOOFRAYS);
//   Domain_Free(ImageP);
//   Result[0] = ZPolyhedronAlloc(Hprime, Image);
//   Basis[0] = Matrix_Copy(B2);

//   /* Free the memory used */
//   Polyhedron_Free(Image);
//   Matrix_Free(B1);
//   Matrix_Free(B2);
//   Matrix_Free(temp);
//   Matrix_Free(T1);
//   Matrix_Free(T2);
//   Matrix_Free(H);
//   Matrix_Free(U);
//   Matrix_Free(Hprime);
//   Matrix_Free(Uprime);
//   value_clear(tmp);
//   return;
// } /* CanonicalForm */

// /*
//  * Given a Z-polyhedron 'A' in which the Lattice is not integral, return the
//  * Z-polyhedron which contains all the integral points in the input lattice.
//  */
// ZPolyhedron *IntegraliseLattice(ZPolyhedron *A) {

//   ZPolyhedron *Result;
//   Lattice *M = NULL, *Id;
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
//     Result = EmptyZPolyhedron(A->Lat->NbRows - 1);
//   else {
//     Preim = DomainPreimage(Im, M, MAXNOOFRAYS);
//     Result = ZPolyhedronAlloc(M, Preim);
//   }
//   Matrix_Free(M);
//   Domain_Free(Im);
//   Domain_Free(Preim);
//   return Result;
// } /* IntegraliseLattice */

/*
 * Return the simplified representation of the Z-domain 'ZDom'. It attempts to
 * convexize unions of polyhedra when they correspond to the same lattices and
 * to simplify union of lattices when they correspond to the same polyhdera.
 */
ZPolyhedron *ZDomainSimplify(ZPolyhedron *ZDom) {

  ZPolyhedron *Ztmp, *Result;
  ForSimplify *Head, *Prev, *Curr;
  ZPolyhedron *ZDomHead, *Emp;

  if (ZDom == NULL) {
    fprintf(stderr, "\nError in ZDomainSimplify - ZDomHead = NULL\n");
    return NULL;
  }
  if (ZDom->next == NULL)
    return (ZPolyhedron_Copy(ZDom));
  Emp = EmptyZPolyhedron(ZDom->Lat->NbRows - 1);
  ZDomHead = ZDomainUnion(ZDom, Emp);
  ZPolyhedron_Free(Emp);
  Head = NULL;
  Ztmp = ZDomHead;
  do {
    Polyhedron *Img;
    Img = DomainImage(Ztmp->P, Ztmp->Lat, MAXNOOFRAYS);
    for (Curr = Head; Curr != NULL; Curr = Curr->next) {
      Polyhedron *Diff1;
      Bool flag = False;

      Diff1 = DomainDifference(Img, Curr->Pol, MAXNOOFRAYS);
      if (emptyQ(Diff1)) {
        Polyhedron *Diff2;

        Diff2 = DomainDifference(Curr->Pol, Img, MAXNOOFRAYS);
        if (emptyQ(Diff2))
          flag = True;
        Domain_Free(Diff2);
      }
      Domain_Free(Diff1);
      if (flag == True) {
        LatticeUnion *temp;

        temp = (LatticeUnion *)malloc(sizeof(LatticeUnion));
        temp->M = (Lattice *)Matrix_Copy((Matrix *)Ztmp->Lat);
        temp->next = Curr->LatUni;
        Curr->LatUni = temp;
        break;
      }
    }
    if (Curr == NULL) {
      Curr = (ForSimplify *)malloc(sizeof(ForSimplify));
      Curr->Pol = Domain_Copy(Img);
      Curr->LatUni = (LatticeUnion *)malloc(sizeof(LatticeUnion));
      Curr->LatUni->M = (Lattice *)Matrix_Copy((Matrix *)Ztmp->Lat);
      Curr->LatUni->next = NULL;
      Curr->next = Head;
      Head = Curr;
    }
    Domain_Free(Img);
    Ztmp = Ztmp->next;
  } while (Ztmp != NULL);

  for (Curr = Head; Curr != NULL; Curr = Curr->next)
    Curr->LatUni = LatticeSimplify(Curr->LatUni);
  Result = NULL;
  for (Curr = Head; Curr != NULL; Curr = Curr->next) {
    LatticeUnion *L;
    for (L = Curr->LatUni; L != NULL; L = L->next) {
      Polyhedron *Preim;
      ZPolyhedron *Zpol;

      Preim = DomainPreimage(Curr->Pol, L->M, MAXNOOFRAYS);
      Zpol = ZPolyhedronAlloc(L->M, Preim);
      Zpol->next = Result;
      Result = Zpol;
      Domain_Free(Preim);
    }
  }
  Curr = Head;
  while (Curr != NULL) {
    Prev = Curr;
    Curr = Curr->next;
    LatticeUnion_Free(Prev->LatUni);
    Domain_Free(Prev->Pol);
    free(Prev);
  }
  return Result;
} /* ZDomainSimplify */

/*
 * 
 *
*/

ZPolyhedron *SplitZpolyhedron(ZPolyhedron *ZPol, Lattice *B) {

  Matrix *H, *U1, *X, *Y;
  ZPolyhedron *zpnew, *Result;
  LatticeUnion *Head = NULL, *tempHead = NULL;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered SplitZpolyhedron \n");
  fclose(fp);
#endif

  if (B->NbRows != B->NbColumns) {
    fprintf(
        stderr,
        "\n SplitZpolyhedron : The Input Matrix B is not a proper Lattice \n");
    return NULL;
  }

  if (ZPol->Lat->NbRows != B->NbRows) {
    fprintf(stderr,
            "\nSplitZpolyhedron : The Lattice in Zpolyhedron and B have ");
    fprintf(stderr, "incompatible dimensions \n");
    return NULL;
  }

  if (isNormalLattice(ZPol->Lat) != True) {
    AffineHermite(ZPol->Lat, &H, &U1);
    X = Matrix_Copy(H);
    Matrix_Free(U1);
    Matrix_Free(H);
  } else
    X = Matrix_Copy(ZPol->Lat);

  if (isNormalLattice(B) != True) {
    AffineHermite(B, &H, &U1);
    Y = Matrix_Copy(H);
    Matrix_Free(H);
    Matrix_Free(U1);
  } else
    Y = Matrix_Copy(B);
  if (isEmptyLattice(X)) {
    return NULL;
  }

  Head = Lattice2LatticeUnion(X, Y);

  /* If the spliting operation can't be done the result is the original
   * Zplyhedron. */

  if (Head == NULL) {
    Matrix_Free(X);
    Matrix_Free(Y);
    return ZPolyhedron_Copy(ZPol);
  }

  Result = NULL;

  while (Head) {
    tempHead = Head;
    Head = Head->next;
    zpnew = ZPolyhedronAlloc(tempHead->M, ZPol->P);
    Result = ZDconcatenate(zpnew, Result);
    tempHead->next = NULL;
    Matrix_Free(tempHead->M);
    free(tempHead);
  }
  Matrix_Free(X);
  Matrix_Free(Y);
  return Result;
} /* SplitZpolyhedron */

/*
 * get the matrix of equalities from a polyhedron
 * (without the first columns of 0's)
 */
static Matrix *get_equalities(Polyhedron *P)
{
  // Eq is the matrix of linear equations of P (including the constant)
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
 * Remove the equalities from (A->Lat, A->P).
 * In place. A->P is a domain.
 */
static void ZP_Remove_Equalities(ZPolyhedron *A, Matrix *Equalities)
{
  // if A->P has equalities, remove them and spread the lattice
  if (A->P->Dimension > 0 && A->P->NbEq != 0) {
    Matrix *ker, *H = NULL, *NewL;

    // remove equalities in domain P and change Lat to spread the original space
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "P has equalities\n");
      fprintf(stderr, "Equality matrix (including constants): ");
      Matrix_Print(stderr, P_VALUE_FMT, Equalities);
    #endif

    // compute Ker(Eq)
    ker = int_ker(Equalities);
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "ker of eq: ");
      Matrix_Print(stderr, P_VALUE_FMT, ker);
    #endif
    
    Matrix_Move_Homogeneous_Dim_First(ker);
    left_hermite(ker, &H, NULL, NULL);
    Matrix_Move_Homogeneous_Dim_Last(H);
    Matrix_Free(ker);

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "Matrix H: ");
      Matrix_Print(stderr, P_VALUE_FMT, H);
      fprintf(stderr, "Lattice of A: ");
      Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
    #endif

    // TODO: CHECK THAT THIS IS CORRECT!

    // if the bottom right value of H is not one, this means that
    // the transformation matrix is not integer
    // just make it integer and do not bother about never taken
    // rational values.
    value_set_si(H->p[H->NbRows-1][H->NbColumns-1], 1);

    // NewL = L . H
    NewL = Matrix_Alloc(A->Lat->NbRows, H->NbColumns);
    Matrix_Product(A->Lat, H, NewL);
    // NewP = H^{-1} . P
    Polyhedron* NewP = DomainPreimage(A->P, H, MAXNOOFRAYS);

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
} /* ZP_Remove_Equalities */

/*
 * Set the lattice to normal form in (A->Lat, A->P).
 * In place. A->P is a domain.
 */
static void ZP_Normalize_Lat(ZPolyhedron *A)
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

    // remove the columns of zero's
    int nbZeros = count_zeroCols(H);
    if(nbZeros) {
      Matrix* Transformation = Matrix_Alloc(A->Lat->NbColumns-nbZeros, A->Lat->NbColumns);
      for (int  i = 0; i < Transformation->NbRows; i++) {
        for (int j = 0; j < Transformation->NbColumns; j++) {
          if(i==j && i!=Transformation->NbRows-1) {
            value_set_si(Transformation->p[i][j],1);
          }else{
            value_set_si(Transformation->p[i][j],0);
          }
        }
      }
      value_set_si(Transformation->p[Transformation->NbRows-1][Transformation->NbColumns-1],1);

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
        ZPolyhedronPrint(stderr, P_VALUE_FMT, A);
      #endif
    }


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
}
/*
 * The function takes a Zpolyhedron 
 * --- containing a domain (list of polyhedra), and a single lattice ---
 * and modifies it in place to be in canonical form as described by Gautam
 * (A->Lat in HNF and no equalities in A->P)
 * IN PLACE: modifies A itself
 * 
 * WARNING: this function transforms a Zpolyhedron into a Zdomain: it may
 * add new ZPolyhedra to the ZDomain (list of ZPolyhedra) just after A
 */
static void Canonical_ZPolyhedron_Gautam(ZPolyhedron* A) {
  
  #ifdef CANONICAL_DEBUG
  fprintf(stderr, "Entering CanonicalZPGautam\n");
  fprintf(stderr, "--------- Input Lat: ");
  Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
  fprintf(stderr, "--------- Input P: ");
  Polyhedron_Print(stderr, P_VALUE_FMT, A->P);
  #endif
  
  if (A->P->Dimension+1 != A->Lat->NbColumns) {
    errormsg1("Polyhedron_Image", "dimincomp", "incompatible dimensions");
    return;
  }

  // if the polyhedron is empty
  if(emptyQ(A->P)) {
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "The polyhedron is empty\n");
    #endif

    if(A->next) {
      #ifdef CANONICAL_DEBUG
        fprintf(stderr, "... But the next one is not\n");
      #endif
      // if there is something else after an empty ZP, need to replace the
      // current ZP with the next ZP.
      // Replace A with next and free A->next
      ZPolyhedron *remove;
      remove = A->next;
      Domain_Free(A->P);
      Matrix_Free(A->Lat);
      A->P = remove->P;
      A->Lat = remove->Lat;
      A->next = remove->next;
      free(remove);

      // now, canonicalize the (new) current ZP itself:
      Canonical_ZPolyhedron_Gautam(A);

      return;
    }

    // A is empty and alone.
    // Verify that it is canonical and return.
    int dimension = A->Lat->NbRows;
    if(A->P->Dimension > 0) {
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
      fprintf(stderr, "We return the empty ZPolyhedron: ");
      ZPolyhedronPrint(stderr, P_VALUE_FMT, A);
    #endif
    return;
  }

  // A->next will be treated separately by the callee

  // change P such that all polyhedra in this list have the same set of equalities,
  // that is, the equalities of the first one.
  // all the other ones are added to a new ZPolyhedron, linked to (ZDomain) A in A->next
  #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Checking for equalites in P\n");
  #endif
  ZPolyhedron *new = NULL;
  Matrix * Equalities = get_equalities(A->P); // get eq from the first one
  Polyhedron *nextpp, *prevpp = A->P;

  for(Polyhedron *pp = A->P->next; pp; prevpp = pp, pp = nextpp) {
    // check that the equalities of pp->Constraints are the same as the ones of matrix Equalities.
    if(!same_equalities(Equalities, pp)) {
      // here, get pp out.
      if(!new) {
        new = malloc(sizeof(*new));
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
    // include new in the ZDomain list A
    new->next = A->next;
    A->next = new;
  }
  
  // NOW REMOVE EQUALITIES FROM A (Lat, P)
  ZP_Remove_Equalities(A, Equalities);
  Matrix_Free(Equalities);


  // this was useful if for some reason the homogeneous dimension spread
  // to something different than 1... (bottom right value of Lat matrix)
  // seems impossible with Hermite(homogeneous_dim_first(Ker(Eq)))
#ifdef I_THINK_THAT_THIS_IS_NOT_NECESSARY
  // Check if the constant part (last column of Lat) is normal (simplified by its gcd)
  // if not, divide it by the gcd, and compute the new polyhedron P' = image(T, P) with
  // T = Id   0
  //      0  gcd
  Value gcd;
  value_init(gcd);
  // calul gcd sur la dernière colonne
  value_absolute(gcd,A->Lat->p[0][A->Lat->NbColumns-1]);
  // value_print(stderr, P_VALUE_FMT, gcd);
  for (int i = 1; i < A->Lat->NbRows; i++){
    Gcd(gcd,A->Lat->p[i][A->Lat->NbColumns-1],&gcd);
  }

  // si gcd != 1 faire la simplification
  if (value_notone_p(gcd)) {
    // diviser la dernière colonne de A->Lat (en place)
    for (int i = 0; i < A->Lat->NbRows; i++)
    {
      value_division(A->Lat->p[i][A->Lat->NbColumns-1], A->Lat->p[i][A->Lat->NbColumns-1], gcd);
    }
    // et construire T, puis faire l'image par T de A->P
    Matrix *T = Identity(A->P->Dimension + 1); // ajouter la constante (gcd) en bas à droite
    for (int i = 0; i < T->NbColumns; i++){
      for (int j = 0; j < T->NbColumns; j++){
        if (i==j){
          value_set_si(T->p[i][j],1);
        }
        else{
          value_set_si(T->p[i][j],0);
        } 
      }  
    }
    value_assign(T->p[T->NbRows-1][T->NbColumns-1],gcd);
    Polyhedron* tmp=Polyhedron_Image(A->P,T,MAXNOOFRAYS);
    Matrix_Free(T);
    Polyhedron_Free(A->P);
    A->P=tmp;
  }
  value_clear(gcd);
#endif


  // NOW NORMALIZE LATTICE A->Lat
  ZP_Normalize_Lat(A);

  return;
} /* Canonical_ZPolyhedron_Gautam */

/*
 * The function takes a ZDomain
 * (single or multiple lattices and single or multiple polyhedra)
 * and transforms it into a canonical form of the ZDomain as described
 * by Gautam:
 *  - all lattices in HNF, and
 *  -  no equalities in all polyhedral domains.
 * Performs the operation IN PLACE (modifies A)
 */
void Canonical_ZDomain(ZPolyhedron *A) {

  // here, just transform every ZPolyhedron of the list individually
  // careful, this may add a new ZPolyhedron to the list A itself
  // (after zp) but they will be scanned by this loop :)
  for(ZPolyhedron *zp = A; zp; zp = zp->next) {
    Canonical_ZPolyhedron_Gautam(zp);
  }

  // check if Ztmp->Lat is present twice in Result, and
  // if it is, add this polyhedron to the existing one
  // and remove the second reference
  for(ZPolyhedron *zp = A; zp; zp = zp->next) {
    ZPolyhedron *ZZ;
    if((ZZ = FindLatticePred(zp->Lat, zp))) {
      ZPolyhedron *remove;
      Polyhedron *nextpp;
      // add all polyhedra of the domain ZZ->next->P to zp
      // consumes ZZ->next->P.
      Polyhedron *pp = ZZ->next->P;
      while(pp) {
        nextpp = pp->next;
        pp->next = NULL;
        // this consumes pp, so need to get next before
        zp->P = AddPolyToDomain(pp, zp->P);
        pp = nextpp;
      }
      // remove ZZ->next by changing the linked list
      remove = ZZ->next;
      ZZ->next = ZZ->next->next;
      Matrix_Free(remove->Lat);
      free(remove);
    }
  }
} /* Canonical_ZDomain */

/*
 * Find if a given lattice is present in a zpolyhedron.
 * Returns the address of the ***previous*** zpolyhedron (such that ZZ->next->Lat == L),
 * NULL if not found
 */
static ZPolyhedron *FindLatticePred(Lattice *L, ZPolyhedron *A) {
  ZPolyhedron* tmp;

  for(tmp = A; tmp->next; tmp=tmp->next) {
    if(sameLattice(L, tmp->next->Lat)) {
      return (tmp);
    }
  }
  return (NULL);
} /* FindLatticePred */

int count_zeroCols (Matrix* M){
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
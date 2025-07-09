#include <polylib/polylib.h>
#include <stdlib.h>

// debug this file:
// #define DEBUG
#ifdef DEBUG
  #define CANONICAL_DEBUG 1
  #define NEWINTERSECTION_DEBUG 1
  #define ADDZPTOZD_DEBUG 1
#endif

static ZPolyhedron *ZPolyhedronIntersection(ZPolyhedron *, ZPolyhedron *);
static ZPolyhedron *ZPolyhedron_Copy(ZPolyhedron *A);
static void ZPolyhedron_Free(ZPolyhedron *Zpol);
static ZPolyhedron *ZPolyhedronDifference(ZPolyhedron *, ZPolyhedron *);
ZPolyhedron *ZPolyhedronDifferenceGautam(ZPolyhedron *, ZPolyhedron *);
static ZPolyhedron *ZPolyhedronImage(ZPolyhedron *, Matrix *);
static ZPolyhedron *ZPolyhedronPreimage(ZPolyhedron *, Matrix *);
static ZPolyhedron *AddZPolytoZDomain(ZPolyhedron *A, ZPolyhedron *Head);
static void ZPolyhedronPrint(FILE *fp, const char *format, ZPolyhedron *A);
static void CanonicalZPGautam(ZPolyhedron* A);

typedef struct forsimplify {
  Polyhedron *Pol;
  LatticeUnion *LatUni;
  struct forsimplify *next;
} ForSimplify;

/*
 * Returns True if 'Zpol' is empty, otherwise returns False
 */
Bool isEmptyZPolyhedron(ZPolyhedron *Zpol) {

  if (Zpol == NULL)
    return True;
  if (emptyQ(Zpol->P))
    return True;
  return False;
} /* isEmptyZPolyhedron */

/*
 * Given Lattice 'Lat' and a Polyhderon 'Poly', allocate space, and return
 * the Z-polyhderon corresponding to the image of the polyhderon 'Poly' by the
 * lattice 'Lat'. If the input lattice 'Lat' is not integeral, it integralises
 * it, i.e. the lattice of the Z-polyhderon returned is integeral.
 */
ZPolyhedron *ZPolyhedron_Alloc(Lattice *Lat, Polyhedron *Poly) {

  ZPolyhedron *A;

  POL_ENSURE_FACETS(Poly);
  POL_ENSURE_VERTICES(Poly);

  if (Lat->NbColumns != Poly->Dimension + 1) {
    fprintf(stderr, "\nInZPolyAlloc - The Lattice  and the Polyhedron");
    fprintf(stderr, " are not compatible to form a ZPolyhedra\n");
    return NULL;
  }
  // if ((!(isEmptyLattice(Lat))) && (!isfulldim(Lat))) {
  //   fprintf(stderr, "\nZPolAlloc: Lattice not Full Dimensional\n");
  //   fprintf(stderr, "\nZPolAlloc: is empty latice: %d \n",((isEmptyLattice(Lat))? 1: 0));
  //   fprintf(stderr, "\nZPolAlloc: is fulldim latice: %d \n",((isfulldim(Lat))? 1: 0));
  //   Matrix_Print(stderr,P_VALUE_FMT,Lat);
  //   return NULL;
  // }
  A = (ZPolyhedron *)malloc(sizeof(ZPolyhedron));
  if (!A) {
    fprintf(stderr, "ZPolAlloc : Out of Memory\n");
    return NULL;
  }
  A->next = NULL;
  A->P = Domain_Copy(Poly);
  A->Lat = Matrix_Copy(Lat);

  CanonicalZPGautam(A);
  return A;
} /* ZPolyhedron_Alloc */

/*
 * Free the memory used by the Z-domain 'Head'
 */
void ZDomain_Free(ZPolyhedron *Head) {

  if (Head == NULL)
    return;
  if (Head->next != NULL)
    ZDomain_Free(Head->next);
  ZPolyhedron_Free(Head);
} /* ZDomain_Free */

/*
 * Free the memory used by the Z-polyhderon 'Zpol'
 */
static void ZPolyhedron_Free(ZPolyhedron *Zpol) {

  if (Zpol == NULL)
    return;
  if(Zpol->Lat)
    Matrix_Free((Matrix *)Zpol->Lat);
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
 * Return a copy of the Z-polyhderon 'A'
 */
static ZPolyhedron *ZPolyhedron_Copy(ZPolyhedron *A) {

  ZPolyhedron *Zpol;

  Zpol = ZPolyhedron_Alloc(A->Lat, A->P);
  return Zpol;
} /* ZPolyhderon_Copy */

/*
 * Add the ZPolyhedron 'Zpol' to the Z-domain 'Result' and return a pointer
 * to the new Z-domain.
 */
static ZPolyhedron *AddZPoly2ZDomain(ZPolyhedron *Zpol, ZPolyhedron *Result) {

  ZPolyhedron *A;

  if (isEmptyZPolyhedron(Zpol))
    return Result;
  A = ZPolyhedron_Copy(Zpol);
  A->next = NULL;

  if (isEmptyZPolyhedron(Result)) {
    ZDomain_Free(Result);
    return A;
  }
  A->next = Result;
  return A;
} /* AddZPoly2ZDomain */

/*
 * Given a Z-polyhderon 'A' and a Z-domain 'Head', return a new Z-domain with
 * 'A' added to it. If the new Z-polyhedron 'A', is already included in the
 * Z-domain 'Head', it is not added in the list. Othewise, the function checks
 * if the new Z-polyhedron 'A' to be added to the Z-domain 'Head' has a common
 * lattice with some other Z-polyhderon already present in the Z-domain. If it
 * is so, it takes the union of the underlying polyhdera; domains and returns.
 * The function tries to make sure that the added Z-polyhedron 'A' is in the
 * canonical form.
 */
static ZPolyhedron *AddZPolytoZDomain(ZPolyhedron *A, ZPolyhedron *Head) {

  ZPolyhedron *Zpol, *temp, *temp1;
  Polyhedron *i;
  Bool Added;

  #ifdef ADDZPTOZD_DEBUG
    fprintf(stderr, "----- ENTERING ADD ZP TO ZD -----\n");
  #endif
  if ((A == NULL) || (isEmptyZPolyhedron(A)))
    return Head;

  /* For each "underlying" Pol, find the Cnf and add Zpol in Cnf*/
  for (i = A->P; i != NULL; i = i->next) {
    ZPolyhedron *Z, *Z1;
    Polyhedron *Image;
    Matrix *H, *U;
    Lattice *Lat;

    Added = False;
    // copy of ZP(Lat, i)
    Polyhedron *tmp = i->next;
    i->next = NULL;
    Zpol = ZPolyhedron_Alloc(A->Lat, i);
    i->next = tmp;
    CanonicalZPGautam(Zpol);
    // CanonicalForm(Z1, &Z, &H);
    // ZDomain_Free(Z1);
    // Lat = (Lattice *)Matrix_Alloc(H->NbRows, Z->Lat->NbColumns);
    // Matrix_Product(H, Z->Lat, (Matrix *)Lat);
    // Matrix_Free(H);
    // AffineHermite(Lat, (Lattice **)&H, &U);
    // Image = DomainImage(Z->P, U, MAXNOOFRAYS);
    // ZDomain_Free(Z);
    #ifdef ADDZPTOZD_DEBUG
      fprintf(stderr, "Head =\n");
      ZDomainPrint(stderr, P_VALUE_FMT, Head);
      fprintf(stderr, "Adding Zpol:\n");
      ZPolyhedronPrint(stderr, P_VALUE_FMT, Zpol);
    #endif

    if ((Head == NULL) || (isEmptyZPolyhedron(Head))) {
      if(Head != NULL)
        ZPolyhedron_Free(Head);
      Head = Zpol;
      continue;
    }

    /* Check if the curr pol is included in the zpol or vice versa. */
    for (temp1 = temp = Head; temp != NULL; temp = temp->next) {
      if (ZPolyhedronIncludes(Zpol, temp)) {
        // already there, go to next
        ZPolyhedron_Free(Zpol);
        Added = True;
        break;
      }
      else if (ZPolyhedronIncludes(temp, Zpol)) {
        if (temp == Head) {
          Zpol->next = temp->next;
          Head = Zpol;
          ZPolyhedron_Free(temp);
          Added = True;
          break;
        }
        temp1->next = Zpol;
        Zpol->next = temp->next;
        ZPolyhedron_Free(temp);
        Added = True;
        break;
      }
      temp1 = temp;
    }
    if (Added == True)
      continue;

    // check if the same lattice is already there, and add P to this one if it is
    for (temp = Head; temp != NULL; temp = temp->next) {
      if (sameLattice(temp->Lat, Zpol->Lat)) {
        Polyhedron *Union;

        // Add Zpol->P to temp->P
        temp->P = AddPolyToDomain(Zpol->P, temp->P);
        Zpol->P = NULL;
        ZPolyhedron_Free(Zpol);
        Added = True;
        break;
      }
      temp1 = temp;
    }
    if (Added == False)
      temp1->next = Zpol;
  }

  #ifdef ADDZPTOZD_DEBUG
    fprintf(stderr, "Final Head =\n");
    ZDomainPrint(stderr, P_VALUE_FMT, Head);
    fprintf(stderr, "----- EXIT ADD ZP TO ZD -----\n");
  #endif

  return Head;
} /* AddZPolytoZDomain */

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
    for(int j=0 ; j<=dimension; j++) {
        value_set_si(E->p[0][j], 0);
    }
  value_set_si(E->p[0][E->NbRows-1],0);

  P = Empty_Polyhedron(0);

  Zpol = ZPolyhedron_Alloc(E, P);
  Matrix_Free((Matrix *)E);
  Domain_Free(P);
  return Zpol;
} /* EmptyZPolyhedron */

/*
 * Given Z-domains 'A' and 'B', return True if A is included in 'B', otherwise
 * return False.
 */
Bool ZDomainIncludes(ZPolyhedron *A, ZPolyhedron *B) {

  ZPolyhedron *Diff;
  Bool ret = False;

  Diff = ZDomainDifference(A, B);
  if (isEmptyZPolyhedron(Diff))
    ret = True;

  ZDomain_Free(Diff);
  return ret;
} /* ZDomainIncludes */

/*
 * Given Z-polyhedra 'A' and 'B', return True if 'A' is included in 'B',
 * otherwise return False
 */
Bool ZPolyhedronIncludes(ZPolyhedron *A, ZPolyhedron *B) {

  Polyhedron *Diff = NULL;
  Bool retval = False;
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZPOLYHEDRONINCLUDES\n");
  fclose(fp);
#endif
  
  if (LatticeIncludes(A->Lat, B->Lat) == True) {
    Polyhedron *ImageA, *ImageB;

    ImageA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
    ImageB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);

    Diff = DomainDifference(ImageA, ImageB, MAXNOOFRAYS);
    if (emptyQ(Diff))
      retval = True;

    Domain_Free(ImageA);
    Domain_Free(ImageB);
    Domain_Free(Diff);
  }
  return retval;
} /* ZPolyhedronIncludes */

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

  for( ; A; A=A->next)
  {
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

  ZPolyhedron *Result = NULL, *temp;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZDOMAINUNION\n");
  fclose(fp);
#endif

  for (temp = A; temp != NULL; temp = temp->next)
    Result = AddZPolytoZDomain(temp, Result);
  for (temp = B; temp != NULL; temp = temp->next)
    Result = AddZPolytoZDomain(temp, Result);
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
      Result = AddZPolytoZDomain(Zpol, Result);
      ZPolyhedron_Free(Zpol);
    }
  if (Result == NULL)
    return EmptyZPolyhedron(A->Lat->NbColumns - 1);
  return Result;
} /* ZDomainIntersection */

/*
 * Return the Z-domain difference of the domains 'A' and 'B'. The dimensions of
 * the Z-domains 'A' and 'B' must be equal. Note that the difference of two
 * Z-polyhedra is a Union of Z-polyhedra. The algorithms is as given below :-
 * Algorithm: (Given Z-domains A and B)
 *           Result <-- NULL
 *           for every Z-polyhderon Zpoly of A {
 *               temp <-- Zpoly;
 *               for every Z-polyhderon Z1 of B
 *                  temp = temp - Z1;
 *               }
 *           Add temp to Result;
 *           return;
 */
ZPolyhedron *ZDomainDifference(ZPolyhedron *A, ZPolyhedron *B) { 

  // need to close the hole where when we have A < B we have as result the empty, already works for when A==B so take a look at that.

  ZPolyhedron *Result = NULL, *tempA = NULL, *tempB = NULL, *test=NULL;
  ZPolyhedron *templist, *res, *i, *j;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZDOMAINDIFFERENCE\n");
  fclose(fp);
#endif

  if (A->Lat->NbRows != B->Lat->NbRows) {
    fprintf(stderr, "In ZDomainDifference : the input ZDomains ");
    fprintf(stderr, "do not have compatible dimensions\n");
    fprintf(stderr, "ZDomainDifference not performed\n");
    return NULL;
  }

  // //the test to see if the Zdomains have an intersection or not
  // Bool EmptyIntersection = True;

  // for (tempA = A; tempA != NULL; tempA = tempA->next) {//test all polyhedrons in A

  //   for ( tempB = B; tempB != NULL; tempB = tempB->next){//test all polyhedrons in B

  //     test=ZPolyhedronIntersection(tempA,tempB);//find their intersection

  //     if ( !isEmptyZPolyhedron(test)){// if we find an intersection
  //       EmptyIntersection=False;//set the boolean to false 
  //       test=NULL;
  //       break;//no need to test anymore
  //     }
  //   }
  // }
  
  // tempA=NULL; tempB=NULL;//reset A and B

  // if ( EmptyIntersection ){ //if there is an empty inteserction
  //   printf("There are no elements in common between the two polyhedrons");
  //   Result=ZDomain_Copy(A);
  //   return Result;//return A unchanged 
  // }

  Result = NULL;
  for (tempA = A; tempA != NULL; tempA = tempA->next) {
    ZPolyhedron *temp = NULL;

    res = ZPolyhedron_Copy(tempA);
    for (tempB = B; tempB != NULL; tempB = tempB->next) {
      ZPolyhedron *tmpres = res;
      res = ZPolyhedronDifferenceGautam(tmpres, tempB);
      ZDomain_Free(tmpres);
    }
    // here: res = tempA - B

    for(i=res ; i != NULL ; i = i->next) {
      Result = AddZPolytoZDomain(i, Result);
    }

    ZDomain_Free(res);
  }
 
  if (Result == NULL)
    return (EmptyZPolyhedron(A->Lat->NbRows - 1));
  return Result;
} /* ZDomainDifference */

/*
 * Return the image of the Z-domain 'A' under the invertible, affine, rational
 * transformation function 'Func'. The matrix representing the function 'Func'
 * must be non-singular and the number of rows of the function must be equal
 * to the number of rows in the matrix representing the lattice of 'A'.
 * Note:: Image((Z1 U Z2),F) = Image(Z1,F) U Image(Z2 U F).
 */
ZPolyhedron *ZDomainImage(ZPolyhedron *A, Matrix *Func) {

  ZPolyhedron *Result = NULL, *temp;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZDOMAINIMAGE\n");
  fclose(fp);
#endif

  for (temp = A; temp != NULL; temp = temp->next) {
    ZPolyhedron *Zpol;
    Zpol = ZPolyhedronImage(temp, Func);
    if(!(isEmptyZPolyhedron(Zpol)))
    {
      Result = AddZPolytoZDomain(Zpol, Result);
    }
    ZPolyhedron_Free(Zpol);
  }
  if (Result == NULL)
    return EmptyZPolyhedron(A->Lat->NbRows - 1);
  return Result;
} /* ZDomainImage */

/*
 * Return the preimage of the Z-domain 'A' under the invertible, affine, ratio-
 * nal transformation 'Func'. The number of rows of the matrix representing
 * the function 'Func' must be equal to the number of rows of the matrix repr-
 * senting the lattice of 'A'.
 */
ZPolyhedron *ZDomainPreimage(ZPolyhedron *A, Matrix *Func) {

  ZPolyhedron *Result = NULL, *temp;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZDOMAINPREIMAGE\n");
  fclose(fp);
#endif

  if (A->Lat->NbRows != Func->NbRows) {
    fprintf(stderr, "\nError : In ZDomainPreimage, ");
    fprintf(stderr, "Incompatible dimensions of ZPolyhedron ");
    fprintf(stderr, "and the Function \n");
    return (EmptyZPolyhedron(Func->NbColumns - 1));
  }
  for (temp = A; temp != NULL; temp = temp->next) {
    ZPolyhedron *Zpol;
    Zpol = ZPolyhedronPreimage(temp, Func);
    Result = AddZPolytoZDomain(Zpol, Result);
    ZPolyhedron_Free(Zpol);
  }
  if (Result == NULL)
    return (EmptyZPolyhedron(Func->NbColumns - 1));
  return Result;
} /* ZDomainPreimage */

/*
 * Return the Z-polyhedron intersection of the Z-polyhedra 'A' and 'B'.
 * We are based on the intersection of the two lattices of the polyhedra, named LInter.
 * If LInter is empty(null), we return the empty Zpolyhedron.
 * Otherwise, we calculate the intersection of the polyhedra on A and B, called PInter.
 * We calculate the Preimage of PInter by LInter and finally we allocate the result,
 * a Zpolyhedron in Gautam Canonical form. 
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

  LInter = NewLatticeIntersection(A->Lat, B->Lat);
  if (isEmptyLattice(LInter)) {
    Matrix_Free(LInter);
    return (EmptyZPolyhedron(A->Lat->NbRows - 1));
  }
  ImageA = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
  ImageB = DomainImage(B->P, B->Lat, MAXNOOFRAYS);
  PInter = DomainIntersection(ImageA, ImageB, MAXNOOFRAYS);
  if (emptyQ(PInter))
    Result = EmptyZPolyhedron(LInter->NbRows - 1);
  else {
    PreImage = DomainPreimage(PInter, LInter, MAXNOOFRAYS);
    Result = ZPolyhedron_Alloc(LInter, PreImage);
    Domain_Free(PreImage);
  }
  Matrix_Free(LInter);
  Domain_Free(PInter);
  Domain_Free(ImageA);
  Domain_Free(ImageB);
  return Result;
} /* ZPolyhedronIntersection */

/*
 * Return the difference of two Z-polyhedra A and B using the method Gautam
 * describes in his thesis.
*/
ZPolyhedron *ZPolyhedronDifferenceGautam(ZPolyhedron* A, ZPolyhedron* B){
  ZPolyhedron *Z1=NULL,*Z2,*Result =NULL, *Ztmp;
  LatticeUnion *LatDiff,*tmp;
  Polyhedron *ap=NULL, *imA, *imB, *PolyInter, *PolyDiff, *temp;
  Lattice *LatInter;

  LatDiff= LatticeDifference(A->Lat,B->Lat); //can simplify here
  LatDiff=LatticeSimplify(LatDiff);

  LatInter= LatticeIntersection(A->Lat,B->Lat);
  
  imA= DomainImage(A->P,A->Lat,MAXNOOFRAYS);
  imB= DomainImage(B->P,B->Lat,MAXNOOFRAYS);

  PolyInter=DomainIntersection(imA,imB,MAXNOOFRAYS);
  PolyDiff=DomainDifference(imA,imB,MAXNOOFRAYS);

  Polyhedron_Free(imA);
  Polyhedron_Free(imB);

  for(tmp=LatDiff;tmp!=NULL;tmp=tmp->next){
    Ztmp=malloc(sizeof(*Ztmp));
    Ztmp->Lat=tmp->M;
    Ztmp->P=Polyhedron_Preimage(PolyInter,tmp->M,MAXNOOFRAYS);
    Ztmp->next=Z1;
    Z1=Ztmp;
  }

  ap=DomainPreimage(PolyDiff,LatInter,MAXNOOFRAYS);

  Polyhedron_Free(PolyInter);
  Polyhedron_Free(PolyDiff);
  
  if(LatInter == NULL)
    return Z1;

  Z2=malloc(sizeof(*Z2));

  Z2->Lat=LatInter;
  Z2->P=ap;
  Z2->next=Z1;

  return Z2;
}/*ZPolyhedronDifferenceGautam*/


/*
 * Return the difference of the two Z-polyhedra 'A' and 'B'. Below is the
 * procedure to find the difference of 'A' and 'B' :-
 * Procedure:
 *     Let A = L1 (intersect) P1' and B = L2 (intersect) P2' where
 *     (P1' = DomImage(P1,L1) and P2' = DomImage(P2,L2)). Then
 *     A-B = L1 (intersect) (P1'-P2') Union
 *           (L1-L2) (intersect) (P1' (intersect) P2')
 */
static ZPolyhedron *ZPolyhedronDifference(ZPolyhedron *A, ZPolyhedron *B) {

  ZPolyhedron *Result = NULL;
  LatticeUnion *LatDiff, *temp;
  Polyhedron *DomDiff, *DomInter, *PreImage, *ImageA, *ImageB;
  Bool flag = False;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZPOLYHEDRONDIFFERENCE\n");
  fclose(fp);
#endif
  //basic tests
  if (isEmptyZPolyhedron(A))
    return NULL;
  if (isEmptyZPolyhedron(B)) {
    Result = ZDomain_Copy(A);
    return Result;
  }


  //placing A and B on another domain using their lattice
  ImageA = DomainImage(A->P, (Matrix *)A->Lat, MAXNOOFRAYS);
  ImageB = DomainImage(B->P, (Matrix *)B->Lat, MAXNOOFRAYS);
  DomDiff = DomainDifference(ImageA, ImageB, MAXNOOFRAYS);
  if (emptyQ(DomDiff))
    flag = True;
  else {
    ZPolyhedron *Z;
    PreImage = DomainPreimage(DomDiff, A->Lat, MAXNOOFRAYS);
    Z = ZPolyhedron_Alloc(A->Lat, PreImage);
    Result = AddZPolytoZDomain(Z, Result);
  }
  if (flag == True) /* DomDiff = NULL; DomInter = A */
    DomInter = Domain_Copy(ImageA);
  else {
    DomInter = DomainIntersection(ImageA, ImageB, MAXNOOFRAYS);
    if (emptyQ(DomInter)) {
      if (flag == True)
        return (EmptyZPolyhedron(A->Lat->NbRows - 1));
      else
        return Result;
    }
  }
  LatDiff = LatticeDifference(A->Lat, B->Lat);
  if (LatDiff == NULL)
    if (flag == True)
      return (EmptyZPolyhedron(A->Lat->NbRows - 1));

  while (LatDiff != NULL) {
    ZPolyhedron *tempZ = NULL;

    PreImage = DomainPreimage(DomInter, LatDiff->M, MAXNOOFRAYS);
    tempZ = ZPolyhedron_Alloc(LatDiff->M, PreImage);
    Domain_Free(PreImage);
    Result = AddZPoly2ZDomain(tempZ, Result);
    ZPolyhedron_Free(tempZ);
    temp = LatDiff;
    LatDiff = LatDiff->next;
    Matrix_Free((Matrix *)temp->M);
    free(temp);
  }
  Domain_Free(DomInter);
  Domain_Free(DomDiff);
  return Result;
} /* ZPolyhedronDifference */

/*
 * Return the image of the Z-polyhedron 'ZPol' under the invertible, affine,
 * rational transformation function 'Func'. The matrix representing the funct-
 * ion must be non-singular and the number of rows of the function must be
 * equal to the number of rows in the matrix representing the lattice of 'ZPol'
 * Algorithm:
 *         1)  Let ZPol = L (intersect) Q
 *         2)  L1 = LatticeImage(L,F)
 *         3)  Q1 = DomainImage(Q,F)
 *         4)  Z1 = L1(Inverse(L1)*Q1)
 *         5)  Return Z1
 */
static ZPolyhedron *ZPolyhedronImage(ZPolyhedron *ZPol, Matrix *Func) {

  ZPolyhedron *Result = NULL;
  Matrix *LatIm;
  Polyhedron *Pol, *PolImage;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZPOLYHEDRONIMAGE\n");
  fclose(fp);
#endif

  if ((Func->NbRows != ZPol->Lat->NbRows) ||
      (Func->NbColumns != ZPol->Lat->NbColumns)) {
    fprintf(stderr, "In ZPolImage - The Function, is not compatible with the "
                    "ZPolyhedron\n");
    return NULL;
  }
  LatIm = LatticeImage(ZPol->Lat, Func);
  if (isEmptyLattice(LatIm)) {
    Matrix_Free(LatIm);
    return NULL;
  }
  Pol = DomainImage(ZPol->P, ZPol->Lat, MAXNOOFRAYS);
  PolImage = DomainImage(Pol, Func, MAXNOOFRAYS);
  Domain_Free(Pol);
  if (emptyQ(PolImage)) {
    Matrix_Free(LatIm);
    Domain_Free(PolImage);
    return NULL;
  }
  Pol = DomainPreimage(PolImage, LatIm, MAXNOOFRAYS);
  Result = ZPolyhedron_Alloc(LatIm, Pol);
  Domain_Free(Pol);
  Domain_Free(PolImage);
  Matrix_Free(LatIm);
  return Result;
} /* ZPolyhedronImage */

/*
 * Return the preimage of the Z-polyhedron 'Zpol' under an affine transformati-
 * on function 'G'. The number of rows of matrix representing the function 'G',
 * must be equal to the number of rows of the matrix representing the lattice
 * of Z1.
 * Algorithm:
 *            1) Let Zpol = L (intersect) Q
 *            2) L1 =LatticePreimage(L,F);
 *            3) Q1 = DomainPreimage(Q,F);
 *            4) Z1 = L1(Inverse(L1)*Q1);
 *            5) Return Z1
 */
static ZPolyhedron *ZPolyhedronPreimage(ZPolyhedron *Zpol, Matrix *G) {

  Lattice *Latpreim;
  Polyhedron *Qprime, *Q, *Polpreim;
  ZPolyhedron *Result;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered ZPOLYHEDRONPREIMAGE\n");
  fclose(fp);
#endif

  if (G->NbRows != Zpol->Lat->NbRows) {
    fprintf(stderr, "\nIn ZPolyhedronPreimage: Error, The dimensions of the ");
    fprintf(stderr, "function are not compatible with that of the Zpolyhedron");
    return EmptyZPolyhedron(G->NbColumns - 1);
  }
  Q = DomainImage(Zpol->P, Zpol->Lat, MAXNOOFRAYS);
  Polpreim = DomainPreimage(Q, G, MAXNOOFRAYS);
  if (emptyQ(Polpreim))
    Result = NULL;
  else {
    Latpreim = LatticePreimage(Zpol->Lat, G);
    if (isEmptyLattice(Latpreim))
      Result = NULL;
    else {
      Qprime = DomainPreimage(Polpreim, Latpreim, MAXNOOFRAYS);
      Result = ZPolyhedron_Alloc(Latpreim, Qprime);
      Domain_Free(Qprime);
    }
    Matrix_Free(Latpreim);
  }
  Domain_Free(Q);
  return Result;
} /* ZPolyhedronPreimage */

/*
 * Return the Z-polyhderon 'Zpol' in canonical form: 'Result' (for the Z-poly-
 * hedron in canonical form) and Basis 'Basis' (for the basis with respect to
 * which 'Result' is in canonical form.
 */
void CanonicalForm(ZPolyhedron *Zpol, ZPolyhedron **Result, Matrix **Basis) {

  Matrix *B1 = NULL, *B2 = NULL, *T1, *B2inv;
  int i, l1, l2;
  Value tmp;
  Polyhedron *Image, *ImageP;
  Matrix *H, *U, *temp, *Hprime, *Uprime, *T2;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered CANONICALFORM\n");
  fclose(fp);
#endif

  if (isEmptyZPolyhedron(Zpol)) {
    Basis[0] = Identity(Zpol->Lat->NbRows);
    Result[0] = ZDomain_Copy(Zpol);
    return;
  }
  value_init(tmp);
  l1 = FindHermiteBasisofDomain(Zpol->P, &B1);
  Image = DomainImage(Zpol->P, (Matrix *)Zpol->Lat, MAXNOOFRAYS);
  l2 = FindHermiteBasisofDomain(Image, &B2);

  if (l1 != l2)
    fprintf(stderr, "In CNF : Something wrong with the Input Zpolyhedra \n");

  B2inv = Matrix_Alloc(B2->NbRows, B2->NbColumns);
  temp = Matrix_Copy(B2);
  Matrix_Inverse(temp, B2inv);
  Matrix_Free(temp);

  temp = Matrix_Alloc(B2inv->NbRows, Zpol->Lat->NbColumns);
  T1 = Matrix_Alloc(temp->NbRows, B1->NbColumns);
  Matrix_Product(B2inv, (Matrix *)Zpol->Lat, temp);
  Matrix_Product(temp, B1, T1);
  Matrix_Free(temp);

  T2 = ChangeLatticeDimension(T1, l1);
  temp = ChangeLatticeDimension(T2, T2->NbRows + 1);

  /* Adding the affine part */
  for (i = 0; i < l1; i++)
    value_assign(temp->p[i][temp->NbColumns - 1], T1->p[i][T1->NbColumns - 1]);

  AffineHermite(temp, &H, &U);
  Hprime = ChangeLatticeDimension(H, Zpol->Lat->NbRows);

  /* Exchanging the Affine part */
  for (i = 0; i < l1; i++) {
    value_assign(tmp, Hprime->p[i][Hprime->NbColumns - 1]);
    value_assign(Hprime->p[i][Hprime->NbColumns - 1],
                 Hprime->p[i][H->NbColumns - 1]);
    value_assign(Hprime->p[i][H->NbColumns - 1], tmp);
  }
  Uprime = ChangeLatticeDimension(U, Zpol->Lat->NbRows);

  /* Exchanging the Affine part */
  for (i = 0; i < l1; i++) {
    value_assign(tmp, Uprime->p[i][Uprime->NbColumns - 1]);
    value_assign(Uprime->p[i][Uprime->NbColumns - 1],
                 Uprime->p[i][U->NbColumns - 1]);
    value_assign(Uprime->p[i][U->NbColumns - 1], tmp);
  }
  Polyhedron_Free(Image);
  Matrix_Free(B2inv);
  B2inv = Matrix_Alloc(B1->NbRows, B1->NbColumns);
  Matrix_Inverse(B1, B2inv);
  ImageP = DomainImage(Zpol->P, B2inv, MAXNOOFRAYS);
  Matrix_Free(B2inv);
  Image = DomainImage(ImageP, Uprime, MAXNOOFRAYS);
  Domain_Free(ImageP);
  Result[0] = ZPolyhedron_Alloc(Hprime, Image);
  Basis[0] = Matrix_Copy(B2);

  /* Free the variables */
  Polyhedron_Free(Image);
  Matrix_Free(B1);
  Matrix_Free(B2);
  Matrix_Free(temp);
  Matrix_Free(T1);
  Matrix_Free(T2);
  Matrix_Free(H);
  Matrix_Free(U);
  Matrix_Free(Hprime);
  Matrix_Free(Uprime);
  value_clear(tmp);
  return;
} /* CanonicalForm */

/*
 * Given a Z-polyhedron 'A' in which the Lattice is not integral, return the
 * Z-polyhedron which contains all the integral points in the input lattice.
 */
ZPolyhedron *IntegraliseLattice(ZPolyhedron *A) {

  ZPolyhedron *Result;
  Lattice *M = NULL, *Id;
  Polyhedron *Im = NULL, *Preim = NULL;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp, "\nEntered INTEGRALISELATTICE\n");
  fclose(fp);
#endif

  Im = DomainImage(A->P, A->Lat, MAXNOOFRAYS);
  Id = Identity(A->Lat->NbRows);
  M = LatticeImage(Id, A->Lat);
  if (isEmptyLattice(M))
    Result = EmptyZPolyhedron(A->Lat->NbRows - 1);
  else {
    Preim = DomainPreimage(Im, M, MAXNOOFRAYS);
    Result = ZPolyhedron_Alloc(M, Preim);
  }
  Matrix_Free(M);
  Domain_Free(Im);
  Domain_Free(Preim);
  return Result;
} /* IntegraliseLattice */

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
      Zpol = ZPolyhedron_Alloc(L->M, Preim);
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

ZPolyhedron *SplitZpolyhedron(ZPolyhedron *ZPol, Lattice *B) {

  Lattice *Intersection = NULL;
  Lattice *B1 = NULL, *B2 = NULL, *newB1 = NULL, *newB2 = NULL;
  Matrix *U = NULL, *M1 = NULL, *M2 = NULL, *M1Inverse = NULL,
         *MtProduct = NULL;
  Matrix *Vinv, *V, *temp, *DiagMatrix;
  Matrix *H, *U1, *X, *Y;
  ZPolyhedron *zpnew, *Result;
  LatticeUnion *Head = NULL, *tempHead = NULL;
  int i;
  Value k;

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

  if (isinHnf(ZPol->Lat) != True) {
    AffineHermite(ZPol->Lat, &H, &U1);
    X = Matrix_Copy(H);
    Matrix_Free(U1);
    Matrix_Free(H);
  } else
    X = Matrix_Copy(ZPol->Lat);

  if (isinHnf(B) != True) {
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
    return ZPol;
  }

  Result = NULL;

  if (Head)
    while (Head) {
      tempHead = Head;
      Head = Head->next;
      zpnew = ZPolyhedron_Alloc(tempHead->M, ZPol->P);
      Result = AddZPoly2ZDomain(zpnew, Result);
      ZPolyhedron_Free(zpnew);
      tempHead->next = NULL;
      free(tempHead);
    }

  return Result;
}

/*
 * Moves the constant part (last line and last row) as first line and row
 * of the matrix.
 * This is useful to perform the HNF and keeping the affine part as top-left
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
}

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
}
/*
* The function takes a polyhedron and modifies it in place
* to be in canonical form as described by Gautam
* (A->Lat in HNF and no equalities in A->P)
*/
void CanonicalZPGautam(ZPolyhedron* A) {
  
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

  // if(emptyQ(A->P)) {
  //   //before return construct zpol empty with dim p=0 if dim a.p is not 0 we need to reallocate it and with a lattice = (0..0 1) do these manually
  //   return;
  // }


    if(emptyQ(A->P)){
      #ifdef CANONICAL_DEBUG
        printf("The polyhedron is empty\n");
      #endif
      Polyhedron* tmp = Empty_Polyhedron(A->P->Dimension);
      Polyhedron_Free(A->P);
      A->P=tmp;
      if(A->Lat->NbColumns>1){
        Matrix_Free(A->Lat);
        Lattice* L=Matrix_Alloc(A->P->Dimension+1,1);
        A->Lat=L;
      }
      for (int i = 0; i < A->Lat->NbRows-1; i++) {
        value_set_si(A->Lat->p[i][0],0);
      }
      value_set_si(A->Lat->p[A->Lat->NbRows-1][0],1);

      #ifdef CANONICAL_DEBUG
        printf("We return:\n");
        ZPolyhedronPrint(stderr,P_VALUE_FMT,A);
      #endif

      return;
    }
      
  

  // Check if P contains equalities
  #ifdef CANONICAL_DEBUG
    printf("Checking for equalites in P\n");
  #endif

  if (A->P->Dimension > 0 && A->P->NbEq != 0) {

    // remove equalities in P and change Lat to spread the original space

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "P has equalities\n");
      fprintf(stderr, "Matrix of P: ");
      Matrix_Print(stderr,P_VALUE_FMT, A->Lat);
    #endif

    // Eq is the matrix of linear equations of P (including the constant)
    Matrix* Eq = Matrix_Alloc(A->P->NbEq, A->P->Dimension+1);
    if(!Eq){
      errormsg1("CanonicalZPGautam", "outofmem", "Not enough memory space!");
      return;
    }

    int eqnum = 0;
    // get equalities
    for(int i=0; i<A->P->NbConstraints; i++) {
      if(A->P->Constraint[i][0] == 0) // Equality
      {
        for(int j=0; j<Eq->NbColumns; j++){
          value_assign(Eq->p[eqnum][j], A->P->Constraint[i][j+1]);
        }
        if(++eqnum >= A->P->NbEq)
          break;
      }
    }

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "Equality matrix (including constants): ");
      Matrix_Print(stderr, P_VALUE_FMT, Eq);
    #endif

    // compute Ker(Eq)
    Matrix* ker;
    ker = int_ker(Eq);
    Matrix_Free(Eq);

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "ker of eq: ");
      Matrix_Print(stderr, P_VALUE_FMT, ker);
    #endif
    
    Matrix *U, *V; 
    
    Matrix* T;
    Matrix_Move_Homogeneous_Dim_First(ker);
    left_hermite(ker,&T,&V,&U);
    Matrix_Move_Homogeneous_Dim_Last(T);
    Matrix_Free(ker);
    Matrix_Free(V);
    Matrix_Free(U);
    if(!T){
      errormsg1("CanonicalZPGautam", "outofmem", "Not enough memory space!");
      return;
    }

    #ifdef CANONICAL_DEBUG
    fprintf(stderr, "Matrix T: ");
    Matrix_Print(stderr, P_VALUE_FMT, T);
    fprintf(stderr, "lat of a: ");
    Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
    #endif

    // NewL = L . T
    Matrix* NewL = Matrix_Alloc(A->Lat->NbRows, T->NbColumns);
    if(!NewL) {
      errormsg1("CanonicalZPGautam", "outofmem", "Not enough memory space!");
      return;
    }
    Matrix_Product(A->Lat, T, NewL);

    Polyhedron* NewP = Polyhedron_Preimage(A->P, T, MAXNOOFRAYS);
    Polyhedron_Free(A->P);
    A->P = NewP;

    // update Lat
    Matrix_Free(A->Lat);
    A->Lat = NewL;
  
    // freeing 
    Matrix_Free(T);

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "----------- New Lat: ");
      Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
      fprintf(stderr, "----------- New P: ");
      Polyhedron_Print(stderr, P_VALUE_FMT, A->P);
      fprintf(stderr, "-----------\n");
    #endif

  } // P contains no equalities
  else {
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "P has no equalities.\n");
    #endif
  }


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


  // check if A->Lat is in Hermite form
  if(!isinHnf(A->Lat)) {
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "A is not HNF\n");
    #endif
    Matrix* U = NULL;
    Matrix* H = NULL;

    // // temporarily remove the homogeneous part of the matrix (dirty)
    // A->Lat->NbColumns--;
    // A->Lat->NbRows--;

    // for left hermite to include the constant:
    Matrix_Move_Homogeneous_Dim_First(A->Lat);
    // to compute HNF of the lattice (constant part left-top)
    // We will use U = Q^{-1}, such that LU = H.
    left_hermite(A->Lat, &H, NULL, &U);
    Matrix_Free(A->Lat);

    // Move the constant back to right-bottom
    Matrix_Move_Homogeneous_Dim_Last(H);
    Matrix_Move_Homogeneous_Dim_Last(U);

    // set the new Lat matrix as H
    A->Lat = H;

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "New Lat (HNF): ");
      Matrix_Print(stderr, P_VALUE_FMT, A->Lat);
    #endif

    // Now update of A->P using the premimage by U.
    Polyhedron *NewP = Polyhedron_Preimage(A->P, U, MAXNOOFRAYS);
    Polyhedron_Free(A->P);
    A->P = NewP;
    Matrix_Free(U);

    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "New P: ");
      Polyhedron_Print(stderr, P_VALUE_FMT, A->P);
    #endif
  } // A->Lat in Hermite form
  else {
    #ifdef CANONICAL_DEBUG
      fprintf(stderr, "A is HNF.\n");
    #endif
  }

  return;
}

/*
 * Takes into parameters two lattices A and B of the form:
 *  A =   A' | a      B =   B' | b
 *      0..0 | 1          0..0 | 1
 * 
 *  Copies them in a matrix (Tmp), used to calculate the left hermite,
 *  the Lattice of size A->Nbrows x min(A->NbCols, B->NbCols) is the
 *  intersection of A and B.
 * 
 */
Lattice* NewLatticeIntersection(Lattice* A, Lattice* B) {
  Lattice *Tmp, *H, *Res;
  if(A->NbRows != B->NbRows){
    errormsg1("NewLatticeIntersection", "dimincomp", "incompatible dimensions!");
    return NULL;
  }
  #ifdef NEWINTERSECTION_DEBUG
  fprintf(stderr,"Matrix A:\n");
  Matrix_Print(stderr, P_VALUE_FMT, A);
  fprintf(stderr,"Matrix B:\n");
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
    errormsg1("NewLatticeIntersection", "outofmem", "Not enough memory space!");
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
  
  #ifdef NEWINTERSECTION_DEBUG
    fprintf(stderr,"\n Tmp init:\n");
    Matrix_Print(stderr,P_VALUE_FMT, Tmp);
  #endif
  

  // left_hermite of the TMP
  // H is the matrix that contains the solution. it is of the form:
  // 
  // H =   D  |   0          D is a square matrix
  //     ------------
  //       X  | 1 0.0
  //          | r  R
  // 
  // with  R    r
  //      0..0  1   being our result
  // if the number above r is not 1 then the intersection is not integer
  // (no solution to the intersection)

  left_hermite(Tmp, &H, NULL, NULL);


  #ifdef NEWINTERSECTION_DEBUG
    fprintf(stderr,"\nH:\n");
    Matrix_Print(stderr,P_VALUE_FMT,H);
  #endif
  Matrix_Free(Tmp);

  // recuperating the result. if the top-left value of R is not 1 then we have an empty solution.

  if(value_notone_p(H->p[A->NbRows][A->NbRows])) {
    #ifdef NEWINTERSECTION_DEBUG
      fprintf(stderr,"\n Empty intersection\n");
    #endif
    Matrix_Free(H);
    return NULL;
  }

  int nbcol = (A->NbColumns<B->NbColumns)?A->NbColumns:B->NbColumns;

  Res=Matrix_Alloc(A->NbRows, nbcol);
  if (!Res) {
    errormsg1("NewLatticeIntersection", "outofmem", "Not enough memory space!");
    return NULL;
  }
  

  for (int i = 0; i < A->NbRows; i++) {
    for (int j = 0; j < nbcol; j++) {
        value_assign(Res->p[i][j], H->p[i + A->NbRows][j + A->NbRows]);
    }
  }
  Matrix_Free(H);
  // put Res in the proper affine form (didn't want to write a loop, had a pre-written function.)
  Matrix_Move_Homogeneous_Dim_Last(Res);

  #ifdef NEWINTERSECTION_DEBUG
    fprintf(stderr, "\n NewLaticceIntersection result: ");
    Matrix_Print(stderr, P_VALUE_FMT, Res);
  #endif

  return Res;
}
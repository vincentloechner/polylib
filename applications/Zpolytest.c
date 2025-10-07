/* zpolytest.c
This is a testbench for the Zpolylib (part of polylib manipulating 
Z-polyhedra. */

#include <stdio.h>
#include <polylib/polylib.h>

#define WS 0

char s[128];

int main() {
  
  Matrix *a=NULL, *b=NULL, *c=NULL, *d=NULL, *e=NULL, *f=NULL, *g;
  LatticeUnion *l1, *l2;
  Polyhedron *A=NULL, *B=NULL, *C=NULL, *D = NULL;
  LBL *ZA=NULL, *ZB=NULL, *ZC=NULL, *ZD=NULL;
  int  nbPol, nbMat, func;

  // The structure of the input file to this program is the following:

  // A line starting with a `#' is considered as a comment.

  // - First a line containing:
  //     M nbMat
  //   Where nbMat is an integer indicating how many Matrices will be given in
  //   the following (max 3). Next the matrices are given in PolyLib format.
  //   For each matrix, the first row is two integers:
  //     nbRows nbColumns
  //   Then the matrix is given row by row.

  // - Then a line containing:
  //     D nbDomain
  //   where nbDomain is an integer indicating how many domains will be given
  //   in the following (max 3). Domains are given in PolyLib format:
  //     the first row is two integers:
  //       nbConstraints dimension
  //     then the constraints line by line in Polylib format.

  // - The last line of the input file contains:
  //     F numTest
  //   which indicates which test will be performed.
  
  // All lines below are ignored.
  
  nbPol = nbMat = 0;
  // read matrices
  do {
    fgets(s, 128, stdin);
  }
  while ((*s=='#') ||
	  ((sscanf(s, "D %d", &nbPol)<1) && (sscanf(s, "M %d", &nbMat)<1)));

  switch (nbMat) {

  case 1:
    a = Matrix_Read();
    break;
  
  case 2: 
    a = Matrix_Read();
    b = Matrix_Read();
    break;
  
  case 3:
    a = Matrix_Read();
    b = Matrix_Read();
    c = Matrix_Read();
    break;
  }

  // read polyhedra
  fgets(s, 128, stdin);
  while ((*s=='#') ||
	 ((sscanf(s, "D %d", &nbPol)<1) && (sscanf(s, "M %d", &nbMat)<1)) )
    fgets(s, 128, stdin);
  
  switch (nbPol) { 
  
  case 1:  
    g = Matrix_Read();
    A = Constraints2Polyhedron(g, WS);
    Matrix_Free(g);
    break;
  
  case 2:         
    g = Matrix_Read();
    A = Constraints2Polyhedron(g, WS);
    Matrix_Free(g);
    g = Matrix_Read();
    B = Constraints2Polyhedron(g, WS);
    Matrix_Free(g);
    break;
  
  case 3:
    g = Matrix_Read();
    A = Constraints2Polyhedron(g, WS);
    Matrix_Free(g);
    g = Matrix_Read();
    B = Constraints2Polyhedron(g, WS);
    Matrix_Free(g);
    g = Matrix_Read();
    C = Constraints2Polyhedron(g, WS);
    Matrix_Free(g);
    break;
  }

  // read function
  do {
    fgets(s, 128, stdin);
  }
  while((*s=='#') || (sscanf(s, "F %d", &func)<1));

  switch (func) {
  case 1:
    
    /* just a test of polylib functions */
    C = DomainUnion(A, B, 200);
    D = DomainConvex(C, 200);
    d = Polyhedron2Constraints(D);
    Matrix_Print(stdout, P_VALUE_FMT, d);
    break;
    
  case 2: /* AffineHermite */
    if(isNormalLattice(a))
      printf("input matrix is normal\n");
    else
      printf("input matrix is not normal\n");
    AffineHermite(a,&b,&c);
    Matrix_Print(stdout, P_VALUE_FMT, b);
    Matrix_Print(stdout, P_VALUE_FMT, c);
    break;
    
  case 3: /* LatticeIntersection */
    
    c = LatticeIntersection(a,b);
    Matrix_Print(stdout, P_VALUE_FMT, c);
    break;
    
  case 4: /* LatticeIncluded */
    AffineHermite(a, &d, NULL);
    AffineHermite(b, &e, NULL);
    AffineHermite(c, &f, NULL);
    printf(" 2 in 1: %d\n", LatticeIncluded(e, d));
    printf(" 1 in 3: %d\n", LatticeIncluded(f, d));
    printf(" 1 in 2: %d\n", LatticeIncluded(d, e));
    break;
  
  case 5: /* LatticeDifference */
    
    l1 = LatticeDifference(a, b);
    l2 = LatticeDifference(b, a);
    printf("L1 - L2:\n");
    PrintLatticeUnion(stdout, P_VALUE_FMT, l1);
    printf("L2 - L1:\n");
    PrintLatticeUnion(stdout, P_VALUE_FMT, l2);
    LatticeUnion_Free(l1);
    LatticeUnion_Free(l2);
    break;
    
  case 6: /* isEmptyLBL */

    ZA = LBLAlloc(a,A);
    printf("is Empty? :%d\n", isEmptyLBL(ZA));
    break;
    
  case 7: /* LBLIntersection */
        
    ZA = LBLAlloc(a, A);
    ZB = LBLAlloc(b, B);
    ZC = LBLIntersection(ZA, ZB);
    LBLPrint(stdout, P_VALUE_FMT, ZC);
    break;
    
  case 8: /* LBLUnion */
    
    ZA = LBLAlloc(a, A);
    ZB = LBLAlloc(b, B);
    ZC = LBLUnion(ZA,ZB);
    LBLPrint(stdout, P_VALUE_FMT, ZC);
    break;
    
  case 9: /* LBLDifference */
    
    ZA = LBLAlloc(a, A);
    ZB = LBLAlloc(b, B);
    // LBLPrint(stdout, P_VALUE_FMT, ZA);
    // LBLPrint(stdout, P_VALUE_FMT, ZB);
    ZC = LBLDifference(ZA, ZB);
    ZD = LBLDifference(ZB, ZA);
    LBLSimplify(ZC);
    LBLSimplify(ZD);
    printf("A - B = ");
    LBLPrint(stdout, P_VALUE_FMT, ZC);
    printf("\n\nB - A = ");
    LBLPrint(stdout, P_VALUE_FMT, ZD);
    break;
    
  case 10: /* LBLImage */
    
    ZA = LBLAlloc(a, A);
    ZC = LBLImage(ZA, b); 
    LBLPrint(stdout, P_VALUE_FMT, ZC);
    break;
    
  case 11: /* LBLPreimage */
    
    ZA = LBLAlloc(a, A);
    ZC = LBLPreimage(ZA, b); 
    LBLPrint(stdout, P_VALUE_FMT, ZC);
    break;
    
  case 12: /* difference between image of preimage and original*/

    ZA = LBLAlloc(a, A);
    ZC = LBLPreimage(ZA, b);
    ZD = LBLImage(ZC, b);

    printf("PreIm(A) = ");
    LBLPrint(stdout, P_VALUE_FMT, ZC);
    printf("Im(PreIm(A)) = ");
    LBLPrint(stdout, P_VALUE_FMT, ZC);
    // ZD should be included in ZA
    printf("The image of the preimage is included in the original LBL ");
    printf("(should always be true)? %d\n", LBLIncluded(ZD, ZA));
    ZB = LBLDifference(ZA, ZD);
    printf("The image of the preimage is exactly the original LBL? %d\n",
	    isEmptyLBL(ZB));
    break;
  
  case 13:  /* LBLSimplify */
    
    ZA = LBLAlloc(a, A);
    printf("A = ");
    LBLPrint(stdout, P_VALUE_FMT, ZA);
    LBLSimplify(ZA);
    printf("LBLSimplify(A) = ");
    LBLPrint(stdout, P_VALUE_FMT, ZA);
    break;
    
  case 14:  /* EmptyLBL */
        
    ZA = EmptyLBL(3);
    printf("is Empty? :%d\n", isEmptyLBL(ZA));
    break;
    
  case 15:  /* LBLIncluded */
  
    ZA = LBLAlloc(a, A);
    ZB = LBLAlloc(b, B);
    printf("A in B  :%d\nB in A  :%d\n", 
    LBLIncluded(ZA, ZB),
    LBLIncluded(ZB, ZA));
    break;
  
  case 16:  /* LBLComplement */
  
    ZA = LBLAlloc(a, A);
    ZB = LBLComplement(ZA);
    printf("A = ");
    LBLPrint(stdout, P_VALUE_FMT, ZA);
    printf("\nComplement(A) = ");
    LBLPrint(stdout, P_VALUE_FMT, ZB);
    ZC = LBLComplement(ZB);
    printf("\nZB AFTER COMPLEMENT = ");
    LBLPrint(stdout, P_VALUE_FMT, ZB);
    printf("\nComplement(Complement(A)) = ");
    LBLPrint(stdout, P_VALUE_FMT, ZC);
    printf("\nZB inter ZC = ");
    LBLPrint(stdout, P_VALUE_FMT, LBLIntersection(ZB, ZC));
    printf("\nIs the the complement of the complement the original?\n");
    printf("  A is included in C(C(A)): %d\n", LBLIncluded(ZA, ZC));
    printf("  C(C(A)) is included in A: %d\n", LBLIncluded(ZC, ZA));
    {
      LBL *diff1, *diff2;
      diff1 = LBLDifference(ZA, ZC);
      LBLSimplify(diff1);
      printf("ZA - ZC = ");
      LBLPrint(stdout, P_VALUE_FMT, diff1);
      LBLFree(diff1);
      diff2 = LBLDifference(ZC, ZA);
      LBLSimplify(diff2);
      printf("ZC - ZA = ");
      LBLPrint(stdout, P_VALUE_FMT, diff2);
      LBLFree(diff2);
    }
    break;

  case 17:  /* LBL2Zdomain */
    ZA = LBLAlloc(a, A);
    printf("A = ");
    LBLPrint(stdout, P_VALUE_FMT, ZA);
    ZB = LBL2ZDomain(ZA);
    printf("LBL2ZDomain(A) = ");
    LBLPrint(stdout, P_VALUE_FMT, ZB);

    // check equality between ZA and ZB
    ZC = LBLDifference(ZA, ZB);
    LBLSimplify(ZC);
    printf("A - ZD(A) = ");
    LBLPrint(stdout, P_VALUE_FMT, ZC);
    ZD = LBLDifference(ZB, ZA);
    LBLSimplify(ZD);
    printf("ZD(A) - A = ");
    LBLPrint(stdout, P_VALUE_FMT, ZD);
    break;


  case 19:  /* check that complement(A) inter A = empty */
            /* and complement(A) union A = universe */
    ZA = LBLAlloc(a, A);
    ZB = LBLComplement(ZA);
    printf("A = ");
    LBLPrint(stdout, P_VALUE_FMT, ZA);
    printf("\nComplement(A) = ");
    LBLPrint(stdout, P_VALUE_FMT, ZB);

    ZC = LBLIntersection(ZA, ZB);
    LBLSimplify(ZC);
    printf("\nA inter Complement(A) (should be empty) = ");
    LBLPrint(stdout, P_VALUE_FMT, ZC);
    LBLFree(ZC);

    ZC = UniverseLBL(ZA->Lat->NbRows - 1);
    ZD = LBLUnion(ZA, ZB);
    // printf("\nA union Complement(A) = ");
    // LBLPrint(stdout, P_VALUE_FMT, ZD);
    LBLFree(ZB);
    ZB = LBLDifference(ZC, ZD);
    LBLSimplify(ZB);
    printf("\nUniverse - (A union Complement(A)) (should be empty) = ");
    LBLPrint(stdout, P_VALUE_FMT, ZB);
    break;

  // no longer in use tests:

  // case 21: /* AffineSmith */
  
  //   AffineSmith(a,&b,&c, &d);
  //   printf("A = U . Diag . V\n");
  //   Matrix_Print(stdout, P_VALUE_FMT, b);
  //   Matrix_Print(stdout, P_VALUE_FMT, d);
  //   Matrix_Print(stdout, P_VALUE_FMT, c);
  //   break;
  
  // case 22: /* SolveDiophantine */

  //   rank=SolveDiophantine(a,&d,&v);
  //   Matrix_Print(stdout, P_VALUE_FMT, a);
  //   printf( "rank: %d\n ", rank);
  //   Matrix_Print(stdout, P_VALUE_FMT, d);
  //   Vector_Print(stdout, P_VALUE_FMT, v);
  //   rank=SolveDiophantine(b, &d, &v);
  //   Matrix_Print(stdout, P_VALUE_FMT, b);
  //   printf( "rank: %d\n ", rank);
  //   Matrix_Print(stdout, P_VALUE_FMT, d);
  //   Vector_Print(stdout, P_VALUE_FMT, v);
  //   rank=SolveDiophantine(c, &d, &v);
  //   Matrix_Print(stdout, P_VALUE_FMT, c);
  //   printf( "rank: %d\n ",rank);
  //   Matrix_Print(stdout, P_VALUE_FMT, d);
  //   Vector_Print(stdout, P_VALUE_FMT, v);
  //   Vector_Free(v);
  //   break;

  
  case 100: /* just alloc and normalize */
    ZA = LBLAlloc(a,A);
    LBLPrint(stdout, P_VALUE_FMT, ZA);
    break;
    
  default:
    printf("? unknown function\n");
  }

  // free 4 allocated matrices, polyhedra, Z-polyhedra.
  if (a)    Matrix_Free(a);
  if (b)    Matrix_Free(b);
  if (c)    Matrix_Free(c);
  if (d)    Matrix_Free(d);

  if (A)    Domain_Free(A);
  if (B)    Domain_Free(B);
  if (C)    Domain_Free(C);
  if (D)    Domain_Free(D);

  if (ZA)    LBLFree(ZA);
  if (ZB)    LBLFree(ZB);
  if (ZC)    LBLFree(ZC);
  if (ZD)    LBLFree(ZD);
  
  // free all remaining cache memory of PolyLib:
  polylib_close();

  return 0;
} /* main */

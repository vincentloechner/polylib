/* zpolytest.c
This is a testbench for the Zpolylib (part of polylib manipulating 
Z-polyhedra. */

#include <stdio.h>
#include <polylib/polylib.h>

#define WS 0

char s[128];

int main() {
  
  Matrix *a=NULL, *b=NULL, *c=NULL, *d=NULL, *g;
  LatticeUnion *l1, *l2, *temp;
  Polyhedron *A=NULL, *B=NULL, *C=NULL, *D;
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
    Domain_Free(D);
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
    
  case 4: /* LatticeDifference */
        
    printf(" 2 in 1 : %d\n", LatticeIncludes(b, a));
    printf(" 1 in 3 : %d\n", LatticeIncludes(c, a));
    printf(" 1 in 2 : %d\n", LatticeIncludes(a, b));
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
    printf(
      "The image of the preimage is included in the original Z-polyhedron (should always be true)? %d\n",
	    LBLIncludes(ZD, ZA));
    ZB = LBLDifference(ZA, ZD);
    printf(
      "The image of the preimage is exactly the original Z-polyhedron? %d\n",
	    isEmptyLBL(ZB));
    break;
  
  // case 13:  /* LBLSimplify */
    
  //   ZA = LBLAlloc(a, A);
  //   ZA->next = LBLAlloc(b, B);
  //   LBLPrint(stdout, P_VALUE_FMT, ZA);
  //   ZD = LBLSimplify(ZA);
  //   LBLPrint(stdout, P_VALUE_FMT, ZD);
  //   break;
    
  case 14:  /* EmptyLBL */
        
    ZA = EmptyLBL(3);
    printf("is Empty? :%d\n", isEmptyLBL(ZA));
    break;
    
  case 15:  /* LBLInclude */
  
    ZA = LBLAlloc(a, A);
    ZB = LBLAlloc(b, B);
    printf("A in B  :%d\nB in A  :%d\n", 
	    LBLIncludes(ZA, ZB),
	    LBLIncludes(ZB, ZA));
    break;
  

  // case 18:  /* EmptyLattice */
        
  //   printf("is Empty? :%d\n", isEmptyLattice(a));
  //   printf("is Empty? :%d\n", isEmptyLattice(EmptyLattice(3)));
  //   break;
  
  // case 19:  /* CanonicalForm */
     
  //   ZA=LBLAlloc(a,A);
  //   ZB=LBLAlloc(a,B);
  //   CanonicalForm(ZA, &ZC, &c);
  //   CanonicalForm(ZB, &ZD, &d);
  //   LBLPrint(stdout, P_VALUE_FMT, ZC);
  //   LBLPrint(stdout, P_VALUE_FMT, ZD);
  //   break;
    
  // case 20: /* LatticeSimplify */
    
  //   l1=LatticeUnion_Alloc();
  //   l2=LatticeUnion_Alloc();
  //   l1->M=Matrix_Copy(a);
  //   l1->next=l2;
  //   l2->M=Matrix_Copy(b);
  //   l1=LatticeSimplify(l1);
  //   PrintLatticeUnion(stdout, P_VALUE_FMT, l1);
  //   LatticeUnion_Free(l1);
  //   break;

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

  // case 23: /* SplitLBL */
        
  //   ZA = LBLAlloc(a, A);
  //   ZC = SplitLBL(ZA, b);
  //   LBLPrint(stdout, P_VALUE_FMT, ZC);
  //   break;

  case 24: /* left_hermite */
    left_hermite(a, &b, &c, NULL);
    printf("A = H . Q\nH = ");
    Matrix_Print(stdout, P_VALUE_FMT, b);
    printf("Q = ");
    Matrix_Print(stdout, P_VALUE_FMT, c);
    break;
  
  case 25:/* move homogenous dimension */

    Matrix_Print(stdout, P_VALUE_FMT, a);
    Matrix_Move_Homogeneous_Dim_First(a);
    b = Matrix_Copy(a);
    Matrix_Move_Homogeneous_Dim_Last(b);
    Matrix_Print(stdout, P_VALUE_FMT, a);
    Matrix_Print(stdout, P_VALUE_FMT, b);
    break;

  case 28:
    b = int_ker(a);
    Matrix_Print(stdout, P_VALUE_FMT, b);
    break;

  case 100: /* just alloc and normalize */
    ZA = LBLAlloc(a,A);
    LBLPrint(stdout, P_VALUE_FMT, ZA);
    break;
    
  default:
    printf("? unknown function\n");
  }

  // free 4 allocated matrices, polyhedra, Z-polyhedra.
  if (a)
    Matrix_Free(a);
  if (b)
    Matrix_Free(b);
  if (c)
    Matrix_Free(c);
  if(d)
    Matrix_Free(d);

  if (A)
    Domain_Free(A);
  if (B)
    Domain_Free(B);
  if (C)
    Domain_Free(C);

  if (ZA)
    LBLFree(ZA);
  if (ZB)
    LBLFree(ZB);
  if (ZC)
    LBLFree(ZC);
  if (ZD)
    LBLFree(ZD);
  
  // free all memory (sanitizer):
  polylib_close();

  return 0;
} /* main */

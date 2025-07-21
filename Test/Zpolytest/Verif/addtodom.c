// gcc -I ../../include/ -L ../.. armin.c -lpolylib32

#include <polylib/polylib32.h>

int main(void) {
    Matrix *ineqs, *lattice;
    Polyhedron *p;
    ZPolyhedron *zp1, *zp2, *zp3, *res, *res2;

    ineqs = Matrix_Alloc(2, 3);
    lattice = Matrix_Alloc(2, 2);

    // 0 <= i <= 10
    ineqs->p[0][0] = 1;
    ineqs->p[0][1] = 1;
    ineqs->p[0][2] = 0;
    ineqs->p[1][0] = 1;
    ineqs->p[1][1] = -1;
    ineqs->p[1][2] = 10;

    // i  -->  2*i
    lattice->p[0][0] = 2;
    lattice->p[0][1] = 0;
    lattice->p[1][0] = 0;
    lattice->p[1][1] = 1;
    
    p = Constraints2Polyhedron(ineqs, 10000);
    
    // zp1 is {0,2,...,20}
    zp1 = ZPolyhedron_Alloc(lattice, p);
    
    // CanonicalZPGautam(zp1);
    // ZPolyhedron* zpe=EmptyZPolyhedron(1);
    
    //test for the empty polyhedrons
    // printf("\n %u \n",isEmptyZPolyhedron(zpe));
    // printf("\n %u \n",isEmptyZPolyhedron(zp1));
    
    // ZDomainPrint(stdout," %d", zpe);
    
    // // i  -->  2*i+1
    lattice->p[0][1] = 1;
    // zp2 is {1,3,...,21}
    zp2 = ZPolyhedron_Alloc(lattice, p);
    
    
    lattice->p[0][0] = 3;
    lattice->p[0][1] = 0;
    lattice->p[1][0] = 0;
    lattice->p[1][1] = 1;
    
    res = ZPolyhedron_Alloc(lattice, p);


    // Set zp3 to the union of {1,3,...,21} and {0,2,...,20}.
    // The call to ZDomainDifference below succeeds, if
    // zp3 is computed as  zp3 = ZDomainUnion(zp1,zp2); .
    zp3 = ZDomainUnion(zp2, zp1);

    printf("\n----------------- The poyhedron before the call: -----------------\n");
    ZDomainPrint(stdout, " %d",zp3);
    
    
    AddZPolytoZDomain1(res,zp3);
    
    
    
    printf("\n----------------- The poyhedron after the call: -----------------\n");
    ZDomainPrint(stdout, " %d",zp3);


    return 0;
}


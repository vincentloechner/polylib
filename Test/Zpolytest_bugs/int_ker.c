// gcc -I ../../include/ -L ../.. armin.c -lpolylib32

#include <polylib/polylib32.h>

int main(void) {
    Matrix *lattice1, *lattice2;

    lattice1 = Matrix_Alloc(2, 4);

    lattice1->p[0][0] = 1;
    lattice1->p[0][1] = -1;
    lattice1->p[0][2] = 0;
    lattice1->p[0][3] = 0;
    lattice1->p[1][0] = 0;
    lattice1->p[1][1] = 0;
    lattice1->p[1][2] = 1;
    lattice1->p[1][3] = -1;

    // lattice2->p[0][0] = 2;
    // lattice2->p[0][1] = 1;
    // lattice2->p[1][0] = 0;
    // lattice2->p[1][1] = 1;

    Matrix_Print(stdout, "%3d", lattice1);
    printf("Ker\n");

    Lattice *ker = int_ker(lattice1);
    Matrix_Print(stdout, "%3d", ker);

    return 0;
}

// rendre P full-dimensional :
// ZP = {L, P}  ->  {L.ker, preimage(ker, P)}

// normaliser L :
// ZP = {L, P}  ->  { H, image(U, P)}  avec L = HU forme normale de Hermite de L.

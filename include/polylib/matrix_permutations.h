//  Permutations on matrices
// Matrices are seen either as transformations (mtransformation) or as polyhedra (mpolyhedron)
// B. Meister

#ifndef __BM_MATRIX_PERMUTATIONS_H__
#define __BM_MATRIX_PERMUTATIONS_H__

#include<polylib/polylib.h>
#include<assert.h>

// Permutations here are vectors that give the future position for each variable

// utility function : bit count (i know, there are faster methods)
unsigned int nb_bits(unsigned long long int x);

// gives the inverse permutation vector of a permutation vector
unsigned int * permutation_inverse(unsigned int * perm, unsigned int nb_elems);

// Given a linear tranformation on initial variables, and a variable permutation, compute the tranformation for the permuted variables.
// perm is a vector giving the new "position of the k^th variable, k \in [1..n]
// we can call it a "permutation vector" if you wish
// transf[x][y] -> permuted[permutation(x)][permutation(y)]
Matrix * mtransformation_permute(Matrix * transf, unsigned int * permutation);

// permutes the variables of a matrix seen as a polyhedron
Matrix * mpolyhedron_permute(Matrix * polyh, unsigned int * permutation);

// find a valid permutation : for a set of m equations, find m variables that will be put at the beginning (to be eliminated)
// it must be possible to eliminate these variables : the submatrix built with their columns must be full-rank.
// brute force method, that tests all the combinations until finding one which works.
// LIMITATIONS : up to x-1 variables, where the long long format is x-1 bits (often 64 in year 2005).
unsigned int * find_a_permutation(Matrix * Eqs, unsigned int nb_parms);

// compute the permutation of variables and parameters, according to some variables to keep.
// put the variables not to be kept at the beginning, then the parameters and finally the variables to be kept.
// strongly related to the function compress_to_full_dim2
unsigned int * permutation_for_full_dim2(unsigned int * vars_to_keep, unsigned int nb_keep, unsigned int nb_vars_parms, unsigned int nb_parms);

#endif //__BM_MATRIX_PERMUTATIONS_H__
#include <stdio.h>
#include <stdlib.h>


typedef struct {
  int count;
  int *fac;
} factor;

// compute the prime factors of n, including n itself if it is prime.
factor prime_factors(int n)
{
  int tabsize = 10;
  int div = 2;
  int rest = n;
  factor result;

  // result init
  result.count = 0;
  result.fac = malloc(tabsize * sizeof(int));

  while(div*div <= rest)
  {
    if(rest % div == 0)
    {
      // double the size of the array if necessary
      if(result.count == tabsize)
      {
        tabsize *= 2;
        result.fac = realloc(result.fac, tabsize * sizeof(int));
      }
      // add div to the result
      result.fac[result.count++] = div;
      rest /= div;
    }
    else
      div += (1 + (1&div)); // 2, 3, 5, 7, 9, ...
  }
  if(rest != 1)
  {
    // add rest
    if(result.count == tabsize)
      {
        tabsize += 1;
        result.fac = realloc(result.fac, tabsize * sizeof(int));
      }
      result.fac[result.count++] = rest;
  }
  return result;
}

// compute all numbers in [1, n[ dividing a given number n
factor all_factors(int n)
{
  // compute the prime decomposition (including n if it is prime)
  factor primes = prime_factors(n);
  factor result;

  int size = (1<<(primes.count));   // max size of result = 2^(prime.count)
  result.count = 1;
  result.fac = malloc(size * sizeof(int));
  result.fac[0] = 1;

  for(int mask=1; mask < size; mask++) {
    // multiply the prime factors for the bits that are 1 in the mask :)
    int val = 1;
    for(int bits=0, rest=mask; rest ; bits++, rest>>=1) {
      if((rest & 1))
        val *= primes.fac[bits];
    }
    // add val to res.fac only if it is not already there
    // (this will suppress duplicates, in case a prime factor is there several times)
    if(val > result.fac[result.count-1]) {
      result.fac[result.count++] = val;
    }
  }
  // always remove the last one, that is n.
  result.count--;

  free(primes.fac);
  return result;
}


// // test
// int main(int argc, char **argv)
// {

//   // factor f = prime_factors(atoi(argv[1]));
//   factor f = all_factors(atoi(argv[1]));

//   for(int i=0; i<f.count; i++)
//   {
//     printf(" %d", f.fac[i]);
//   }
//   printf("\n");

//   free(f.fac);
//   return 0;
// }

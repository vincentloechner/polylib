#ifndef _SolveDio_h_
#define _SolveDio_h_

#if defined(__cplusplus)
extern "C" {
#endif

// used by Barvinok
extern int SolveDiophantine(Matrix *M, Matrix **U, Vector **X);

#if defined(__cplusplus)
}
#endif

#endif /* _SolveDio_h_ */

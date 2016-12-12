#include "mkl_types.h"
#include "mkl_lapacke.h"


/*Complex.h Header Here To Be Compact********************************************/
#ifndef _NR_COMPLEX_H_
#define _NR_COMPLEX_H_

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

MKL_Complex16 Cadd(MKL_Complex16 a, MKL_Complex16 b);
MKL_Complex16 Csub(MKL_Complex16 a, MKL_Complex16 b);
MKL_Complex16 Cmul(MKL_Complex16 a, MKL_Complex16 b);
MKL_Complex16 Complex(double re, double im);
MKL_Complex16 Conjg(MKL_Complex16 z);
MKL_Complex16 Cdiv(MKL_Complex16 a, MKL_Complex16 b);
double Cabs(MKL_Complex16 z);
MKL_Complex16 Csqrt(MKL_Complex16 z);
MKL_Complex16 RCmul(double x, MKL_Complex16 a);

#else /* ANSI */
/* traditional - K&R */

MKL_Complex16 Cadd();
MKL_Complex16 Csub();
MKL_Complex16 Cmul();
MKL_Complex16 Complex();
MKL_Complex16 Conjg();
MKL_Complex16 Cdiv();
double Cabs();
MKL_Complex16 Csqrt();
MKL_Complex16 RCmul();

#endif /* ANSI */

#endif /* _NR_COMPLEX_H_ */
/*******************************************************************************/

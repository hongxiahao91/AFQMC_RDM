#include "complex-mkl.h"
#include "mkl_lapacke.h" 
#include "mkl_types.h" 
#include <math.h>

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

MKL_Complex16 Cadd(MKL_Complex16 a, MKL_Complex16 b)
{
        MKL_Complex16 c;
        c.real=a.real+b.real;
        c.imag=a.imag+b.imag;
        return c;
}

MKL_Complex16 Csub(MKL_Complex16 a, MKL_Complex16 b)
{
        MKL_Complex16 c;
        c.real=a.real-b.real;
        c.imag=a.imag-b.imag;
        return c;
}


MKL_Complex16 Cmul(MKL_Complex16 a, MKL_Complex16 b)
{
        MKL_Complex16 c;
        c.real=a.real*b.real-a.imag*b.imag;
        c.imag=a.imag*b.real+a.real*b.imag;
        return c;
}

MKL_Complex16 Complex(double re, double im)
{
        MKL_Complex16 c;
        c.real=re;
        c.imag=im;
        return c;
}

double Real(MKL_Complex16 number) 
{
       return number.real; 
}

MKL_Complex16 Conjg(MKL_Complex16 z)
{
        MKL_Complex16 c;
        c.real=z.real;
        c.imag = -z.imag;
        return c;
}

MKL_Complex16 Cdiv(MKL_Complex16 a, MKL_Complex16 b)
{
        MKL_Complex16 c;
        double r,den;
        if (fabs(b.real) >= fabs(b.imag)) {
                r=b.imag/b.real;
                den=b.real+r*b.imag;
                c.real=(a.real+r*a.imag)/den;
                c.imag=(a.imag-r*a.real)/den;
        } else {
                r=b.real/b.imag;
                den=b.imag+r*b.real;
                c.real=(a.real*r+a.imag)/den;
                c.imag=(a.imag*r-a.real)/den;
        }
        return c;
}
MKL_Complex16 Cexp(MKL_Complex16 a)
{

        MKL_Complex16 c;
        c.real=exp(a.real)*cos(a.imag);
        c.imag=exp(a.real)*sin(a.imag);
return c;
}

double Cabs(MKL_Complex16 z)
{
        double x,y,ans,temp;
        x=fabs(z.real);
        y=fabs(z.imag);
        if (x == 0.0)
                ans=y;
        else if (y == 0.0)
                ans=x;
        else if (x > y) {
                temp=y/x;
                ans=x*sqrt(1.0+temp*temp);
        } else {
                temp=x/y;
                ans=y*sqrt(1.0+temp*temp);
        }
        return ans;
}

MKL_Complex16 Csqrt(MKL_Complex16 z)
{
        MKL_Complex16 c;
        double x,y,w,r;
        if ((z.real == 0.0) && (z.imag == 0.0)) {
                c.real=0.0;
                c.imag=0.0;
                return c;
        } else {
                x=fabs(z.real);
                y=fabs(z.imag);
if (x >= y) {
                        r=y/x;
                        w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
                } else {
                        r=x/y;
                        w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
                }
                if (z.real >= 0.0) {
                        c.real=w;
                        c.imag=z.imag/(2.0*w);
                } else {
                        c.imag=(z.imag >= 0) ? w : -w;
                        c.real=z.imag/(2.0*c.imag);
                }
                return c;
        }
}

MKL_Complex16 RCmul(double x, MKL_Complex16 a)
{
        MKL_Complex16 c;
        c.real=x*a.real;
        c.imag=x*a.imag;
        return c;
}

#else /* ANSI */
/* traditional - K&R */
MKL_Complex16 Cadd(a,b)
MKL_Complex16 a,b;
{
        MKL_Complex16 c;
        c.real=a.real+b.real;
        c.imag=a.imag+b.imag;
        return c;
}

MKL_Complex16 Csub(a,b)
MKL_Complex16 a,b;
{
        MKL_Complex16 c;
        c.real=a.real-b.real;
        c.imag=a.imag-b.imag;
        return c;
}


MKL_Complex16 Cmul(a,b)
MKL_Complex16 a,b;
{
        MKL_Complex16 c;
        c.real=a.real*b.real-a.imag*b.imag;
        c.imag=a.imag*b.real+a.real*b.imag;
        return c;
}

MKL_Complex16 Complex(re,im)
double im,re;
{
        MKL_Complex16 c;
        c.real=re;
c.imag=im;
        return c;
}

MKL_Complex16 Conjg(z)
MKL_Complex16 z;
{
        MKL_Complex16 c;
        c.real=z.real;
        c.imag = -z.imag;
        return c;
}

MKL_Complex16 Cexp(MKL_Complex16 a)
MKL_Complex16 a;
{

        MKL_Complex16 c;
        c.real=exp(a.real)*cos(a.imag);
        c.imag=exp(a.real)*sin(a.imag);
        return c;
}

MKL_Complex16 Cdiv(a,b)
MKL_Complex16 a,b;
{
        MKL_Complex16 c;
        double r,den;
        if (fabs(b.real) >= fabs(b.imag)) {
                r=b.imag/b.real;
                den=b.real+r*b.imag;
                c.real=(a.real+r*a.imag)/den;
                c.imag=(a.imag-r*a.real)/den;
} else {
                r=b.real/b.imag;
                den=b.imag+r*b.real;
                c.real=(a.real*r+a.imag)/den;
                c.imag=(a.imag*r-a.real)/den;
        }
        return c;
}

double Cabs(z)
MKL_Complex16 z;
{
        double x,y,ans,temp;
        x=fabs(z.real);
        y=fabs(z.imag);
        if (x == 0.0)
                ans=y;
        else if (y == 0.0)
                ans=x;
        else if (x > y) {
                temp=y/x;
                ans=x*sqrt(1.0+temp*temp);
        } else {
                temp=x/y;
                ans=y*sqrt(1.0+temp*temp);
        }
        return ans;
}

MKL_Complex16 Csqrt(z)
MKL_Complex16 z;
{
        MKL_Complex16 c;
double x,y,w,r;
        if ((z.real == 0.0) && (z.imag == 0.0)) {
                c.real=0.0;
                c.imag=0.0;
                return c;
        } else {
                x=fabs(z.real);
                y=fabs(z.imag);
                if (x >= y) {
                        r=y/x;
                        w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
                } else {
                        r=x/y;
                        w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
                }
                if (z.real >= 0.0) {
                        c.real=w;
                        c.imag=z.imag/(2.0*w);
                } else {
                        c.imag=(z.imag >= 0) ? w : -w;
                        c.real=z.imag/(2.0*c.imag);
                }
                return c;
        }
}

MKL_Complex16 RCmul(x,a)
MKL_Complex16 a;
double x;
{
        MKL_Complex16 c;
        c.real=x*a.real;
        c.iimag=x*a.imag;
        return c;
}

#endif /* ANSI */

#include "afqmc.h"

/*******************************************************************************/

double ran_nrc(long *idum)
{
  int    j;
  long   k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;
  
  if (*idum <= 0){  /* Initialize */
    if (-(*idum) < 1) *idum = 1;    /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2 = (*idum);

    for (j = NTAB+7; j >= 0; j--){/* Load the shuffle table after 8 warm-ups */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k * IQ1) - k * IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;      /* Start Here when not initializing */
  *idum = IA1 * (*idum - k * IQ1) - k * IR1; /* Compute idum = IA1*idum % IM1*/
  if (*idum < 0) *idum += IM1;               /* without overflow by Scharge  */

  k = idum2 / IQ2;                           /* method.                      */
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2; /* Compute idum2=IA2*idum % IM2 */
  if (idum2 < 0) idum2 += IM2;               /* likewise.                    */
  
  j = iy / NDIV;          /* Will be in the range 0..NTAB-1 */
  iy = iv[j] - idum2;     /* Here idum is shuffled, idum and idum2 are */
  iv[j] = *idum;          /* combined to generate output               */

  if (iy < 1) iy += IMM1;
  if ((temp = AM * iy) > RNMX) return RNMX; /* Because users don't expect */
  else return temp;                         /* endpoint value             */
}

/****************************************************************************/

float gasdev(long *idum) {
   /*Returns a Random Number from a Gaussian Distribution*/

static int iset=0;
static float gset;
float fac, rsq, v1, v2;

if (iset ==0) {
   do {
    v1=2.0*ran1(idum)-1.0;
    v2=2.0*ran1(idum)-1.0;
    rsq=v1*v1+v2*v2;
   }   while (rsq >=1.0 || rsq==0.0);
   fac=sqrt(-2.0*log(rsq)/rsq);
   gset=v1*fac;
   iset=1;
   return v2*fac;
 }
 else {
   iset=0;
   return gset;
 }

}

/**************************************************************************/

void Randomize()
{
  double tt = (double)time(0);
  long   rn = (long)tt % 1000;
  int    i;
  
  for (i = 0; i < rn; i++) random();
  return;
}

/****************************************************************************/

double ran()
{
  double     rnd;
  int static once = 0;
  
  if(!once) {
    struct timeval tv;
    struct timezone tpz;

    gettimeofday(&tv,&tpz);
    srandom(tv.tv_sec);
    once = 1;
  }
  rnd = (double)random() / 2147483647.0;

  return rnd;
}

/****************************************************************************/

float ran1(long *idum)
{
int j;
long k;
float temp;
static long iy=0;
static long iv[NTAB];

if (*idum <= 0 || !iy) {
   if (-(*idum) < 1) *idum=1;
   else *idum = -(*idum);
  for (j=NTAB+7;j>=0;j--) {
     k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
     if (*idum < 0) *idum += IM;
     if (j < NTAB) iv[j] = *idum;
   }
   iy=iv[0];
  }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;

}


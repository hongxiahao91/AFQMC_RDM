#include "afqmc.h"

/*************************************************/

void dzero_vec(double *vec, int size) {

  int i; 

  for ( i=0; i<size; i++) {
     vec[i] = 0.0; 
  } 

return; 
}

/*************************************************/

void izero_vec(int *vec, int size) {

  int i;

  for ( i=0; i<size; i++) {
     vec[i] = 0.0;
  }

return;
}

/*************************************************/ 

void czero_vec(MKL_Complex16*vec, int size) {

  int i;
  MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0; 

  for ( i=0; i<size; i++) {
     vec[i] = Zero;  
  }

return;
}

/*************************************************/

void cunity_vec(MKL_Complex16 *vec, int size) {

   int i, j; 
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0; 
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0; 

   for ( i=0; i<size; i++) {
     for ( j=0; j<size; j++) {
        vec[i*size+j] = Zero; 
     }
     vec[i*size+i] = One; 
   }

return; 
}

/***************************************************/ 

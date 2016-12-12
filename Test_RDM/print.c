#include "afqmc.h"

/********************************************************************************/

void print_dvec(double *vec, int size, char *str) {

   FILE *pf = fopen(str, "a+"); 
   int i; 

   for (i=0; i<size; i++) {
     fprintf(pf, "%f\n", vec[i]); 
   }
   fprintf(pf, "\n\n"); 
   fflush(pf); 

fclose(pf); 
return; 
}

/********************************************************************************/

void print_cvec(MKL_Complex16 *vec, int size, char *str) {

   FILE *pf = fopen(str, "a+");
   int i;

   for (i=0; i<size; i++) {
     fprintf(pf, "%f+%fi\n", vec[i].real, vec[i].imag);
   }
   fprintf(pf, "\n\n");
   fflush(pf);

fclose(pf);
return;
}

/********************************************************************************/

void print_dmat(double *mat, int size1, int size2, char *str) {

    FILE *pf = fopen(str, "a+"); 
    int i, j; 

    for (i=0; i<size1; i++) {
      for (j=0; j<size2; j++) {
        fprintf(pf, "%f\t", mat[i*size2+j]); 
      }
      fprintf(pf, "\n"); 
    }
    fprintf(pf, "\n\n"); 
    fflush(pf); 

fclose(pf); 
return; 
}

/*********************************************************************************/

void print_cmat(MKL_Complex16 *mat, int size1, int size2, char *str) {

    FILE *pf = fopen(str, "a+");
    int i, j;

    for (i=0; i<size1; i++) {
      for (j=0; j<size2; j++) {
        fprintf(pf, "%f+%fi\t", mat[i*size2+j].real, mat[i*size2+j].imag);
      }
      fprintf(pf, "\n");
    }
    fprintf(pf, "\n\n");
    fflush(pf);

fclose(pf);
return;
}

/*********************************************************************************/
  

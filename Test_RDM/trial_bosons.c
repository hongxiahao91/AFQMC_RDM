#include "afqmc.h"

void trial_identity_bosons(MKL_Complex16 *trial_wf_bosons,int_st ist) {

  /*Sets the Trial Wavefunction to the Identity Matrix*/
  
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;  

  czero_vec(trial_wf_bosons,ist.n_sites); 

  trial_wf_bosons[0] = One; 
 
return; 
}

/***********************************************************************************/

void trial_free_bosons(MKL_Complex16 *trial_wf_bosons,double *kinetic_eigs_bosons,double *kinetic_eigvecs_bosons,int_st ist,cns_st cns) {

   /*Sets the Trial Wavefunction to the Trial Kinetic WF*/        
   int i, j, k;
   double p;  
   FILE *pf2 = fopen("trial_wf.dat", "a+");

   /*Order Eigenvalues and Eigenvectors in Descending Order First*/
   for (i=1;i<ist.n_sites;i++) {
     k=i;
     p=kinetic_eigs_bosons[k-1];
     for (j=i+1;j<=ist.n_sites;j++) {
       if (kinetic_eigs_bosons[j-1] >= p) {k=j; p=kinetic_eigs_bosons[k-1]; };
     }
     if (k != i) {
        kinetic_eigs_bosons[k-1]=kinetic_eigs_bosons[i-1];
        kinetic_eigs_bosons[i-1]=p;
        for (j=1;j<=ist.n_sites;j++) {
           p=kinetic_eigvecs_bosons[(j-1)*ist.n_sites+i-1];
           kinetic_eigvecs_bosons[(j-1)*ist.n_sites+i-1]=kinetic_eigvecs_bosons[(j-1)*ist.n_sites+k-1];
           kinetic_eigvecs_bosons[(j-1)*ist.n_sites+k-1]=p;
         }
      }
    }

   /*Copy Lowest Up Down Eigenvectors Into Matrices*/
   for (j=0; j<ist.n_sites; j++) {
      trial_wf_bosons[j].real = kinetic_eigvecs_bosons[j*ist.n_sites+(ist.n_sites-1)]; 
      trial_wf_bosons[j].imag = 0.0; 
   }

   /*Save Trial Energy*/
   cns.trial_energy = kinetic_eigs_bosons[ist.n_sites-1]; 

   /*Print Trial Matrices*/
   print_cvec(trial_wf_bosons,ist.n_sites,"trial_wf.dat");
   fprintf(pf2, "Trial Energy: %g\n", cns.trial_energy); fflush(pf2);   

fclose(pf2); 
return; 
}  

/************************************************************************************/

void trial_random_bosons(MKL_Complex16 *trial_wf_bosons,int_st ist) {

  int ip; 
  double value; 
  long tidum;
  FILE *pf = fopen("trial_wf.dat", "a+"); 

  Randomize(); tidum = -random();

  czero_vec(trial_wf_bosons,ist.n_sites);

  value = 0.0; 
  for (ip=0; ip < ist.n_sites; ip++ ) {
    trial_wf_bosons[ip].real = ran1(&tidum);  
    trial_wf_bosons[ip].imag = 0.0; 
    value += trial_wf_bosons[ip].real * trial_wf_bosons[ip].real; 
  }
  value = 1.0/sqrt(value); 

  for (ip = 0; ip<ist.n_sites; ip++) {
    trial_wf_bosons[ip].real *= value; 
  } 

  print_cvec(trial_wf_bosons,ist.n_sites,"trial_wf.dat"); 

fclose(pf); 
return;
}

/***********************************************************************************/

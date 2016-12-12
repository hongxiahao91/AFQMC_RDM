#include "afqmc.h"

void trial_identity_fermions_spin_polarized(MKL_Complex16 *trial_wf_up,int_st ist) {

  /*Sets the Trial Wavefunction to the Identity Matrix*/
  
  int ip;
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;  

  czero_vec(trial_wf_up,ist.n_sites*ist.n_up); 

  for ( ip = 0; ip < ist.n_up ; ip++ ) {
     trial_wf_up[ip * ist.n_up + ip] = One; 
  }
 
return; 
}

/***********************************************************************************/

void trial_free_fermions_spin_polarized(MKL_Complex16 *trial_wf_up,double *kinetic_eigs,double *kinetic_eigvecs,int_st ist,cns_st cns) {

   /*Sets the Trial Wavefunction to the Trial Kinetic WF*/        
   int i, j, k;
   double p;  
   FILE *pf2 = fopen("trial_wf.par", "f");

   /*Order Eigenvalues and Eigenvectors in Descending Order First*/
   for (i=1;i<ist.n_sites;i++) {
     k=i;
     p=kinetic_eigs[k-1];
     for (j=i+1;j<=ist.n_sites;j++) {
       if (kinetic_eigs[j-1] >= p) {k=j; p=kinetic_eigs[k-1]; };
     }
     if (k != i) {
        kinetic_eigs[k-1]=kinetic_eigs[i-1];
        kinetic_eigs[i-1]=p;
        for (j=1;j<=ist.n_sites;j++) {
           p=kinetic_eigvecs[(j-1)*ist.n_sites+i-1];
           kinetic_eigvecs[(j-1)*ist.n_sites+i-1]=kinetic_eigvecs[(j-1)*ist.n_sites+k-1];
           kinetic_eigvecs[(j-1)*ist.n_sites+k-1]=p;
         }
      }
    }

   /*Copy Lowest Up Down Eigenvectors Into Matrices*/
   for (i=0; i<ist.n_up; i++) {
      for (j=0; j<ist.n_sites; j++) {
        trial_wf_up[j*ist.n_up+i].real = kinetic_eigvecs[j*ist.n_sites+(ist.n_sites-i-1)]; 
        trial_wf_up[j*ist.n_up+i].imag = 0.0; 
      }
   }

   /*Save Trial Energy*/
   cns.trial_energy = kinetic_eigs[ist.n_sites-1]; 

   /*Print Trial Matrices*/
   print_cmat(trial_wf_up,ist.n_sites,ist.n_up,"trial_wf.dat");
   fprintf(pf2, "Trial Energy: %g\n", cns.trial_energy); fflush(pf2);   

fclose(pf2); 
return; 
}  

/************************************************************************************/

void trial_rhf_fermions_spin_polarized(MKL_Complex16 *trial_wf_up,int_st ist) {

   /*Sets the Trial Wavefunction to the Restricted Hartree-Fock Wavefunction Using Coefficients from MolPro*/
   int in, ip;
   FILE *pf = fopen("trial_rhf_wf.par", "r");

   czero_vec(trial_wf_up,ist.n_sites*ist.n_up);

   /*Get Up Components First*/
   for ( ip = 0; ip < ist.n_up; ip++ ) {
    for ( in = 0; in < ist.n_sites; in++ ) {
      fscanf(pf, "%lf", &trial_wf_up[ip * ist.n_sites + in].real);
    }
   }

fclose(pf);
return;
}

/********************************************************************************/

void trial_random_fermions_spin_polarized(MKL_Complex16 *trial_wf_up,int_st ist) {

  int ip, jp;
  long tidum;
  FILE *pf = fopen("trial_wf.dat", "a+"); 

  Randomize(); tidum = -random();

  czero_vec(trial_wf_up,ist.n_sites*ist.n_up);

  for (ip=0; ip < ist.n_sites; ip++ ) {
    for (jp=0; jp < ist.n_up; jp++ ) {
       trial_wf_up[ip*ist.n_up + jp].real = ran_nrc(&tidum);  
       trial_wf_up[ip*ist.n_up + jp].imag = 0.0; 

    }
  }

  print_cmat(trial_wf_up,ist.n_sites,ist.n_up,"trial_wf.dat"); 

fclose(pf); 
return;
}

/***********************************************************************************/

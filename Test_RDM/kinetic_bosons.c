#include "afqmc.h"

/*******************************************************************/

void propagate_forwards_kinetic_bosons(double *kinetic_full_bosons,MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i; 
   MKL_Complex16 previous_overlap_bosons; 
   MKL_Complex16 *stored_product_bosons; 

   stored_product_bosons=(MKL_Complex16*)calloc(ist.n_sites,sizeof(MKL_Complex16)); 

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) { 
  
      /*Store Previous Overlaps*/
      previous_overlap_bosons = overlap_bosons[i]; 

      dmat_cmat(kinetic_full_bosons,&wf_bosons[i*ist.n_sites],stored_product_bosons,ist.n_sites,ist.n_sites,1);

      cblas_zcopy(ist.n_sites,stored_product_bosons,1,&wf_bosons[i*ist.n_sites],1);
      //copy_cmat(&wf_bosons[i*ist.n_sites],stored_product_bosons,ist.n_sites); 

      /*Get Product of Trial and Actual*/
      update_overlaps_bosons(&wf_bosons[i*ist.n_sites],trial_wf_bosons,&overlap_bosons[i],ist);

      /*Get Overall Weight*/
      weights[i] = Cmul(weights[i], Cdiv(overlap_bosons[i], previous_overlap_bosons));  
 
    } /*weights walkers*/

   } /*i*/

free(stored_product_bosons); 
return; 
}

/******************************************************************/

void propagate_half_backwards_kinetic_bosons(double *kinetic_backwards_half_bosons,MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *weights,int_st ist){ 

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 previous_overlap_bosons; 
    MKL_Complex16 *stored_product_bosons; 
    FILE *pf = fopen("checkpropagatehalfbackwards.dat", "a+"); 

    stored_product_bosons=(MKL_Complex16*)calloc(ist.n_sites,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) { 

       /*Previous Overlaps*/
       previous_overlap_bosons = overlap_bosons[i]; 

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_bosons,&wf_bosons[i*ist.n_sites],stored_product_bosons,ist.n_sites,ist.n_sites,1);

       cblas_zcopy(ist.n_sites,stored_product_bosons,1,&wf_bosons[i*ist.n_sites],1);
       //copy_cmat(&wf_bosons[i*ist.n_sites],stored_product_bosons,ist.n_sites);

       /*Get Product of Trial and Actual*/
       update_overlaps_bosons(&wf_bosons[i*ist.n_sites],trial_wf_bosons,&overlap_bosons[i],ist); 

       if ( i == 0 ) {
         fprintf(pf, "trial wf %f+%fi %f+%fi\n", trial_wf_bosons[0].real, trial_wf_bosons[0].imag, trial_wf_bosons[1].real, trial_wf_bosons[1].imag); fflush(pf); 
         fprintf(pf, "overlap %f\n", overlap_bosons[0].real); 
       }
 
       /*Get Overall Weights*/
       weights[i] = Cmul(weights[i], Cdiv(overlap_bosons[i], previous_overlap_bosons));  
  
      } /*Non Zero Weights*/ 

     }

free(stored_product_bosons);
}

/***********************************************************************/

void propagate_half_forwards_kinetic_bosons(double *kinetic_forwards_half_bosons,MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i; 
    MKL_Complex16 previous_overlap_bosons; 
    MKL_Complex16 *stored_product_bosons; 

    stored_product_bosons=(MKL_Complex16*)calloc(ist.n_sites,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_bosons = overlap_bosons[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_bosons,&wf_bosons[i*ist.n_sites],stored_product_bosons,ist.n_sites,ist.n_sites,1);

       cblas_zcopy(ist.n_sites,stored_product_bosons,1,&wf_bosons[i*ist.n_sites],1);
       //copy_cmat(&wf_bosons[i*ist.n_sites],stored_product_bosons,ist.n_sites);

       /*Get Product of Trial and Actual*/
       update_overlaps_bosons(&wf_bosons[i*ist.n_sites],trial_wf_bosons,&overlap_bosons[i],ist);

       /*Get Overall Weights*/
       weights[i] = Cmul(weights[i], Cdiv(overlap_bosons[i], previous_overlap_bosons)); 
  
     }

   }

free(stored_product_bosons);
return;
} 

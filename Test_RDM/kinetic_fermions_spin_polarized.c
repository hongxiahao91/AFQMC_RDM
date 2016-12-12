#include "afqmc.h"

/*******************************************************************/

void propagate_forwards_kinetic_fermions_spin_polarized(double *kinetic_full_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;  
   int n_sites_n_up = ist.n_sites * ist.n_up; 
   int n_up_sq = ist.n_up * ist.n_up; 
   MKL_Complex16 previous_overlap_up; 
   MKL_Complex16 weights_numerator, weights_denominator; 
   MKL_Complex16 *stored_product_up; 
   MKL_Complex16 *stored_product_up_2; 
   MKL_Complex16 *One, *Zero;

   One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   One[0].real = 1.0; One[0].imag = 0.0;
   Zero[0].real = 0.0; Zero[0].imag = 0.0;

   stored_product_up=(MKL_Complex16*)calloc(n_sites_n_up,sizeof(MKL_Complex16)); 
   stored_product_up_2=(MKL_Complex16*)calloc(n_up_sq,sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) { 

      /*Store Previous Overlaps*/
      previous_overlap_up = overlap_up[i]; 

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_up[i*n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);

      cblas_zcopy(ist.n_sites_n_up,stored_product_up,1,&wf_up[i*ist.n_sites_n_up],1);
      //copy_cmat(&wf_up[i*n_sites_n_up],stored_product_up,n_sites_n_up); 

      /*Get Product of Trial and Actual*/
      //transpose_cmat_cmat(trial_wf_up,&wf_up[i*n_sites_n_up],stored_product_up_2,ist.n_sites,ist.n_up,ist.n_up); 
      cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,&wf_up[i*ist.n_sites_n_up],ist.n_up,Zero,stored_product_up_2,ist.n_up);
      overlap_up[i]=complex_inverse_det_fermions(stored_product_up_2,&overlap_inverse_up[i*n_up_sq],ist.n_up); 

      /*Get Overall Weight*/
      weights_numerator = overlap_up[i]; 
      weights_denominator = previous_overlap_up; 
      weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator)); 

    } /*weights walkers*/

   } /*i*/

free(stored_product_up); 
free(stored_product_up_2); 
free(Zero); 
free(One); 
return; 
}

/******************************************************************/

void propagate_half_backwards_kinetic_fermions_spin_polarized(double *kinetic_backwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int_st ist){ 

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    int n_sites_n_up = ist.n_sites * ist.n_up; 
    int n_up_sq = ist.n_up * ist.n_up; 
    MKL_Complex16 previous_overlap_up; 
    MKL_Complex16 weights_numerator, weights_denominator; 
    MKL_Complex16 *stored_product_up; 
    MKL_Complex16 *stored_product_up_2; 
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_up=(MKL_Complex16*)calloc(n_sites_n_up,sizeof(MKL_Complex16));
    stored_product_up_2=(MKL_Complex16*)calloc(n_up_sq,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) { 

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i]; 

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_up[i*n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);

       cblas_zcopy(ist.n_sites_n_up,stored_product_up,1,&wf_up[i*ist.n_sites_n_up],1);
       //copy_cmat(&wf_up[i*n_sites_n_up],stored_product_up,n_sites_n_up);

       /*Get Product of Trial and Actual*/
       cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,&wf_up[i*ist.n_sites_n_up],ist.n_up,Zero,stored_product_up_2,ist.n_up);
       //transpose_cmat_cmat(trial_wf_up,&wf_up[i*n_sites_n_up],stored_product_up_2,ist.n_sites,ist.n_up,ist.n_up);
       overlap_up[i]=complex_inverse_det_fermions(stored_product_up_2,&overlap_inverse_up[i*n_up_sq],ist.n_up);
 
       /*Get Overall Weights*/
       weights_numerator = overlap_up[i]; 
       weights_denominator = previous_overlap_up; 
       weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));  
  
      } /*Non Zero Weights*/ 

     }

free(One); 
free(Zero); 
free(stored_product_up);
free(stored_product_up_2); 
return;   
}

/***********************************************************************/

void propagate_half_forwards_kinetic_fermions_spin_polarized(double *kinetic_forwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i; 
    int n_sites_n_up = ist.n_sites * ist.n_up;
    int n_up_sq = ist.n_up * ist.n_up; 
    MKL_Complex16 weight_numerator, weight_denominator; 
    MKL_Complex16 previous_overlap_up; 
    MKL_Complex16 *stored_product_up; 
    MKL_Complex16 *stored_product_up_2; 
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_up=(MKL_Complex16*)calloc(n_sites_n_up,sizeof(MKL_Complex16));
    stored_product_up_2=(MKL_Complex16*)calloc(n_up_sq,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_up[i*n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);

       cblas_zcopy(ist.n_sites_n_up,stored_product_up,1,&wf_up[i*ist.n_sites_n_up],1);
       //copy_cmat(&wf_up[i*n_sites_n_up],stored_product_up,n_sites_n_up);

       /*Get Product of Trial and Actual*/
       cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,&wf_up[i*ist.n_sites_n_up],ist.n_up,Zero,stored_product_up_2,ist.n_up); 
       //transpose_cmat_cmat(trial_wf_up,&wf_up[i*n_sites_n_up],stored_product_up_2,ist.n_sites,ist.n_up,ist.n_up);
       overlap_up[i]=complex_inverse_det_fermions(stored_product_up_2,&overlap_inverse_up[i*n_up_sq],ist.n_up);
          
       /*Get Overall Weights*/
       weight_numerator = overlap_up[i]; 
       weight_denominator = previous_overlap_up; 
       weights[i] = Cmul(weights[i], Cdiv(weight_numerator, weight_denominator)); 
  
      }

     }

free(stored_product_up);
free(stored_product_up_2);
free(One); 
free(Zero); 
return;
} 

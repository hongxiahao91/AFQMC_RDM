#include "afqmc.h"

/*******************************************************************/

void propagate_forwards_kinetic_fermions(double *kinetic_full_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;  
   MKL_Complex16 previous_overlap_up, previous_overlap_down; 
   MKL_Complex16 weights_numerator, weights_denominator; 
   MKL_Complex16 *stored_product_up, *stored_product_down;    
   MKL_Complex16 *stored_product_up_2, *stored_product_down_2; 
   MKL_Complex16 *One, *Zero;

   One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   One[0].real = 1.0; One[0].imag = 0.0;
   Zero[0].real = 0.0; Zero[0].imag = 0.0;

   stored_product_up=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16)); 
   stored_product_down=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));

   stored_product_up_2=(MKL_Complex16*)calloc(ist.n_up_sq,sizeof(MKL_Complex16));
   stored_product_down_2=(MKL_Complex16*)calloc(ist.n_down_sq,sizeof(MKL_Complex16)); 

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) { 
  
      /*Store Previous Overlaps*/
      previous_overlap_up = overlap_up[i]; 
      previous_overlap_down = overlap_down[i]; 

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_up[i*ist.n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);
      cblas_zcopy(ist.n_sites_n_up,stored_product_up,1,&wf_up[i*ist.n_sites_n_up],1);

      /*Get Product of Trial and Actual*/
      cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,&wf_up[i*ist.n_sites_n_up],ist.n_up,Zero,stored_product_up_2,ist.n_up); 
      overlap_up[i]=complex_inverse_det_fermions(stored_product_up_2,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up); 

      /*Get Down Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_down[i*ist.n_sites_n_down],stored_product_down,ist.n_sites,ist.n_sites,ist.n_down);
      cblas_zcopy(ist.n_sites_n_down,stored_product_down,1,&wf_down[i*ist.n_sites_n_down],1);
   
      /*Get Product of Trial and Down Actual*/
      cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_sites,One,trial_wf_down,ist.n_down,&wf_down[i*ist.n_sites_n_down],ist.n_down,Zero,stored_product_down_2,ist.n_down); 
      overlap_down[i]=complex_inverse_det_fermions(stored_product_down_2,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down); 

      /*Get Overall Weight*/
      weights_numerator = Cmul(overlap_up[i], overlap_down[i]); 
      weights_denominator = Cmul(previous_overlap_up, previous_overlap_down); 
      weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator)); 
 
    } /*weights walkers*/

   } /*i*/

free(One); 
free(Zero); 
free(stored_product_up); 
free(stored_product_down); 
free(stored_product_up_2); 
free(stored_product_down_2); 
return; 
}

/*******************************************************************/

void propagate_forwards_kinetic_fermions_up(double *kinetic_full_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;
   MKL_Complex16 previous_overlap_up; 
   MKL_Complex16 weights_numerator, weights_denominator;
   MKL_Complex16 *stored_product_up; 
   MKL_Complex16 *stored_product_up_2; 
   MKL_Complex16 *One, *Zero;

   One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   One[0].real = 1.0; One[0].imag = 0.0;
   Zero[0].real = 0.0; Zero[0].imag = 0.0;

   stored_product_up=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));
   stored_product_up_2=(MKL_Complex16*)calloc(ist.n_up_sq,sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) {

      /*Store Previous Overlaps*/
      previous_overlap_up = overlap_up[i];

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_up[i*ist.n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);
      cblas_zcopy(ist.n_sites_n_up,stored_product_up,1,&wf_up[i*ist.n_sites_n_up],1);

      /*Get Product of Trial and Actual*/
      cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,&wf_up[i*ist.n_sites_n_up],ist.n_up,Zero,stored_product_up_2,ist.n_up);
      overlap_up[i]=complex_inverse_det_fermions(stored_product_up_2,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up);

      /*Get Overall Weight*/
      weights_numerator = overlap_up[i]; 
      weights_denominator = previous_overlap_up; 
      weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));

    } /*weights walkers*/

   } /*i*/

free(One);
free(Zero);
free(stored_product_up);
free(stored_product_up_2);
return;
}

/*******************************************************************/

void propagate_forwards_kinetic_fermions_down(double *kinetic_full_fermions,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;
   MKL_Complex16 previous_overlap_down;  
   MKL_Complex16 weights_numerator, weights_denominator;
   MKL_Complex16 *stored_product_down;  
   MKL_Complex16 *stored_product_down_2;  
   MKL_Complex16 *One, *Zero;
    
   One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   One[0].real = 1.0; One[0].imag = 0.0;
   Zero[0].real = 0.0; Zero[0].imag = 0.0;
    
   stored_product_down=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));
   stored_product_down_2=(MKL_Complex16*)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) {

      /*Store Previous Overlaps*/
      previous_overlap_down = overlap_down[i];

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_down[i*ist.n_sites_n_down],stored_product_down,ist.n_sites,ist.n_sites,ist.n_down);
      cblas_zcopy(ist.n_sites_n_down,stored_product_down,1,&wf_down[i*ist.n_sites_n_down],1);

      /*Get Product of Trial and Actual*/
      cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_sites,One,trial_wf_down,ist.n_down,&wf_down[i*ist.n_sites_n_down],ist.n_down,Zero,stored_product_down_2,ist.n_down);
      overlap_down[i]=complex_inverse_det_fermions(stored_product_down_2,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down);

      /*Get Overall Weight*/
      weights_numerator = overlap_down[i]; 
      weights_denominator = previous_overlap_down; 
      weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));

    } /*weights walkers*/

   } /*i*/

free(One);
free(Zero);
free(stored_product_down);
free(stored_product_down_2);
return;
}

/******************************************************************/

void propagate_half_backwards_kinetic_fermions(double *kinetic_backwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){ 

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i; 
    MKL_Complex16 previous_overlap_up, previous_overlap_down; 
    MKL_Complex16 weights_numerator, weights_denominator; 
    MKL_Complex16 *stored_product_up, *stored_product_down;
    MKL_Complex16 *stored_product_up_2, *stored_product_down_2; 
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_up=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));

    stored_product_up_2=(MKL_Complex16*)calloc(ist.n_up_sq,sizeof(MKL_Complex16));
    stored_product_down_2=(MKL_Complex16*)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) { 

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i]; 
       previous_overlap_down = overlap_down[i]; 

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_up[i*ist.n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);
       cblas_zcopy(ist.n_sites_n_up,stored_product_up,1,&wf_up[i*ist.n_sites_n_up],1);

       /*Get Product of Trial and Actual*/
       cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,&wf_up[i*ist.n_sites_n_up],ist.n_up,Zero,stored_product_up_2,ist.n_up);
       overlap_up[i]=complex_inverse_det_fermions(stored_product_up_2,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up);
 
       /*Get Down Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_down[i*ist.n_sites_n_down],stored_product_down,ist.n_sites,ist.n_sites,ist.n_down);
       cblas_zcopy(ist.n_sites_n_down,stored_product_down,1,&wf_down[i*ist.n_sites_n_down],1);

       /*Get Product of Trial and Actual*/
       cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_sites,One,trial_wf_down,ist.n_down,&wf_down[i*ist.n_sites_n_down],ist.n_down,Zero,stored_product_down_2,ist.n_down);
       overlap_down[i]=complex_inverse_det_fermions(stored_product_down_2,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down);

       /*Get Overall Weights*/
       weights_numerator = Cmul(overlap_up[i], overlap_down[i]); 
       weights_denominator = Cmul(previous_overlap_up, previous_overlap_down); 
       weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));  
  
      } /*Non Zero Weights*/ 

     }

free(One); 
free(Zero); 
free(stored_product_up);
free(stored_product_down);
free(stored_product_up_2); 
free(stored_product_down_2);
return;   
}

/******************************************************************/

void propagate_half_backwards_kinetic_fermions_up(double *kinetic_backwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 previous_overlap_up; 
    MKL_Complex16 weights_numerator, weights_denominator;
    MKL_Complex16 *stored_product_up; 
    MKL_Complex16 *stored_product_up_2; 
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_up=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));
    stored_product_up_2=(MKL_Complex16*)calloc(ist.n_up_sq,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_up[i*ist.n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);
       cblas_zcopy(ist.n_sites_n_up,stored_product_up,1,&wf_up[i*ist.n_sites_n_up],1);

       /*Get Product of Trial and Actual*/
       cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,&wf_up[i*ist.n_sites_n_up],ist.n_up,Zero,stored_product_up_2,ist.n_up);
       overlap_up[i]=complex_inverse_det_fermions(stored_product_up_2,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up);

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

/******************************************************************/

void propagate_half_backwards_kinetic_fermions_down(double *kinetic_backwards_half_fermions,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 previous_overlap_down;
    MKL_Complex16 weights_numerator, weights_denominator;
    MKL_Complex16 *stored_product_down;
    MKL_Complex16 *stored_product_down_2;
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_down=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));
    stored_product_down_2=(MKL_Complex16*)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_down = overlap_down[i];

       /*Get Down Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_down[i*ist.n_sites_n_down],stored_product_down,ist.n_sites,ist.n_sites,ist.n_down);
       cblas_zcopy(ist.n_sites_n_down,stored_product_down,1,&wf_down[i*ist.n_sites_n_down],1);

       /*Get Product of Trial and Actual*/
       cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_sites,One,trial_wf_down,ist.n_down,&wf_down[i*ist.n_sites_n_down],ist.n_down,Zero,stored_product_down_2,ist.n_down);
       overlap_down[i]=complex_inverse_det_fermions(stored_product_down_2,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down);

       /*Get Overall Weights*/
       weights_numerator = overlap_down[i]; 
       weights_denominator = previous_overlap_down; 
       weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));

      } /*Non Zero Weights*/

     }

free(One);
free(Zero);
free(stored_product_down);
free(stored_product_down_2);
return;
}

/******************************************************************/

void propagate_half_forwards_kinetic_fermions(double *kinetic_forwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i; 
    MKL_Complex16 weight_numerator, weight_denominator; 
    MKL_Complex16 previous_overlap_up, previous_overlap_down;
    MKL_Complex16 *stored_product_up, *stored_product_down;
    MKL_Complex16 *stored_product_up_2, *stored_product_down_2;
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_up=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));

    stored_product_up_2=(MKL_Complex16*)calloc(ist.n_up_sq,sizeof(MKL_Complex16));
    stored_product_down_2=(MKL_Complex16*)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i];
       previous_overlap_down = overlap_down[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_up[i*ist.n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);
       cblas_zcopy(ist.n_sites_n_up,stored_product_up,1,&wf_up[i*ist.n_sites_n_up],1);

       /*Get Product of Trial and Actual*/
       cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,&wf_up[i*ist.n_sites_n_up],ist.n_up,Zero,stored_product_up_2,ist.n_up);
       overlap_up[i]=complex_inverse_det_fermions(stored_product_up_2,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up);
          
       /*Get Down Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_down[i*ist.n_sites_n_down],stored_product_down,ist.n_sites,ist.n_sites,ist.n_down);
       cblas_zcopy(ist.n_sites_n_down,stored_product_down,1,&wf_down[i*ist.n_sites_n_down],1);

       /*Get Product of Trial and Actual*/
       cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_sites,One,trial_wf_down,ist.n_down,&wf_down[i*ist.n_sites_n_down],ist.n_down,Zero,stored_product_down_2,ist.n_down);
       overlap_down[i]=complex_inverse_det_fermions(stored_product_down_2,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down);

       /*Get Overall Weights*/
       weight_numerator = Cmul(overlap_up[i], overlap_down[i]); 
       weight_denominator = Cmul(previous_overlap_up, previous_overlap_down); 
       weights[i] = Cmul(weights[i], Cdiv(weight_numerator, weight_denominator)); 
  
      }

     }

free(One); 
free(Zero); 
free(stored_product_up);
free(stored_product_down);
free(stored_product_up_2);
free(stored_product_down_2);
return;
}

/******************************************************************/

void propagate_half_forwards_kinetic_fermions_up(double *kinetic_forwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 weight_numerator, weight_denominator;
    MKL_Complex16 previous_overlap_up; 
    MKL_Complex16 *stored_product_up; 
    MKL_Complex16 *stored_product_up_2; 
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_up=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));
    stored_product_up_2=(MKL_Complex16*)calloc(ist.n_up_sq,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_up[i*ist.n_sites_n_up],stored_product_up,ist.n_sites,ist.n_sites,ist.n_up);
       cblas_zcopy(ist.n_sites_n_up,stored_product_up,1,&wf_up[i*ist.n_sites_n_up],1);

       /*Get Product of Trial and Actual*/
       cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,&wf_up[i*ist.n_sites_n_up],ist.n_up,Zero,stored_product_up_2,ist.n_up);
       overlap_up[i]=complex_inverse_det_fermions(stored_product_up_2,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up);

       /*Get Overall Weights*/
       weight_numerator = overlap_up[i]; 
       weight_denominator = previous_overlap_up; 
       weights[i] = Cmul(weights[i], Cdiv(weight_numerator, weight_denominator));

      }

     }

free(One);
free(Zero);
free(stored_product_up);
free(stored_product_up_2);
return;
}

/******************************************************************/

void propagate_half_forwards_kinetic_fermions_down(double *kinetic_forwards_half_fermions,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 weight_numerator, weight_denominator;
    MKL_Complex16 previous_overlap_down;
    MKL_Complex16 *stored_product_down;
    MKL_Complex16 *stored_product_down_2;
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_down=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));
    stored_product_down_2=(MKL_Complex16*)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_down = overlap_down[i];

       /*Get Down Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_down[i*ist.n_sites_n_down],stored_product_down,ist.n_sites,ist.n_sites,ist.n_down);
       cblas_zcopy(ist.n_sites_n_down,stored_product_down,1,&wf_down[i*ist.n_sites_n_down],1);

       /*Get Product of Trial and Actual*/
       cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_sites,One,trial_wf_down,ist.n_down,&wf_down[i*ist.n_sites_n_down],ist.n_down,Zero,stored_product_down_2,ist.n_down);
       overlap_down[i]=complex_inverse_det_fermions(stored_product_down_2,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down);

       /*Get Overall Weights*/
       weight_numerator = overlap_down[i]; 
       weight_denominator = previous_overlap_down; 
       weights[i] = Cmul(weights[i], Cdiv(weight_numerator, weight_denominator));

      }

     }

free(One);
free(Zero);
free(stored_product_down);
free(stored_product_down_2);
return;
} 

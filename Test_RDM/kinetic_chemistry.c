#include "afqmc.h"

/*******************************************************************/

void propagate_forwards_kinetic_chemistry_real_restricted(double *kinetic_full_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/
   int i;
   MKL_Complex16 previous_overlap_up, previous_overlap_down;
   MKL_Complex16 weights_numerator, weights_denominator;
   MKL_Complex16 *stored_product_up, *stored_product_down;

   stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
   stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) {

      /*Store Previous Overlaps*/
      previous_overlap_up = overlap_up[i];
      previous_overlap_down = overlap_down[i];

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
      cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

      /*Get Product of Trial and Actual*/
      overlap_up[i] = compute_overlap_inverse(trial_wf_up,&wf_up[i*ist.n_spatial_orbitals_n_up],&overlap_inverse_up[i*ist.n_up_sq],ist,0);

      /*Get Down Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
      cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

      /*Get Product of Trial and Down Actual*/
      overlap_down[i] = compute_overlap_inverse(trial_wf_down,&wf_down[i*ist.n_spatial_orbitals_n_down],&overlap_inverse_down[i*ist.n_down_sq],ist,1);

      /*Get Overall Weight*/
      weights_numerator = Cmul(overlap_up[i], overlap_down[i]);
      weights_denominator = Cmul(previous_overlap_up, previous_overlap_down);
      weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));

    } /*weights walkers*/

   } /*i*/


free(stored_product_up);
free(stored_product_down);
return;
}

/***********************************************************************************/

void propagate_forwards_kinetic_chemistry_real_multi_restricted(double *kinetic_full_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   /*Determine What the Spin Situation Is*/
   if ( ist.n_up > 0 ) {
    if ( ist.n_down > 0 ) { 
      propagate_forwards_kinetic_chemistry_real_multi_restricted_both(kinetic_full_fermions,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,ist); 
    }
    else {
      propagate_forwards_kinetic_chemistry_real_multi_up(kinetic_full_fermions,wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_total,det_overlap_up,weights,ist);  
    }
  }
  else {
    if ( ist.n_down > 0 ) {
      propagate_forwards_kinetic_chemistry_real_multi_down(kinetic_full_fermions,wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_down,overlap_total,det_overlap_down,weights,ist);  
    }
  }   


return; 
}

/***********************************************************************************/

void propagate_forwards_kinetic_chemistry_real_multi_unrestricted(double *kinetic_full_fermions_up,double *kinetic_full_fermions_down,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *phases,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   /*Determine What the Spin Situation Is*/
   if ( ist.n_up > 0 ) {
    if ( ist.n_down > 0 ) {
      propagate_forwards_kinetic_chemistry_real_multi_unrestricted_both(kinetic_full_fermions_up,kinetic_full_fermions_down,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,phases,ist);
    }
    else {
      propagate_forwards_kinetic_chemistry_real_multi_up(kinetic_full_fermions_up,wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_total,det_overlap_up,weights,ist);
    }
  }
  else {
    if ( ist.n_down > 0 ) {
      propagate_forwards_kinetic_chemistry_real_multi_down(kinetic_full_fermions_down,wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_down,overlap_total,det_overlap_down,weights,ist);
    }
  }


return;
}

/***********************************************************************************/

void propagate_forwards_kinetic_chemistry_real_multi_restricted_both(double *kinetic_full_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;
   MKL_Complex16 previous_overlap; 
   MKL_Complex16 *stored_product_up, *stored_product_down;

   stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
   stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) {

      /*Store Previous Overlaps*/
      previous_overlap = overlap_total[i];

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
      cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

      /*Get Down Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
      cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

      /*Get Product of Trial and Down Actual*/
      overlap_total[i] = compute_overlap_inverse_multi_phaseless_both(trial_wf_up_phaseless,trial_wf_down_phaseless,&wf_up[i*ist.n_spatial_orbitals_n_up],&wf_down[i*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_up[i*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[i*ist.n_determinants_phaseless_n_down_sq],&det_overlap_up[i*ist.n_determinants_trial_phaseless],&det_overlap_down[i*ist.n_determinants_trial_phaseless],ist); 

      /*Get Overall Weight*/
      weights[i] = Cmul(weights[i], Cdiv(overlap_total[i], previous_overlap)); 

    } /*weights walkers*/

   } /*i*/

free(stored_product_up);
free(stored_product_down);
return;
}

/***********************************************************************************/

void propagate_forwards_kinetic_chemistry_real_multi_unrestricted_both(double *kinetic_full_fermions_up,double *kinetic_full_fermions_down,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *phases,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;
   double phase, theta_1, theta_2, projected_ratio_determinant; 
   MKL_Complex16 previous_overlap, determinant_after; 
   MKL_Complex16 *stored_product_up, *stored_product_down;

   stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
   stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) {

      /*Store Previous Overlaps*/
      previous_overlap = overlap_total[i];

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions_up,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
      cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

      /*Get Down Matrices*/
      dmat_cmat(kinetic_full_fermions_down,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
      cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

      /*Get Product of Trial and Down Actual*/
      overlap_total[i] = compute_overlap_inverse_multi_phaseless_both(trial_wf_up_phaseless,trial_wf_down_phaseless,&wf_up[i*ist.n_spatial_orbitals_n_up],&wf_down[i*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_up[i*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[i*ist.n_determinants_phaseless_n_down_sq],&det_overlap_up[i*ist.n_determinants_trial_phaseless],&det_overlap_down[i*ist.n_determinants_trial_phaseless],ist);

      /*Multiply By New Overlap*/
      determinant_after = overlap_total[i];

      /*Project Down Last Part of Phase Just In Case Near 0 in Axis*/
      theta_1 = acos(determinant_after.real/Cabs(determinant_after));
      theta_2 = acos(previous_overlap.real/Cabs(previous_overlap));
      phase = theta_1 - theta_2;
      phases[i].real = phase;   
     /*
      if ( theta_1 < .5 * 3.1415 && theta_1 >= 0  ) {
       projected_ratio_determinant = 1; //cos(phase);
      }
      else {
        //projected_ratio_determinant = 1.0; 
        projected_ratio_determinant = 0.0;
      }
      weights[i] = RCmul(projected_ratio_determinant, weights[i]);
      */


    } /*weights walkers*/

   } /*i*/

free(stored_product_up); 
free(stored_product_down); 
return; 
}

/***********************************************************************************/

void propagate_forwards_kinetic_chemistry_real_multi_up(double *kinetic_full_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,MKL_Complex16 *det_overlap_up,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;
   MKL_Complex16 previous_overlap_up; 
   MKL_Complex16 weights_numerator, weights_denominator;
   MKL_Complex16 *stored_product_up; 

   stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) {

      /*Store Previous Overlaps*/
      previous_overlap_up = overlap_up[i];

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
      cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

      /*Get Product of Trial and Actual*/
      overlap_up[i] = compute_overlap_inverse_multi_phaseless_up(trial_wf_up_phaseless,&wf_up[i*ist.n_spatial_orbitals_n_up],trial_determinant_coefficients_phaseless,&overlap_inverse_up[i*ist.n_determinants_phaseless_n_up_sq],&det_overlap_up[i*ist.n_determinants_trial_phaseless],ist);

      /*Get Overall Weight*/
      weights_numerator = overlap_up[i]; 
      weights_denominator = previous_overlap_up; 
      weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));

    } /*weights walkers*/

   } /*i*/


free(stored_product_up);
return;
}

/*****************************************************************************************************************/

void propagate_forwards_kinetic_chemistry_real_multi_down(double *kinetic_full_fermions,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_down,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;
   MKL_Complex16 previous_overlap_down;  
   MKL_Complex16 weights_numerator, weights_denominator;
   MKL_Complex16 *stored_product_down;  

   stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) {

      /*Store Previous Overlaps*/
      previous_overlap_down = overlap_down[i];

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
      cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

      /*Get Product of Trial and Actual*/
      overlap_down[i] = compute_overlap_inverse_multi_phaseless_down(trial_wf_down_phaseless,&wf_down[i*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_down[i*ist.n_determinants_phaseless_n_down_sq],&det_overlap_down[i*ist.n_determinants_trial_phaseless],ist);

      /*Get Overall Weight*/
      weights_numerator = overlap_down[i]; 
      weights_denominator = previous_overlap_down; 
      weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));

    } /*weights walkers*/

   } /*i*/


free(stored_product_down);
return;
}

/***************************************************************************************/

void propagate_forwards_kinetic_chemistry_real_unrestricted(double *kinetic_full_fermions_up,double *kinetic_full_fermions_down,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;
   MKL_Complex16 previous_overlap_up, previous_overlap_down;
   MKL_Complex16 weights_numerator, weights_denominator;
   MKL_Complex16 *stored_product_up, *stored_product_down;

   stored_product_up=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
   stored_product_down=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

    /*Ensure that the Weights of the Walkers Is Not Zero*/
    if ( Cabs(weights[i]) != 0 ) {

      /*Store Previous Overlaps*/
      previous_overlap_up = overlap_up[i];
      previous_overlap_down = overlap_down[i];

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions_up,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
      cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

      /*Get Product of Trial and Actual*/
      overlap_up[i] = compute_overlap_inverse(trial_wf_up,&wf_up[i*ist.n_spatial_orbitals_n_up],&overlap_inverse_up[i*ist.n_up_sq],ist,0);

      /*Get Down Matrices*/
      dmat_cmat(kinetic_full_fermions_down,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
      cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

      /*Get Product of Trial and Down Actual*/
      overlap_down[i] = compute_overlap_inverse(trial_wf_down,&wf_down[i*ist.n_spatial_orbitals_n_down],&overlap_inverse_down[i*ist.n_down_sq],ist,1);

      /*Get Overall Weight*/
      weights_numerator = Cmul(overlap_up[i], overlap_down[i]);
      weights_denominator = Cmul(previous_overlap_up, previous_overlap_down);
      weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));

    } /*weights walkers*/

   } /*i*/

free(stored_product_up);
free(stored_product_down);
return;
}

/******************************************************************/

void propagate_half_backwards_kinetic_chemistry_real_multi_restricted(double *kinetic_backwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice - Determine Which Wavefunctions Must Be Propagated*/
    if ( ist.n_up > 0 ) {
      if ( ist.n_down > 0 ) {
        propagate_half_backwards_kinetic_chemistry_real_multi_restricted_both(kinetic_backwards_half_fermions,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist); 
      }
      else { 
        propagate_half_backwards_kinetic_chemistry_real_multi_up(kinetic_backwards_half_fermions,wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,overlap_inverse_up,weights,ist); 
      }
    }
    else {
      if ( ist.n_down > 0 ) {
       propagate_half_backwards_kinetic_chemistry_real_multi_down(kinetic_backwards_half_fermions,wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_down,overlap_inverse_down,weights,ist); 
      }
    } 
  
return;
}

/******************************************************************/

void propagate_half_backwards_kinetic_chemistry_real_multi_unrestricted(double *kinetic_backwards_half_fermions_up,double *kinetic_backwards_half_fermions_down,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice - Determine Which Wavefunctions Must Be Propagated*/
    if ( ist.n_up > 0 ) {
      if ( ist.n_down > 0 ) {
        propagate_half_backwards_kinetic_chemistry_real_multi_unrestricted_both(kinetic_backwards_half_fermions_up,kinetic_backwards_half_fermions_down,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);
      }
      else {
        propagate_half_backwards_kinetic_chemistry_real_multi_up(kinetic_backwards_half_fermions_up,wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,overlap_inverse_up,weights,ist);
      }
    }
    else {
      if ( ist.n_down > 0 ) {
       propagate_half_backwards_kinetic_chemistry_real_multi_down(kinetic_backwards_half_fermions_down,wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_down,overlap_inverse_down,weights,ist);
      }
    }

return;
}

/******************************************************************/

void propagate_half_backwards_kinetic_chemistry_real_multi_restricted_both(double *kinetic_backwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i, j, k;
    MKL_Complex16 previous_overlap;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap= overlap_total[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

       /*Get Product of Trial and Actual*/
       overlap_total[i] = compute_overlap_inverse_multi_phaseless_both(trial_wf_up_phaseless,trial_wf_down_phaseless,&wf_up[i*ist.n_spatial_orbitals_n_up],&wf_down[i*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_up[i*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[i*ist.n_determinants_phaseless_n_down_sq],&det_overlap_up[i*ist.n_determinants_trial_phaseless],&det_overlap_down[i*ist.n_determinants_trial_phaseless],ist);

       /*Get Overall Weights*/
       weights[i] = Cmul(weights[i], Cdiv(overlap_total[i], previous_overlap));


      } /*Non Zero Weights*/
     }

free(stored_product_up);
free(stored_product_down);
return;
}

/**************************************************************************************************/

void propagate_half_backwards_kinetic_chemistry_real_multi_unrestricted_both(double *kinetic_backwards_half_fermions_up,double *kinetic_backwards_half_fermions_down,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i, j, k;
    MKL_Complex16 previous_overlap; 
    MKL_Complex16 *stored_product_up, *stored_product_down;
    /*HERE*/

    stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap= overlap_total[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions_up,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions_down,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

       /*Get Product of Trial and Actual*/
       overlap_total[i] = compute_overlap_inverse_multi_phaseless_both(trial_wf_up_phaseless,trial_wf_down_phaseless,&wf_up[i*ist.n_spatial_orbitals_n_up],&wf_down[i*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_up[i*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[i*ist.n_determinants_phaseless_n_down_sq],&det_overlap_up[i*ist.n_determinants_trial_phaseless],&det_overlap_down[i*ist.n_determinants_trial_phaseless],ist); 

       /*Get Overall Weights*/
       weights[i] = Cmul(weights[i], Cdiv(overlap_total[i], previous_overlap)); 


      } /*Non Zero Weights*/
     }

free(stored_product_up);
free(stored_product_down);
return;
}

/******************************************************************/

void propagate_half_backwards_kinetic_chemistry_real_multi_up(double *kinetic_backwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i, j, k;
    MKL_Complex16 previous_overlap; 
    MKL_Complex16 *stored_product_up; 

    stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap = overlap_total[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Product of Trial and Actual*/
       overlap_total[i] = compute_overlap_inverse_multi_phaseless_up(trial_wf_up_phaseless,&wf_up[i*ist.n_spatial_orbitals_n_up],trial_determinant_coefficients_phaseless,&overlap_inverse_up[i*ist.n_determinants_phaseless_n_up_sq],&det_overlap_up[i*ist.n_determinants_trial_phaseless],ist);

       /*Get Overall Weights*/
       weights[i] = Cmul(weights[i], Cdiv(overlap_total[i], previous_overlap)); 

      } /*Non Zero Weights*/
     }

free(stored_product_up);
return;
}

/******************************************************************/

void propagate_half_backwards_kinetic_chemistry_real_multi_down(double *kinetic_backwards_half_fermions,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i, j, k;
    MKL_Complex16 previous_overlap; 
    MKL_Complex16 *stored_product_down;  

    stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap = overlap_total[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

       /*Get Product of Trial and Actual*/
       overlap_total[i] = compute_overlap_inverse_multi_phaseless_down(trial_wf_down_phaseless,&wf_down[i*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_down[i*ist.n_determinants_phaseless_n_down_sq],&det_overlap_down[i*ist.n_determinants_trial_phaseless],ist);

       /*Get Overall Weights*/
       weights[i] = Cmul(weights[i], Cdiv(overlap_total[i], previous_overlap)); 

      } /*Non Zero Weights*/
     }

free(stored_product_down);
return;
}

/***********************************************************************/

void propagate_half_backwards_kinetic_chemistry_real_restricted(double *kinetic_backwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 previous_overlap_up, previous_overlap_down;
    MKL_Complex16 weights_numerator, weights_denominator;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i];
       previous_overlap_down = overlap_down[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Product of Trial and Actual*/
       overlap_up[i] = compute_overlap_inverse(trial_wf_up,&wf_up[i*ist.n_spatial_orbitals_n_up],&overlap_inverse_up[i*ist.n_up_sq],ist,0);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

       /*Get Product of Trial and Actual*/
       overlap_down[i] = compute_overlap_inverse(trial_wf_down,&wf_down[i*ist.n_spatial_orbitals_n_down],&overlap_inverse_down[i*ist.n_down_sq],ist,1);

       /*Get Overall Weights*/
       weights_numerator = Cmul(overlap_up[i], overlap_down[i]);
       weights_denominator = Cmul(previous_overlap_up, previous_overlap_down);

       weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));

      } /*Non Zero Weights*/
     }

free(stored_product_up);
free(stored_product_down);
return;
}

/***********************************************************************/

void propagate_half_backwards_kinetic_chemistry_real_unrestricted(double *kinetic_backwards_half_fermions_up,double *kinetic_backwards_half_fermions_down,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 previous_overlap_up, previous_overlap_down;
    MKL_Complex16 weights_numerator, weights_denominator;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check Weights Non Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i];
       previous_overlap_down = overlap_down[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions_up,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Product of Trial and Actual*/
       overlap_up[i] = compute_overlap_inverse(trial_wf_up,&wf_up[i*ist.n_spatial_orbitals_n_up],&overlap_inverse_up[i*ist.n_up_sq],ist,0);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions_down,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

       /*Get Product of Trial and Actual*/
       overlap_down[i] = compute_overlap_inverse(trial_wf_down,&wf_down[i*ist.n_spatial_orbitals_n_down],&overlap_inverse_down[i*ist.n_down_sq],ist,1); 

       /*Get Overall Weights*/
       weights_numerator = Cmul(overlap_up[i], overlap_down[i]);
       weights_denominator = Cmul(previous_overlap_up, previous_overlap_down);
       weights[i] = Cmul(weights[i], Cdiv(weights_numerator, weights_denominator));

      } /*Non Zero Weights*/

     }

free(stored_product_up);
free(stored_product_down);
return;
}

/***********************************************************************/

void propagate_half_forwards_kinetic_chemistry_real_multi_unrestricted(double *kinetic_forwards_half_fermions_up,double *kinetic_forwards_half_fermions_down,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,int_st ist){

     /*Propagate a Half Step Forwards Using Kinetic Operator*/

     /*Choose Correct Spin*/
     if ( ist.n_up > 0 ) {
      if ( ist.n_down > 0 ) {
         propagate_half_forwards_kinetic_chemistry_real_multi_unrestricted_both(kinetic_forwards_half_fermions_up,kinetic_forwards_half_fermions_down,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,weights,ist); 
      }
      else {
         propagate_half_forwards_kinetic_chemistry_real_multi_up(kinetic_forwards_half_fermions_up,wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,det_overlap_up,weights,ist); 
      }
     }
     else {
         propagate_half_forwards_kinetic_chemistry_real_multi_down(kinetic_forwards_half_fermions_down,wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_down,det_overlap_down,weights,ist); 
     }   

return;
} 

/***********************************************************************/

void propagate_half_forwards_kinetic_chemistry_real_multi_restricted(double *kinetic_forwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,int_st ist){

     /*Propagate a Half Step Forwards Using Kinetic Operator*/

     /*Choose Correct Spin*/
     if ( ist.n_up > 0 ) {
      if ( ist.n_down > 0 ) {
         propagate_half_forwards_kinetic_chemistry_real_multi_restricted_both(kinetic_forwards_half_fermions,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,weights,ist);
      }
      else {
         propagate_half_forwards_kinetic_chemistry_real_multi_up(kinetic_forwards_half_fermions,wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,det_overlap_up,weights,ist);
      }
     }
     else {
         propagate_half_forwards_kinetic_chemistry_real_multi_down(kinetic_forwards_half_fermions,wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_down,det_overlap_down,weights,ist);
     }

return;
}

/***********************************************************************/

void propagate_half_forwards_kinetic_chemistry_real_multi_restricted_both(double *kinetic_forwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,int_st ist){


    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/
    int i;
    MKL_Complex16 previous_overlap_up; 
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_total[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

       /*Get Product of Trial and Actual*/
       overlap_total[i] = compute_overlap_inverse_multi_phaseless_both(trial_wf_up_phaseless,trial_wf_down_phaseless,&wf_up[i*ist.n_spatial_orbitals_n_up],&wf_down[i*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_up[i*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[i*ist.n_determinants_phaseless_n_down_sq],&det_overlap_up[i*ist.n_determinants_trial_phaseless],&det_overlap_down[i*ist.n_determinants_trial_phaseless],ist); 

       /*Get Overall Weights*/
       weights[i] = Cmul(weights[i], Cdiv(overlap_total[i], previous_overlap_up)); 

      }
     }

free(stored_product_up);
free(stored_product_down);
return;
}

/****************************************************************************************/

void propagate_half_forwards_kinetic_chemistry_real_multi_unrestricted_both(double *kinetic_forwards_half_fermions_up,double *kinetic_forwards_half_fermions_down,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,int_st ist){


    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/
    int i;
    MKL_Complex16 previous_overlap_up;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_total[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions_up,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions_down,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

       /*Get Product of Trial and Actual*/
       overlap_total[i] = compute_overlap_inverse_multi_phaseless_both(trial_wf_up_phaseless,trial_wf_down_phaseless,&wf_up[i*ist.n_spatial_orbitals_n_up],&wf_down[i*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_up[i*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[i*ist.n_determinants_phaseless_n_down_sq],&det_overlap_up[i*ist.n_determinants_trial_phaseless],&det_overlap_down[i*ist.n_determinants_trial_phaseless],ist);

       /*Get Overall Weights*/
       weights[i] = Cmul(weights[i], Cdiv(overlap_total[i], previous_overlap_up));

      }
     }

free(stored_product_up);
free(stored_product_down);
return; 
}

/******************************************************************/

void propagate_half_forwards_kinetic_chemistry_real_multi_up(double *kinetic_forwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *det_overlap_up,MKL_Complex16 *weights,int_st ist){


    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/
    int i;
    MKL_Complex16 weight_numerator, weight_denominator;
    MKL_Complex16 previous_overlap_up; 
    MKL_Complex16 *stored_product_up; 

    stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Product of Trial and Actual*/
       overlap_up[i] = compute_overlap_inverse_multi_phaseless_up(trial_wf_up_phaseless,&wf_up[i*ist.n_spatial_orbitals_n_up],trial_determinant_coefficients_phaseless,&overlap_inverse_up[i*ist.n_determinants_phaseless_n_up_sq],&det_overlap_up[i*ist.n_determinants_trial_phaseless],ist);

       /*Get Overall Weights*/
       weight_numerator = overlap_up[i]; 
       weight_denominator = previous_overlap_up; 

       weights[i] = Cmul(weights[i], Cdiv(weight_numerator, weight_denominator));

      }

     }

free(stored_product_up);
return;
}
/**********************************************************************************/

void propagate_half_forwards_kinetic_chemistry_real_multi_down(double *kinetic_forwards_half_fermions,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/
    int i;
    MKL_Complex16 weight_numerator, weight_denominator;
    MKL_Complex16 previous_overlap_down;  
    MKL_Complex16 *stored_product_down;  

    stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_down = overlap_down[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

       /*Get Product of Trial and Actual*/
       overlap_down[i] = compute_overlap_inverse_multi_phaseless_down(trial_wf_down_phaseless,&wf_down[i*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_down[i*ist.n_determinants_phaseless_n_down_sq],&det_overlap_down[i*ist.n_determinants_trial_phaseless],ist);

       /*Get Overall Weights*/
       weight_numerator = overlap_down[i]; 
       weight_denominator = previous_overlap_down; 

       weights[i] = Cmul(weights[i], Cdiv(weight_numerator, weight_denominator));

      }
     }

free(stored_product_down);
return;
}
/**********************************************************************************/

void propagate_half_forwards_kinetic_chemistry_real_restricted(double *kinetic_forwards_half_fermions,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 weight_numerator, weight_denominator;
    MKL_Complex16 previous_overlap_up, previous_overlap_down;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {

       /*Previous Overlaps*/
       previous_overlap_up = overlap_up[i];
       previous_overlap_down = overlap_down[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Product of Trial and Actual*/
       overlap_up[i] = compute_overlap_inverse(trial_wf_up,&wf_up[i*ist.n_spatial_orbitals_n_up],&overlap_inverse_up[i*ist.n_up_sq],ist,0);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

       /*Get Product of Trial and Actual*/
       overlap_down[i] = compute_overlap_inverse(trial_wf_down,&wf_down[i*ist.n_spatial_orbitals_n_down],&overlap_inverse_down[i*ist.n_down_sq],ist,1);

       /*Get Overall Weights*/
       weight_numerator = Cmul(overlap_up[i], overlap_down[i]);
       weight_denominator = Cmul(previous_overlap_up, previous_overlap_down);

       weights[i] = Cmul(weights[i], Cdiv(weight_numerator, weight_denominator));

      }

     }

free(stored_product_up);
free(stored_product_down);
return;
}

/*********************************************************************************/

void propagate_half_forwards_kinetic_chemistry_real_unrestricted(double *kinetic_forwards_half_fermions_up,double *kinetic_forwards_half_fermions_down,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int_st ist){

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 weight_numerator, weight_denominator;
    MKL_Complex16 previous_overlap_up, previous_overlap_down;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

     /*Check if Weights Are Zero*/
     if ( Cabs(weights[i]) != 0 ) {/*Previous Overlaps*/
       previous_overlap_up = overlap_up[i];
       previous_overlap_down = overlap_down[i];

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions_up,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Product of Trial and Actual*/
       overlap_up[i] = compute_overlap_inverse(trial_wf_up,&wf_up[i*ist.n_spatial_orbitals_n_up],&overlap_inverse_up[i*ist.n_up_sq],ist,0);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions_down,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&wf_down[i*ist.n_spatial_orbitals_n_down],1);

       /*Get Product of Trial and Actual*/
       overlap_down[i] = compute_overlap_inverse(trial_wf_down,&wf_down[i*ist.n_spatial_orbitals_n_down],&overlap_inverse_down[i*ist.n_down_sq],ist,1);

       /*Get Overall Weights*/
       weight_numerator = Cmul(overlap_up[i], overlap_down[i]);
       weight_denominator = Cmul(previous_overlap_up, previous_overlap_down);
       weights[i] = Cmul(weights[i], Cdiv(weight_numerator, weight_denominator));

      }

     }

free(stored_product_up);
free(stored_product_down);
return;
}

#include "afqmc.h"

/*****************************************************************************/

void orthogonalize_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,MKL_Complex16 *R_up,MKL_Complex16 *R_down,int_st ist){

    /*Performs a Gram-Schmidt Orthogonalization on WF To Prevent OverFlow*/
    int walker;
    int n_sites_n_up = ist.n_sites * ist.n_up; 
    int n_sites_n_down = ist.n_sites * ist.n_down; 
    MKL_Complex16 determinant; 

    for (walker = 0; walker < ist.n_walkers; walker++) {

      /*Check Non-Zero Weight*/
      if ( Cabs(weights[walker]) != 0 ) {

         complex_modified_gram_schmidt(&wf_up[walker*n_sites_n_up],R_up,ist.n_sites,ist.n_up); 
         complex_det(R_up,ist.n_up,&determinant);
         overlap_up[walker] = Cdiv(overlap_up[walker], determinant); 

         complex_modified_gram_schmidt(&wf_down[walker*n_sites_n_down],R_down,ist.n_sites,ist.n_down);
         complex_det(R_down,ist.n_down,&determinant); 
         overlap_down[walker] = Cdiv(overlap_down[walker], determinant);

      } 

    }  
return; 
}

/*****************************************************************************/

void orthogonalize_fermions_up(MKL_Complex16 *wf_up,MKL_Complex16 *overlap_up,MKL_Complex16 *weights,MKL_Complex16 *R_up,int_st ist){

    /*Performs a Gram-Schmidt Orthogonalization on WF To Prevent OverFlow*/
    int walker;
    int n_sites_n_up = ist.n_sites * ist.n_up;
    MKL_Complex16 determinant;

    for (walker = 0; walker < ist.n_walkers; walker++) {

      /*Check Non-Zero Weight*/
      if ( Cabs(weights[walker]) != 0 ) {

         complex_modified_gram_schmidt(&wf_up[walker*n_sites_n_up],R_up,ist.n_sites,ist.n_up);
         complex_det(R_up,ist.n_up,&determinant);
         overlap_up[walker] = Cdiv(overlap_up[walker], determinant);

      }

    }
return;
}

/*****************************************************************************/

void orthogonalize_fermions_down(MKL_Complex16 *wf_down,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,MKL_Complex16 *R_down,int_st ist){

    /*Performs a Gram-Schmidt Orthogonalization on WF To Prevent OverFlow*/
    int walker;
    int n_sites_n_down = ist.n_sites * ist.n_down;
    MKL_Complex16 determinant;

    for (walker = 0; walker < ist.n_walkers; walker++) {

      /*Check Non-Zero Weight*/
      if ( Cabs(weights[walker]) != 0 ) {

         complex_modified_gram_schmidt(&wf_down[walker*n_sites_n_down],R_down,ist.n_sites,ist.n_down);
         complex_det(R_down,ist.n_down,&determinant);
         overlap_down[walker] = Cdiv(overlap_down[walker], determinant);

      }

    }
return;
}

/*****************************************************************************/

void orthogonalize_chemistry(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,MKL_Complex16 *R_up,MKL_Complex16 *R_down,int_st ist){

    /*Performs a Gram-Schmidt Orthogonalization on WF To Prevent OverFlow*/
    int walker;
    int n_spatial_orbitals_n_up = ist.n_spatial_orbitals * ist.n_up;
    int n_spatial_orbitals_n_down = ist.n_spatial_orbitals * ist.n_down;
    MKL_Complex16 determinant;

    for (walker = 0; walker < ist.n_walkers; walker++) {

      /*Check Non-Zero Weight*/
      if ( Cabs(weights[walker]) != 0 ) {

         complex_modified_gram_schmidt(&wf_up[walker*n_spatial_orbitals_n_up],R_up,ist.n_spatial_orbitals,ist.n_up);
         complex_det(R_up,ist.n_up,&determinant);
         overlap_up[walker] = Cdiv(overlap_up[walker], determinant);

         complex_modified_gram_schmidt(&wf_down[walker*n_spatial_orbitals_n_down],R_down,ist.n_spatial_orbitals,ist.n_down);
         complex_det(R_down,ist.n_down,&determinant);
         overlap_down[walker] = Cdiv(overlap_down[walker], determinant);

      }

    }
return;
}

/*****************************************************************************/

void orthogonalize_chemistry_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *R_up,MKL_Complex16*R_down,int_st ist){

    /*Performs a Gram-Schmidt Orthogonalization on WF To Prevent OverFlow*/
    /*Check Spin Situation*/
    if ( ist.n_up > 0 ) {
      if ( ist.n_down > 0 ) {
         orthogonalize_chemistry_multi_both(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,weights,R_up,R_down,ist); 
      }
      else {
          orthogonalize_chemistry_multi_up(wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,det_overlap_up,weights,R_up,ist); 
      }
    }
    else {
          orthogonalize_chemistry_multi_down(wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_down,det_overlap_down,weights,R_down,ist); 
    }

return;
}

/****************************************************************************/

void orthogonalize_chemistry_multi_both(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *R_up,MKL_Complex16*R_down,int_st ist){

    /*Performs a Gram-Schmidt Orthogonalization on WF To Prevent OverFlow*/
    int walker, i, j;
    int n_spatial_orbitals_n_up = ist.n_spatial_orbitals * ist.n_up;
    int n_spatial_orbitals_n_down = ist.n_spatial_orbitals * ist.n_down;
    MKL_Complex16 determinant;

    for (walker = 0; walker < ist.n_walkers; walker++) {

      /*Check Non-Zero Weight*/
      if ( Cabs(weights[walker]) != 0 ) {

         complex_modified_gram_schmidt(&wf_up[walker*n_spatial_orbitals_n_up],R_up,ist.n_spatial_orbitals,ist.n_up);
         complex_modified_gram_schmidt(&wf_down[walker*n_spatial_orbitals_n_down],R_down,ist.n_spatial_orbitals,ist.n_down);

         overlap_total[walker] = compute_overlap_inverse_multi_phaseless_both(trial_wf_up_phaseless,trial_wf_down_phaseless,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

      }

    }
return;
}

/****************************************************************************/

void orthogonalize_chemistry_multi_up(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *det_overlap_up,MKL_Complex16 *weights,MKL_Complex16 *R_up,int_st ist){

    /*Performs a Gram-Schmidt Orthogonalization on WF To Prevent OverFlow*/
    int walker;
    int n_spatial_orbitals_n_up = ist.n_spatial_orbitals * ist.n_up;
    MKL_Complex16 determinant;

    for (walker = 0; walker < ist.n_walkers; walker++) {

      /*Check Non-Zero Weight*/
      if ( Cabs(weights[walker]) != 0 ) {

         complex_modified_gram_schmidt(&wf_up[walker*n_spatial_orbitals_n_up],R_up,ist.n_spatial_orbitals,ist.n_up);
         overlap_up[walker] = compute_overlap_inverse_multi_phaseless_up(trial_wf_up_phaseless,&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],ist);

      }

    }
return;
}

/****************************************************************************/

void orthogonalize_chemistry_multi_down(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *R_down,int_st ist){

    /*Performs a Gram-Schmidt Orthogonalization on WF To Prevent OverFlow*/
    int walker;
    int n_spatial_orbitals_n_down = ist.n_spatial_orbitals * ist.n_down;
    MKL_Complex16 determinant;

    for (walker = 0; walker < ist.n_walkers; walker++) {

      /*Check Non-Zero Weight*/
      if ( Cabs(weights[walker]) != 0 ) {

         complex_modified_gram_schmidt(&wf_down[walker*n_spatial_orbitals_n_down],R_down,ist.n_spatial_orbitals,ist.n_down);
         overlap_down[walker] = compute_overlap_inverse_multi_phaseless_down(trial_wf_down_phaseless,&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_determinant_coefficients_phaseless,&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

      }

    }
return;
}

/****************************************************************************/

void orthogonalize_fermions_spin_polarized(MKL_Complex16 *wf_up,MKL_Complex16 *overlap_up,MKL_Complex16 *weights,MKL_Complex16 *R_up,int_st ist){

    /*Performs a Gram-Schmidt Orthogonalization on WF To Prevent OverFlow*/
    int walker;
    int n_sites_n_up = ist.n_sites * ist.n_up;
    MKL_Complex16 determinant;

    for (walker = 0; walker < ist.n_walkers; walker++) {

      /*Check Non-Zero Weight*/
      if ( Cabs(weights[walker]) != 0 ) {

         complex_modified_gram_schmidt(&wf_up[walker*n_sites_n_up],R_up,ist.n_sites,ist.n_up);
         complex_det(R_up,ist.n_up,&determinant);
         overlap_up[walker] = Cdiv(overlap_up[walker], determinant);

      }

    }
return;
}

/*********************************************************************************/

void normalize_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *weights,int_st ist) {
 
     /*Normalize the Single Orbital Wavefunctions*/
     int walker; 
     int i, iw; 
     MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0; 
     MKL_Complex16 normalization; 

     for (walker = 0; walker < ist.n_walkers; walker++) {

        iw = walker * ist.n_sites; 
 
        /*Check for non-Zero Weight*/
        if ( Cabs(weights[walker]) != 0 ) {

            /*Find Normalization Constant*/
            normalization = Zero; 
            for (i=0; i<ist.n_sites; i++) {
              normalization = Cadd(normalization, Cmul(conjugate(wf_bosons[iw+i]), wf_bosons[iw+i]));   
            }
            normalization = Csqrt(normalization); 
 
            /*Divide by Constant*/
            for (i=0; i<ist.n_sites; i++) {
             wf_bosons[iw+i] = Cdiv(wf_bosons[iw+i], normalization); 
           }          
           weights[walker] = Cmul(weights[walker], normalization); 

        } /*weights check*/
     } /*walkers*/     

return; 
}

/*********************************************************************************/ 

void normalize_fermions_spin_polarized(MKL_Complex16 *wf_fermions,MKL_Complex16 *weights,int_st ist) {

     /*Normalize the Single Orbital Wavefunctions*/
     int walker;
     int i, j, iw; 
     MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
     MKL_Complex16 normalization;

     for (walker = 0; walker < ist.n_walkers; walker++) {

        iw = walker * ist.n_sites;

        /*Check for non-Zero Weight*/
        if ( Cabs(weights[walker]) != 0 ) {

            for (j=0; j<ist.n_up; j++) {

               /*Find Normalization Constant*/
               normalization = Zero;
               for (i=0; i<ist.n_sites; i++) {
                 normalization = Cadd(normalization, Cmul(conjugate(wf_fermions[iw+i*ist.n_up+j]), wf_fermions[iw+i]));
               }
               normalization = Csqrt(normalization);

               /*Divide by Constant*/
               for (i=0; i<ist.n_sites; i++) {
                wf_fermions[iw+i*ist.n_up+j] = Cdiv(wf_fermions[iw+i*ist.n_up+j], normalization);
              }
              weights[walker] = Cmul(weights[walker], normalization);

            }

        } /*weights check*/
     } /*walkers*/

return; 
}

/*********************************************************************************************/

void normalize_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *weights,int_st ist) {

   /*Normalize Fermions for Spin Up and Down*/
   int walker; 
   int i, j; 
   int iw;
   MKL_Complex16 Zero; Zero.real = 0.0;  Zero.imag = 0.0;  
   MKL_Complex16 normalization; 

   for (walker = 0; walker < ist.n_walkers; walker++) {

        iw = walker * ist.n_sites;

        /*Check for non-Zero Weight*/
        if ( Cabs(weights[walker]) != 0 ) {

            for (j=0; j<ist.n_up; j++) {

               /*Find Normalization Constant*/
               normalization = Zero;
               for (i=0; i<ist.n_sites; i++) {
                 normalization = Cadd(normalization, Cmul(conjugate(wf_up[iw+i*ist.n_up+j]), wf_up[iw+i*ist.n_up+j]));
               }
               normalization = Csqrt(normalization);

               /*Divide by Constant*/
               for (i=0; i<ist.n_sites; i++) {
                wf_up[iw+i*ist.n_up+j] = Cdiv(wf_up[iw+i*ist.n_up+j], normalization);
              }
              weights[walker] = Cmul(weights[walker], normalization);

            }

            /*Spin Down*/
            for (j=0; j<ist.n_down; j++) {

               /*Find Normalization Constant*/
               normalization = Zero;
               for (i=0; i<ist.n_sites; i++) {
                 normalization = Cadd(normalization, Cmul(conjugate(wf_down[iw+i*ist.n_down+j]), wf_down[iw+i*ist.n_down+j]));
               }
               normalization = Csqrt(normalization);

               /*Divide by Constant*/
               for (i=0; i<ist.n_sites; i++) {
                wf_down[iw+i*ist.n_down+j] = Cdiv(wf_down[iw+i*ist.n_down+j], normalization);
              }
              weights[walker] = Cmul(weights[walker], normalization);

           }

        } /*weights check*/
     } /*walkers*/ 

return; 
}

/*********************************************************************************************/

void normalize_chemistry(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *weights,int_st ist) {

   /*Normalize Fermions for Spin Up and Down*/
   int walker;
   int i, j;
   int iw;
   MKL_Complex16 Zero; Zero.real = 0.0;  Zero.imag = 0.0;
   MKL_Complex16 normalization;

   for (walker = 0; walker < ist.n_walkers; walker++) {

        iw = walker * ist.n_spatial_orbitals;

        /*Check for non-Zero Weight*/
        if ( Cabs(weights[walker]) != 0 ) {

            for (j=0; j<ist.n_up; j++) {

               /*Find Normalization Constant*/
               normalization = Zero;
               for (i=0; i<ist.n_spatial_orbitals; i++) {
                 normalization = Cadd(normalization, Cmul(conjugate(wf_up[iw+i*ist.n_up+j]), wf_up[iw+i*ist.n_up+j]));
               }
               normalization = Csqrt(normalization);

               /*Divide by Constant*/
               for (i=0; i<ist.n_spatial_orbitals; i++) {
                wf_up[iw+i*ist.n_up+j] = Cdiv(wf_up[iw+i*ist.n_up+j], normalization);
              }
              weights[walker] = Cmul(weights[walker], normalization);

            }

            /*Spin Down*/
            for (j=0; j<ist.n_down; j++) {

               /*Find Normalization Constant*/
               normalization = Zero;
               for (i=0; i<ist.n_spatial_orbitals; i++) {
                 normalization = Cadd(normalization, Cmul(conjugate(wf_down[iw+i*ist.n_down+j]), wf_down[iw+i*ist.n_down+j]));
               }
               normalization = Csqrt(normalization);

               /*Divide by Constant*/
               for (i=0; i<ist.n_spatial_orbitals; i++) {
                wf_down[iw+i*ist.n_down+j] = Cdiv(wf_down[iw+i*ist.n_down+j], normalization);
              }
              weights[walker] = Cmul(weights[walker], normalization);

           }

        } /*weights check*/
     } /*walkers*/

return;
}

/***************************************************************************************************************/

void get_normalized_wf_up(MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *wf_up,int_st ist){

    /*Find the Up Wavefunction Normalized Through Columns*/
    int i, j; 
    MKL_Complex16 normalization; 
    MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;    
 
    for (j=0; j<ist.n_up; j++) {

      /*Find Normalization Constant*/
      normalization = Zero;
      for (i=0; i<ist.n_sites; i++) {
         normalization = Cadd(normalization, Cmul(conjugate(wf_up[i*ist.n_up+j]), wf_up[i*ist.n_up+j]));
      }
      normalization = Csqrt(normalization);

      /*Divide by Constant*/
      for (i=0; i<ist.n_sites; i++) {
        normalized_wave_function_up[i*ist.n_up+j] = Cdiv(wf_up[i*ist.n_up+j], normalization);
      }   
   }

return; 
}

/******************************************************************************************************************/

void get_normalized_wf_up_chemistry(MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *wf_up,int_st ist){

    /*Find the Up Wavefunction Normalized Through Columns*/
    int i, j;
    MKL_Complex16 normalization;
    MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;

    for (j=0; j<ist.n_up; j++) {

      /*Find Normalization Constant*/
      normalization = Zero;
      for (i=0; i<ist.n_spatial_orbitals; i++) {
         normalization = Cadd(normalization, Cmul(conjugate(wf_up[i*ist.n_up+j]), wf_up[i*ist.n_up+j]));
      }
      normalization = Csqrt(normalization);

      /*Divide by Constant*/
      for (i=0; i<ist.n_spatial_orbitals; i++) {
        normalized_wave_function_up[i*ist.n_up+j] = Cdiv(wf_up[i*ist.n_up+j], normalization);
      }
   }

return;
}

/*******************************************************************************************************************/

void get_normalized_wf_down(MKL_Complex16 *normalized_wave_function_down,MKL_Complex16 *wf_down,int_st ist){

   /*Find the Up Wavefunction Normalized Through Columns*/  
    int i, j;    
    MKL_Complex16 normalization;
    MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0; 

    for (j=0; j<ist.n_down; j++) {

      /*Find Normalization Constant*/
      normalization = Zero;
      for (i=0; i<ist.n_sites; i++) {
         normalization = Cadd(normalization, Cmul(conjugate(wf_down[i*ist.n_down+j]), wf_down[i*ist.n_down+j]));
      }
      normalization = Csqrt(normalization);

      /*Divide by Constant*/
      for (i=0; i<ist.n_sites; i++) {
        normalized_wave_function_down[i*ist.n_down+j] = Cdiv(wf_down[i*ist.n_down+j], normalization);
      }        
   }

return; 
}

/*******************************************************************************************************************/

void get_normalized_wf_down_chemistry(MKL_Complex16 *normalized_wave_function_down,MKL_Complex16 *wf_down,int_st ist){

   /*Find the Up Wavefunction Normalized Through Columns*/
    int i, j;
    MKL_Complex16 normalization;
    MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;

    for (j=0; j<ist.n_down; j++) {

      /*Find Normalization Constant*/
      normalization = Zero;
      for (i=0; i<ist.n_spatial_orbitals; i++) {
         normalization = Cadd(normalization, Cmul(conjugate(wf_down[i*ist.n_down+j]), wf_down[i*ist.n_down+j]));
      }
      normalization = Csqrt(normalization);

      /*Divide by Constant*/
      for (i=0; i<ist.n_spatial_orbitals; i++) {
        normalized_wave_function_down[i*ist.n_down+j] = Cdiv(wf_down[i*ist.n_down+j], normalization);
      }
   }

return;
}


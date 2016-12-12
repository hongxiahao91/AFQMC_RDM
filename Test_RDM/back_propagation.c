/*Includes All Back Propagation Functions*/ 

#include "afqmc.h"

/***************************************************************************************************************/

void back_propagation_chemistry_restricted(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,MKL_Complex16 *back_propagation_matrices_up,MKL_Complex16 *back_propagation_matrices_down,double *kinetic_backwards_half_fermions_real_up,double *kinetic_full_fermions_real_up,double *kinetic_forwards_half_fermions_real_up,int_st ist,cns_st cns) {

     /*Get the Back Propagated Wavefunctions for Cases of Chemistry*/
     int steps; 

     /*First Initialize Back Propagated WFs*/
     init_bp_wf_restricted(trial_wf_up,trial_wf_down,bp_wf_up,bp_wf_down,ist); 

     /*Propagate Half Step Backwards So That You Can use Full Propagator Going Forward*/
     propagate_half_backwards_kinetic_chemistry_real_restricted_bp(kinetic_backwards_half_fermions_real_up,bp_wf_up,bp_wf_down,ist);

     /*Now Walk Matrices Backwards According to Weights*/
     for (steps = 0; steps < ist.n_steps_back_propagation; steps++ ) { 

         /*Propagate Forwards*/
         propagate_forwards_kinetic_chemistry_real_restricted_bp(kinetic_full_fermions_real_up,bp_wf_up,bp_wf_down,ist);

         /*Propagate Forwards By Potential Matrices*/
         propagate_forwards_potential_chemistry_bp(&back_propagation_matrices_up[steps*cns.max_number_walkers*ist.n_spatial_orbitals_sq],&back_propagation_matrices_down[steps*cns.max_number_walkers*ist.n_spatial_orbitals_sq],bp_wf_up,bp_wf_down,ist,cns);  

     }

    /*Propagate Forwards One More Half Step*/
    propagate_half_forwards_kinetic_chemistry_real_restricted_bp(kinetic_forwards_half_fermions_real_up,bp_wf_up,bp_wf_down,ist);

return; 
}

/*****************************************************************************************************************/

void back_propagation_chemistry_multi(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *trial_determinant_coefficients,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,MKL_Complex16 *back_propagation_matrices_up,MKL_Complex16 *back_propagation_matrices_down,double *kinetic_backwards_half_fermions_real_up,double *kinetic_full_fermions_real_up,double *kinetic_forwards_half_fermions_real_up,int_st ist,cns_st cns) {

     /*Get the Back Propagated Wavefunctions for Cases of Chemistry*/
     int steps;

     /*First Initialize Back Propagated WFs*/
     init_bp_wf_multi(trial_wf_up,trial_wf_down,trial_determinant_coefficients,bp_wf_up,bp_wf_down,ist);

     /*Propagate Half Step Backwards So That You Can use Full Propagator Going Forward*/
     propagate_half_backwards_kinetic_chemistry_real_restricted_bp(kinetic_backwards_half_fermions_real_up,bp_wf_up,bp_wf_down,ist);

     /*Now Walk Matrices Backwards According to Weights*/
     for (steps = 0; steps < ist.n_steps_back_propagation; steps++ ) {

         /*Propagate Forwards*/
         propagate_forwards_kinetic_chemistry_real_restricted_bp(kinetic_full_fermions_real_up,bp_wf_up,bp_wf_down,ist);

         /*Propagate Forwards By Potential Matrices*/
         propagate_forwards_potential_chemistry_bp(&back_propagation_matrices_up[steps*cns.max_number_walkers*ist.n_spatial_orbitals_sq],&back_propagation_matrices_down[steps*cns.max_number_walkers*ist.n_spatial_orbitals_sq],bp_wf_up,bp_wf_down,ist,cns);

     }

    /*Propagate Forwards One More Half Step*/
    propagate_half_forwards_kinetic_chemistry_real_restricted_bp(kinetic_forwards_half_fermions_real_up,bp_wf_up,bp_wf_down,ist);

return;
}

/*****************************************************************************************************************/

void back_propagation_chemistry_unrestricted(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,MKL_Complex16 *back_propagation_matrices_up,MKL_Complex16 *back_propagation_matrices_down,double *kinetic_backwards_half_fermions_real_up,double *kinetic_full_fermions_real_up,double *kinetic_forwards_half_fermions_real_up,double *kinetic_backwards_half_fermions_real_down,double *kinetic_full_fermions_real_down,double *kinetic_forwards_half_fermions_real_down,int_st ist,cns_st cns) {

     /*Get the Back Propagated Wavefunctions for Cases of Chemistry with UHF Wave functions*/
     int steps;

     /*First Initialize Back Propagated WFs*/
     init_bp_wf_restricted(trial_wf_up,trial_wf_down,bp_wf_up,bp_wf_down,ist);

     /*Propagate Half Step Backwards So That You Can use Full Propagator Going Forward*/
     propagate_half_backwards_kinetic_chemistry_real_unrestricted_bp(kinetic_backwards_half_fermions_real_up,kinetic_backwards_half_fermions_real_down,bp_wf_up,bp_wf_down,ist);

     /*Now Walk Matrices Backwards According to Weights*/
     for (steps = 0; steps < ist.n_steps_back_propagation; steps++ ) {

         /*Propagate Forwards*/
         propagate_forwards_kinetic_chemistry_real_unrestricted_bp(kinetic_full_fermions_real_up,kinetic_full_fermions_real_down,bp_wf_up,bp_wf_down,ist);

         /*Propagate Forwards By Potential Matrices*/
         propagate_forwards_potential_chemistry_bp(&back_propagation_matrices_up[steps*cns.max_number_walkers*ist.n_spatial_orbitals_sq],&back_propagation_matrices_down[steps*cns.max_number_walkers*ist.n_spatial_orbitals_sq],bp_wf_up,bp_wf_down,ist,cns);

     }

    /*Propagate Forwards One More Half Step*/
    propagate_half_forwards_kinetic_chemistry_real_unrestricted_bp(kinetic_forwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_down,bp_wf_up,bp_wf_down,ist);

return;
}

/***************************************************************************************************************/

void init_bp_wf_restricted(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,int_st ist) {

     /*Set the Initial BP  WF to the Trial WF*/ 
     int walker; 

     for (walker = 0; walker < ist.n_walkers; walker++) {
        cblas_zcopy(ist.n_spatial_orbitals_n_up,trial_wf_up,1,&bp_wf_up[walker*ist.n_spatial_orbitals_n_up],1); 
        cblas_zcopy(ist.n_spatial_orbitals_n_down,trial_wf_down,1,&bp_wf_down[walker*ist.n_spatial_orbitals_n_down],1);
     } 
    
return; 
}

/***************************************************************************************************************/

void init_bp_wf_multi(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *trial_determinant_coefficients,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,int_st ist) {

     /*Set the Initial BP  WF to the Trial WF*/
     int walker, i; 
     int largest_hold, largest_hold_up, largest_hold_down;  
     double largest_value = 0; 

     /*Find Largest Determinant*/
     for (i=0; i<ist.n_determinants_trial_energy; i++) {
       if ( Cabs(trial_determinant_coefficients[i]) > largest_value) {
          largest_value = Cabs(trial_determinant_coefficients[i]);
          largest_hold = i; 
       }
     }    

     /*Copy Largest Determinant Into Trial WF*/
     largest_hold_up = largest_hold * ist.n_spatial_orbitals_n_up; 
     largest_hold_down = largest_hold * ist.n_spatial_orbitals_n_down; 
     for (walker = 0; walker < ist.n_walkers; walker++) {
        cblas_zcopy(ist.n_spatial_orbitals_n_up,&trial_wf_up[largest_hold_up],1,&bp_wf_up[walker*ist.n_spatial_orbitals_n_up],1);
        cblas_zcopy(ist.n_spatial_orbitals_n_down,&trial_wf_down[largest_hold_down],1,&bp_wf_down[walker*ist.n_spatial_orbitals_n_down],1);
     }

return;
}

/***************************************************************************************************************/

void propagate_half_backwards_kinetic_chemistry_real_restricted_bp(double *kinetic_backwards_half_fermions_real_up,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,int_st ist) {

    /*Propagate a Half of a Step Backwards The BP WF*/
    int i;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions_real_up,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions_real_up,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],1);

    }

free(stored_product_up);
free(stored_product_down);
return;
}

/*************************************************************************************************************/

void propagate_forwards_kinetic_chemistry_real_unrestricted_bp(double *kinetic_full_fermions_up,double *kinetic_full_fermions_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;
   MKL_Complex16 *stored_product_up, *stored_product_down;

   stored_product_up=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
   stored_product_down=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions_up,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
      cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],1);

      /*Get Down Matrices*/
      dmat_cmat(kinetic_full_fermions_down,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
      cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],1);

   } /*i*/

free(stored_product_up);
free(stored_product_down);
return;
}

/************************************************************************************************************/

void propagate_forwards_kinetic_chemistry_real_restricted_bp(double *kinetic_full_fermions,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,int_st ist) {

   /*Moves the Trial Wavefunction One Half Kinetic Step Backwards*/

   int i;
   MKL_Complex16 *stored_product_up, *stored_product_down;

   stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
   stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

   /*Run Through Walkers*/
   for (i=0; i<ist.n_walkers; i++) {

      /*First Get Up Matrices*/
      dmat_cmat(kinetic_full_fermions,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
      cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],1);

      /*Get Down Matrices*/
      dmat_cmat(kinetic_full_fermions,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
      cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],1);

   } /*i*/


free(stored_product_up);
free(stored_product_down);
return;
}

/**********************************************************************************************************/

void propagate_half_forwards_kinetic_chemistry_real_restricted_bp(double *kinetic_forwards_half_fermions,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,int_st ist) {

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],1);

    }

free(stored_product_up);
free(stored_product_down);
return;
}

/**********************************************************************************************************/

void propagate_half_backwards_kinetic_chemistry_real_unrestricted_bp(double *kinetic_backwards_half_fermions_up,double *kinetic_backwards_half_fermions_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,int_st ist) {

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions_up,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_backwards_half_fermions_down,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],1);

    }

free(stored_product_up);
free(stored_product_down);
return;
}

/**********************************************************************************************************/

void propagate_half_forwards_kinetic_chemistry_real_unrestricted_bp(double *kinetic_forwards_half_fermions_up,double *kinetic_forwards_half_fermions_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,int_st ist) {

    /*Propagate the Wavefunction Forward By One Full Kinetic Time Slice*/

    int i;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    stored_product_up=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    stored_product_down=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

    /*Run Through Walkers*/
    for (i=0; i<ist.n_walkers; i++) {

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions_up,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions_down,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],1);

    }

free(stored_product_up);
free(stored_product_down);
}

/**********************************************************************************************************/

void copy_back_propagation_matrices_forward(MKL_Complex16* back_propagation_matrix_up, MKL_Complex16 *back_propagation_matrix_down,int_st ist,cns_st cns) {

       /*Copy Back Propagation Matrices One Step Back So That You Can Advance Accoding to the Time Slices*/
       MKL_Complex16 *new_back_propagation_matrix_up, *new_back_propagation_matrix_down; 
       int back_prop_size = ist.n_spatial_orbitals_sq*cns.max_number_walkers*ist.n_steps_back_propagation; 
       int back_prop_size_one_less = ist.n_spatial_orbitals_sq*cns.max_number_walkers*(ist.n_steps_back_propagation-1);  

       new_back_propagation_matrix_up=(MKL_Complex16 *)calloc(back_prop_size,sizeof(MKL_Complex16)); 
       new_back_propagation_matrix_down=(MKL_Complex16 *)calloc(back_prop_size,sizeof(MKL_Complex16)); 

       czero_vec(new_back_propagation_matrix_up,back_prop_size); 
       czero_vec(new_back_propagation_matrix_down,back_prop_size); 

       /*Copy New Entries Down*/
       cblas_zcopy(back_prop_size_one_less,back_propagation_matrix_up,1,&new_back_propagation_matrix_up[ist.n_spatial_orbitals_sq*cns.max_number_walkers],1); 
       cblas_zcopy(back_prop_size_one_less,back_propagation_matrix_down,1,&new_back_propagation_matrix_down[ist.n_spatial_orbitals_sq*cns.max_number_walkers],1);

       /*Copy Back Into Original Matrix*/ 
       cblas_zcopy(back_prop_size,new_back_propagation_matrix_up,1,back_propagation_matrix_up,1); 
       cblas_zcopy(back_prop_size,new_back_propagation_matrix_down,1,back_propagation_matrix_down,1); 

free(new_back_propagation_matrix_up); free(new_back_propagation_matrix_down); 
return; 
}

/************************************************************************************************************/

void copy_prev_wf_forward(MKL_Complex16 *prev_wf_up,MKL_Complex16 *prev_wf_down,int_st ist,cns_st cns) {

   /*Move the Stored WFs One WF Forward So That You Have the Correct WF For Back Propagation*/

   MKL_Complex16 *new_prev_wf_up, *new_prev_wf_down; 
   int prev_size_up = ist.n_spatial_orbitals_n_up*cns.max_number_walkers*ist.n_steps_back_propagation; 
   int prev_size_down = ist.n_spatial_orbitals_n_down*cns.max_number_walkers*ist.n_steps_back_propagation; 
   int prev_size_up_one_less = ist.n_spatial_orbitals_n_up*cns.max_number_walkers*(ist.n_steps_back_propagation-1); 
   int prev_size_down_one_less = ist.n_spatial_orbitals_n_down*cns.max_number_walkers*(ist.n_steps_back_propagation-1);   

   new_prev_wf_up=(MKL_Complex16 *)calloc(prev_size_up,sizeof(MKL_Complex16)); 
   new_prev_wf_down=(MKL_Complex16 *)calloc(prev_size_down,sizeof(MKL_Complex16)); 
   
   czero_vec(new_prev_wf_up,prev_size_up); 
   czero_vec(new_prev_wf_down,prev_size_down);   

   /*Copy New Entries Down*/
   cblas_zcopy(prev_size_up_one_less,prev_wf_up,1,&new_prev_wf_up[ist.n_spatial_orbitals_n_up*cns.max_number_walkers],1); 
   cblas_zcopy(prev_size_down_one_less,prev_wf_down,1,&new_prev_wf_down[ist.n_spatial_orbitals_n_down*cns.max_number_walkers],1); 

   /*Copy Back Into Original Matrix*/
   cblas_zcopy(prev_size_up,new_prev_wf_up,1,prev_wf_up,1);
   cblas_zcopy(prev_size_down,new_prev_wf_down,1,prev_wf_down,1);

free(new_prev_wf_up); free(new_prev_wf_down); 
return; 
}

/************************************************************************************************************/ 

void propagate_forwards_potential_chemistry_bp(MKL_Complex16 *back_propagation_matrices_up,MKL_Complex16 *back_propagation_matrices_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,int_st ist,cns_st cns) {

    /*Do Back Propagation Use Back Propagation Matrices for the Potential*/
    int i; 
    MKL_Complex16 *One, *Zero; 
    MKL_Complex16 *new_bp_wf_up, *new_bp_wf_down;   
    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    /*New BP WF*/
    new_bp_wf_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16)); 
    new_bp_wf_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16)); 

    /*Run Through Walkers*/
    for (i=0; i<cns.max_number_walkers; i++) {
  
       /*Zero the New Wave Functions First*/ 
       czero_vec(new_bp_wf_up,ist.n_spatial_orbitals_n_up); 
       czero_vec(new_bp_wf_down,ist.n_spatial_orbitals_n_down);   
       
       /*Multiply Matrix Into Wave function*/
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_up,ist.n_spatial_orbitals,One,&back_propagation_matrices_up[i*ist.n_spatial_orbitals_sq],ist.n_spatial_orbitals,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],ist.n_up,Zero,new_bp_wf_up,ist.n_up); 
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_down,ist.n_spatial_orbitals,One,&back_propagation_matrices_down[i*ist.n_spatial_orbitals_sq],ist.n_spatial_orbitals,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],ist.n_down,Zero,new_bp_wf_down,ist.n_down);      
 
       cblas_zcopy(ist.n_spatial_orbitals_n_up,new_bp_wf_up,1,&bp_wf_up[i*ist.n_spatial_orbitals_n_up],1); 
       cblas_zcopy(ist.n_spatial_orbitals_n_down,new_bp_wf_down,1,&bp_wf_down[i*ist.n_spatial_orbitals_n_down],1); 

    }

free(One); free(Zero);  
free(new_bp_wf_up); 
free(new_bp_wf_down); 
return; 
}

/**************************************************************************************************************/

void store_prev_wf(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *prev_wf_up,MKL_Complex16 *prev_wf_down,double *kinetic_forwards_half_fermions_real_up,int_st ist) {

     /*Store Wave Functions for Back Propagation and Advance Them One Half Step KE Forward To Compare With WFs in Computation of Energy*/
     int i; 
     MKL_Complex16 *stored_product_up, *stored_product_down;

     stored_product_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
     stored_product_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));
 
     /*Run Through Walkers Storing Previous Wave Function*/
     for (i=0; i<ist.n_walkers; i++) { 

       /*First Get Up Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions_real_up,&wf_up[i*ist.n_spatial_orbitals_n_up],stored_product_up,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up);
       cblas_zcopy(ist.n_spatial_orbitals_n_up,stored_product_up,1,&prev_wf_up[i*ist.n_spatial_orbitals_n_up],1);

       /*Get Down Matrices*/
       dmat_cmat(kinetic_forwards_half_fermions_real_up,&wf_down[i*ist.n_spatial_orbitals_n_down],stored_product_down,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down);
       cblas_zcopy(ist.n_spatial_orbitals_n_down,stored_product_down,1,&prev_wf_down[i*ist.n_spatial_orbitals_n_down],1);

     }

free(stored_product_up); 
free(stored_product_down); 
return; 
}



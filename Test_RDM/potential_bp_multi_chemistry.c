#include "afqmc.h"

/*******************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside; 
   int gamma, i, l; 
   int i_spin, l_spin, i_orbital, l_orbital;  
   double field_one, field_two;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 propagation_exponent_2, propagation_exponent_conjugate_2; 
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down; 
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;    
   MKL_Complex16 ratio_determinant; 
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16)); 
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));    
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16)); 
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16)); 

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16)); 
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16)); 

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Run Through Sites*/
        flag_inside = 0; 
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {
 
          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);
           field_two = gasdev(&tidum);

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2; 
            i_orbital = (int)(i/2.0);   

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2; 
             l_orbital = (int)(l/2.0);  

             /*Only Works If Two Spins are the Same*/ 
             if ( i_spin == l_spin ) { 

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]); 
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma])); 

                  propagation_exponent_2 = Cmul(potential_eigs_fermions_second_transform[gamma], potential_eigvecs_fermions[i * ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]); 
                  propagation_exponent_conjugate_2 = Cmul(potential_eigs_fermions_second_transform[gamma], conjugate(potential_eigvecs_fermions[l * ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma])); 

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 || Cabs(propagation_exponent_2)>0 || Cabs(propagation_exponent_conjugate_2)>0 ) {
                    flag_inside++; 

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(RCmul(-1.0*field_one,propagation_exponent_1), RCmul(-1.0*field_two,propagation_exponent_2)));
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(RCmul(-1.0*field_one,propagation_exponent_conjugate_1), RCmul(1.0*field_two, propagation_exponent_conjugate_2)));
                    }
                    else { 
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(RCmul(-1.0*field_one,propagation_exponent_1), RCmul(-1.0*field_two,propagation_exponent_2)));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(RCmul(-1.0*field_one,propagation_exponent_conjugate_1), RCmul(1.0*field_two, propagation_exponent_conjugate_2)));
                   } 

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/  

            } /*Check Eig*/

           } /*Gamma*/ 

           /*If Matrices Non-Zero*/
           if ( flag_inside != 0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals); 
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);  

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);   
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);            

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);  

            /*Find Ratios of Determinants*/
            ratio_determinant = determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_up);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_down)); 

            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1); 

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            weights[walker] = Cmul(weights[walker],ratio_determinant);       
           }/*sum inside*/ 

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up); 
free(new_wavefunction_down);
 
free(eigs_propagation_matrix_up); 
free(eigvecs_propagation_matrix_up);  
free(exponential_propagation_matrix_up); 
free(propagation_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
return;
}

/*********************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_real_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside; 
   int gamma, i, l; 
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one; 
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 ratio_determinant;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        flag_inside = 0;  
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 ) {
                    flag_inside++;  
      
                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], RCmul(-1.0*field_one,propagation_exponent_1)); 
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], RCmul(-1.0*field_one,propagation_exponent_conjugate_1));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], RCmul(-1.0*field_one,propagation_exponent_1));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], RCmul(-1.0*field_one,propagation_exponent_conjugate_1));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) { 

            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

            ratio_determinant = determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_up);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_down));

            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            weights[walker] = Cmul(weights[walker],ratio_determinant);

           }

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
return;
}

/**************************************************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_multi_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix, MKL_Complex16 *potential_original_matrix,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside; 
   int gamma, i, l; 
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one, field_two;
   double local_energy; 
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0; 
   MKL_Complex16 field_shift_one, field_shift_two;
   MKL_Complex16 field_final_one, field_final_two;   
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 propagation_exponent_2, propagation_exponent_conjugate_2;
   MKL_Complex16 local_ke, local_pe; 
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_up, *density_matrix_down; 
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16)); 
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16)); 

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);       
    
        flag_inside = 0;  
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);
           field_two = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,gamma,ist,0); 
           field_shift_two = get_shift(density_matrix_up,density_matrix_down,potential_eigs_fermions_second_transform[gamma],potential_eigvecs_fermions,gamma,ist,1);
 
           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real; 
           field_final_one.imag = -field_shift_one.imag; 

           field_final_two.real = field_two - field_shift_two.real; 
           field_final_two.imag = -field_shift_two.imag;   

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                  propagation_exponent_2 = Cmul(potential_eigs_fermions_second_transform[gamma], potential_eigvecs_fermions[i * ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_2 = Cmul(potential_eigs_fermions_second_transform[gamma], conjugate(potential_eigvecs_fermions[l * ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 || Cabs(propagation_exponent_2)>0 || Cabs(propagation_exponent_conjugate_2)>0 ) {

                    flag_inside++; 
                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2)))); 
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)),Cmul(field_final_two, propagation_exponent_conjugate_2)));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1)), Cmul(field_final_two, propagation_exponent_conjugate_2)));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           /*Ensure Matrix Is Nonzero*/
           if ( flag_inside != 0 ) { 
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

            /*Obtain the Local Energy*/
            local_ke = compute_kinetic_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,kinetic_original_matrix,ist);
            local_pe = compute_potential_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,potential_original_matrix,ist);
            local_energy = local_ke.real + local_pe.real; 

            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]); 
           } /*sum inside*/

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up); 

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down); 
return;
}

/*********************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_multi_unrestricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *kinetic_original_matrix_down,MKL_Complex16 *potential_original_matrix_up,MKL_Complex16 *potential_original_matrix_down,MKL_Complex16 *potential_original_matrix_updown,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;
   int gamma, i, l;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one, field_two;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one, field_shift_two;
   MKL_Complex16 field_final_one, field_final_two;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 propagation_exponent_2, propagation_exponent_conjugate_2;
   MKL_Complex16 local_ke, local_pe;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);

        flag_inside = 0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);
           field_two = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,gamma,ist,0);
           field_shift_two = get_shift(density_matrix_up,density_matrix_down,potential_eigs_fermions_second_transform[gamma],potential_eigvecs_fermions,gamma,ist,1);

           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
           field_final_one.imag = -field_shift_one.imag;

           field_final_two.real = field_two - field_shift_two.real;
           field_final_two.imag = -field_shift_two.imag;

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));propagation_exponent_2 = Cmul(potential_eigs_fermions_second_transform[gamma], potential_eigvecs_fermions[i * ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_2 = Cmul(potential_eigs_fermions_second_transform[gamma], conjugate(potential_eigvecs_fermions[l * ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 || Cabs(propagation_exponent_2)>0 || Cabs(propagation_exponent_conjugate_2)>0 ) {

                    flag_inside++;
                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)),Cmul(field_final_two, propagation_exponent_conjugate_2)));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1)), Cmul(field_final_two, propagation_exponent_conjugate_2)));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           /*Ensure Matrix Is Nonzero*/
           if ( flag_inside != 0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

            /*Obtain the Local Energy*/
            local_ke = compute_kinetic_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,kinetic_original_matrix_up,kinetic_original_matrix_down,ist);
            local_pe = compute_potential_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,ist);
            local_energy = local_ke.real + local_pe.real;

            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);
           } /*sum inside*/

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}



/*********************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_real_multi_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix, MKL_Complex16 *potential_original_matrix,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int gamma, i, l;
   int flag_inside;  
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one; 
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one; 
   MKL_Complex16 field_final_one;  
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 local_ke, local_pe;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);

        flag_inside = 0;  
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,gamma,ist,0);

           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
           field_final_one.imag = -field_shift_one.imag;

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 ) {
                    flag_inside++; 

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1))); 
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1))); 
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1))); 
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1))); 
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);


            /*Obtain the Local Energy*/
            local_ke = compute_kinetic_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,kinetic_original_matrix,ist);
            local_pe = compute_potential_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,potential_original_matrix,ist);
            local_energy = local_ke.real + local_pe.real;

            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);
           }

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/***********************************************************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_real_multi_unrestricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *kinetic_original_matrix_down, MKL_Complex16 *potential_original_matrix_up,MKL_Complex16 *potential_original_matrix_down,MKL_Complex16 *potential_original_matrix_updown,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int gamma, i, l;
   int flag_inside;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one;
   MKL_Complex16 field_final_one;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 local_ke, local_pe;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);

        flag_inside = 0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,gamma,ist,0);

           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
           field_final_one.imag = -field_shift_one.imag;

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 ) {
                    flag_inside++;

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)));
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1)));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1)));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);


            /*Obtain the Local Energy*/
            local_ke = compute_kinetic_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,kinetic_original_matrix_up,kinetic_original_matrix_down,ist);
            local_pe = compute_potential_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,ist);
            local_energy = local_ke.real + local_pe.real;

            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);
           }

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/*********************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_meanfield_multi_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix, MKL_Complex16 *potential_original_matrix,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int gamma, i, l;
   int flag_inside; 
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one, field_two;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one, field_shift_two;
   MKL_Complex16 field_final_one, field_final_two;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 propagation_exponent_2, propagation_exponent_conjugate_2;
   MKL_Complex16 local_ke, local_pe;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals*sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals*sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_up*sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_down*sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);

        flag_inside = 0; 
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);
           field_two = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_both(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);
           field_shift_two = get_shift_meanfield_both(density_matrix_up,density_matrix_down,potential_eigs_fermions_second_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,1);

           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
           field_final_one.imag = -field_shift_one.imag;

           field_final_two.real = field_two - field_shift_two.real;
           field_final_two.imag = -field_shift_two.imag;

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                  propagation_exponent_2 = Cmul(potential_eigs_fermions_second_transform[gamma], potential_eigvecs_fermions[i * ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_2 = Cmul(potential_eigs_fermions_second_transform[gamma], conjugate(potential_eigvecs_fermions[l * ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 || Cabs(propagation_exponent_2)>0 || Cabs(propagation_exponent_conjugate_2)>0 ) {
                    flag_inside++;  
          
                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)),Cmul(field_final_two, propagation_exponent_conjugate_2)));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1)), Cmul(field_final_two, propagation_exponent_conjugate_2)));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

            /*Obtain the Local Energy*/
            local_ke = compute_kinetic_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,kinetic_original_matrix,ist);
            local_pe = compute_potential_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,potential_original_matrix,ist);
            local_energy = local_ke.real + local_pe.real;

            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

             weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);
            }

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/***********************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_meanfield_multi_unrestricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *kinetic_original_matrix_down, MKL_Complex16 *potential_original_matrix_up,MKL_Complex16 *potential_original_matrix_down,MKL_Complex16 *potential_original_matrix_updown,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int gamma, i, l;
   int flag_inside;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one, field_two;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one, field_shift_two;
   MKL_Complex16 field_final_one, field_final_two;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 propagation_exponent_2, propagation_exponent_conjugate_2;
   MKL_Complex16 local_ke, local_pe;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals*sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals*sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_up*sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_down*sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);

        flag_inside = 0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);
           field_two = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_both(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);
           field_shift_two = get_shift_meanfield_both(density_matrix_up,density_matrix_down,potential_eigs_fermions_second_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,1);

           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
           field_final_one.imag = -field_shift_one.imag;

           field_final_two.real = field_two - field_shift_two.real;
           field_final_two.imag = -field_shift_two.imag;

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                  propagation_exponent_2 = Cmul(potential_eigs_fermions_second_transform[gamma], potential_eigvecs_fermions[i * ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_2 = Cmul(potential_eigs_fermions_second_transform[gamma], conjugate(potential_eigvecs_fermions[l * ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 || Cabs(propagation_exponent_2)>0 || Cabs(propagation_exponent_conjugate_2)>0 ) {
                    flag_inside++;

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)),Cmul(field_final_two, propagation_exponent_conjugate_2)));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1)), Cmul(field_final_two, propagation_exponent_conjugate_2)));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

            /*Obtain the Local Energy*/
            local_ke = compute_kinetic_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,kinetic_original_matrix_up,kinetic_original_matrix_down,ist);
            local_pe = compute_potential_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,ist);
            local_energy = local_ke.real + local_pe.real;

            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

             weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);
            }

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/***********************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_meanfield_real_multi_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix, MKL_Complex16 *potential_original_matrix,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside; 
   int gamma, i, j, l, m, n;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one; 
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one; 
   MKL_Complex16 field_final_one; 
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 local_ke, local_pe;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);
      
        flag_inside = 0; 
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {
           
           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_both(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);


           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
           field_final_one.imag = -field_shift_one.imag;

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 ) {
                    flag_inside++;   

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1))); 
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1))); 
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1))); 
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1))); 
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {

             /*Need to Diagonalize and Exponentiate Propagation Matrix*/
             complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
             complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

             /*Exponentiate Propagation Matrix*/
             exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
             exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

             /*Store Matrices for Back Propagation*/
             cblas_zcopy(ist.n_spatial_orbitals_sq,exponential_propagation_matrix_up,1,&potential_matrices_back_propagation_up[walker*ist.n_spatial_orbitals_sq],1); 
             cblas_zcopy(ist.n_spatial_orbitals_sq,exponential_propagation_matrix_down,1,&potential_matrices_back_propagation_down[walker*ist.n_spatial_orbitals_sq],1); 

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,kinetic_original_matrix,ist);
             local_pe = compute_potential_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,potential_original_matrix,ist);
             local_energy = local_ke.real + local_pe.real;
 
             cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
             cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

             update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

             /*Multiply By Trial Energy Outside of This Function*/ 
             weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);

           } /*If There Is a Sum Inside and The Matrix Is NonZero*/

    } /*Check Weights*/
  
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/*********************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *kinetic_original_matrix_down,MKL_Complex16 *potential_original_matrix_up,MKL_Complex16 *potential_original_matrix_down,MKL_Complex16 *potential_original_matrix_updown,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;
   int gamma, i, j, l, m, n;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one;
   MKL_Complex16 field_final_one;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 local_ke, local_pe;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);

        flag_inside = 0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_both(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);


           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
           field_final_one.imag = -field_shift_one.imag;

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 ) {
                    flag_inside++;

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)));
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1)));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1)));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {

             /*Need to Diagonalize and Exponentiate Propagation Matrix*/
             complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
             complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

             /*Exponentiate Propagation Matrix*/
             exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
             exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

             /*Store Matrices for Back Propagation*/
             cblas_zcopy(ist.n_spatial_orbitals_sq,exponential_propagation_matrix_up,1,&potential_matrices_back_propagation_up[walker*ist.n_spatial_orbitals_sq],1);
             cblas_zcopy(ist.n_spatial_orbitals_sq,exponential_propagation_matrix_down,1,&potential_matrices_back_propagation_down[walker*ist.n_spatial_orbitals_sq],1);

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,kinetic_original_matrix_up,kinetic_original_matrix_down,ist);
             local_pe = compute_potential_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,ist);
             local_energy = local_ke.real + local_pe.real;

             cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
             cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

             update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

             /*Multiply By Trial Energy Outside of This Function*/
             weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);

           } /*If There Is a Sum Inside and The Matrix Is NonZero*/

    } /*Check Weights*/

   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/*********************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside; 
   int gamma, i, l; 
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one, field_two;
   double phase, projected_ratio_determinant, abs_value;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one, field_shift_two;
   MKL_Complex16 field_final_one, field_final_two;
   MKL_Complex16 shift_product_one, shift_product_two;
   MKL_Complex16 shift_exponential_one, shift_exponential_two;
   MKL_Complex16 total_shift_exponential;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 propagation_exponent_2, propagation_exponent_conjugate_2;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 ratio_determinant;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);

        flag_inside = 0; 
        total_shift_exponential.real = 1.0; total_shift_exponential.imag = 0.0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);
           field_two = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,gamma,ist,0);
           field_shift_two = get_shift(density_matrix_up,density_matrix_down,potential_eigs_fermions_second_transform[gamma],potential_eigvecs_fermions,gamma,ist,1);

           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
field_final_one.imag = -field_shift_one.imag;

           field_final_two.real = field_two - field_shift_two.real;
           field_final_two.imag = -field_shift_two.imag;

           /*Find Exponentials of Shift Terms*/
           shift_product_one = Csub(RCmul(field_one, field_shift_one), RCmul(.5, Cmul(field_shift_one,field_shift_one)));
           shift_exponential_one.real = exp(shift_product_one.real) * cos(shift_product_one.imag);
           shift_exponential_one.imag = exp(shift_product_one.real) * sin(shift_product_one.imag);

           shift_product_two = Csub(RCmul(field_two, field_shift_two), RCmul(.5, Cmul(field_shift_two,field_shift_two)));
           shift_exponential_two.real = exp(shift_product_two.real) * cos(shift_product_two.imag);
           shift_exponential_two.imag = exp(shift_product_two.real) * sin(shift_product_two.imag);

           total_shift_exponential = Cmul(total_shift_exponential, Cmul(shift_exponential_one, shift_exponential_two));

           /*HHEERE*/
           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                  propagation_exponent_2 = Cmul(potential_eigs_fermions_second_transform[gamma], potential_eigvecs_fermions[i * ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_2 = Cmul(potential_eigs_fermions_second_transform[gamma], conjugate(potential_eigvecs_fermions[l * ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 || Cabs(propagation_exponent_2)>0 || Cabs(propagation_exponent_conjugate_2)>0 ) {
                     flag_inside++; 
  
                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)),Cmul(field_final_two, propagation_exponent_conjugate_2)));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1)), Cmul(field_final_two, propagation_exponent_conjugate_2)));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

            /*Ratio Determinants*/
            ratio_determinant = determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_up);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_down));
            phase = atan(ratio_determinant.imag/ratio_determinant.real);

            /*Multiply By Shifts*/
            ratio_determinant = RCmul(cns.exp_trial_energy, Cmul(ratio_determinant, total_shift_exponential));

            /*Update Walkers and Overlap Matrices*/
            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            /*Project Down Last Part of Phase Just In Case Near 0 in Axis*/
            abs_value = Cabs(ratio_determinant);
            if ( fabs(phase) < .5 * 3.1415 && abs_value > 0.0001  ) {
               projected_ratio_determinant = abs_value * cos(phase);
            }
            else {
               projected_ratio_determinant = 0.0;
            }
            weights[walker] = Cmul(weights[walker], ratio_determinant); 
            //weights[walker] = RCmul(projected_ratio_determinant, weights[walker]);

          } /*sum inside*/

    } /*Check Weights*/
   } /*Run Through Walkers*/


free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/*****************************************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_meanfield_real_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_meanfield_first_constant, MKL_Complex16 *potential_meanfield_second_constant,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside; 
   int gamma, i, l; 
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one; 
   double phase, projected_ratio_determinant, abs_value;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one; 
   MKL_Complex16 field_final_one; 
   MKL_Complex16 shift_product_one; 
   MKL_Complex16 shift_exponential_one; 
   MKL_Complex16 total_shift_exponential; 
   MKL_Complex16 total_meanfield_exponential; 
   MKL_Complex16 exponent_mean_field, exponential_mean_field; 
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1; 
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 ratio_determinant;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;
 
   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);

        total_shift_exponential.real = 1.0; total_shift_exponential.imag = 0.0;  
        total_meanfield_exponential.real = 1.0; total_meanfield_exponential.imag = 0.0;

        flag_inside = 0;  
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_both(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);
           field_shift_one.real = field_shift_one.imag = 0.0; 

           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
           field_final_one.imag = -field_shift_one.imag;

           /*Find Exponentials of Shift Terms*/
           shift_product_one = Csub(RCmul(field_one, field_shift_one), RCmul(.5, Cmul(field_shift_one,field_shift_one))); 
           shift_exponential_one.real = exp(shift_product_one.real) * cos(shift_product_one.imag); 
           shift_exponential_one.imag = exp(shift_product_one.real) * sin(shift_product_one.imag);  

           total_shift_exponential = shift_exponential_one; 
 
           /*Find the Mean Field Exponential*/
           if ( Cabs(potential_meanfield_first_constant[gamma]) > .0001 ) { 
             exponent_mean_field = Cmul(RCmul(-1, field_final_one), potential_meanfield_first_constant[gamma]); 
             exponential_mean_field.real = exp(exponent_mean_field.real) * cos(exponent_mean_field.imag); 
             exponential_mean_field.imag = exp(exponent_mean_field.real) * sin(exponent_mean_field.imag);
             total_meanfield_exponential = Cmul(total_meanfield_exponential, exponential_mean_field);   
           }  

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 ) {
                    flag_inside++; 

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1))); 
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1))); 
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1)));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside !=0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

            /*Ratio Determinants*/
            ratio_determinant = determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_up);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_down));
            phase = atan(ratio_determinant.imag/ratio_determinant.real); 

            /*Multiply By Shifts*/
            ratio_determinant = RCmul(cns.exp_trial_energy, Cmul(ratio_determinant, Cmul(total_shift_exponential, total_meanfield_exponential))); 

            /*Update Walkers and Overlap Matrices*/
            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            /*Project Down Last Part of Phase Just In Case Near 0 in Axis*/
            abs_value = Cabs(ratio_determinant);
            if ( fabs(phase) < .5 * 3.1415 && abs_value > 0.0001  ) {
               projected_ratio_determinant = abs_value * cos(phase);
            }
            else {
               projected_ratio_determinant = 0.0;
            } 

            weights[walker] = RCmul(projected_ratio_determinant, weights[walker]); 
           }

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/*********************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_meanfield_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_meanfield_first_constant, MKL_Complex16 *potential_meanfield_second_constant,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside; 
   int gamma, i, l; 
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one, field_two; 
   double phase, projected_ratio_determinant, abs_value;
   double phase_two; 
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one, field_shift_two; 
   MKL_Complex16 field_final_one, field_final_two; 
   MKL_Complex16 shift_product_one, shift_product_two; 
   MKL_Complex16 shift_exponential_one, shift_exponential_two; 
   MKL_Complex16 total_shift_exponential;
   MKL_Complex16 total_meanfield_exponential;
   MKL_Complex16 product_exponentials; 
   MKL_Complex16 exponent_mean_field, exponential_mean_field;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1; 
   MKL_Complex16 propagation_exponent_2, propagation_exponent_conjugate_2;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 ratio_determinant;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));
 density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);

        flag_inside = 0; 
        total_shift_exponential.real = 1.0; total_shift_exponential.imag = 0.0;
        total_meanfield_exponential.real = 1.0; total_meanfield_exponential.imag = 0.0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);
           field_two = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_both(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);
           field_shift_two = get_shift_meanfield_both(density_matrix_up,density_matrix_down,potential_eigs_fermions_second_transform[gamma],potential_eigvecs_fermions,potential_meanfield_second_constant[gamma],gamma,ist,1);

           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
           field_final_one.imag = -field_shift_one.imag;

           field_final_two.real = field_two - field_shift_two.real;
           field_final_two.imag = -field_shift_two.imag;

           /*Find Exponentials of Shift Terms*/
           shift_product_one = Csub(RCmul(field_one, field_shift_one), RCmul(.5, Cmul(field_shift_one,field_shift_one)));
           shift_exponential_one.real = exp(shift_product_one.real) * cos(shift_product_one.imag);
           shift_exponential_one.imag = exp(shift_product_one.real) * sin(shift_product_one.imag);

           shift_product_two = Csub(RCmul(field_two, field_shift_two), RCmul(.5, Cmul(field_shift_two,field_shift_two)));
           shift_exponential_two.real = exp(shift_product_two.real) * cos(shift_product_two.imag);
           shift_exponential_two.imag = exp(shift_product_two.real) * sin(shift_product_two.imag);

           total_shift_exponential = Cmul(total_shift_exponential, Cmul(shift_exponential_one, shift_exponential_two));

           /*Find the Mean Field Exponential*/
           exponent_mean_field = Cadd(Cmul(RCmul(-1, field_final_one), potential_meanfield_first_constant[gamma]), Cmul(RCmul(-1.0,field_final_two), potential_meanfield_second_constant[gamma]));
           exponential_mean_field.real = exp(exponent_mean_field.real) * cos(exponent_mean_field.imag);
           exponential_mean_field.imag = exp(exponent_mean_field.real) * sin(exponent_mean_field.imag);
           total_meanfield_exponential = Cmul(total_meanfield_exponential, exponential_mean_field);

           /*HHEERE*/
           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                  propagation_exponent_2 = Cmul(potential_eigs_fermions_second_transform[gamma], potential_eigvecs_fermions[i * ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_2 = Cmul(potential_eigs_fermions_second_transform[gamma], conjugate(potential_eigvecs_fermions[l * ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 || Cabs(propagation_exponent_2)>0 || Cabs(propagation_exponent_conjugate_2)>0 ) {
                    flag_inside++;  

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)),Cmul(field_final_two, propagation_exponent_conjugate_2)));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1)), Cmul(Negative_One,Cmul(field_final_two,propagation_exponent_2))));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1)), Cmul(field_final_two, propagation_exponent_conjugate_2)));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

            /*Ratio Determinants*/
            ratio_determinant = determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_up);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_down));
            phase = atan(ratio_determinant.imag/ratio_determinant.real);

            /*Multiply By Shifts*/
            product_exponentials=Cmul(total_shift_exponential,total_meanfield_exponential); 
            phase_two = atan(product_exponentials.imag/product_exponentials.real);

            ratio_determinant = RCmul(cns.exp_trial_energy, Cmul(ratio_determinant, Cmul(total_shift_exponential, total_meanfield_exponential)));


            /*Update Walkers and Overlap Matrices*/
            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            /*Project Down Last Part of Phase Just In Case Near 0 in Axis*/
            abs_value = Cabs(ratio_determinant);
            if ( fabs(phase) < .5 * 3.1415 && abs_value > 0.0001  ) {
               projected_ratio_determinant = abs_value * cos(phase);
            }
            else {
               projected_ratio_determinant = 0.0;
            }

            weights[walker] = Cmul(ratio_determinant, weights[walker]);   
            //weights[walker] = RCmul(projected_ratio_determinant, weights[walker]);

           } /*sum inside*/

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/*********************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_real_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside; 
   int gamma, i, l; 
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one; 
   double phase, projected_ratio_determinant, abs_value;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one; 
   MKL_Complex16 field_final_one; 
   MKL_Complex16 shift_product_one; 
   MKL_Complex16 shift_exponential_one; 
   MKL_Complex16 total_shift_exponential;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 ratio_determinant;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_up,density_matrix_down,ist);

        flag_inside = 0; 
        total_shift_exponential.real = 1.0; total_shift_exponential.imag = 0.0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift(density_matrix_up,density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,gamma,ist,0);

           /*Get Total Shift Field*/
           field_final_one.real = field_one - field_shift_one.real;
           field_final_one.imag = -field_shift_one.imag;

           /*Find Exponentials of Shift Terms*/
           shift_product_one = Csub(RCmul(field_one, field_shift_one), RCmul(.5, Cmul(field_shift_one,field_shift_one)));
           shift_exponential_one.real = exp(shift_product_one.real) * cos(shift_product_one.imag);
           shift_exponential_one.imag = exp(shift_product_one.real) * sin(shift_product_one.imag);

           total_shift_exponential = Cmul(total_shift_exponential, shift_exponential_one); 

           /*HHEERE*/
           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 ) {
                    flag_inside++; 

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1))); 
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1))); 
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_1))); 
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One, Cmul(field_final_one,propagation_exponent_conjugate_1))); 
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {
              /*Need to Diagonalize and Exponentiate Propagation Matrix*/
              complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
              complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

              /*Exponentiate Propagation Matrix*/
              exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

              propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

              /*Ratio Determinants*/
              ratio_determinant = determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_up);
              ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_down));
              phase = atan(ratio_determinant.imag/ratio_determinant.real); 

              /*Multiply By Shifts*/
              ratio_determinant = RCmul(cns.exp_trial_energy, Cmul(ratio_determinant, total_shift_exponential));

              /*Update Walkers and Overlap Matrices*/
              cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
              cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

              update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);
 
             /*Project Down Last Part of Phase Just In Case Near 0 in Axis*/
             abs_value = Cabs(ratio_determinant);
             if ( fabs(phase) < .5 * 3.1415 && abs_value > .00001 ) {
                projected_ratio_determinant = abs_value * cos(phase);
             }
             else {
                projected_ratio_determinant = 0.0;
             }

             weights[walker] = RCmul(projected_ratio_determinant, weights[walker]); 
            
            } 

           } /*Check Weights*/
          } /*Run Through Walkers*/
            
free(new_wavefunction_up);
free(new_wavefunction_down);
            
free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);
            
free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/***************************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_meanfield_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;  
   int gamma, i, l; 
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one, field_two;
   MKL_Complex16 exponent_mean_field, exponential_gamma; 
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 propagation_exponent_2, propagation_exponent_conjugate_2;
   MKL_Complex16 total_gamma; 
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 ratio_determinant;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Run Through Sites*/
        flag_inside = 0; 
        total_gamma.real = 1.0; total_gamma.imag = 0.0;   
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);
           field_two = gasdev(&tidum);

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                  propagation_exponent_2 = Cmul(potential_eigs_fermions_second_transform[gamma], potential_eigvecs_fermions[i * ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
                  propagation_exponent_conjugate_2 = Cmul(potential_eigs_fermions_second_transform[gamma], conjugate(potential_eigvecs_fermions[l * ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 || Cabs(propagation_exponent_2)>0 || Cabs(propagation_exponent_conjugate_2)>0 ) {
                    flag_inside++; 

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(RCmul(-1.0*field_one,propagation_exponent_1), RCmul(-1.0*field_two,propagation_exponent_2)));
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(RCmul(-1.0*field_one,propagation_exponent_conjugate_1), RCmul(1.0*field_two, propagation_exponent_conjugate_2)));
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cadd(RCmul(-1.0*field_one,propagation_exponent_1), RCmul(-1.0*field_two,propagation_exponent_2)));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cadd(RCmul(-1.0*field_one,propagation_exponent_conjugate_1), RCmul(1.0*field_two, propagation_exponent_conjugate_2)));
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

             /*Find the Correct Mean Field Shift Taking Fields Into Account*/
             if ( Cabs(potential_meanfield_first_constant[gamma]) > .00001) {
               exponent_mean_field = Cadd(RCmul(-1*field_one,potential_meanfield_first_constant[gamma]), RCmul(-1*field_two,potential_meanfield_second_constant[gamma]));
               exponential_gamma.real = exp(exponent_mean_field.real) * cos(exponent_mean_field.imag);
               exponential_gamma.imag = exp(exponent_mean_field.real) * sin(exponent_mean_field.imag);
               total_gamma = Cmul(exponential_gamma, total_gamma);
             }


            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

            ratio_determinant = determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_up);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_down));

            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            /*Multiply Weights In From Gamma*/
            ratio_determinant = Cmul(ratio_determinant, total_gamma);  

            weights[walker] = Cmul(weights[walker],ratio_determinant);
           }

    } /*Check Weights*/
   } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
return;
}

/**********************************************************************************************************************************************/

void propagate_forwards_potential_continuous_bp_chemistry_meanfield_real_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_matrices_back_propagation_up,MKL_Complex16 *potential_matrices_back_propagation_down,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside; 
   int gamma, i, l;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one; 
   MKL_Complex16 exponent_mean_field, exponential_gamma;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 total_gamma;
   MKL_Complex16 *new_wavefunction_up, *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 ratio_determinant;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);
        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Run Through Sites*/
        flag_inside = 0; 
        total_gamma.real = 1.0; total_gamma.imag = 0.0;
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                  /*Get Propagation Exponent for Given Combination of Gamma, Sites*/
                  propagation_exponent_1 = Cmul(potential_eigs_fermions_first_transform[gamma], potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma]);
propagation_exponent_conjugate_1 = Cmul(potential_eigs_fermions_first_transform[gamma], conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]));

                   if ( Cabs(propagation_exponent_1)>0 || Cabs(propagation_exponent_conjugate_1)>0 ) {
                    flag_inside++;   

                    /*Form Exponential Matrix*/
                    if ( i_spin == 0 ) {
                      propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_up[i_orbital*ist.n_spatial_orbitals+l_orbital], RCmul(-1.0*field_one,propagation_exponent_1)); 
                      propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_up[l_orbital*ist.n_spatial_orbitals+i_orbital], RCmul(-1.0*field_one,propagation_exponent_conjugate_1)); 
                    }
                    else {
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], RCmul(-1.0*field_one,propagation_exponent_1)); 
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], RCmul(-1.0*field_one,propagation_exponent_conjugate_1)); 
                   }

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

             /*Find the Correct Mean Field Shift Taking Fields Into Account*/
             exponent_mean_field = RCmul(-1*field_one,potential_meanfield_first_constant[gamma]); 
             exponential_gamma.real = exp(exponent_mean_field.real) * cos(exponent_mean_field.imag);
             exponential_gamma.imag = exp(exponent_mean_field.real) * sin(exponent_mean_field.imag);

             /*Multiply Each Factor From Gamma Due to Mean Fields Into the Total Gamma*/
             total_gamma = Cmul(exponential_gamma, total_gamma);


            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {
            /*Need to Diagonalize and Exponentiate Propagation Matrix*/
            complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            /*Exponentiate Propagation Matrix*/
            exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);
            exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

            propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

            ratio_determinant = determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_up);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_chemistry_multi_phaseless(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,ist,ist.n_spatial_orbitals,ist.n_down));

            cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
            cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

            update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

            /*Multiply Weights In From Gamma*/
            ratio_determinant = Cmul(ratio_determinant, total_gamma);

            weights[walker] = Cmul(weights[walker],ratio_determinant);
           }

    } /*Check Weights*/
 } /*Run Through Walkers*/

free(new_wavefunction_up);
free(new_wavefunction_down);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);

free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
return;
}

/*********************************************************************************************************/


#include "afqmc.h"

/*All Potential Propagators that Are Comparible with the Mixed Estimator*/

/*******************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,int_st ist,cns_st cns,long *idum){

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

void propagate_forwards_potential_continuous_mixed_chemistry_real_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,int_st ist,cns_st cns,long *idum){

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

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_multi_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix, MKL_Complex16 *potential_original_matrix,int_st ist,cns_st cns,long *idum){

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

/**********************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_multi_unrestricted_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,int_st ist,cns_st cns,long *idum){

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
            local_ke = compute_kinetic_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,kinetic_ij_sparse_up,kinetic_ij_sparse_down,ist);
            local_pe = compute_potential_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,ist);
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

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_real_multi_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix, MKL_Complex16 *potential_original_matrix,int_st ist,cns_st cns,long *idum){

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

/*********************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_real_multi_unrestricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,int_st ist,cns_st cns,long *idum){

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
            local_ke = compute_kinetic_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,kinetic_ij_sparse_up,kinetic_ij_sparse_down,ist);
            local_pe = compute_potential_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,ist);
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

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_multi_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix, MKL_Complex16 *potential_original_matrix,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

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

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_multi_unrestricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_matrix_original_up,MKL_Complex16 *kinetic_matrix_original_down,MKL_Complex16 *potential_matrix_original_up,MKL_Complex16 *potential_matrix_original_down,MKL_Complex16 *potential_matrix_original_updown,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

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
            local_ke = compute_kinetic_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,kinetic_matrix_original_up,kinetic_matrix_original_down,ist);
            local_pe = compute_potential_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,potential_matrix_original_up,potential_matrix_original_down,potential_matrix_original_updown,ist);
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

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *potential_original_matrix_up,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   if ( ist.n_up > 0 ) {
    if ( ist.n_down > 0 ) {
      propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_restricted_both(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_original_matrix_up,potential_original_matrix_up,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum); 

    }
    else {
      propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_up(wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_total,det_overlap_up,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_original_matrix_up,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum); 
    }
   }
   else {
     propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_down(wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_down,overlap_total,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_original_matrix_up,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);
   } 

return;
}

/***********************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_restricted_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse,double *potential_matrix_sparse,int *kinetic_ij_sparse,int *potential_ij_sparse,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   if ( ist.n_up > 0 ) {
    if ( ist.n_down > 0 ) {
      propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_restricted_both_sparse(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse,potential_matrix_sparse,kinetic_ij_sparse,potential_ij_sparse,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);

    }
    else {
      propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_up_sparse_equilibration(wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_total,det_overlap_up,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse,kinetic_ij_sparse,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);
    }
   }
   else {
     propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_down_sparse_equilibration(wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_down,overlap_total,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse,kinetic_ij_sparse,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);
   }

return;
}

/***********************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *kinetic_original_matrix_down,MKL_Complex16 *potential_original_matrix_up,MKL_Complex16 *potential_original_matrix_down,MKL_Complex16 *potential_original_matrix_updown,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/

   if ( ist.n_up > 0 ) {
    if ( ist.n_down > 0 ) {
        propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_multi_unrestricted_both(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_original_matrix_up,kinetic_original_matrix_down,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,ist,cns,idum); 
    }
    else {
      propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_up(wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_total,det_overlap_up,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_original_matrix_up,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);
    }
   }
   else {
     propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_down(wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_down,overlap_total,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_original_matrix_down,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);
   }

return;
}

/*********************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_multi_unrestricted_both(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *kinetic_original_matrix_down,MKL_Complex16 *potential_original_matrix_up,MKL_Complex16 *potential_original_matrix_down,MKL_Complex16 *potential_original_matrix_updown,int_st ist,cns_st cns,long *idum){

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

                } /*Set Spins Equal*/

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

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted_sparse_equilibration(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *phases,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   if ( ist.n_up > 0 ) {
    if ( ist.n_down > 0 ) {
      propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted_sparse_both_equilibration(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,phases,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);

    }
    else {
      propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_up_sparse_equilibration(wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_total,det_overlap_up,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse_up,kinetic_ij_sparse_up,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);
    }
   }
   else {
     propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_down_sparse_equilibration(wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_down,overlap_total,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse_down,kinetic_ij_sparse_down,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);
   }

return;
}

/***********************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted_sparse_production(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *phases,int *local_energy_flag,int *field_flag,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   if ( ist.n_up > 0 ) {
    if ( ist.n_down > 0 ) {
      propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted_sparse_both_production(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,phases,local_energy_flag,field_flag,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);

    }
    else {
      propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_up_sparse_production(wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_total,det_overlap_up,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse_up,kinetic_ij_sparse_up,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);
    }
   }
   else {
     propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_down_sparse_production(wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_down,overlap_total,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse_down,kinetic_ij_sparse_down,potential_meanfield_first_constant,potential_meanfield_second_constant,ist,cns,idum);
   }

return;
}

/*********************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_restricted_both(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix,MKL_Complex16 *potential_original_matrix,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

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

/**********************************************************************************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_restricted_both_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse,double *potential_matrix_sparse,int *kinetic_ij_sparse,int *potential_ij_sparse,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

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

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_restricted_sparse(density_matrix_up,density_matrix_down,kinetic_matrix_sparse,kinetic_ij_sparse,ist);
             local_pe = compute_potential_energy_density_matrix_restricted_sparse(density_matrix_up,density_matrix_down,potential_matrix_sparse,potential_ij_sparse,ist);
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

/***************************************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted_sparse_both_equilibration(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *phases,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*HERE*/

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;
   int gamma, i, j, l, m, n;
   int i_spin, l_spin, i_orbital, l_orbital;
   double phase, theta_1, theta_2;  
   double projected_ratio_determinant, abs_value; 
   double field_one;
   double local_energy;
   MKL_Complex16 determinant_before, determinant_after; 
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one;
   MKL_Complex16 ratio_determinant; 
   MKL_Complex16 field_final_one;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 exponent_mean_field, exponential_mean_field; 
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

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,kinetic_ij_sparse_up,kinetic_ij_sparse_down,ist); 
             local_pe = compute_potential_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,ist); 
             local_energy = local_ke.real + local_pe.real;

             cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
             cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

             /*Obtain Constant from Mean Field*/ 
             if ( Cabs(potential_meanfield_first_constant[gamma]) > .0000001 ) {
               exponent_mean_field = Cmul(RCmul(-1, field_final_one), potential_meanfield_first_constant[gamma]);
               exponential_mean_field.real = exp(exponent_mean_field.real) * cos(exponent_mean_field.imag);
               exponential_mean_field.imag = exp(exponent_mean_field.real) * sin(exponent_mean_field.imag);
             }
             else {  
               exponential_mean_field.real = 1.0; exponential_mean_field.imag = 0.0; 
             }

             /*Store The Old Overlap*/
             determinant_before = overlap_total[walker];

           update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

             /*Multiply By New Overlap*/
             determinant_after = overlap_total[walker];

             /*Project Down Last Part of Phase Just In Case Near 0 in Axis*/
             theta_1 = acos(determinant_after.real/Cabs(determinant_after));
             theta_2 = acos(determinant_before.real/Cabs(determinant_before));
             phase = theta_1 - theta_2;
             phases[walker*3].real = phase; 
             phases[walker*3+1] = determinant_before; 
             phases[walker*3+2] = determinant_after; 

             /*if ( theta_1 < .5 * 3.1415 && theta_1 >= 0 ) {
                projected_ratio_determinant = 1; //cos(phase);
             }
             else {
                projected_ratio_determinant = 0.0;
             }*/
             projected_ratio_determinant = 1.0; 
             projected_ratio_determinant *= fabs(exp(-1*cns.dtau*local_energy));
             weights[walker] = RCmul(projected_ratio_determinant, weights[walker]);
             weights[walker] = Cmul(weights[walker], exponential_mean_field);          
 

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

/******************************************************************************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted_sparse_both_production(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *phases,int *local_energy_flag,int *field_flag,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;
   int gamma, i, l;
   int i_spin, l_spin, i_orbital, l_orbital;
   double phase, theta_1, theta_2;
   double projected_ratio_determinant; 
   double field_one;
   double local_energy;
   MKL_Complex16 determinant_before, determinant_after;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one;
   MKL_Complex16 field_final_one;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 exponent_mean_field, exponential_mean_field;
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

           field_flag[walker] = 0; 
           /*Place Cap*/
           if ( Cabs(field_shift_one) > 1.0 ) {
             field_shift_one = RCmul(1.0/Cabs(field_shift_one), field_shift_one); 
             field_flag[walker] = 1;  
           }


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

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_both(new_wavefunction_up,new_wavefunction_down,&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_up,exponential_propagation_matrix_down,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,kinetic_ij_sparse_up,kinetic_ij_sparse_down,ist);
             local_pe = compute_potential_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,ist);
             local_energy = local_ke.real + local_pe.real;

             local_energy_flag[walker] = 0;  
             /*Place Caps on Local Energy*/
             if ( local_energy > cns.energy_max ) {
               local_energy = cns.energy_max; 
               local_energy_flag[walker] = 1; 
             }
             else if (local_energy < cns.energy_min ) {
               local_energy = cns.energy_min;
               local_energy_flag[walker] = 1;  
             }

             cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);
             cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

             /*Obtain Constant from Mean Field*/
             if ( Cabs(potential_meanfield_first_constant[gamma]) > .0000001 ) {
               exponent_mean_field = Cmul(RCmul(-1, field_final_one), potential_meanfield_first_constant[gamma]);
               exponential_mean_field.real = exp(exponent_mean_field.real) * cos(exponent_mean_field.imag);
               exponential_mean_field.imag = exp(exponent_mean_field.real) * sin(exponent_mean_field.imag);
             }
             else {
               exponential_mean_field.real = 1.0; exponential_mean_field.imag = 0.0;
             }

             /*Store The Old Overlap*/
             determinant_before = overlap_total[walker];
             update_overlaps_chemistry_multi_phaseless_both(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

             /*Multiply By New Overlap*/
             determinant_after = overlap_total[walker];

             /*Project Down Last Part of Phase Just In Case Near 0 in Axis*/
             theta_1 = acos(determinant_after.real/Cabs(determinant_after));
             theta_2 = acos(determinant_before.real/Cabs(determinant_before));
             phase = theta_1 - theta_2;
             phases[walker*3].real = phase;
             phases[walker*3+1] = determinant_before;
             phases[walker*3+2] = determinant_after;

             projected_ratio_determinant = 1.0;
             projected_ratio_determinant *= fabs(exp(-1*cns.dtau*local_energy));
             weights[walker] = RCmul(projected_ratio_determinant, weights[walker]);
             weights[walker] = Cmul(weights[walker], exponential_mean_field);

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

/******************************************************************************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_up_sparse_equilibration(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,MKL_Complex16 *det_overlap_up,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse_up,int *kinetic_ij_sparse_up,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*HERE*/

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;
   int gamma, i, l;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one;
   MKL_Complex16 field_final_one;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 local_ke; 
   MKL_Complex16 *new_wavefunction_up; 
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *density_matrix_up; 
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi_up(&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_up[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],density_matrix_up,ist);

        flag_inside = 0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_up(density_matrix_up,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);

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

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {

             /*Need to Diagonalize and Exponentiate Propagation Matrix*/
             complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);

             /*Exponentiate Propagation Matrix*/
             exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_up(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],exponential_propagation_matrix_up,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_restricted_up_sparse(density_matrix_up,kinetic_matrix_sparse_up,kinetic_ij_sparse_up,ist);
             local_energy = local_ke.real; 

             cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);

             update_overlaps_chemistry_multi_phaseless_up(&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_up[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],ist);

             /*Multiply By Trial Energy Outside of This Function*/
             weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);

           } /*If There Is a Sum Inside and The Matrix Is NonZero*/

    } /*Check Weights*/

   } /*Run Through Walkers*/

free(new_wavefunction_up);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);
return;
}

/**************************************************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_up_sparse_production(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,MKL_Complex16 *det_overlap_up,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse_up,int *kinetic_ij_sparse_up,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*HERE*/

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;
   int gamma, i, l;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one;
   MKL_Complex16 field_final_one;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 local_ke, local_pe;
   MKL_Complex16 *new_wavefunction_up;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *density_matrix_up;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {
 czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi_up(&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_up[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],density_matrix_up,ist);

        flag_inside = 0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_up(density_matrix_up,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);

           /*Put Cap on Field Shift*/
           if ( Cabs(field_shift_one) > 1.0 ) {
              field_shift_one = RCmul(1.0/Cabs(field_shift_one), field_shift_one); 
           }

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

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {

             /*Need to Diagonalize and Exponentiate Propagation Matrix*/
             complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);

             /*Exponentiate Propagation Matrix*/
             exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_up(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],exponential_propagation_matrix_up,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_restricted_up_sparse(density_matrix_up,kinetic_matrix_sparse_up,kinetic_ij_sparse_up,ist);
             local_energy = local_ke.real + local_pe.real;

             /*Cap the Energies*/
             if ( local_energy > cns.energy_max ) {
               local_energy = cns.energy_max; 
             }
             else if ( local_energy < cns.energy_min ) {
               local_energy = cns.energy_min; 
             }

             cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);

             update_overlaps_chemistry_multi_phaseless_up(&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_up[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],ist);

             /*Multiply By Trial Energy Outside of This Function*/
             weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);

           } /*If There Is a Sum Inside and The Matrix Is NonZero*/
  } /*Check Weights*/

   } /*Run Through Walkers*/

free(new_wavefunction_up);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);
return;
}

/*************************************************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_up(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,MKL_Complex16 *det_overlap_up,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;
   int gamma, i, l;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one;
   MKL_Complex16 field_final_one;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 local_ke; 
   MKL_Complex16 *new_wavefunction_up;
   MKL_Complex16 *exponential_propagation_matrix_up, *propagation_matrix_up, *eigs_propagation_matrix_up, *eigvecs_propagation_matrix_up;
   MKL_Complex16 *density_matrix_up;
   long tidum;

   exponential_propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_up = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {
	  czero_vec(propagation_matrix_up,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi_up(&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_up[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],density_matrix_up,ist);

        flag_inside = 0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_up(density_matrix_up,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);

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

                 }

               } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {

             /*Need to Diagonalize and Exponentiate Propagation Matrix*/
             complex_eigvecs_2(propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);

             /*Exponentiate Propagation Matrix*/
             exponentiate_matrix_complex(exponential_propagation_matrix_up,eigs_propagation_matrix_up,eigvecs_propagation_matrix_up,ist.n_spatial_orbitals);

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_up(new_wavefunction_up,&wf_up[walker*ist.n_spatial_orbitals_n_up],exponential_propagation_matrix_up,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_restricted_up(density_matrix_up,kinetic_original_matrix_up,ist);
             local_energy = local_ke.real; 

             cblas_zcopy(ist.n_spatial_orbitals_n_up,new_wavefunction_up,1,&wf_up[walker*ist.n_spatial_orbitals_n_up],1);

             update_overlaps_chemistry_multi_phaseless_up(&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_up[walker*ist.n_determinants_phaseless_n_up_sq],&overlap_up[walker],&det_overlap_up[walker*ist.n_determinants_trial_phaseless],ist);

             /*Multiply By Trial Energy Outside of This Function*/
             weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);

           } /*If There Is a Sum Inside and The Matrix Is NonZero*/

    } /*Check Weights*/

   } /*Run Through Walkers*/

free(new_wavefunction_up);

free(eigs_propagation_matrix_up);
free(eigvecs_propagation_matrix_up);
free(exponential_propagation_matrix_up);
free(propagation_matrix_up);
free(density_matrix_up);
return;
}

/*********************************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_down_sparse_equilibration(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_down,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse_up,int *kinetic_ij_sparse_up,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;
   int gamma, i, l; 
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one;
   MKL_Complex16 field_final_one;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 local_ke; 
   MKL_Complex16 *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi_down(&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_down[walker],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_down,ist);

        flag_inside = 0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_down(density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);

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
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)));
                    }

                 }
                 } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {

             /*Need to Diagonalize and Exponentiate Propagation Matrix*/
             complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

             /*Exponentiate Propagation Matrix*/
             exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_down(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_down,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_restricted_down_sparse(density_matrix_down,kinetic_matrix_sparse_up,kinetic_ij_sparse_up,ist); 
             local_energy = local_ke.real; 

             cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

              update_overlaps_chemistry_multi_phaseless_down(&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_down[walker],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

             /*Multiply By Trial Energy Outside of This Function*/
             weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);

           } /*If There Is a Sum Inside and The Matrix Is NonZero*/

    } /*Check Weights*/

   } /*Run Through Walkers*/

free(new_wavefunction_down);
free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/***************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_down_sparse_production(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_down,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,double *kinetic_matrix_sparse_up,int *kinetic_ij_sparse_up,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;
   int gamma, i, l;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one;
   MKL_Complex16 field_final_one;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 local_ke;
   MKL_Complex16 *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
     if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

 /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi_down(&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_down[walker],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_down,ist);

        flag_inside = 0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_down(density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);

           /*Cap The Shift*/
           if ( Cabs(field_shift_one) > 1.0 ) {
             field_shift_one = RCmul(1.0/Cabs(field_shift_one), field_shift_one); 
           }

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
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)));
                    }

                 }
                 } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {

             /*Need to Diagonalize and Exponentiate Propagation Matrix*/
             complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

             /*Exponentiate Propagation Matrix*/
             exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_down(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_down,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_restricted_down_sparse(density_matrix_down,kinetic_matrix_sparse_up,kinetic_ij_sparse_up,ist);
             local_energy = local_ke.real;

             /*Cap Energy*/
             if ( local_energy > cns.energy_max ) {
               local_energy = cns.energy_max; 
             }
             else if ( local_energy < cns.energy_min ) {
               local_energy = cns.energy_min; 
             }  

             cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

              update_overlaps_chemistry_multi_phaseless_down(&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_down[walker],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

             /*Multiply By Trial Energy Outside of This Function*/
             weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);

           } /*If There Is a Sum Inside and The Matrix Is NonZero*/

    } /*Check Weights*/

   } /*Run Through Walkers*/

free(new_wavefunction_down);
free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
} 

/***************************************************************************************/
void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_down(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_down,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Uses the Local Energy Version of the Phaseless Approximation*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int flag_inside;
   int gamma, i, l;
   int i_spin, l_spin, i_orbital, l_orbital;
   double field_one;
   double local_energy;
   MKL_Complex16 Negative_One; Negative_One.real = -1.0; Negative_One.imag = 0.0;
   MKL_Complex16 field_shift_one;
   MKL_Complex16 field_final_one;
   MKL_Complex16 propagation_exponent_1, propagation_exponent_conjugate_1;
   MKL_Complex16 local_ke; 
   MKL_Complex16 *new_wavefunction_down;
   MKL_Complex16 *exponential_propagation_matrix_down, *propagation_matrix_down, *eigs_propagation_matrix_down, *eigvecs_propagation_matrix_down;
   MKL_Complex16 *density_matrix_down;
   long tidum;

   exponential_propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   propagation_matrix_down = (MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   eigs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals,sizeof(MKL_Complex16));
   eigvecs_propagation_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   new_wavefunction_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

     /*Check WEights*/
 if ( Cabs(weights[walker]) != 0 ) {

        czero_vec(propagation_matrix_down,ist.n_spatial_orbitals_sq);

        /*Compute Density Matrix From Wavefunction and Trial Wavefunction*/
        compute_density_matrix_chemistry_mixed_multi_down(&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_down[walker],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],density_matrix_down,ist);

        flag_inside = 0;
        /*Run Through Sites*/
        for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

          if ( Cabs(potential_eigs_fermions_first_transform[gamma]) > 0.000001 ) {

           /*Select Random Field*/
           field_one = gasdev(&tidum);

           /*Get Field Shift*/
           field_shift_one = get_shift_meanfield_down(density_matrix_down,potential_eigs_fermions_first_transform[gamma],potential_eigvecs_fermions,potential_meanfield_first_constant[gamma],gamma,ist,0);

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
                      propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital] = Cadd(propagation_matrix_down[i_orbital*ist.n_spatial_orbitals+l_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_1)));
                      propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital] = Cadd(propagation_matrix_down[l_orbital*ist.n_spatial_orbitals+i_orbital], Cmul(Negative_One,Cmul(field_final_one,propagation_exponent_conjugate_1)));
                    }

                 }
                 } /*Check Spins*/

              } /*l*/
             } /*i*/

            } /*Check Eig*/

           } /*Gamma*/

           if ( flag_inside != 0 ) {

             /*Need to Diagonalize and Exponentiate Propagation Matrix*/
             complex_eigvecs_2(propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

             /*Exponentiate Propagation Matrix*/
             exponentiate_matrix_complex(exponential_propagation_matrix_down,eigs_propagation_matrix_down,eigvecs_propagation_matrix_down,ist.n_spatial_orbitals);

             /*Propagate Wave Function*/
             propagate_wave_functions_continuous_chemistry_2_down(new_wavefunction_down,&wf_down[walker*ist.n_spatial_orbitals_n_down],exponential_propagation_matrix_down,ist);

             /*Obtain the Local Energy*/
             local_ke = compute_kinetic_energy_density_matrix_restricted_down(density_matrix_down,kinetic_original_matrix,ist);
             local_energy = local_ke.real;

             cblas_zcopy(ist.n_spatial_orbitals_n_down,new_wavefunction_down,1,&wf_down[walker*ist.n_spatial_orbitals_n_down],1);

              update_overlaps_chemistry_multi_phaseless_down(&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,&overlap_inverse_down[walker*ist.n_determinants_phaseless_n_down_sq],&overlap_down[walker],&det_overlap_down[walker*ist.n_determinants_trial_phaseless],ist);

             /*Multiply By Trial Energy Outside of This Function*/
             weights[walker] = RCmul(exp(-1*cns.dtau*local_energy), weights[walker]);

           } /*If There Is a Sum Inside and The Matrix Is NonZero*/

    } /*Check Weights*/

   } /*Run Through Walkers*/

free(new_wavefunction_down);
free(eigs_propagation_matrix_down);
free(eigvecs_propagation_matrix_down);
free(exponential_propagation_matrix_down);
free(propagation_matrix_down);
free(density_matrix_down);
return;
}

/******************************************************************************************************************************************/


void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,int_st ist,cns_st cns,long *idum){

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

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_meanfield_real_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_meanfield_first_constant, MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

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
           else {
             exponential_mean_field.real = 1.0; exponential_mean_field.imag = 0.0; 
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

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_meanfield_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix, MKL_Complex16 *potential_original_matrix,MKL_Complex16 *potential_meanfield_first_constant, MKL_Complex16 *potential_meanfield_second_constant,int_st ist,cns_st cns,long *idum){

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

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_real_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,int_st ist,cns_st cns,long *idum){

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

MKL_Complex16 determine_overlap_ratio_continuous_chemistry_multi_phaseless(MKL_Complex16 *new_wf,MKL_Complex16 *wf,MKL_Complex16 *trial_wf_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,int_st ist,int size1,int size2) {

     int i;  
     int size2_sq = size2*size2; 
     MKL_Complex16 determinant_first, determinant_second; 
     MKL_Complex16 total_overlap_first, total_overlap_second, ratio; 
     MKL_Complex16 *overlap_matrix_first, *overlap_matrix_second; 
     MKL_Complex16 *One, *Zero;  

     Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
     Zero[0].real = 0.0; Zero[0].imag = 0.0;
     One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
     One[0].real = 1.0; One[0].imag = 0.0;

     overlap_matrix_first = (MKL_Complex16 *)calloc(size2_sq,sizeof(MKL_Complex16)); 
     overlap_matrix_second = (MKL_Complex16 *)calloc(size2_sq,sizeof(MKL_Complex16));  

     total_overlap_first = total_overlap_second = Zero[0]; 
     /*Run Through All Determinants in a Multi Determinant Expansion*/
     for (i=0; i<ist.n_determinants_trial_phaseless; i++) {
         cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,size2,size2,size1,One,new_wf,size2,&trial_wf_phaseless[i*size1*size2],size2,Zero,overlap_matrix_first,size2);
         cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,size2,size2,size1,One,wf,size2,&trial_wf_phaseless[i*size1*size2],size2,Zero,overlap_matrix_second,size2);

         complex_det(overlap_matrix_first,size2,&determinant_first); 
         total_overlap_first = Cadd(total_overlap_first, Cmul(conjugate(trial_determinant_coefficients_phaseless[i]), determinant_first)); 
    
         complex_det(overlap_matrix_second,size2,&determinant_second); 
         total_overlap_second = Cadd(total_overlap_second, Cmul(conjugate(trial_determinant_coefficients_phaseless[i]), determinant_second));  
     }
     ratio = Cdiv(total_overlap_first, total_overlap_second); 


free(overlap_matrix_first); 
free(overlap_matrix_second); 
free(One);
free(Zero);  
return(ratio); 
}

/***************************************************************************************************************/

MKL_Complex16 determine_overlap_ratio_continuous_chemistry_multi_energy(MKL_Complex16 *new_wf,MKL_Complex16 *wf,MKL_Complex16 *trial_wf_energy,MKL_Complex16 *trial_determinant_coefficients_energy,int_st ist,int size1,int size2) {

     int i;
     int size2_sq = size2*size2;
     MKL_Complex16 determinant_first, determinant_second;
     MKL_Complex16 total_overlap_first, total_overlap_second, ratio;
     MKL_Complex16 *overlap_matrix_first, *overlap_matrix_second;
     MKL_Complex16 *One, *Zero;

     Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
     Zero[0].real = 0.0; Zero[0].imag = 0.0;
     One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
     One[0].real = 1.0; One[0].imag = 0.0;

     overlap_matrix_first = (MKL_Complex16 *)calloc(size2_sq,sizeof(MKL_Complex16));
     overlap_matrix_second = (MKL_Complex16 *)calloc(size2_sq,sizeof(MKL_Complex16));

     total_overlap_first = total_overlap_second = Zero[0];
     /*Run Through All Determinants in a Multi Determinant Expansion*/
     for (i=0; i<ist.n_determinants_trial_energy; i++) {
         cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,size2,size2,size1,One,new_wf,size2,&trial_wf_energy[i*size1*size2],size2,Zero,overlap_matrix_first,size2);
         cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,size2,size2,size1,One,wf,size2,&trial_wf_energy[i*size1*size2],size2,Zero,overlap_matrix_second,size2);

         complex_det(overlap_matrix_first,size2,&determinant_first);
         total_overlap_first = Cadd(total_overlap_first, Cmul(conjugate(trial_determinant_coefficients_energy[i]), determinant_first));

         complex_det(overlap_matrix_second,size2,&determinant_second);
         total_overlap_second = Cadd(total_overlap_second, Cmul(conjugate(trial_determinant_coefficients_energy[i]), determinant_second));
     }
     ratio = Cdiv(total_overlap_first, total_overlap_second);


free(overlap_matrix_first);
free(overlap_matrix_second);
free(One);
free(Zero);
return(ratio);
}

/**************************************************************************************/

void update_overlaps_chemistry_multi_phaseless_both(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,int_st ist) {

    /*Get Product of Trial and Actual*/
    (*overlap_total) = compute_overlap_inverse_multi_phaseless_both(trial_wf_up_phaseless,trial_wf_down_phaseless,wf_up,wf_down,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,ist);

return;
}

/**************************************************************************************/

void update_overlaps_chemistry_multi_phaseless_up(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,MKL_Complex16 *det_overlap_up,int_st ist) {

    /*Get Product of Trial and Actual*/
    (*overlap_up) = compute_overlap_inverse_multi_phaseless_up(trial_wf_up_phaseless,wf_up,trial_determinant_coefficients_phaseless,overlap_inverse_up,det_overlap_up,ist);

return;
}

/**************************************************************************************/

void update_overlaps_chemistry_multi_phaseless_down(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_down,MKL_Complex16 *det_overlap_down,int_st ist) {

    /*Get Product of Trial and Actual*/
    (*overlap_down) = compute_overlap_inverse_multi_phaseless_down(trial_wf_down_phaseless,wf_down,trial_determinant_coefficients_phaseless,overlap_inverse_down,det_overlap_down,ist);

return;
}

/**************************************************************************************/

void update_overlaps_chemistry_multi_energy_both(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,int_st ist) {

    /*Get Product of Trial and Actual*/
    (*overlap_total) = compute_overlap_inverse_multi_energy_both(trial_wf_up_energy,trial_wf_down_energy,wf_up,wf_down,trial_determinant_coefficients_energy,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,ist);

return;
}

/**************************************************************************************/

void update_overlaps_chemistry_multi_energy_up(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,MKL_Complex16 *det_overlap_up,int_st ist) {

    /*Get Product of Trial and Actual*/
    (*overlap_up) = compute_overlap_inverse_multi_energy_up(trial_wf_up_energy,wf_up,trial_determinant_coefficients_energy,overlap_inverse_up,det_overlap_up,ist);

return;
}

/**************************************************************************************/

void update_overlaps_chemistry_multi_energy_down(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_down,MKL_Complex16 *det_overlap_down,int_st ist) {

    /*Get Product of Trial and Actual*/
    (*overlap_down) = compute_overlap_inverse_multi_energy_down(trial_wf_down_energy,wf_down,trial_determinant_coefficients_energy,overlap_inverse_down,det_overlap_down,ist);

return;
}

/*************************************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_meanfield_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,MKL_Complex16 *potential_eigvecs_fermions,int_st ist,cns_st cns,long *idum){

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

void propagate_forwards_potential_continuous_mixed_chemistry_meanfield_real_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_meanfield_first_constant,MKL_Complex16 *potential_meanfield_second_constant,MKL_Complex16 *potential_eigvecs_fermions,int_st ist,cns_st cns,long *idum){

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

void compute_density_matrix_chemistry_mixed_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *density_matrix_up,MKL_Complex16 *density_matrix_down,int_st ist) {

   /*Computes the Spin Up and Down Density Matrices for a MultiDeterminant Trial WF*/
   /*Sort Based Upon Number of electrons*/
   if ( ist.n_up > 0 ) {
    if (ist.n_down > 0 ) {
      compute_density_matrix_chemistry_mixed_multi_both(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,density_matrix_up,density_matrix_down,ist); 
    }
    else {
      compute_density_matrix_chemistry_mixed_multi_up(wf_up,trial_wf_up_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_total,det_overlap_up,density_matrix_up,ist);
    }
   }
   else {
    compute_density_matrix_chemistry_mixed_multi_down(wf_down,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_down,overlap_total,det_overlap_down,density_matrix_down,ist);
   }       

return;
}

/*********************************************************************************************************/

void compute_density_matrix_chemistry_mixed_multi_both(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *density_matrix_up,MKL_Complex16 *density_matrix_down,int_st ist) {

   /*Computes the Spin Up and Down Density Matrices for a MultiDeterminant Trial WF*/
   int i; 
   MKL_Complex16 coefficient_up, coefficient_down;
   MKL_Complex16 total_overlap_up;
   MKL_Complex16 total_overlap_down;
   MKL_Complex16 *local_density_matrix_up, *local_density_matrix_down;
   MKL_Complex16 *stored_product1_up, *stored_product1_down;
   MKL_Complex16 *One, *Zero;

   stored_product1_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
   stored_product1_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

   local_density_matrix_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   local_density_matrix_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   One[0].real = 1.0; One[0].imag = 0.0;
   Zero[0].real = 0.0; Zero[0].imag = 0.0;

    /*Determine Full Overlap*/
    czero_vec(density_matrix_up,ist.n_spatial_orbitals_sq);
    czero_vec(density_matrix_down,ist.n_spatial_orbitals_sq);

    total_overlap_up.real = total_overlap_up.imag = 0.0; 
    total_overlap_down.real = total_overlap_down.imag = 0.0; 
    for (i=0; i<ist.n_determinants_trial_phaseless; i++) {

      /*Overlap Inverse Time Trial WF*/
      cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_up,ist.n_spatial_orbitals,ist.n_up,One,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up,&trial_wf_up_phaseless[i*ist.n_spatial_orbitals_n_up],ist.n_up,Zero,stored_product1_up,ist.n_spatial_orbitals);
      cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_down,ist.n_spatial_orbitals,ist.n_down,One,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down,&trial_wf_down_phaseless[i*ist.n_spatial_orbitals_n_down],ist.n_down,Zero,stored_product1_down,ist.n_spatial_orbitals);

      /*Now Get Total Density Matrix*/
      cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up,One,wf_up,ist.n_up,stored_product1_up,ist.n_spatial_orbitals,Zero,local_density_matrix_up,ist.n_spatial_orbitals);
      cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down,One,wf_down,ist.n_down,stored_product1_down,ist.n_spatial_orbitals,Zero,local_density_matrix_down,ist.n_spatial_orbitals);

      /*Get Weighting Coefficient for Determinant*/
      coefficient_up = Cmul(det_overlap_up[i], conjugate(trial_determinant_coefficients_phaseless[i]));
      total_overlap_up = Cadd(total_overlap_up, coefficient_up); 

      coefficient_down = Cmul(det_overlap_down[i], conjugate(trial_determinant_coefficients_phaseless[i]));
      total_overlap_down = Cadd(total_overlap_down, coefficient_down); 

      /*Add Appropriately Weighted Density Matrix To Total*/
      /*Check If This Works*/
      cblas_zaxpy(ist.n_spatial_orbitals_sq,&coefficient_up,local_density_matrix_up,1,density_matrix_up,1);
      cblas_zaxpy(ist.n_spatial_orbitals_sq,&coefficient_down,local_density_matrix_down,1,density_matrix_down,1);
    }
    /*Take Inverse*/
    total_overlap_up = Cdiv(One[0], total_overlap_up); 
    total_overlap_down = Cdiv(One[0], total_overlap_down); 

    /*Divide By Overall Weighting*/
    cblas_zscal(ist.n_spatial_orbitals_sq,&total_overlap_up,density_matrix_up,1);
    cblas_zscal(ist.n_spatial_orbitals_sq,&total_overlap_down,density_matrix_down,1);


free(stored_product1_up);
free(stored_product1_down);
free(local_density_matrix_up);
free(local_density_matrix_down);
return;
}

/*********************************************************************************************************/

void compute_density_matrix_chemistry_mixed_multi_up(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,MKL_Complex16 *det_overlap_up,MKL_Complex16 *density_matrix_up,int_st ist) {

   /*Computes the Spin Up and Down Density Matrices for a MultiDeterminant Trial WF*/
   int i;
   MKL_Complex16 coefficient_up; 
   MKL_Complex16 total_overlap_up;
   MKL_Complex16 *local_density_matrix_up; 
   MKL_Complex16 *stored_product1_up; 
   MKL_Complex16 *One, *Zero;

   stored_product1_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
   local_density_matrix_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   One[0].real = 1.0; One[0].imag = 0.0;
   Zero[0].real = 0.0; Zero[0].imag = 0.0;

    /*Determine Full Overlap*/
    czero_vec(density_matrix_up,ist.n_spatial_orbitals_sq);

    for (i=0; i<ist.n_determinants_trial_phaseless; i++) {

      /*Overlap Inverse Time Trial WF*/
      cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_up,ist.n_spatial_orbitals,ist.n_up,One,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up,&trial_wf_up_phaseless[i*ist.n_spatial_orbitals_n_up],ist.n_up,Zero,stored_product1_up,ist.n_spatial_orbitals);

      /*Now Get Total Density Matrix*/
      cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up,One,wf_up,ist.n_up,stored_product1_up,ist.n_spatial_orbitals,Zero,local_density_matrix_up,ist.n_spatial_orbitals);

      /*Get Weighting Coefficient for Determinant*/
      coefficient_up = Cmul(det_overlap_up[i], conjugate(trial_determinant_coefficients_phaseless[i]));

      /*Add Appropriately Weighted Density Matrix To Total*/
      /*Check If This Works*/
      cblas_zaxpy(ist.n_spatial_orbitals_sq,&coefficient_up,local_density_matrix_up,1,density_matrix_up,1);
    }
    /*Take Inverse*/
    total_overlap_up = Cdiv(One[0], overlap_up[0]);

    /*Divide By Overall Weighting*/
    cblas_zscal(ist.n_spatial_orbitals_sq,&total_overlap_up,density_matrix_up,1);


free(stored_product1_up);
free(local_density_matrix_up);
return;
}

/*********************************************************************************************************/

void compute_density_matrix_chemistry_mixed_multi_down(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_down,MKL_Complex16 *det_overlap_down,MKL_Complex16 *density_matrix_down,int_st ist) {

   /*Computes the Spin Up and Down Density Matrices for a MultiDeterminant Trial WF*/
   int i;
   MKL_Complex16 coefficient_down;  
   MKL_Complex16 total_overlap_down;
   MKL_Complex16 *local_density_matrix_down;  
   MKL_Complex16 *stored_product1_down;  
   MKL_Complex16 *One, *Zero;

   stored_product1_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));
   local_density_matrix_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   One[0].real = 1.0; One[0].imag = 0.0;
   Zero[0].real = 0.0; Zero[0].imag = 0.0;

    /*Determine Full Overlap*/
    czero_vec(density_matrix_down,ist.n_spatial_orbitals_sq);

    for (i=0; i<ist.n_determinants_trial_phaseless; i++) {

      /*Overlap Inverse Time Trial WF*/
      cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_down,ist.n_spatial_orbitals,ist.n_down,One,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down,&trial_wf_down_phaseless[i*ist.n_spatial_orbitals_n_down],ist.n_down,Zero,stored_product1_down,ist.n_spatial_orbitals);

      /*Now Get Total Density Matrix*/
      cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down,One,wf_down,ist.n_down,stored_product1_down,ist.n_spatial_orbitals,Zero,local_density_matrix_down,ist.n_spatial_orbitals);
      
      /*Get Weighting Coefficient for Determinant*/
      coefficient_down = Cmul(det_overlap_down[i], conjugate(trial_determinant_coefficients_phaseless[i]));

      /*Add Appropriately Weighted Density Matrix To Total*/
      /*Check If This Works*/
      cblas_zaxpy(ist.n_spatial_orbitals_sq,&coefficient_down,local_density_matrix_down,1,density_matrix_down,1);
    }
    /*Take Inverse*/
    total_overlap_down = Cdiv(One[0], overlap_down[0]);

    /*Divide By Overall Weighting*/
    cblas_zscal(ist.n_spatial_orbitals_sq,&total_overlap_down,density_matrix_down,1);


free(stored_product1_down);
free(local_density_matrix_down);
return;
}

/***************************************************************************************************************************************************/

void propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_multi_unrestricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *weights,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *kinetic_original_matrix_down,MKL_Complex16 *potential_original_matrix_up,MKL_Complex16 *potential_original_matrix_down,MKL_Complex16 *potential_original_matrix_updown,int_st ist,cns_st cns,long *idum){

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

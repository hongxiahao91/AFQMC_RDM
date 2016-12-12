#include "afqmc.h"

/***********************************************************************************************/

void compute_energy_chemistry_mixed_restricted_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,double *kinetic_matrix_sparse,double *potential_matrix_sparse,int *kinetic_ij_sparse,int *potential_ij_sparse,MKL_Complex16 *energy_shifts,MKL_Complex16 *total_potential_energy,MKL_Complex16 *total_kinetic_energy,MKL_Complex16 *total_energy_2,MKL_Complex16 *total_weight,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *pe_error,MKL_Complex16 *ke,MKL_Complex16 *ke_error,MKL_Complex16 *te,MKL_Complex16 *te_error,MKL_Complex16 *tw,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   //Computes the Energy of the System Using Mixed Estimator for a Restricted Situation
   int walker, i, j, k;
   int number_of_accumulated_steps;  
   int number_walkers_active = 0; 
   double weight_factor = 1.0;
   double energy_error, kinetic_energy_error_real, potential_energy_error_real;
   MKL_Complex16 kinetic_energy_error_complex, potential_energy_error_complex, total_energy_error_complex; 
   MKL_Complex16 *average_wave_function_up, *average_wave_function_down;
   MKL_Complex16 local_pe, local_ke, local_pe_squared, local_ke_squared;
   MKL_Complex16 Zero; 
   MKL_Complex16 potential_energy, kinetic_energy, total_energy;
   MKL_Complex16 potential_energy_squared, kinetic_energy_squared; 
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator;
   MKL_Complex16 energy_average; 

   Zero.real = Zero.imag = 0.0; 
   potential_energy.real = potential_energy.imag = 0.0; 
   kinetic_energy.real = kinetic_energy.imag = 0.0; 
   total_energy.real = total_energy.imag = 0.0; 
   potential_energy_squared.real = potential_energy_squared.imag = 0.0; 
   kinetic_energy_squared.real = kinetic_energy_squared.imag; 
   weight_denominator.real = weight_denominator.imag = 0.0; 

   average_wave_function_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_up*sizeof(MKL_Complex16));
   average_wave_function_down=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_down*sizeof(MKL_Complex16));

   czero_vec(average_wave_function_up,ist.n_spatial_orbitals*ist.n_up);
   czero_vec(average_wave_function_down,ist.n_spatial_orbitals*ist.n_down);

   //Obtain Population Control Rescaling Factor
   for (i = 0; i < 10 ; i++) {
    weight_factor *= weight_rescaling[i];
   }

   for (walker = 0; walker < ist.n_walkers ; walker++ ) {

     //Check That Walker Weights Are Non-Zero
     if ( Cabs(weights[walker]) != 0 ) {
        number_walkers_active++; 

        //Accumulate Weight in Denominator/
        current_weight = RCmul(weight_factor, weights[walker]);
        weight_denominator = Cadd(weight_denominator, current_weight);

        /*Compute The Energy If Non Multi*/
        compute_walker_energy_chemistry_mixed_restricted_sparse(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],kinetic_matrix_sparse,potential_matrix_sparse,kinetic_ij_sparse,potential_ij_sparse,&local_ke,&local_pe,ist); 

        //Obtain Kinetic Energy
        local_ke_squared = Cmul(local_ke, conjugate(local_ke)); 
        local_ke = Cmul(local_ke, current_weight);
        local_ke_squared = Cmul(local_ke_squared, current_weight); 

        //Get Potential Energy
        local_pe_squared = Cmul(local_pe, conjugate(local_pe)); 
        local_pe = Cmul(local_pe, current_weight); 
        local_pe_squared = Cmul(local_pe_squared, current_weight); 

        //Add to Total
        potential_energy = Cadd(local_pe, potential_energy);
        potential_energy_squared = Cadd(local_pe_squared, potential_energy_squared); 
        kinetic_energy = Cadd(local_ke, kinetic_energy);
        kinetic_energy_squared = Cadd(kinetic_energy_squared, local_ke_squared); 

        //Obtain Average Wavefunction As Well
        for (i=0; i<ist.n_spatial_orbitals; i++) {
         for (j=0; j<ist.n_up; j++) {
           average_wave_function_up[i*ist.n_up+j] = Cadd(average_wave_function_up[i*ist.n_up+j], Cmul(current_weight, wf_up[walker*ist.n_spatial_orbitals_n_up+i*ist.n_up+j]));
         }
         for (j=0; j<ist.n_down; j++) {
           average_wave_function_down[i*ist.n_down+j] = Cadd(average_wave_function_down[i*ist.n_down+j], Cmul(current_weight, wf_down[walker*ist.n_spatial_orbitals_n_down+i*ist.n_down+j]));
         }
        }

      } //Non-Zero Walker Weights
     }  //Walkers

    //Accumulate the Energy of Each Block To Get Overall Statistics
    if ( step >= ist.n_steps_equilibration ) {
     // accumulate_energy_blocks(potential_energy, kinetic_energy, Cadd(potential_energy, kinetic_energy), weight_denominator, total_potential_energy, total_kinetic_energy, total_energy_2, total_weight); 
    }

   //Divide By Total Weight
   potential_energy = Cdiv(potential_energy, weight_denominator);
   potential_energy_squared = Cdiv(potential_energy_squared, weight_denominator); 
   kinetic_energy = Cdiv(kinetic_energy, weight_denominator);
   kinetic_energy_squared = Cdiv(kinetic_energy_squared, weight_denominator); 
   total_energy = Cadd(potential_energy, kinetic_energy);

   //Also Add Initial Energy Shift from MolPro
   total_energy.real += energy_shifts[0].real;

   //Print Desired Information*******************************************
   //Find Errors
   potential_energy_error_complex = Csub(potential_energy_squared, Cmul(potential_energy, conjugate(potential_energy)));
   if ( potential_energy_error_complex.real >  0 ) {
     potential_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * potential_energy_error_complex.real);
   }
   else {
     potential_energy_error_complex.real = 0.0;
   }

   kinetic_energy_error_complex = Csub(kinetic_energy_squared, Cmul(kinetic_energy, conjugate(kinetic_energy)));
   if ( kinetic_energy_error_complex.real > 0 ) {
     kinetic_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * kinetic_energy_error_complex.real);
   }
   else {
    kinetic_energy_error_complex.real = 0.0;
   }
   total_energy_error_complex.real = sqrt(potential_energy_error_complex.real*potential_energy_error_complex.real +   kinetic_energy_error_complex.real*kinetic_energy_error_complex.real);

   //Print Desired Information*******************************************
   *ke = kinetic_energy; 
   *pe = potential_energy; 
   *te = total_energy; 
   *ke_error = kinetic_energy_error_complex; 
   *pe_error = potential_energy_error_complex; 
   *te_error = total_energy_error_complex; 
   *tw = weight_denominator; 

  //Obtain and Print Average WFs
  for (i=0; i<ist.n_spatial_orbitals; i++) {
    for (j=0; j<ist.n_up; j++) {
      average_wave_function_up[i*ist.n_up+j] = Cdiv(average_wave_function_up[i*ist.n_up+j], weight_denominator);
   } 
 } 
  for (i=0; i<ist.n_spatial_orbitals; i++) {
    for (j=0; j<ist.n_down; j++) {
      average_wave_function_down[i*ist.n_down+j] = Cdiv(average_wave_function_down[i*ist.n_down+j], weight_denominator);
    } 
  }

 //Find Normalized WF
 get_normalized_wf_up_chemistry(normalized_wave_function_up,average_wave_function_up,ist);
 get_normalized_wf_down_chemistry(normalized_wave_function_down,average_wave_function_down,ist); 

  //Find Statistics So Far*******************************************
  /*if ( step >= ist.n_steps_equilibration ) {
    number_of_accumulated_steps = 1+(int)((step-ist.n_steps_equilibration)/(ist.n_steps_energy)); 
    fprintf(pf7, "%d\t", step); 

    //Get PE
    energy_average = Cdiv(total_potential_energy[0], total_weight[0]);  
    potential_energy_error_real = error_real(total_potential_energy[1].real, total_potential_energy[2].real, number_of_accumulated_steps); 
    fprintf(pf7, "%f   %f\t", energy_average.real, potential_energy_error_real);   

    //Get KE
    energy_average = Cdiv(total_kinetic_energy[0], total_weight[0]); 
    kinetic_energy_error_real = error_real(total_kinetic_energy[1].real, total_kinetic_energy[2].real, number_of_accumulated_steps);
    fprintf(pf7, "%f   %f\t", energy_average.real, kinetic_energy_error_real);

    //Show Shift
    fprintf(pf7, "%f   %f\t", energy_shifts[0].real, 0.0); 

    //Get Total
    energy_average = Cdiv(total_energy_2[0], total_weight[0]);
    energy_error = sqrt(potential_energy_error_real*potential_energy_error_real+kinetic_energy_error_real*kinetic_energy_error_real);  
    //energy_error = error_real(total_energy_2[1].r, total_energy_2[2].r, number_of_accumulated_steps);
    fprintf(pf7, "%f   %f\n", energy_average.real+energy_shifts[0].real, energy_error);
    fflush(pf7); 
   }
*/

free(average_wave_function_up); 
free(average_wave_function_down); 
return; 
}

/*************************************************************************************************/

void compute_energy_chemistry_mixed_multi_unrestricted_sparse_equilibration(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,MKL_Complex16 *energies,MKL_Complex16 *energies2,MKL_Complex16 *energies3,MKL_Complex16 *overlaps,MKL_Complex16 *density_matrices_up,MKL_Complex16 *density_matrices_down,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *energy_shifts,MKL_Complex16 *total_potential_energy,MKL_Complex16 *total_kinetic_energy,MKL_Complex16 *total_energy_2,MKL_Complex16 *total_weight,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *pe_error,MKL_Complex16 *ke,MKL_Complex16 *ke_error,MKL_Complex16 *te,MKL_Complex16 *te_error,MKL_Complex16 *tw,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   /*Computes the Energy of the System Using Mixed Estimator for a Restricted Situation*/
   int walker, i, j, k;
   int number_of_accumulated_steps;
   int number_walkers_active = 0;
   double weight_factor = 1.0;
   double energy_error, kinetic_energy_error_real, potential_energy_error_real;
   MKL_Complex16 kinetic_energy_error_complex, potential_energy_error_complex, total_energy_error_complex;
   MKL_Complex16 *average_wave_function_up, *average_wave_function_down;
   MKL_Complex16 *density_matrix_up, *density_matrix_down; 
   MKL_Complex16 local_pe, local_ke, local_pe_squared, local_ke_squared;
   MKL_Complex16 local_total; 
   MKL_Complex16 Zero;
   MKL_Complex16 potential_energy, kinetic_energy, total_energy;
   MKL_Complex16 potential_energy_squared, kinetic_energy_squared;
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator;
   MKL_Complex16 energy_average;

   /*HERE*/
   Zero.real = Zero.imag = 0.0;
   potential_energy.real = potential_energy.imag = 0.0;
   kinetic_energy.real = kinetic_energy.imag = 0.0;
   total_energy.real = total_energy.imag = 0.0;
   potential_energy_squared.real = potential_energy_squared.imag = 0.0;
   kinetic_energy_squared.real = kinetic_energy_squared.imag;
   weight_denominator.real = weight_denominator.imag = 0.0;

   average_wave_function_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_up*sizeof(MKL_Complex16));
   average_wave_function_down=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_down*sizeof(MKL_Complex16));

   density_matrix_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16)); 
   density_matrix_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16)); 

   czero_vec(average_wave_function_up,ist.n_spatial_orbitals*ist.n_up);
   czero_vec(average_wave_function_down,ist.n_spatial_orbitals*ist.n_down);

   /*Obtain Population Control Rescaling Factor*/
   for (i = 0; i < 10 ; i++) {
    weight_factor *= weight_rescaling[i];
   }

   for (walker = 0; walker < ist.n_walkers ; walker++ ) {

     /*Check That Walker Weights Are Non-Zero*/
     if ( Cabs(weights[walker]) != 0 ) {
        number_walkers_active++;

        /*Accumulate Weight in Denominator*/
        current_weight = RCmul(weight_factor, weights[walker]);

        /*Compute the Energies Using the Mixed Estimator, but for a Multideterminant WF*/
        /*Using the Full Trial Energy Wavefunction To Get Exact Energies To Some Degree*/
        if ( ist.n_up > 0 ) {
         if ( ist.n_down > 0 ) {
              compute_walker_energy_chemistry_mixed_multi_energy_unrestricted_sparse(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_energy,trial_wf_down_energy,&energies2[3*walker],&energies3[3*walker],trial_determinant_coefficients_energy,density_matrix_up,density_matrix_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,&local_ke,&local_pe,ist);
         }
         else {
              compute_walker_energy_chemistry_mixed_multi_energy_up_sparse(&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_energy,trial_determinant_coefficients_energy,kinetic_matrix_sparse_up,kinetic_ij_sparse_up,&local_ke,&local_pe,ist);
         }
        }
        else {
          if ( ist.n_down > 0 ) {
            compute_walker_energy_chemistry_mixed_multi_energy_down_sparse(&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_energy,trial_determinant_coefficients_energy,kinetic_matrix_sparse_down,kinetic_ij_sparse_down,&local_ke,&local_pe,ist);
          }
        }

        local_total = Cadd(local_ke, local_pe); 
        local_total.real += energy_shifts[0].real + energy_shifts[1].real + energy_shifts[2].real + energy_shifts[3].real + energy_shifts[4].real;    

        /*Obtain Kinetic Energy*/
        local_ke_squared = Cmul(local_ke, conjugate(local_ke));
        local_ke = Cmul(local_ke, current_weight);
        local_ke_squared = Cmul(local_ke_squared, current_weight);

        /*Get Potential Energy*/
        local_pe_squared = Cmul(local_pe, conjugate(local_pe));
        local_pe = Cmul(local_pe, current_weight);
        local_pe_squared = Cmul(local_pe_squared, current_weight);

        /*Add to Total*/
        potential_energy = Cadd(local_pe, potential_energy);
        potential_energy_squared = Cadd(local_pe_squared, potential_energy_squared);
        kinetic_energy = Cadd(local_ke, kinetic_energy);
        kinetic_energy_squared = Cadd(kinetic_energy_squared, local_ke_squared);

          /*Obtain Average Wavefunction As Well*/
          for (i=0; i<ist.n_spatial_orbitals; i++) {
           for (j=0; j<ist.n_up; j++) {
             average_wave_function_up[i*ist.n_up+j] = Cadd(average_wave_function_up[i*ist.n_up+j], Cmul(current_weight, wf_up[walker*ist.n_spatial_orbitals_n_up+i*ist.n_up+j]));
           }
           for (j=0; j<ist.n_down; j++) {
             average_wave_function_down[i*ist.n_down+j] = Cadd(average_wave_function_down[i*ist.n_down+j], Cmul(current_weight, wf_down[walker*ist.n_spatial_orbitals_n_down+i*ist.n_down+j]));
           }
          }

        /*Weight Denominator*/ 
        weight_denominator = Cadd(weight_denominator, current_weight);

        energies[walker*3] = local_ke; 
        energies[walker*3+1] = local_pe; 
        energies[walker*3+2] = Cadd(local_pe, local_ke); 
        overlaps[walker] = overlap_total[walker];  

        for (i=0; i<ist.n_spatial_orbitals_sq; i++) {
          density_matrices_up[walker*ist.n_spatial_orbitals_sq+i] = density_matrix_up[i]; 
          density_matrices_down[walker*ist.n_spatial_orbitals_sq+i] = density_matrix_down[i]; 
        }


      } /*Non-Zero Walker Weights*/
     }  /*Walkers*/


    /*Accumulate the Energy of Each Block To Get Overall Statistics*/
    if ( step >= ist.n_steps_equilibration ) {
      accumulate_energy_blocks(potential_energy, kinetic_energy, Cadd(potential_energy, kinetic_energy), weight_denominator, total_potential_energy, total_kinetic_energy, total_energy_2, total_weight);
    }

   /*Divide By Total Weight*/
   potential_energy = Cdiv(potential_energy, weight_denominator);
   potential_energy_squared = Cdiv(potential_energy_squared, weight_denominator);
   kinetic_energy = Cdiv(kinetic_energy, weight_denominator);
   kinetic_energy_squared = Cdiv(kinetic_energy_squared, weight_denominator);
   total_energy = Cadd(potential_energy, kinetic_energy);

   /*Also Add Initial Energy Shift from MolPro*/
   total_energy.real += energy_shifts[0].real;

   /*Print Desired Information*******************************************/
   /*Find Errors*/
   potential_energy_error_complex = Csub(potential_energy_squared, Cmul(potential_energy, conjugate(potential_energy)));
   if ( potential_energy_error_complex.real >  0 ) {
     potential_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * potential_energy_error_complex.real);
   }
   else {
     potential_energy_error_complex.real = 0.0;
   }

   kinetic_energy_error_complex = Csub(kinetic_energy_squared, Cmul(kinetic_energy, conjugate(kinetic_energy)));
   if ( kinetic_energy_error_complex.real > 0 ) {
     kinetic_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * kinetic_energy_error_complex.real);
   }
   else {
    kinetic_energy_error_complex.real = 0.0;
   }
   total_energy_error_complex.real = sqrt(potential_energy_error_complex.real*potential_energy_error_complex.real +   kinetic_energy_error_complex.real*kinetic_energy_error_complex.real);

   /*Print Desired Information*******************************************/
   *ke = kinetic_energy;
   *pe = potential_energy;
   *te = total_energy;
   *ke_error = kinetic_energy_error_complex;
   *pe_error = potential_energy_error_complex;
   *te_error = total_energy_error_complex;
   *tw = weight_denominator;

  /*Obtain and Print Average WFs*/
  for (i=0; i<ist.n_spatial_orbitals; i++) {
    for (j=0; j<ist.n_up; j++) {
      average_wave_function_up[i*ist.n_up+j] = Cdiv(average_wave_function_up[i*ist.n_up+j], weight_denominator);
   }
 }
  for (i=0; i<ist.n_spatial_orbitals; i++) {
    for (j=0; j<ist.n_down; j++) {
      average_wave_function_down[i*ist.n_down+j] = Cdiv(average_wave_function_down[i*ist.n_down+j], weight_denominator);
    }
  }

 /*Find Normalized WF*/
 get_normalized_wf_up_chemistry(normalized_wave_function_up,average_wave_function_up,ist);
 get_normalized_wf_down_chemistry(normalized_wave_function_down,average_wave_function_down,ist);

free(average_wave_function_up);
free(average_wave_function_down);
return;
}

/*********************************************************************************************************************************************************/

void compute_energy_chemistry_mixed_multi_unrestricted_sparse_production(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,MKL_Complex16 *energies,MKL_Complex16 *energies2,MKL_Complex16 *energies3,MKL_Complex16 *overlaps,int *total_energy_flag,MKL_Complex16 *density_matrices_up,MKL_Complex16 *density_matrices_down,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *energy_shifts,MKL_Complex16 *total_potential_energy,MKL_Complex16 *total_kinetic_energy,MKL_Complex16 *total_energy_2,MKL_Complex16 *total_weight,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *pe_error,MKL_Complex16 *ke,MKL_Complex16 *ke_error,MKL_Complex16 *te,MKL_Complex16 *te_error,MKL_Complex16 *tw,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   /*Computes the Energy of the System Using Mixed Estimator for a Restricted Situation*/
   int walker, i, j, k;
   int number_of_accumulated_steps;
   int number_walkers_active = 0;
   double weight_factor = 1.0;
   double energy_error, kinetic_energy_error_real, potential_energy_error_real;
   MKL_Complex16 kinetic_energy_error_complex, potential_energy_error_complex, total_energy_error_complex;
   MKL_Complex16 *average_wave_function_up, *average_wave_function_down;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   MKL_Complex16 local_pe, local_ke, local_pe_squared, local_ke_squared;
   MKL_Complex16 local_total;
   MKL_Complex16 Zero;
   MKL_Complex16 potential_energy, kinetic_energy, total_energy;
   MKL_Complex16 potential_energy_squared, kinetic_energy_squared;
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator;
   MKL_Complex16 energy_average;

   /*HERE*/
   Zero.real = Zero.imag = 0.0;
   potential_energy.real = potential_energy.imag = 0.0;
   kinetic_energy.real = kinetic_energy.imag = 0.0;
   total_energy.real = total_energy.imag = 0.0;
   potential_energy_squared.real = potential_energy_squared.imag = 0.0;
   kinetic_energy_squared.real = kinetic_energy_squared.imag;
   weight_denominator.real = weight_denominator.imag = 0.0;

   average_wave_function_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_up*sizeof(MKL_Complex16));
   average_wave_function_down=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_down*sizeof(MKL_Complex16));

   density_matrix_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
   density_matrix_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

 czero_vec(average_wave_function_up,ist.n_spatial_orbitals*ist.n_up);
   czero_vec(average_wave_function_down,ist.n_spatial_orbitals*ist.n_down);

   /*Obtain Population Control Rescaling Factor*/
   for (i = 0; i < 10 ; i++) {
    weight_factor *= weight_rescaling[i];
   }

   for (walker = 0; walker < ist.n_walkers ; walker++ ) {

     /*Check That Walker Weights Are Non-Zero*/
     if ( Cabs(weights[walker]) != 0 ) {
        number_walkers_active++;

        /*Accumulate Weight in Denominator*/
        current_weight = RCmul(weight_factor, weights[walker]);

        /*Compute the Energies Using the Mixed Estimator, but for a Multideterminant WF*/
        /*Using the Full Trial Energy Wavefunction To Get Exact Energies To Some Degree*/
        if ( ist.n_up > 0 ) {
         if ( ist.n_down > 0 ) {
              compute_walker_energy_chemistry_mixed_multi_energy_unrestricted_sparse(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_energy,trial_wf_down_energy,&energies2[3*walker],&energies3[3*walker],trial_determinant_coefficients_energy,density_matrix_up,density_matrix_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,&local_ke,&local_pe,ist);
         }
         else {
              compute_walker_energy_chemistry_mixed_multi_energy_up_sparse(&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_energy,trial_determinant_coefficients_energy,kinetic_matrix_sparse_up,kinetic_ij_sparse_up,&local_ke,&local_pe,ist);
         }
        }
        else {
          if ( ist.n_down > 0 ) {
            compute_walker_energy_chemistry_mixed_multi_energy_down_sparse(&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_energy,trial_determinant_coefficients_energy,kinetic_matrix_sparse_down,kinetic_ij_sparse_down,&local_ke,&local_pe,ist);
          }
        }

        local_total = Cadd(local_ke, local_pe);

        /*Using The Energy Cap Causes Problems If Energy Not Set Right*/ 
        total_energy_flag[walker] = 0;
        printf("cns energy max %f cns energy min %f\n", cns.energy_max, cns.energy_min); fflush(NULL); 
 
        /*Put Cap on Value*/
        /*if ( local_total.real > cns.energy_max ) {
           current_weight.real = current_weight.imag = 0.0; 
           total_energy_flag[walker] = 1; 
        }
        else if ( local_total.real < cns.energy_min ) {
           current_weight.real = current_weight.imag = 0.0; 
           total_energy_flag[walker] = 1; 
        } */ 

        local_total.real += energy_shifts[0].real + energy_shifts[1].real + energy_shifts[2].real + energy_shifts[3].real + energy_shifts[4].real;

        /*Obtain Kinetic Energy*/
        local_ke_squared = Cmul(local_ke, conjugate(local_ke));
        local_ke = Cmul(local_ke, current_weight);
        local_ke_squared = Cmul(local_ke_squared, current_weight);

        /*Get Potential Energy*/
        local_pe_squared = Cmul(local_pe, conjugate(local_pe));
        local_pe = Cmul(local_pe, current_weight);
        local_pe_squared = Cmul(local_pe_squared, current_weight);

        /*Add to Total*/
        potential_energy = Cadd(local_pe, potential_energy);
        potential_energy_squared = Cadd(local_pe_squared, potential_energy_squared);
        kinetic_energy = Cadd(local_ke, kinetic_energy);
        kinetic_energy_squared = Cadd(kinetic_energy_squared, local_ke_squared);

          /*Obtain Average Wavefunction As Well*/
          for (i=0; i<ist.n_spatial_orbitals; i++) {
           for (j=0; j<ist.n_up; j++) {
             average_wave_function_up[i*ist.n_up+j] = Cadd(average_wave_function_up[i*ist.n_up+j], Cmul(current_weight, wf_up[walker*ist.n_spatial_orbitals_n_up+i*ist.n_up+j]));
           }
           for (j=0; j<ist.n_down; j++) {
             average_wave_function_down[i*ist.n_down+j] = Cadd(average_wave_function_down[i*ist.n_down+j], Cmul(current_weight, wf_down[walker*ist.n_spatial_orbitals_n_down+i*ist.n_down+j]));
           }
          }

        /*Weight Denominator*/
        weight_denominator = Cadd(weight_denominator, current_weight);

        energies[walker*3] = local_ke;
        energies[walker*3+1] = local_pe;
        energies[walker*3+2] = Cadd(local_pe, local_ke);
        overlaps[walker] = overlap_total[walker];

        for (i=0; i<ist.n_spatial_orbitals_sq; i++) {
          density_matrices_up[walker*ist.n_spatial_orbitals_sq+i] = density_matrix_up[i];
          density_matrices_down[walker*ist.n_spatial_orbitals_sq+i] = density_matrix_down[i];
        }


      } /*Non-Zero Walker Weights*/
 }  /*Walkers*/


    /*Accumulate the Energy of Each Block To Get Overall Statistics*/
    if ( step >= ist.n_steps_equilibration ) {
      accumulate_energy_blocks(potential_energy, kinetic_energy, Cadd(potential_energy, kinetic_energy), weight_denominator, total_potential_energy, total_kinetic_energy, total_energy_2, total_weight);
    }

   /*Divide By Total Weight*/
   potential_energy = Cdiv(potential_energy, weight_denominator);
   potential_energy_squared = Cdiv(potential_energy_squared, weight_denominator);
   kinetic_energy = Cdiv(kinetic_energy, weight_denominator);
   kinetic_energy_squared = Cdiv(kinetic_energy_squared, weight_denominator);
   total_energy = Cadd(potential_energy, kinetic_energy);

   /*Also Add Initial Energy Shift from MolPro*/
   total_energy.real += energy_shifts[0].real;

   /*Print Desired Information*******************************************/
   /*Find Errors*/
   potential_energy_error_complex = Csub(potential_energy_squared, Cmul(potential_energy, conjugate(potential_energy)));
   if ( potential_energy_error_complex.real >  0 ) {
     potential_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * potential_energy_error_complex.real);
   }
   else {
     potential_energy_error_complex.real = 0.0;
   }

   kinetic_energy_error_complex = Csub(kinetic_energy_squared, Cmul(kinetic_energy, conjugate(kinetic_energy)));
   if ( kinetic_energy_error_complex.real > 0 ) {
     kinetic_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * kinetic_energy_error_complex.real);
   }
   else {
    kinetic_energy_error_complex.real = 0.0;
   }
   total_energy_error_complex.real = sqrt(potential_energy_error_complex.real*potential_energy_error_complex.real +   kinetic_energy_error_complex.real*kinetic_energy_error_complex.real);

   /*Print Desired Information*******************************************/
   *ke = kinetic_energy;
   *pe = potential_energy;
   *te = total_energy;
   *ke_error = kinetic_energy_error_complex;
   *pe_error = potential_energy_error_complex;
   *te_error = total_energy_error_complex;
   *tw = weight_denominator;

  /*Obtain and Print Average WFs*/
  for (i=0; i<ist.n_spatial_orbitals; i++) {
    for (j=0; j<ist.n_up; j++) {
      average_wave_function_up[i*ist.n_up+j] = Cdiv(average_wave_function_up[i*ist.n_up+j], weight_denominator);
   }
 }
  for (i=0; i<ist.n_spatial_orbitals; i++) {
    for (j=0; j<ist.n_down; j++) {
      average_wave_function_down[i*ist.n_down+j] = Cdiv(average_wave_function_down[i*ist.n_down+j], weight_denominator);
    }
  }

 /*Find Normalized WF*/
 get_normalized_wf_up_chemistry(normalized_wave_function_up,average_wave_function_up,ist);
 get_normalized_wf_down_chemistry(normalized_wave_function_down,average_wave_function_down,ist);

free(average_wave_function_up);
free(average_wave_function_down);
return;
}

/************************************************************************************************************************************************/

void compute_energy_chemistry_mixed_multi_restricted_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,double *kinetic_matrix_sparse,double *potential_matrix_sparse,int *kinetic_ij_sparse,int *potential_ij_sparse,MKL_Complex16 *energy_shifts,MKL_Complex16 *total_potential_energy,MKL_Complex16 *total_kinetic_energy,MKL_Complex16 *total_energy_2,MKL_Complex16 *total_weight,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *pe_error,MKL_Complex16 *ke,MKL_Complex16 *ke_error,MKL_Complex16 *te,MKL_Complex16 *te_error,MKL_Complex16 *tw,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   /*Computes the Energy of the System Using Mixed Estimator for a Restricted Situation*/
   int walker, i, j, k;
   int number_of_accumulated_steps;
   int number_walkers_active = 0;
   double weight_factor = 1.0;
   double energy_error, kinetic_energy_error_real, potential_energy_error_real;
   MKL_Complex16 kinetic_energy_error_complex, potential_energy_error_complex, total_energy_error_complex;
   MKL_Complex16 *average_wave_function_up, *average_wave_function_down;
   MKL_Complex16 local_pe, local_ke, local_pe_squared, local_ke_squared;
   MKL_Complex16 Zero;
   MKL_Complex16 potential_energy, kinetic_energy, total_energy;
   MKL_Complex16 potential_energy_squared, kinetic_energy_squared;
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator;
   MKL_Complex16 energy_average;

   /*HERE*/
   Zero.real = Zero.imag = 0.0;
   potential_energy.real = potential_energy.imag = 0.0;
   kinetic_energy.real = kinetic_energy.imag = 0.0;
   total_energy.real = total_energy.imag = 0.0;
   potential_energy_squared.real = potential_energy_squared.imag = 0.0;
   kinetic_energy_squared.real = kinetic_energy_squared.imag;
   weight_denominator.real = weight_denominator.imag = 0.0;

   average_wave_function_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_up*sizeof(MKL_Complex16));
   average_wave_function_down=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_down*sizeof(MKL_Complex16));

   czero_vec(average_wave_function_up,ist.n_spatial_orbitals*ist.n_up);
   czero_vec(average_wave_function_down,ist.n_spatial_orbitals*ist.n_down);

   /*Obtain Population Control Rescaling Factor*/
   for (i = 0; i < 10 ; i++) {
    weight_factor *= weight_rescaling[i];
   }

   for (walker = 0; walker < ist.n_walkers ; walker++ ) {

     /*Check That Walker Weights Are Non-Zero*/
     if ( Cabs(weights[walker]) != 0 ) {
        number_walkers_active++;

        /*Accumulate Weight in Denominator*/
        current_weight = RCmul(weight_factor, weights[walker]);
        weight_denominator = Cadd(weight_denominator, current_weight);

        /*Compute the Energies Using the Mixed Estimator, but for a Multideterminant WF*/ 
        /*Using the Full Trial Energy Wavefunction To Get Exact Energies To Some Degree*/
        if ( ist.n_up > 0 ) {
         if ( ist.n_down > 0 ) {
              compute_walker_energy_chemistry_mixed_multi_energy_restricted_sparse(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,kinetic_matrix_sparse,potential_matrix_sparse,kinetic_ij_sparse,potential_ij_sparse,&local_ke,&local_pe,ist);
         }
         else {
              compute_walker_energy_chemistry_mixed_multi_energy_up_sparse(&wf_up[walker*ist.n_spatial_orbitals_n_up],trial_wf_up_energy,trial_determinant_coefficients_energy,kinetic_matrix_sparse,kinetic_ij_sparse,&local_ke,&local_pe,ist);
         }
        }
        else { 
          if ( ist.n_down > 0 ) {
            compute_walker_energy_chemistry_mixed_multi_energy_down_sparse(&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_down_energy,trial_determinant_coefficients_energy,kinetic_matrix_sparse,kinetic_ij_sparse,&local_ke,&local_pe,ist); 
          }
        }

        /*Obtain Kinetic Energy*/
        local_ke_squared = Cmul(local_ke, conjugate(local_ke));
        local_ke = Cmul(local_ke, current_weight);
        local_ke_squared = Cmul(local_ke_squared, current_weight);

        /*Get Potential Energy*/
        local_pe_squared = Cmul(local_pe, conjugate(local_pe));
        local_pe = Cmul(local_pe, current_weight);
        local_pe_squared = Cmul(local_pe_squared, current_weight);

        /*Add to Total*/
        potential_energy = Cadd(local_pe, potential_energy);
        potential_energy_squared = Cadd(local_pe_squared, potential_energy_squared);
        kinetic_energy = Cadd(local_ke, kinetic_energy);
        kinetic_energy_squared = Cadd(kinetic_energy_squared, local_ke_squared);

        /*Obtain Average Wavefunction As Well*/
        for (i=0; i<ist.n_spatial_orbitals; i++) {
         for (j=0; j<ist.n_up; j++) {
           average_wave_function_up[i*ist.n_up+j] = Cadd(average_wave_function_up[i*ist.n_up+j], Cmul(current_weight, wf_up[walker*ist.n_spatial_orbitals_n_up+i*ist.n_up+j]));
         }
         for (j=0; j<ist.n_down; j++) {
           average_wave_function_down[i*ist.n_down+j] = Cadd(average_wave_function_down[i*ist.n_down+j], Cmul(current_weight, wf_down[walker*ist.n_spatial_orbitals_n_down+i*ist.n_down+j]));
         }
        }

      } /*Non-Zero Walker Weights*/
     }  /*Walkers*/


    /*Accumulate the Energy of Each Block To Get Overall Statistics*/
    if ( step >= ist.n_steps_equilibration ) {
      accumulate_energy_blocks(potential_energy, kinetic_energy, Cadd(potential_energy, kinetic_energy), weight_denominator, total_potential_energy, total_kinetic_energy, total_energy_2, total_weight);
    }

   /*Divide By Total Weight*/
   potential_energy = Cdiv(potential_energy, weight_denominator);
   potential_energy_squared = Cdiv(potential_energy_squared, weight_denominator);
   kinetic_energy = Cdiv(kinetic_energy, weight_denominator);
   kinetic_energy_squared = Cdiv(kinetic_energy_squared, weight_denominator);
   total_energy = Cadd(potential_energy, kinetic_energy);

   /*Also Add Initial Energy Shift from MolPro*/
   total_energy.real += energy_shifts[0].real;

   /*Print Desired Information*******************************************/
   /*Find Errors*/
   potential_energy_error_complex = Csub(potential_energy_squared, Cmul(potential_energy, conjugate(potential_energy)));
   if ( potential_energy_error_complex.real >  0 ) {
     potential_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * potential_energy_error_complex.real);
   }
   else {
     potential_energy_error_complex.real = 0.0;
   }

   kinetic_energy_error_complex = Csub(kinetic_energy_squared, Cmul(kinetic_energy, conjugate(kinetic_energy)));
   if ( kinetic_energy_error_complex.real > 0 ) {
     kinetic_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * kinetic_energy_error_complex.real);
   }
   else {
    kinetic_energy_error_complex.real = 0.0;
   }
   total_energy_error_complex.real = sqrt(potential_energy_error_complex.real*potential_energy_error_complex.real +   kinetic_energy_error_complex.real*kinetic_energy_error_complex.real);

   /*Print Desired Information*******************************************/
   *ke = kinetic_energy;
   *pe = potential_energy;
   *te = total_energy;
   *ke_error = kinetic_energy_error_complex;
   *pe_error = potential_energy_error_complex;
   *te_error = total_energy_error_complex;
   *tw = weight_denominator;

  /*Obtain and Print Average WFs*/
  for (i=0; i<ist.n_spatial_orbitals; i++) {
    for (j=0; j<ist.n_up; j++) {
      average_wave_function_up[i*ist.n_up+j] = Cdiv(average_wave_function_up[i*ist.n_up+j], weight_denominator);
   }
 }
  for (i=0; i<ist.n_spatial_orbitals; i++) {
    for (j=0; j<ist.n_down; j++) {
      average_wave_function_down[i*ist.n_down+j] = Cdiv(average_wave_function_down[i*ist.n_down+j], weight_denominator);
    }
  }

 /*Find Normalized WF*/
 get_normalized_wf_up_chemistry(normalized_wave_function_up,average_wave_function_up,ist);
 get_normalized_wf_down_chemistry(normalized_wave_function_down,average_wave_function_down,ist);

  /*Find Statistics So Far*******************************************/
  /*if ( step >= ist.n_steps_equilibration ) {
      number_of_accumulated_steps = 1+(int)((step-ist.n_steps_equilibration)/(ist.n_steps_energy)); 
      fprintf(pf7, "%d\t", step); 
 
    //Get PE
    energy_average = Cdiv(total_potential_energy[0], total_weight[0]);  
    potential_energy_error_real = error_real(total_potential_energy[1].real, total_potential_energy[2].real, number_of_accumulated_steps); 
    fprintf(pf7, "%f   %f\t", energy_average.real, potential_energy_error_real);   

    //Get KE
    energy_average = Cdiv(total_kinetic_energy[0], total_weight[0]); 
    kinetic_energy_error_real = error_real(total_kinetic_energy[1].real, total_kinetic_energy[2].real, number_of_accumulated_steps);
    fprintf(pf7, "%f   %f\t", energy_average.real, kinetic_energy_error_real);

    //Show Shift
    fprintf(pf7, "%f   %f\t", energy_shifts[0].real, 0.0); 

    //Get Total
    energy_average = Cdiv(total_energy_2[0], total_weight[0]);
    energy_error = sqrt(potential_energy_error_real*potential_energy_error_real+kinetic_energy_error_real*kinetic_energy_error_real);  
    //energy_error = error_real(total_energy_2[1].r, total_energy_2[2].r, number_of_accumulated_steps);
    fprintf(pf7, "%f   %f\n", energy_average.real+energy_shifts[0].real, energy_error);
    fflush(pf7); 
   }
   */

free(average_wave_function_up);
free(average_wave_function_down);
return;
}

/************************************************************************************************/

void compute_energy_chemistry_mixed_unrestricted_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *energy_shifts,MKL_Complex16 *total_potential_energy,MKL_Complex16 *total_kinetic_energy,MKL_Complex16 *total_energy_2,MKL_Complex16 *total_weight,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *pe_error,MKL_Complex16 *ke,MKL_Complex16 *ke_error,MKL_Complex16 *te,MKL_Complex16 *te_error,MKL_Complex16 *tw,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   /*Computes the Energy of the System For an Unrestricted Situation Using the Mixed Estimator*/
   int walker, i, j, k;
   int number_of_accumulated_steps, number_walkers_active;
   double weight_factor = 1.0;
   double energy_error;
   double kinetic_energy_error_real, potential_energy_error_real;
   MKL_Complex16 kinetic_energy_error_complex, potential_energy_error_complex, total_energy_error_complex;
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
   MKL_Complex16 *average_wave_function_up, *average_wave_function_down;
   MKL_Complex16 local_pe, local_ke, local_pe_squared, local_ke_squared;
   MKL_Complex16 potential_energy = Zero, kinetic_energy = Zero, total_energy = Zero;
   MKL_Complex16 potential_energy_squared = Zero, kinetic_energy_squared = Zero;
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator = Zero;
   MKL_Complex16 energy_average;

   average_wave_function_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   average_wave_function_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   czero_vec(average_wave_function_up,ist.n_spatial_orbitals*ist.n_up);
   czero_vec(average_wave_function_down,ist.n_spatial_orbitals*ist.n_down);

   /*Obtain Population Control Rescaling Factor*/
   for (i = 0; i < 10 ; i++) {
    weight_factor *= weight_rescaling[i];
   }

   number_walkers_active = 0;
   for (walker = 0; walker < ist.n_walkers ; walker++ ) {

     /*Check That Walker Weights Are Non-Zero*/
     if ( Cabs(weights[walker]) != 0 ) {
        number_walkers_active++;

        /*Accumulate Weight in Denominator*/
        current_weight = RCmul(weight_factor, weights[walker]);
        weight_denominator = Cadd(weight_denominator, current_weight);

        /*First Compute Kinetic Energy*/
        compute_walker_energy_chemistry_mixed_unrestricted_sparse(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,&local_ke,&local_pe,ist); 

        /*Obtain Kinetic Energy*/
        local_ke_squared = Cmul(local_ke, conjugate(local_ke));
        local_ke = Cmul(local_ke, current_weight);
        local_ke_squared = Cmul(local_ke_squared, current_weight);

        /*Get Potential Energy*/
        local_pe_squared = Cmul(local_pe, conjugate(local_pe));
        local_pe = Cmul(local_pe, current_weight);
        local_pe_squared = Cmul(local_pe_squared, current_weight);

        /*Add to Total*/
        potential_energy = Cadd(local_pe, potential_energy);
        potential_energy_squared = Cadd(local_pe_squared, potential_energy_squared);
        kinetic_energy = Cadd(local_ke, kinetic_energy);
        kinetic_energy_squared = Cadd(local_ke_squared, kinetic_energy_squared);

        /*Obtain Average Wavefunction As Well*/
        for (i=0; i<ist.n_spatial_orbitals; i++) {
         for (j=0; j<ist.n_up; j++) {
           average_wave_function_up[i*ist.n_up+j] = Cadd(average_wave_function_up[i*ist.n_up+j], Cmul(current_weight, wf_up[walker*ist.n_spatial_orbitals_n_up+i*ist.n_up+j]));
         }
         for (j=0; j<ist.n_down; j++) {
           average_wave_function_down[i*ist.n_down+j] = Cadd(average_wave_function_down[i*ist.n_down+j], Cmul(current_weight, wf_down[walker*ist.n_spatial_orbitals_n_down+i*ist.n_down+j]));
         }
        }

      } /*Non-Zero Walker Weights*/
   }  /*Walkers*/

   /*Accumulate Energies and Weights for Statistics*/
   if ( step > ist.n_steps_equilibration ) {
     accumulate_energy_blocks(potential_energy, kinetic_energy, Cadd(potential_energy, kinetic_energy), weight_denominator, total_potential_energy, total_kinetic_energy, total_energy_2, total_weight);
   }

   /*Divide By Total Weight*/
   potential_energy = Cdiv(potential_energy, weight_denominator);

   potential_energy_squared = Cdiv(potential_energy_squared, weight_denominator);
   kinetic_energy = Cdiv(kinetic_energy, weight_denominator);
   kinetic_energy_squared = Cdiv(kinetic_energy_squared, weight_denominator);
   total_energy = Cadd(potential_energy, kinetic_energy);

   /*Also Add Initial Energy Shift from MolPro*/
   total_energy.real += energy_shifts[0].real+energy_shifts[1].real+energy_shifts[2].real+energy_shifts[3].real;

   /*Find Errors*/
   potential_energy_error_complex = Csub(potential_energy_squared, Cmul(potential_energy, conjugate(potential_energy)));
   if ( potential_energy_error_complex.real >  0 ) {
     potential_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * potential_energy_error_complex.real);
   }
   else {
     potential_energy_error_complex.real = 0.0;
   }

   kinetic_energy_error_complex = Csub(kinetic_energy_squared, Cmul(kinetic_energy, conjugate(kinetic_energy)));
   if ( kinetic_energy_error_complex.real > 0 ) {
     kinetic_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * kinetic_energy_error_complex.real);
   }
   else {
    kinetic_energy_error_complex.real = 0.0;
   }
   total_energy_error_complex.real = sqrt(potential_energy_error_complex.real*potential_energy_error_complex.real + kinetic_energy_error_complex.real*kinetic_energy_error_complex.real);


   /*Print Desired Information*******************************************/
   *ke = kinetic_energy;
   *pe = potential_energy;
   *te = total_energy;
   *ke_error = kinetic_energy_error_complex;
   *pe_error = potential_energy_error_complex;
   *te_error = total_energy_error_complex;
   *tw = weight_denominator;

  /*Obtain and Print Average WFs*/
  for (i=0; i<ist.n_spatial_orbitals; i++) {
    for (j=0; j<ist.n_up; j++) {
      average_wave_function_up[i*ist.n_up+j] = Cdiv(average_wave_function_up[i*ist.n_up+j], weight_denominator);
   }
 }
  for (i=0; i<ist.n_spatial_orbitals; i++) {
    for (j=0; j<ist.n_down; j++) {
      average_wave_function_down[i*ist.n_down+j] = Cdiv(average_wave_function_down[i*ist.n_down+j], weight_denominator);
    }
  }

 /*Find Normalized WF*/
 get_normalized_wf_up_chemistry(normalized_wave_function_up,average_wave_function_up,ist);
 get_normalized_wf_down_chemistry(normalized_wave_function_down,average_wave_function_down,ist);

  /*Restart Files*****************************************************************/
  /* for (i=0; i<ist.n_walkers; i++) {
     fprintf(pf8, "%f\t %f\t", weight_factor * weights[i].real, weight_factor * weights[i].imag);
      for (j=0; j<ist.n_spatial_orbitals; j++) {
        for (k=0; k<ist.n_up; k++) {
          fprintf(pf8, "%f\t %f\t", average_wave_function_up[j*ist.n_up+k].real, average_wave_function_up[j*ist.n_up+k].imag);
        }
      }
 
     for (j=0; j<ist.n_spatial_orbitals; j++) {
       for (k=0; k<ist.n_down; k++) {
         fprintf(pf8, "%f\t %f\t", average_wave_function_down[j*ist.n_down+k].real, average_wave_function_down[j*ist.n_down+k].imag);
       }
     }
  }
 */

  /*Find Statistics So Far*******************************************/
  /*if ( step >= ist.n_steps_equilibration ) {
       number_of_accumulated_steps = 1+(int)((step-ist.n_steps_equilibration)/(ist.n_steps_energy));
       fprintf(pf7, "%d\t", step);
 
    //Get PE
    energy_average = Cdiv(total_potential_energy[0], total_weight[0]);
    potential_energy_error_real = error_real(total_potential_energy[1].real, total_potential_energy[2].real, number_of_accumulated_steps);
    fprintf(pf7, "%f   %f\t", energy_average.real, potential_energy_error_real);

    //Get KE
    energy_average = Cdiv(total_kinetic_energy[0], total_weight[0]);
    kinetic_energy_error_real = error_real(total_kinetic_energy[1].real, total_kinetic_energy[2].real, number_of_accumulated_steps);
    fprintf(pf7, "%f   %f\t", energy_average.real, kinetic_energy_error_real);

    //Print Shifts
    fprintf(pf7, "%f   %f\t", energy_shifts[0].real+energy_shifts[1].real+energy_shifts[2].real+energy_shifts[3].real, 0.0);

    //Get Total
    energy_average = Cdiv(total_energy_2[0], total_weight[0]);
    energy_error = sqrt(potential_energy_error_real*potential_energy_error_real+kinetic_energy_error_real*kinetic_energy_error_real);
    energy_error = error_real(total_energy_2[1].real, total_energy_2[2].real, number_of_accumulated_steps);
    fprintf(pf7, "%f   %f\n", energy_average.real+energy_shifts[0].real+energy_shifts[1].real+energy_shifts[2].real+energy_shifts[3].real,energy_error);
    fflush(pf7); 

   } 
 */

free(average_wave_function_up);
free(average_wave_function_down);
return;
}

/***********************************************************************************************/

MKL_Complex16 compute_kinetic_energy_density_matrix_restricted_sparse(MKL_Complex16 *density_matrix_up,MKL_Complex16 *density_matrix_down,double *kinetic_matrix_sparse,int *kinetic_ij_sparse,int_st ist) {

    /*Computes the Kinetic Energy*/
    int i;
    MKL_Complex16 kinetic_energy; kinetic_energy.real = kinetic_energy.imag = 0;

    /*Multiply Sparse Up and Down Matrices*/
    for (i=0; i< ist.n_kinetic_sparse_up; i++) {
      kinetic_energy = Cadd(kinetic_energy, RCmul(kinetic_matrix_sparse[i], Cadd(density_matrix_up[kinetic_ij_sparse[i]], density_matrix_down[kinetic_ij_sparse[i]]))); 
    }  

return(kinetic_energy);
}

/*******************************************************************************************/

MKL_Complex16 compute_kinetic_energy_density_matrix_restricted_up_sparse(MKL_Complex16 *density_matrix_up,double *kinetic_matrix_sparse,int *kinetic_ij_sparse,int_st ist) {

    /*Computes the Kinetic Energy*/
    int i;
    MKL_Complex16 kinetic_energy; kinetic_energy.real = kinetic_energy.imag = 0;

    /*Multiply Sparse Up and Down Matrices*/
    for (i=0; i< ist.n_kinetic_sparse_up; i++) {
      kinetic_energy = Cadd(kinetic_energy, RCmul(kinetic_matrix_sparse[i], density_matrix_up[kinetic_ij_sparse[i]])); 
    }

return(kinetic_energy);
}

/*******************************************************************************************/

MKL_Complex16 compute_kinetic_energy_density_matrix_restricted_down_sparse(MKL_Complex16 *density_matrix_down,double *kinetic_matrix_sparse,int *kinetic_ij_sparse,int_st ist) {

    /*Computes the Kinetic Energy*/
    int i;
    MKL_Complex16 kinetic_energy; kinetic_energy.real = kinetic_energy.imag = 0;

    /*Multiply Sparse Up and Down Matrices*/
    for (i=0; i< ist.n_kinetic_sparse_down; i++) {
      kinetic_energy = Cadd(kinetic_energy, RCmul(kinetic_matrix_sparse[i], density_matrix_down[kinetic_ij_sparse[i]]));       
    }

return(kinetic_energy);
}

/*******************************************************************************************/

MKL_Complex16 compute_kinetic_energy_density_matrix_unrestricted_sparse(MKL_Complex16 *density_matrix_up,MKL_Complex16 *density_matrix_down,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int_st ist) {

    /*Computes the Kinetic Energy*/
    int i; 
    MKL_Complex16 kinetic_energy; kinetic_energy.real = kinetic_energy.imag = 0;

    /*Multiply Sparse Up and Down Matrices*/
    if ( ist.n_kinetic_sparse_up == ist.n_kinetic_sparse_down ) { 

       for (i=0; i< ist.n_kinetic_sparse_up; i++) {
        kinetic_energy = Cadd(kinetic_energy, RCmul(kinetic_matrix_sparse_up[i], density_matrix_up[kinetic_ij_sparse_up[i]])); 
        kinetic_energy = Cadd(kinetic_energy, RCmul(kinetic_matrix_sparse_down[i], density_matrix_down[kinetic_ij_sparse_down[i]])); 
       }
  
    }
    else if ( ist.n_kinetic_sparse_up > ist.n_kinetic_sparse_down ) {
    
       for (i=0; i<ist.n_kinetic_sparse_up; i++) {
         kinetic_energy = Cadd(kinetic_energy, RCmul(kinetic_matrix_sparse_up[i], density_matrix_up[kinetic_ij_sparse_up[i]]));

         if ( i < ist.n_kinetic_sparse_down ) {
           kinetic_energy = Cadd(kinetic_energy, RCmul(kinetic_matrix_sparse_down[i], density_matrix_down[kinetic_ij_sparse_down[i]]));
         }
       } 

    }
    else {
 
       for (i=0; i<ist.n_kinetic_sparse_down; i++) {
         kinetic_energy = Cadd(kinetic_energy, RCmul(kinetic_matrix_sparse_down[i], density_matrix_down[kinetic_ij_sparse_down[i]]));
         
         if ( i < ist.n_kinetic_sparse_up ) {
           kinetic_energy = Cadd(kinetic_energy, RCmul(kinetic_matrix_sparse_up[i], density_matrix_up[kinetic_ij_sparse_up[i]]));
         }
       }

    }
   
return(kinetic_energy);
}

/********************************************************************************************************/

MKL_Complex16 compute_potential_energy_density_matrix_restricted_sparse(MKL_Complex16 *density_matrix_up,MKL_Complex16 *density_matrix_down,double *potential_matrix_sparse, int *potential_ij_sparse, int_st ist) {

     /*Compute PE Using Density Matrix Using Sparsity*/
     int i; 
     int ik, il, jk, jl; 
     MKL_Complex16 pot, prod; 
 
     /*Zero Potential*/
     pot.real = pot.imag = 0.0; 

     /*Run Through All Non-Zero Elements To Get Potential*/
     for (i=0; i<ist.n_potential_sparse_up; i++) {

       ik = potential_ij_sparse[4*i]; 
       jl = potential_ij_sparse[4*i+1]; 
       il = potential_ij_sparse[4*i+2]; 
       jk = potential_ij_sparse[4*i+3]; 

       /*Coulomb Terms*/
       prod = Cmul(density_matrix_up[ik], density_matrix_up[jl]); 
       prod = Cadd(prod, Cmul(density_matrix_down[ik], density_matrix_down[jl])); 

       prod = Cadd(prod, Cmul(density_matrix_up[ik], density_matrix_down[jl])); 
       prod = Cadd(prod, Cmul(density_matrix_down[ik], density_matrix_up[jl])); 

       /*Exchange Terms*/
       prod = Csub(prod, Cmul(density_matrix_up[il], density_matrix_up[jk])); 
       prod = Csub(prod, Cmul(density_matrix_down[il], density_matrix_down[jk])); 

       /*Multiply Into Potential*/
       prod = RCmul(potential_matrix_sparse[i], prod); 
       pot = Cadd(pot, prod); 
    }   

    /*Half Sum*/
    pot = RCmul(.5, pot); 

return(pot); 
}

/************************************************************************************************************/

MKL_Complex16 compute_potential_energy_density_matrix_unrestricted_sparse(MKL_Complex16 *density_matrix_up,MKL_Complex16 *density_matrix_down, double *potential_matrix_sparse_up, double *potential_matrix_sparse_down, double *potential_matrix_sparse_updown, int *potential_ij_sparse_up, int *potential_ij_sparse_down, int *potential_ij_sparse_updown, int_st ist) {

   /*Computes the PE Using Density Matrix and Sparse Routines*/
   int i;
   int ik, il, jl, jk;  
   MKL_Complex16 pot; 
   MKL_Complex16 prod; 

   /*Zero Pot*/
   pot.real = pot.imag = 0.0;  

   /*May Be a Faster Way of Doing This If Entries Are Zero and Non-Zero at the Same Time*/
   /*Run Through Non-Zero Entries of Up-Up Matrix*/
   for (i=0; i<ist.n_potential_sparse_up; i++) {    

     ik = potential_ij_sparse_up[4*i]; 
     jl = potential_ij_sparse_up[4*i+1]; 
     il = potential_ij_sparse_up[4*i+2]; 
     jk = potential_ij_sparse_up[4*i+3]; 

     /*Get Exchange and Coulomb Terms*/
     prod = Cmul(density_matrix_up[ik], density_matrix_up[jl]); 
     prod = Csub(prod, Cmul(density_matrix_up[il], density_matrix_up[jk]));    

     /*Multiply By Constant*/
     prod = RCmul(potential_matrix_sparse_up[i], prod); 
     pot = Cadd(pot, prod); 
  }

  /*Do Same for Down-Down Matrix*/
  for (i=0; i<ist.n_potential_sparse_down; i++) {

     ik = potential_ij_sparse_down[4*i]; 
     jl = potential_ij_sparse_down[4*i+1];
     il = potential_ij_sparse_down[4*i+2];
     jk = potential_ij_sparse_down[4*i+3];

     /*Get Exchange and Coulomb Terms*/
     prod = Cmul(density_matrix_down[ik], density_matrix_down[jl]);
     prod = Csub(prod, Cmul(density_matrix_down[il], density_matrix_down[jk]));   

     /*Multiply By Constant*/
     prod = RCmul(potential_matrix_sparse_down[i], prod);
     pot = Cadd(pot, prod);
  
  }

  /*Do Same Thing for Up-Down Matrix*/
  for (i=0; i<ist.n_potential_sparse_updown; i++) {

     ik = potential_ij_sparse_updown[2*i];
     jl = potential_ij_sparse_updown[2*i+1];

     /*Get Exchange and Coulomb Terms*/
     prod = Cmul(density_matrix_up[ik], density_matrix_down[jl]);
     prod = RCmul(potential_matrix_sparse_updown[2*i], prod); 
     pot = Cadd(pot, prod); 

     /*Multiply By Constant*/     
     prod = Cmul(density_matrix_up[jl], density_matrix_down[ik]);  
     prod = RCmul(potential_matrix_sparse_updown[2*i+1], prod);
     pot = Cadd(pot, prod);

  }

  /*Half Sum*/
  pot = RCmul(.5, pot); 

return(pot); 
}

/**************************************************************************************************************************************************************/

void compute_walker_energy_chemistry_mixed_multi_phaseless_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,double *kinetic_matrix_sparse,double *potential_matrix_sparse,int *kinetic_ij_sparse, int *potential_ij_sparse, MKL_Complex16 *potential_original_matrix,MKL_Complex16 *total_ke,MKL_Complex16 *total_pe,int_st ist) {

    /*Computes the Kinetic and Potential Energies of a Wavefunction If It Is Crossed with a Multideterminant Trial Wave Function*/
    int i, j, k;
    MKL_Complex16 *One, *Zero;
    MKL_Complex16 *stored_product1_up, *stored_product1_down;
    MKL_Complex16 *density_matrix_up, *density_matrix_down;
    MKL_Complex16 product_dets; 
    MKL_Complex16 local_ke;
    MKL_Complex16 local_pe;
    MKL_Complex16 ke, pe;
    /*HERE*/
    
    /*Computes the Spin Up and Down Greens Functions*/
    stored_product1_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    stored_product1_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

    density_matrix_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
    density_matrix_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    /*Run Through Determinants*/
    ke.real = ke.imag = 0.0;
    pe.real = pe.imag = 0.0;

    /*Determine Full Overlap*/
    for (i=0; i<ist.n_determinants_trial_phaseless; i++) {

     /*Check Whether Determinant Is Non Zero So That Singularities Don't Occur*/
     if ( Cabs(det_overlap_up[i]) > 0 || Cabs(det_overlap_down[i]) > 0 ) {


       /*Zero All Total Vectors*/
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_up,ist.n_spatial_orbitals,ist.n_up,One,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up,&trial_wf_up_phaseless[i*ist.n_spatial_orbitals_n_up],ist.n_up,Zero,stored_product1_up,ist.n_spatial_orbitals);
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_down,ist.n_spatial_orbitals,ist.n_down,One,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down,&trial_wf_down_phaseless[i*ist.n_spatial_orbitals_n_down],ist.n_down,Zero,stored_product1_down,ist.n_spatial_orbitals);

       /*Now Get Total Density Matrix*/
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up,One,wf_up,ist.n_up,stored_product1_up,ist.n_spatial_orbitals,Zero,density_matrix_up,ist.n_spatial_orbitals);
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down,One,wf_down,ist.n_down,stored_product1_down,ist.n_spatial_orbitals,Zero,density_matrix_down,ist.n_spatial_orbitals);

       /*Get Ke and PE*/
       local_ke = compute_kinetic_energy_density_matrix_restricted_sparse(density_matrix_up,density_matrix_down,kinetic_matrix_sparse,kinetic_ij_sparse,ist); 
       local_pe = compute_potential_energy_density_matrix_restricted_sparse(density_matrix_up,density_matrix_down,potential_matrix_sparse,potential_ij_sparse,ist); 

       /*Get Product of dets*/
       product_dets = Cmul(conjugate(trial_determinant_coefficients_phaseless[i]), Cmul(det_overlap_up[i], det_overlap_down[i]));

       /*Scale By Overlap and Coefficient*/
       local_ke = Cmul(local_ke, product_dets);
       local_pe = Cmul(local_pe, product_dets);

       /*Add Up Local Pe and Ke*/
       ke = Cadd(ke, local_ke);
       pe = Cadd(pe, local_pe);

     } /*Check Determinant*/

    }
    *total_ke = Cdiv(ke, overlap_total[0]); 
    *total_pe = Cdiv(pe, overlap_total[0]); 

free(stored_product1_up);
free(stored_product1_down);
free(density_matrix_up);
free(density_matrix_down);
return;
}

/***************************************************************************************************************************************************************************/

void compute_walker_energy_chemistry_mixed_multi_energy_unrestricted_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *energies2,MKL_Complex16 *energies3,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *density_matrix_up_2,MKL_Complex16 *density_matrix_down_2, double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *total_ke,MKL_Complex16 *total_pe,int_st ist) {

    /*Computes the Kinetic and Potential Energies of a Wavefunction If It Is Crossed with a Multideterminant Trial Wave Function*/
    int i, j, k;
    MKL_Complex16 *One, *Zero;
    MKL_Complex16 *stored_product1_up, *stored_product1_down;
    MKL_Complex16 *density_matrix_up, *density_matrix_down;
    MKL_Complex16 *det_overlap_up, *det_overlap_down;
    MKL_Complex16 *overlap_inverse_up, *overlap_inverse_down;
    MKL_Complex16 *overlap_total;
    MKL_Complex16 product_dets;
    MKL_Complex16 local_ke;
    MKL_Complex16 local_pe;
    MKL_Complex16 pe, ke;
   /*HERE*/

    /*Computes the Spin Up and Down Greens Functions*/
    stored_product1_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    stored_product1_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

    density_matrix_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
    density_matrix_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

    det_overlap_up=(MKL_Complex16 *)calloc(ist.n_determinants_trial_energy,sizeof(MKL_Complex16));
    det_overlap_down=(MKL_Complex16 *)calloc(ist.n_determinants_trial_energy,sizeof(MKL_Complex16));

    overlap_inverse_up=(MKL_Complex16 *)calloc(ist.n_determinants_energy_n_up_sq,sizeof(MKL_Complex16));
    overlap_inverse_down=(MKL_Complex16 *)calloc(ist.n_determinants_energy_n_down_sq,sizeof(MKL_Complex16));

    overlap_total=(MKL_Complex16 *)calloc(1,sizeof(MKL_Complex16));

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    /*Determine Overlaps With Full Trial Energy WF Instead of Just Phaseless One*/
    update_overlaps_chemistry_multi_energy_both(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_inverse_up,overlap_inverse_down,&overlap_total[0],det_overlap_up,det_overlap_down,ist);

    /*Run Through Determinants*/
    ke.real = ke.imag = 0.0;
    pe.real = pe.imag = 0.0;

    /*Determine Full Overlap*/
    for (i=0; i<ist.n_determinants_trial_energy; i++) {

     /*Check Whether Determinant Is Non Zero So That Singularities Don't Occur*/
     if ( Cabs(det_overlap_up[i]) > 0 || Cabs(det_overlap_down[i]) > 0 ) {

       /*Zero All Total Vectors*/
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_up,ist.n_spatial_orbitals,ist.n_up,One,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up,&trial_wf_up_energy[i*ist.n_spatial_orbitals_n_up],ist.n_up,Zero,stored_product1_up,ist.n_spatial_orbitals);
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_down,ist.n_spatial_orbitals,ist.n_down,One,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down,&trial_wf_down_energy[i*ist.n_spatial_orbitals_n_down],ist.n_down,Zero,stored_product1_down,ist.n_spatial_orbitals);

       /*Now Get Total Density Matrix*/
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up,One,wf_up,ist.n_up,stored_product1_up,ist.n_spatial_orbitals,Zero,density_matrix_up,ist.n_spatial_orbitals);
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down,One,wf_down,ist.n_down,stored_product1_down,ist.n_spatial_orbitals,Zero,density_matrix_down,ist.n_spatial_orbitals);

       for (j=0; j<ist.n_spatial_orbitals_sq; j++) {
         density_matrix_up_2[j] = density_matrix_up[j]; 
         density_matrix_down_2[j] = density_matrix_down[j]; 
       } 

       /*Get Ke and PE*/
       local_ke = compute_kinetic_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,kinetic_ij_sparse_up,kinetic_ij_sparse_down,ist); 
       local_pe = compute_potential_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,ist); 
       energies2[0] = local_ke; 
       energies2[1] = local_pe; 
       energies2[2] = Cadd(local_ke, local_pe); 


       /*Get Product of dets*/
       product_dets = Cmul(conjugate(trial_determinant_coefficients_energy[i]), Cmul(det_overlap_up[i], det_overlap_down[i]));

       /*Scale By Overlap and Coefficient*/
       local_ke = Cmul(local_ke, product_dets);
       local_pe = Cmul(local_pe, product_dets);

       /*Add Up Local Pe and Ke*/
       ke = Cadd(ke, local_ke);
       pe = Cadd(pe, local_pe);

      }

    } /*Check Determinant*/
    *total_ke = Cdiv(ke, overlap_total[0]);
    *total_pe = Cdiv(pe, overlap_total[0]);

    energies3[0] = *total_ke; 
    energies3[1] = *total_pe; 
    energies3[2] = Cdiv(Cadd(ke,pe), overlap_total[0]);      


free(stored_product1_up);
free(stored_product1_down);
free(density_matrix_up);
free(density_matrix_down);
free(det_overlap_up);
free(det_overlap_down);
free(overlap_inverse_up);
free(overlap_inverse_down);
return; 
}

/***************************************************************************************************************************************/

void compute_walker_energy_chemistry_mixed_multi_energy_restricted_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,double *kinetic_matrix_sparse,double *potential_matrix_sparse,int *kinetic_ij_sparse,int *potential_ij_sparse,MKL_Complex16 *total_ke,MKL_Complex16 *total_pe,int_st ist) {

    /*Computes the Kinetic and Potential Energies of a Wavefunction If It Is Crossed with a Multideterminant Trial Wave Function*/
    int i, j, k;
    MKL_Complex16 *One, *Zero;
    MKL_Complex16 *stored_product1_up, *stored_product1_down;
    MKL_Complex16 *density_matrix_up, *density_matrix_down;
    MKL_Complex16 *det_overlap_up, *det_overlap_down;
    MKL_Complex16 *overlap_inverse_up, *overlap_inverse_down;
    MKL_Complex16 *overlap_total;
    MKL_Complex16 product_dets;
    MKL_Complex16 local_ke;
    MKL_Complex16 local_pe;
    MKL_Complex16 pe, ke;

    /*Computes the Spin Up and Down Greens Functions*/
    stored_product1_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    stored_product1_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

    density_matrix_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
    density_matrix_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

    det_overlap_up=(MKL_Complex16 *)calloc(ist.n_determinants_trial_energy,sizeof(MKL_Complex16));
    det_overlap_down=(MKL_Complex16 *)calloc(ist.n_determinants_trial_energy,sizeof(MKL_Complex16));

    overlap_inverse_up=(MKL_Complex16 *)calloc(ist.n_determinants_energy_n_up_sq,sizeof(MKL_Complex16));
    overlap_inverse_down=(MKL_Complex16 *)calloc(ist.n_determinants_energy_n_down_sq,sizeof(MKL_Complex16));

    overlap_total=(MKL_Complex16 *)calloc(1,sizeof(MKL_Complex16));

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    /*Determine Overlaps With Full Trial Energy WF Instead of Just Phaseless One*/
    update_overlaps_chemistry_multi_energy_both(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_inverse_up,overlap_inverse_down,&overlap_total[0],det_overlap_up,det_overlap_down,ist);

    /*Run Through Determinants*/
    ke.real = ke.imag = 0.0;
    pe.real = pe.imag = 0.0;

    /*Determine Full Overlap*/
    for (i=0; i<ist.n_determinants_trial_energy; i++) {

     /*Check Whether Determinant Is Non Zero So That Singularities Don't Occur*/
     if ( Cabs(det_overlap_up[i]) > 0 || Cabs(det_overlap_down[i]) > 0 ) {

       /*Zero All Total Vectors*/
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_up,ist.n_spatial_orbitals,ist.n_up,One,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up,&trial_wf_up_energy[i*ist.n_spatial_orbitals_n_up],ist.n_up,Zero,stored_product1_up,ist.n_spatial_orbitals);
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_down,ist.n_spatial_orbitals,ist.n_down,One,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down,&trial_wf_down_energy[i*ist.n_spatial_orbitals_n_down],ist.n_down,Zero,stored_product1_down,ist.n_spatial_orbitals);

       /*Now Get Total Density Matrix*/
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up,One,wf_up,ist.n_up,stored_product1_up,ist.n_spatial_orbitals,Zero,density_matrix_up,ist.n_spatial_orbitals);
       cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down,One,wf_down,ist.n_down,stored_product1_down,ist.n_spatial_orbitals,Zero,density_matrix_down,ist.n_spatial_orbitals);

       /*Get Ke and PE*/
       local_ke = compute_kinetic_energy_density_matrix_restricted_sparse(density_matrix_up,density_matrix_down,kinetic_matrix_sparse,kinetic_ij_sparse,ist); 
       local_pe = compute_potential_energy_density_matrix_restricted_sparse(density_matrix_up,density_matrix_down,potential_matrix_sparse,potential_ij_sparse,ist); 

       /*Get Product of dets*/
       product_dets = Cmul(conjugate(trial_determinant_coefficients_energy[i]), Cmul(det_overlap_up[i], det_overlap_down[i]));

       /*Scale By Overlap and Coefficient*/
       local_ke = Cmul(local_ke, product_dets);
       local_pe = Cmul(local_pe, product_dets);

       /*Add Up Local Pe and Ke*/
       ke = Cadd(ke, local_ke);
       pe = Cadd(pe, local_pe);

      }

    } /*Check Determinant*/
    *total_ke = Cdiv(ke, overlap_total[0]);
    *total_pe = Cdiv(pe, overlap_total[0]);

free(stored_product1_up);
free(stored_product1_down);
free(density_matrix_up);
free(density_matrix_down);
free(det_overlap_up);
free(det_overlap_down);
free(overlap_inverse_up);
free(overlap_inverse_down);
return;
}

/**********************************************************************************************************************************************************/

void compute_walker_energy_chemistry_mixed_multi_phaseless_up_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_up,MKL_Complex16 *det_overlap_up,MKL_Complex16 *overlap_inverse_up,double *kinetic_matrix_sparse,int *kinetic_ij_sparse,MKL_Complex16 *total_ke,MKL_Complex16 *total_pe,int_st ist) {

  /*Computes the Kinetic and Potential Energies of a Wavefunction If It Is Crossed with a Multideterminant Trial Wave Function*/
    int i, j, k;
    MKL_Complex16 *One, *Zero;
    MKL_Complex16 *stored_product1_up; 
    MKL_Complex16 *density_matrix_up; 
    MKL_Complex16 product_dets; 
    MKL_Complex16 local_ke;
    MKL_Complex16 local_pe;
    MKL_Complex16 ke, pe;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product1_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    density_matrix_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

    /*Run Through Determinants*/
    ke.real = ke.imag = 0.0;
    pe.real = pe.imag = 0.0;
    /*Computes the Spin Up and Down Density Matrices*/
    for (i=0; i<ist.n_determinants_trial_phaseless; i++) {

       /*Check Whether Determinant Is Non Zero So That Singularities Don't Occur*/
       if ( Cabs(det_overlap_up[i]) > 0 ) {

        /*Zero All Total Vectors*/
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_up,ist.n_spatial_orbitals,ist.n_up,One,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up,&trial_wf_up_phaseless[i*ist.n_spatial_orbitals_n_up],ist.n_up,Zero,stored_product1_up,ist.n_spatial_orbitals);

        /*Now Get Total Density Matrix*/
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up,One,wf_up,ist.n_up,stored_product1_up,ist.n_spatial_orbitals,Zero,density_matrix_up,ist.n_spatial_orbitals);

        /*Get Ke and PE*/
        local_ke = compute_kinetic_energy_density_matrix_restricted_up_sparse(density_matrix_up,kinetic_matrix_sparse,kinetic_ij_sparse,ist);

        /*Get Product of dets*/
        product_dets = Cmul(conjugate(trial_determinant_coefficients_phaseless[i]), det_overlap_up[i]);

        /*Scale By Overlap and Coefficient*/
        local_ke = Cmul(local_ke, product_dets);

        /*Add Up Local Pe and Ke*/
        ke = Cadd(ke, local_ke);

       } /*Check Determinant*/

    }/*for*/
    *total_ke = Cdiv(ke, overlap_up[0]); 
    *total_pe = Cdiv(pe, overlap_up[0]); 


free(density_matrix_up); 
free(stored_product1_up); 
free(Zero); 
free(One); 
return; 
}

/*********************************************************************************************************************************************************************************/

void compute_walker_energy_chemistry_mixed_multi_energy_up_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_determinant_coefficients_energy,double *kinetic_matrix_sparse, int *kinetic_ij_sparse,MKL_Complex16 *total_ke,MKL_Complex16 *total_pe,int_st ist) {

  /*Computes the Kinetic and Potential Energies of a Wavefunction If It Is Crossed with a Multideterminant Trial Wave Function*/
    int i, j, k;
    MKL_Complex16 *One, *Zero;
    MKL_Complex16 *stored_product1_up;
    MKL_Complex16 *density_matrix_up;
    MKL_Complex16 *overlap_inverse_up; 
    MKL_Complex16 *det_overlap_up; 
    MKL_Complex16 overlap_total; 
    MKL_Complex16 product_dets;
    MKL_Complex16 local_ke;
    MKL_Complex16 local_pe;
    MKL_Complex16 ke, pe;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product1_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    density_matrix_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
    overlap_inverse_up=(MKL_Complex16 *)calloc(ist.n_determinants_energy_n_up_sq,sizeof(MKL_Complex16)); 
    det_overlap_up=(MKL_Complex16 *)calloc(ist.n_determinants_trial_energy,sizeof(MKL_Complex16)); 

    /*Update Overlaps*/
    update_overlaps_chemistry_multi_energy_up(wf_up,trial_wf_up_energy,trial_determinant_coefficients_energy,overlap_inverse_up,&overlap_total,det_overlap_up,ist);

    /*Run Through Determinants*/
    ke.real = ke.imag = 0.0;
    pe.real = pe.imag = 0.0;
    /*Computes the Spin Up and Down Density Matrices*/
    for (i=0; i<ist.n_determinants_trial_energy; i++) {

       /*Check Whether Determinant Is Non Zero So That Singularities Don't Occur*/
       if ( Cabs(det_overlap_up[i]) > 0 ) {

        /*Zero All Total Vectors*/
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_up,ist.n_spatial_orbitals,ist.n_up,One,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up,&trial_wf_up_energy[i*ist.n_spatial_orbitals_n_up],ist.n_up,Zero,stored_product1_up,ist.n_spatial_orbitals);

        /*Now Get Total Density Matrix*/
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up,One,wf_up,ist.n_up,stored_product1_up,ist.n_spatial_orbitals,Zero,density_matrix_up,ist.n_spatial_orbitals);

        /*Get Ke and PE*/
        local_ke = compute_kinetic_energy_density_matrix_restricted_up_sparse(density_matrix_up,kinetic_matrix_sparse,kinetic_ij_sparse,ist);

        /*Get Product of dets*/
        product_dets = Cmul(conjugate(trial_determinant_coefficients_energy[i]), det_overlap_up[i]);

        /*Scale By Overlap and Coefficient*/
        local_ke = Cmul(local_ke, product_dets);

        /*Add Up Local Pe and Ke*/
        ke = Cadd(ke, local_ke);

       } /*Check Determinant*/

    }/*for*/
    *total_ke = Cdiv(ke, overlap_total); 
    *total_pe = Cdiv(pe, overlap_total); 


free(density_matrix_up);
free(stored_product1_up);
free(overlap_inverse_up); 
free(det_overlap_up); 
free(Zero);
free(One);
return;
}

/***************************************************************************************************************************************************************************/

void compute_walker_energy_chemistry_mixed_multi_phaseless_down_sparse(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_down,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_down,double *kinetic_matrix_sparse,int *kinetic_ij_sparse,MKL_Complex16 *total_ke,MKL_Complex16 *total_pe,int_st ist) {

  /*Computes the Kinetic and Potential Energies of a Wavefunction If It Is Crossed with a Multideterminant Trial Wave Function*/
    int i, j, k;
    MKL_Complex16 *One, *Zero;
    MKL_Complex16 *stored_product1_down;  
    MKL_Complex16 *density_matrix_down;  
    MKL_Complex16 product_dets; 
    MKL_Complex16 local_ke;
    MKL_Complex16 local_pe;
    MKL_Complex16 ke, pe;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product1_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));
    density_matrix_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

    /*Run Through Determinants*/
    ke.real = ke.imag = 0.0;
    pe.real = pe.imag = 0.0;

    /*Computes the Spin Up and Down Density Matrices*/
    for (i=0; i<ist.n_determinants_trial_phaseless; i++) {

       /*Check Whether Determinant Is Non Zero So That Singularities Don't Occur*/
       if ( Cabs(det_overlap_down[i]) > 0 ) {

        /*Zero All Total Vectors*/
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_down,ist.n_spatial_orbitals,ist.n_down,One,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down,&trial_wf_down_phaseless[i*ist.n_spatial_orbitals_n_down],ist.n_down,Zero,stored_product1_down,ist.n_spatial_orbitals);

        /*Now Get Total Density Matrix*/
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down,One,wf_down,ist.n_down,stored_product1_down,ist.n_spatial_orbitals,Zero,density_matrix_down,ist.n_spatial_orbitals);

        /*Get Ke and PE*/
        local_ke = compute_kinetic_energy_density_matrix_restricted_down_sparse(density_matrix_down,kinetic_matrix_sparse,kinetic_ij_sparse,ist);

        /*Get Product of dets*/
        product_dets = Cmul(conjugate(trial_determinant_coefficients_phaseless[i]), det_overlap_down[i]);

        /*Scale By Overlap and Coefficient*/
        local_ke = Cmul(local_ke, product_dets);

        /*Add Up Local Pe and Ke*/
        ke = Cadd(ke, local_ke);

       } /*Check Determinant*/

    }/*for*/
    *total_ke = Cdiv(ke, overlap_down[0]); 
    *total_pe = Cdiv(pe, overlap_down[0]); 

free(density_matrix_down);
free(stored_product1_down);
return;
}

/**********************************************************************************************************************************************************************/

void compute_walker_energy_chemistry_mixed_multi_energy_down_sparse(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,double *kinetic_matrix_sparse,int *kinetic_ij_sparse,MKL_Complex16 *total_ke,MKL_Complex16 *total_pe,int_st ist) {

  /*Computes the Kinetic and Potential Energies of a Wavefunction If It Is Crossed with a Multideterminant Trial Wave Function*/
    int i, j, k;
    MKL_Complex16 *One, *Zero;
    MKL_Complex16 *stored_product1_down;
    MKL_Complex16 *density_matrix_down;
    MKL_Complex16 *overlap_inverse_down;
    MKL_Complex16 *det_overlap_down;
    MKL_Complex16 overlap_total;
    MKL_Complex16 product_dets;
    MKL_Complex16 local_ke;
    MKL_Complex16 local_pe;
    MKL_Complex16 ke, pe;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product1_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));
    density_matrix_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
    overlap_inverse_down=(MKL_Complex16 *)calloc(ist.n_determinants_energy_n_down_sq,sizeof(MKL_Complex16));
    det_overlap_down=(MKL_Complex16 *)calloc(ist.n_determinants_trial_energy,sizeof(MKL_Complex16));

    /*Update Overlaps*/
    update_overlaps_chemistry_multi_energy_down(wf_down,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_inverse_down,&overlap_total,det_overlap_down,ist);

    /*Run Through Determinants*/
    ke.real = ke.imag = 0.0;
    pe.real = pe.imag = 0.0;
    /*Computes the Spin Up and Down Density Matrices*/
    for (i=0; i<ist.n_determinants_trial_energy; i++) {

       /*Check Whether Determinant Is Non Zero So That Singularities Don't Occur*/
       if ( Cabs(det_overlap_down[i]) > 0 ) {

        /*Zero All Total Vectors*/
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_down,ist.n_spatial_orbitals,ist.n_down,One,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down,&trial_wf_down_energy[i*ist.n_spatial_orbitals_n_down],ist.n_down,Zero,stored_product1_down,ist.n_spatial_orbitals);

        /*Now Get Total Density Matrix*/
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down,One,wf_down,ist.n_down,stored_product1_down,ist.n_spatial_orbitals,Zero,density_matrix_down,ist.n_spatial_orbitals);

        /*Get Ke and PE*/
        local_ke = compute_kinetic_energy_density_matrix_restricted_down_sparse(density_matrix_down,kinetic_matrix_sparse,kinetic_ij_sparse,ist);

        /*Get Product of dets*/
        product_dets = Cmul(conjugate(trial_determinant_coefficients_energy[i]), det_overlap_down[i]);

        /*Scale By Overlap and Coefficient*/
        local_ke = Cmul(local_ke, product_dets);

        /*Add Up Local Pe and Ke*/
        ke = Cadd(ke, local_ke);

       } /*Check Determinant*/

    }/*for*/
    *total_ke = Cdiv(ke, overlap_total);
    *total_pe = Cdiv(pe, overlap_total);


free(density_matrix_down);
free(stored_product1_down);
free(overlap_inverse_down);
free(det_overlap_down);
free(Zero);
free(One);
return;
}

/***************************************************************************************************************************************************************************/

void compute_walker_energy_chemistry_mixed_restricted_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,double *kinetic_matrix_sparse,double *potential_matrix_sparse,int *kinetic_ij_sparse,int *potential_ij_sparse,MKL_Complex16 *ke,MKL_Complex16 *pe,int_st ist) {

      /*Computes the Kinetic and Potential Energies of a Wavefunction If Crossed with a single Determinmant*/
      MKL_Complex16 *density_matrix_up, *density_matrix_down;   

      /*Declare Density Matrices*/
      density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16)); 
      density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

      /*Compute Density Matrices*/
      compute_density_matrix_chemistry_mixed_restricted(wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_inverse_up,overlap_inverse_down,density_matrix_up,density_matrix_down,ist);
     
     /*Compute Ke and Pe*/
     *ke = compute_kinetic_energy_density_matrix_restricted_sparse(density_matrix_up,density_matrix_down,kinetic_matrix_sparse,kinetic_ij_sparse,ist);
     *pe = compute_potential_energy_density_matrix_restricted_sparse(density_matrix_up,density_matrix_down,potential_matrix_sparse,potential_ij_sparse,ist);

free(density_matrix_up); 
free(density_matrix_down); 
return; 
}

/***************************************************************************************************************************************************************/

void compute_walker_energy_chemistry_mixed_unrestricted_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *ke,MKL_Complex16 *pe,int_st ist) {

      /*Computes the Kinetic and Potential Energies of a Wavefunction If Crossed with a single Determinmant*/
      MKL_Complex16 *density_matrix_up, *density_matrix_down;

      /*Declare Density Matrices*/
      density_matrix_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
      density_matrix_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

      /*Compute Density Matrices*/
      compute_density_matrix_chemistry_mixed_restricted(wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_inverse_up,overlap_inverse_down,density_matrix_up,density_matrix_down,ist);

     /*Compute Ke and Pe*/
     *ke = compute_kinetic_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,kinetic_ij_sparse_up,kinetic_ij_sparse_down,ist);
     *pe = compute_potential_energy_density_matrix_unrestricted_sparse(density_matrix_up,density_matrix_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,ist);

free(density_matrix_up);
free(density_matrix_down);
return;
}

/*************************************************************************************************************************************************************/

void compute_energy_chemistry_mixed_two_body_density_matrix_multi_unrestricted_sparse(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,double *kinetic_matrix_sparse_up,double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up,int *kinetic_ij_sparse_down,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *energy_shifts,MKL_Complex16 *total_potential_energy,MKL_Complex16 *total_kinetic_energy,MKL_Complex16 *total_energy_2,MKL_Complex16 *total_weight,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *pe_error,MKL_Complex16 *ke,MKL_Complex16 *ke_error,MKL_Complex16 *te,MKL_Complex16 *te_error,MKL_Complex16 *tw,MKL_Complex16 *average_density_matrix_up,MKL_Complex16 *average_density_matrix_down,MKL_Complex16 *average_two_body_density_matrix,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   /*Computes the Energy of the System for a Mixed Estimator and Compute Two-Body Density Matrix*/
   int walker, i, j, k;
   int number_of_accumulated_steps;
   int number_walkers_active = 0;
   double weight_factor = 1.0;
   double energy_error, kinetic_energy_error_real, potential_energy_error_real;
   MKL_Complex16 kinetic_energy_error_complex, potential_energy_error_complex, total_energy_error_complex;
   MKL_Complex16 *average_wave_function_up, *average_wave_function_down;
   MKL_Complex16 local_pe, local_ke, local_pe_squared, local_ke_squared;
   MKL_Complex16 Zero, One;
   MKL_Complex16 potential_energy, kinetic_energy, total_energy;
   MKL_Complex16 potential_energy_squared, kinetic_energy_squared;
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator, inverse_weight_denominator;
   MKL_Complex16 energy_average;

   Zero.real = Zero.imag = 0.0;
   One.real = 1.0; One.imag = 0.0;
   potential_energy.real = potential_energy.imag = 0.0;
   kinetic_energy.real = kinetic_energy.imag = 0.0;
   total_energy.real = total_energy.imag = 0.0;
   potential_energy_squared.real = potential_energy_squared.imag = 0.0;
   kinetic_energy_squared.real = kinetic_energy_squared.imag;
   weight_denominator.real = weight_denominator.imag = 0.0;

   /*FILE *pf = fopen("energy.dat", "a");*/
   /*FILE *pf9 = fopen("energy_only.dat", "a"); */
   /*FILE *pf8 = fopen("restart.par", "w+");  */
   /*
 *  *  *      FILE *pf7 = fopen("overall_energy.dat", "a+");
 *   *   *        */

   average_wave_function_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_up*sizeof(MKL_Complex16));
 average_wave_function_down=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_down*sizeof(MKL_Complex16));

   czero_vec(average_two_body_density_matrix,ist.n_spin_orbitals_fourth);
   czero_vec(average_wave_function_up,ist.n_spatial_orbitals*ist.n_up);
   czero_vec(average_wave_function_down,ist.n_spatial_orbitals*ist.n_down);
   czero_vec(average_density_matrix_up,ist.n_spatial_orbitals_sq);
   czero_vec(average_density_matrix_down,ist.n_spatial_orbitals_sq);

   /*Obtain Population Control Rescaling Factor*/
   for (i = 0; i < 10 ; i++) {
    weight_factor *= weight_rescaling[i];
   }

   for (walker = 0; walker < ist.n_walkers ; walker++ ) {

     /*Check That Walker Weights Are Non-Zero*/
     if ( Cabs(weights[walker]) != 0 ) {
        number_walkers_active++;

        /*Accumulate Weight in Denominator*/
        current_weight = RCmul(weight_factor, weights[walker]);
        weight_denominator = Cadd(weight_denominator, current_weight);

        /*Get Energies and Density Matrices*/
        compute_walker_energy_chemistry_mixed_twobody_multi_energy_unrestricted_sparse(&wf_up[walker*ist.n_spatial_orbitals_n_up],&wf_down[walker*ist.n_spatial_orbitals_n_down],trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,&overlap_total[walker],&det_overlap_up[walker*ist.n_determinants_trial_energy],&det_overlap_down[walker*ist.n_determinants_trial_energy],&overlap_inverse_up[walker*ist.n_determinants_n_up_sq],&overlap_inverse_down[walker*ist.n_determinants_n_down_sq],kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_updown,average_two_body_density_matrix,current_weight,&local_ke,&local_pe,ist);

        /*Obtain Kinetic Energy*/
        local_ke_squared = Cmul(local_ke, conjugate(local_ke));
        local_ke = Cmul(local_ke, current_weight);
        local_ke_squared = Cmul(local_ke_squared, current_weight);

        /*Get Potential Energy*/
        local_pe_squared = Cmul(local_pe, conjugate(local_pe));
        local_pe = Cmul(local_pe, current_weight);
        local_pe_squared = Cmul(local_pe_squared, current_weight);

        /*Add to Total*/
        potential_energy = Cadd(local_pe, potential_energy);
 potential_energy_squared = Cadd(local_pe_squared, potential_energy_squared);
        kinetic_energy = Cadd(local_ke, kinetic_energy);
        kinetic_energy_squared = Cadd(kinetic_energy_squared, local_ke_squared);

        /*Obtain Average Wavefunction As Well*/
        for (i=0; i<ist.n_spatial_orbitals; i++) {
         for (j=0; j<ist.n_up; j++) {
           average_wave_function_up[i*ist.n_up+j] = Cadd(average_wave_function_up[i*ist.n_up+j], Cmul(current_weight, wf_up[walker*ist.n_spatial_orbitals_n_up+i*ist.n_up+j]));
         }
         for (j=0; j<ist.n_down; j++) {
           average_wave_function_down[i*ist.n_down+j] = Cadd(average_wave_function_down[i*ist.n_down+j], Cmul(current_weight, wf_down[walker*ist.n_spatial_orbitals_n_down+i*ist.n_down+j]));
         }
        }
       } /*Non-Zero Walker Weights*/
     }  /*Walkers*/

    /*Accumulate the Energy of Each Block To Get Overall Statistics*/
    if ( step >= ist.n_steps_equilibration ) {
      accumulate_energy_blocks(potential_energy, kinetic_energy, Cadd(potential_energy, kinetic_energy), weight_denominator, total_potential_energy, total_kinetic_energy, total_energy_2, total_weight);
    }

   /*Divide By Total Weight*/
   potential_energy = Cdiv(potential_energy, weight_denominator);
   potential_energy_squared = Cdiv(potential_energy_squared, weight_denominator);
   kinetic_energy = Cdiv(kinetic_energy, weight_denominator);
   kinetic_energy_squared = Cdiv(kinetic_energy_squared, weight_denominator);
   total_energy = Cadd(potential_energy, kinetic_energy);

   /*Also Add Initial Energy Shift from MolPro*/
   total_energy.real += energy_shifts[0].real;

   /*Print Desired Information*******************************************/
   /*Find Errors*/
   potential_energy_error_complex = Csub(potential_energy_squared, Cmul(potential_energy, conjugate(potential_energy)));
   if ( potential_energy_error_complex.real >  0 ) {
 potential_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * potential_energy_error_complex.real);
   }
   else {
     potential_energy_error_complex.real = 0.0;
   }

   kinetic_energy_error_complex = Csub(kinetic_energy_squared, Cmul(kinetic_energy, conjugate(kinetic_energy)));
   if ( kinetic_energy_error_complex.real > 0 ) {
     kinetic_energy_error_complex.real = sqrt(1.0/((double) number_walkers_active) * kinetic_energy_error_complex.real);
   }
   else {
    kinetic_energy_error_complex.real = 0.0;
   }
   total_energy_error_complex.real = sqrt(potential_energy_error_complex.real*potential_energy_error_complex.real +   kinetic_energy_error_complex.real*kinetic_energy_error_complex.real);

   /*Print Desired Information*******************************************/
   *ke = kinetic_energy;
   *pe = potential_energy;
   *te = total_energy;
   *ke_error = kinetic_energy_error_complex;
   *pe_error = potential_energy_error_complex;
   *te_error = total_energy_error_complex;
   *tw = weight_denominator;

   /*Obtain and Print Average WFs*/
   for (i=0; i<ist.n_spatial_orbitals; i++) {
     for (j=0; j<ist.n_up; j++) {
       average_wave_function_up[i*ist.n_up+j] = Cdiv(average_wave_function_up[i*ist.n_up+j], weight_denominator);
    }
  }
   for (i=0; i<ist.n_spatial_orbitals; i++) {
     for (j=0; j<ist.n_down; j++) {
       average_wave_function_down[i*ist.n_down+j] = Cdiv(average_wave_function_down[i*ist.n_down+j], weight_denominator);
     }
   }

  /*Find Normalized WF*/
  get_normalized_wf_up_chemistry(normalized_wave_function_up,average_wave_function_up,ist);
  get_normalized_wf_down_chemistry(normalized_wave_function_down,average_wave_function_down,ist);
  /*Find Average Reduced Density Matrices*/
  inverse_weight_denominator = Cdiv(One, weight_denominator);
  cblas_zscal(ist.n_spatial_orbitals_sq, &inverse_weight_denominator, average_density_matrix_up, 1);
  cblas_zscal(ist.n_spatial_orbitals_sq, &inverse_weight_denominator, average_density_matrix_down, 1);
  cblas_zscal(ist.n_spin_orbitals_fourth, &inverse_weight_denominator, average_two_body_density_matrix, 1);


  /*Find Statistics So Far*******************************************
 *  *     if ( step >= ist.n_steps_equilibration ) {
 *   *               number_of_accumulated_steps = 1+(int)((step-ist.n_steps_equilibration)/(ist.n_steps_energy));
 *    *                               fprintf(pf7, "%d\t", step);
 *     *
 *      *                                   /*Get PE*/
    /*energy_average = Cdiv(total_potential_energy[0], total_weight[0]);
 *     potential_energy_error_real = error_real(total_potential_energy[1].real, total_potential_energy[2].real, number_of_accumulated_steps);
 *         fprintf(pf7, "%f   %f\t", energy_average.real, potential_energy_error_real);
 *
 *             /*Get KE*/
    /*energy_average = Cdiv(total_kinetic_energy[0], total_weight[0]);
 *     kinetic_energy_error_real = error_real(total_kinetic_energy[1].real, total_kinetic_energy[2].real, number_of_accumulated_steps);
 *         fprintf(pf7, "%f   %f\t", energy_average.real, kinetic_energy_error_real);
 *
 *             /*Show Shift*/
    /*fprintf(pf7, "%f   %f\t", energy_shifts[0].real, 0.0);
 *
 *     /*Get Total*/
    /*energy_average = Cdiv(total_energy_2[0], total_weight[0]);
 *     energy_error = sqrt(potential_energy_error_real*potential_energy_error_real+kinetic_energy_error_real*kinetic_energy_error_real);
 *         energy_error = error_real(total_energy_2[1].r, total_energy_2[2].r, number_of_accumulated_steps);
 *             fprintf(pf7, "%f   %f\n", energy_average.real+energy_shifts[0].real, energy_error);
 *                 fflush(pf7);
 *                     }
 *                         */


free(average_wave_function_up);
free(average_wave_function_down);
return;
}

/**************************************************************************************************************************************************************/

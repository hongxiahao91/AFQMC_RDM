#include "afqmc.h"

/***********************************************************************************************/

void compute_energy_chemistry_bp_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,MKL_Complex16 *prev_wf_up,MKL_Complex16 *prev_wf_down,MKL_Complex16 *weights,MKL_Complex16 *kinetic_original_matrix, MKL_Complex16 *potential_original_matrix,MKL_Complex16 *energy_shifts,MKL_Complex16 *total_potential_energy,MKL_Complex16 *total_kinetic_energy,MKL_Complex16 *total_energy_2,MKL_Complex16 *total_weight,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *pe_error,MKL_Complex16 *ke,MKL_Complex16 *ke_error,MKL_Complex16 *te,MKL_Complex16 *te_error,MKL_Complex16 *tw,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   //Computes the Energy of the System Using Mixed Estimator for a Restricted Situation
   int walker, i, j, k;
   int number_of_accumulated_steps;  
   int number_walkers_active = 0; 
   double weight_factor = 1.0;
   double energy_error, kinetic_energy_error_real, potential_energy_error_real;
   MKL_Complex16 kinetic_energy_error_complex, potential_energy_error_complex, total_energy_error_complex; 
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   MKL_Complex16 *average_wave_function_up, *average_wave_function_down; 
   MKL_Complex16 local_pe, local_ke, local_pe_squared, local_ke_squared;
   MKL_Complex16 Zero; 
   MKL_Complex16 potential_energy, kinetic_energy, total_energy;
   MKL_Complex16 potential_energy_squared, kinetic_energy_squared; 
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator;
   MKL_Complex16 energy_average; 
   FILE *pf = fopen("errors.dat", "a+");  

   Zero.real = Zero.imag = 0.0; 
   potential_energy.real = potential_energy.imag = 0.0; 
   kinetic_energy.real = kinetic_energy.imag = 0.0; 
   total_energy.real = total_energy.imag = 0.0; 
   potential_energy_squared.real = potential_energy_squared.imag = 0.0; 
   kinetic_energy_squared.real = kinetic_energy_squared.imag; 
   weight_denominator.real = weight_denominator.imag = 0.0; 

  // FILE *pf8 = fopen("restart.par", "w+");  
   /*FILE *pf7 = fopen("overall_energy.dat", "a+");  
  */

   density_matrix_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_spatial_orbitals*sizeof(MKL_Complex16));
   density_matrix_down=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_spatial_orbitals*sizeof(MKL_Complex16));

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

        //First Compute the Density Matrix
        compute_density_matrix_chemistry_bp(&prev_wf_up[walker*ist.n_spatial_orbitals_n_up],&prev_wf_down[walker*ist.n_spatial_orbitals_n_down],&bp_wf_up[walker*ist.n_spatial_orbitals_n_up],&bp_wf_down[walker*ist.n_spatial_orbitals_n_down],density_matrix_up,density_matrix_down,ist); 

        //Obtain Kinetic Energy
        local_ke = compute_kinetic_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,kinetic_original_matrix,ist);
        local_ke_squared = Cmul(local_ke, conjugate(local_ke)); 
        local_ke = Cmul(local_ke, current_weight);
        local_ke_squared = Cmul(local_ke_squared, current_weight); 

        //Get Potential Energy
        local_pe = compute_potential_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,potential_original_matrix,ist);
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
      accumulate_energy_blocks(potential_energy, kinetic_energy, Cadd(potential_energy, kinetic_energy), weight_denominator, total_potential_energy, total_kinetic_energy, total_energy_2, total_weight); 
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
/*  if ( step >= ist.n_steps_equilibration ) {
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

free(density_matrix_up);
free(density_matrix_down);
free(average_wave_function_up); 
free(average_wave_function_down); 
return; 
}

/************************************************************************************************/

void compute_energy_chemistry_bp_two_body_density_matrix_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,MKL_Complex16 *prev_wf_up,MKL_Complex16 *prev_wf_down,MKL_Complex16 *weights,MKL_Complex16 *kinetic_original_matrix, MKL_Complex16 *potential_original_matrix,MKL_Complex16 *energy_shifts,MKL_Complex16 *total_potential_energy,MKL_Complex16 *total_kinetic_energy,MKL_Complex16 *total_energy_2,MKL_Complex16 *total_weight,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *pe_error,MKL_Complex16 *ke,MKL_Complex16 *ke_error,MKL_Complex16 *te,MKL_Complex16 *te_error,MKL_Complex16 *tw,MKL_Complex16 *average_density_matrix_up,MKL_Complex16 *average_density_matrix_down,MKL_Complex16 *average_two_body_density_matrix,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   /*Computes the Energy of the System for a Mixed Estimator and Compute Two-Body Density Matrix*/
   int walker, i, j, k;
   int number_of_accumulated_steps;
   int number_walkers_active = 0;
   double weight_factor = 1.0;
   double energy_error, kinetic_energy_error_real, potential_energy_error_real;
   MKL_Complex16 kinetic_energy_error_complex, potential_energy_error_complex, total_energy_error_complex;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
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

   /*FILE *pf8 = fopen("restart.par", "w+");  */
   /*  FILE *pf7 = fopen("overall_energy.dat", "a+");  
  */

   density_matrix_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_spatial_orbitals*sizeof(MKL_Complex16));
   density_matrix_down=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals*ist.n_spatial_orbitals*sizeof(MKL_Complex16));

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

        /*First Compute the Density Matrix*/
        compute_density_matrix_chemistry_bp(&prev_wf_up[walker*ist.n_spatial_orbitals_n_up],&prev_wf_down[walker*ist.n_spatial_orbitals_n_down],&bp_wf_up[walker*ist.n_spatial_orbitals_n_up],&bp_wf_down[walker*ist.n_spatial_orbitals_n_down],density_matrix_up,density_matrix_down,ist); 

        /*Compute 2RDM If Desired*/
        compute_two_body_reduced_density_matrix(density_matrix_up,density_matrix_down,average_two_body_density_matrix,current_weight,ist);                                                                               

        /*Obtain Kinetic Energy*/
        local_ke = compute_kinetic_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,kinetic_original_matrix,ist);
        local_ke_squared = Cmul(local_ke, conjugate(local_ke));
        local_ke = Cmul(local_ke, current_weight);
        local_ke_squared = Cmul(local_ke_squared, current_weight);

        /*Get Potential Energy*/
        local_pe = compute_potential_energy_density_matrix_restricted(density_matrix_up,density_matrix_down,potential_original_matrix,ist);
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
         for (j=0; j<ist.n_spatial_orbitals; j++) {
           average_density_matrix_up[i*ist.n_spatial_orbitals+j] = Cadd(average_density_matrix_up[i*ist.n_spatial_orbitals+j], Cmul(density_matrix_up[i*ist.n_spatial_orbitals+j], current_weight));
           average_density_matrix_down[i*ist.n_spatial_orbitals+j] = Cadd(average_density_matrix_down[i*ist.n_spatial_orbitals+j], Cmul(density_matrix_down[i*ist.n_spatial_orbitals+j], current_weight));
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
  if ( step >= ist.n_steps_equilibration ) {
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

free(density_matrix_up);
free(density_matrix_down);
free(average_wave_function_up);
free(average_wave_function_down);
return;
}


/***********************************************************************************************/

void compute_energy_chemistry_bp_unrestricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,MKL_Complex16 *prev_wf_up,MKL_Complex16 *prev_wf_down,MKL_Complex16 *weights,MKL_Complex16 *kinetic_original_matrix_up, MKL_Complex16 *kinetic_original_matrix_down, MKL_Complex16 *potential_original_matrix_up, MKL_Complex16 *potential_original_matrix_down, MKL_Complex16 *potential_original_matrix_updown, MKL_Complex16 *energy_shifts,MKL_Complex16 *total_potential_energy,MKL_Complex16 *total_kinetic_energy,MKL_Complex16 *total_energy_2,MKL_Complex16 *total_weight,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *pe_error,MKL_Complex16 *ke,MKL_Complex16 *ke_error,MKL_Complex16 *te,MKL_Complex16 *te_error,MKL_Complex16 *tw,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   /*Computes the Energy of the System For an Unrestricted Situation Using the Mixed Estimator*/
   int walker, i, j, k;
   int number_of_accumulated_steps, number_walkers_active;
   double weight_factor = 1.0;
   double energy_error;
   double kinetic_energy_error_real, potential_energy_error_real;
   MKL_Complex16 kinetic_energy_error_complex, potential_energy_error_complex, total_energy_error_complex;
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   MKL_Complex16 *average_wave_function_up, *average_wave_function_down;
   MKL_Complex16 local_pe, local_ke, local_pe_squared, local_ke_squared;
   MKL_Complex16 potential_energy = Zero, kinetic_energy = Zero, total_energy = Zero;
   MKL_Complex16 potential_energy_squared = Zero, kinetic_energy_squared = Zero;
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator = Zero;
   MKL_Complex16 energy_average;

   density_matrix_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

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

        /*First Compute the Density Matrix*/
        compute_density_matrix_chemistry_bp(&prev_wf_up[walker*ist.n_spatial_orbitals_n_up],&prev_wf_down[walker*ist.n_spatial_orbitals_n_down],&bp_wf_up[walker*ist.n_spatial_orbitals_n_up],&bp_wf_down[walker*ist.n_spatial_orbitals_n_down],density_matrix_up,density_matrix_down,ist);

        /*Obtain Kinetic Energy*/
        local_ke = compute_kinetic_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,kinetic_original_matrix_up,kinetic_original_matrix_down,ist);
        local_ke_squared = Cmul(local_ke, conjugate(local_ke));
        local_ke = Cmul(local_ke, current_weight);
        local_ke_squared = Cmul(local_ke_squared, current_weight);

        /*Get Potential Energy*/
        local_pe = compute_potential_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,ist);
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

free(density_matrix_up);
free(density_matrix_down);
free(average_wave_function_up);
free(average_wave_function_down);
return;
}



/***********************************************************************************************/

void compute_energy_chemistry_unrestricted_bp_two_body_density_matrix(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,MKL_Complex16 *prev_wf_up,MKL_Complex16 *prev_wf_down,MKL_Complex16 *weights,MKL_Complex16 *kinetic_original_matrix_up, MKL_Complex16 *kinetic_original_matrix_down, MKL_Complex16 *potential_original_matrix_up, MKL_Complex16 *potential_original_matrix_down, MKL_Complex16 *potential_original_matrix_updown, MKL_Complex16 *energy_shifts,MKL_Complex16 *total_potential_energy,MKL_Complex16 *total_kinetic_energy,MKL_Complex16 *total_energy_2,MKL_Complex16 *total_weight,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *pe_error,MKL_Complex16 *ke,MKL_Complex16 *ke_error,MKL_Complex16 *te,MKL_Complex16 *te_error,MKL_Complex16 *tw,MKL_Complex16 *average_density_matrix_up,MKL_Complex16 *average_density_matrix_down,MKL_Complex16 *average_two_body_density_matrix,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   /*Computes the Energy of the System and Outputs Two-Body Density Matrices Using Mixed Estimator for an Unrestricted Situation*/
   int walker, i, j, k;  
   int number_of_accumulated_steps, number_walkers_active; 
   double weight_factor = 1.0;
   double energy_error; 
   double kinetic_energy_error_real, potential_energy_error_real; 
   MKL_Complex16 kinetic_energy_error_complex, potential_energy_error_complex, total_energy_error_complex;
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0; 
   MKL_Complex16 *density_matrix_up, *density_matrix_down;
   MKL_Complex16 *average_wave_function_up, *average_wave_function_down;
   MKL_Complex16 local_pe, local_ke, local_pe_squared, local_ke_squared; 
   MKL_Complex16 potential_energy = Zero, kinetic_energy = Zero, total_energy = Zero;
   MKL_Complex16 potential_energy_squared = Zero, kinetic_energy_squared = Zero; 
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator = Zero, inverse_weight_denominator;
   MKL_Complex16 energy_average; 

   density_matrix_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   density_matrix_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   average_wave_function_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_up,sizeof(MKL_Complex16));
   average_wave_function_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals*ist.n_down,sizeof(MKL_Complex16));

   czero_vec(average_two_body_density_matrix,ist.n_spin_orbitals_fourth); 
   czero_vec(average_wave_function_up,ist.n_spatial_orbitals*ist.n_up);
   czero_vec(average_wave_function_down,ist.n_spatial_orbitals*ist.n_down);
   czero_vec(average_density_matrix_up,ist.n_spatial_orbitals_sq);
   czero_vec(average_density_matrix_down,ist.n_spatial_orbitals_sq);

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

        /*First Compute the Density Matrix*/
        compute_density_matrix_chemistry_bp(&prev_wf_up[walker*ist.n_spatial_orbitals_n_up],&prev_wf_down[walker*ist.n_spatial_orbitals_n_down],&bp_wf_up[walker*ist.n_spatial_orbitals_n_up],&bp_wf_down[walker*ist.n_spatial_orbitals_n_down],density_matrix_up,density_matrix_down,ist);

        /*Compute Two-Body Density Matrix*/ 
        compute_two_body_reduced_density_matrix(density_matrix_up,density_matrix_down,average_two_body_density_matrix,current_weight,ist);

        /*Obtain Kinetic Energy*/
        local_ke = compute_kinetic_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,kinetic_original_matrix_up,kinetic_original_matrix_down,ist);
        local_ke_squared = Cmul(local_ke, conjugate(local_ke)); 
        local_ke = Cmul(local_ke, current_weight);
        local_ke_squared = Cmul(local_ke_squared, current_weight); 

        /*Get Potential Energy*/
        local_pe = compute_potential_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,ist);
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
         for (j=0; j<ist.n_spatial_orbitals; j++) {
           average_density_matrix_up[i*ist.n_spatial_orbitals+j] = Cadd(average_density_matrix_up[i*ist.n_spatial_orbitals+j], Cmul(density_matrix_up[i*ist.n_spatial_orbitals+j], current_weight));
           average_density_matrix_down[i*ist.n_spatial_orbitals+j] = Cadd(average_density_matrix_down[i*ist.n_spatial_orbitals+j], Cmul(density_matrix_down[i*ist.n_spatial_orbitals+j], current_weight));
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
  

   //Print Desired Information*******************************************
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
/*  for (i=0; i<ist.n_walkers; i++) {
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

   /*Find Average Reduced Density Matrices*/
  inverse_weight_denominator = Cdiv(One, weight_denominator);
  cblas_zscal(ist.n_spatial_orbitals_sq, &inverse_weight_denominator, average_density_matrix_up, 1);
  cblas_zscal(ist.n_spatial_orbitals_sq, &inverse_weight_denominator, average_density_matrix_down, 1);
  cblas_zscal(ist.n_spin_orbitals_fourth, &inverse_weight_denominator, average_two_body_density_matrix, 1);

  /*Find Statistics So Far*******************************************/
/*  if ( step >= ist.n_steps_equilibration ) {
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

   } //
 */

//fclose(pf);
//fclose(pf3);
//fclose(pf4);
//fclose(pf5);
//fclose(pf6);
//fclose(pf7);   
//fclose(pf8); 
free(density_matrix_up);
free(density_matrix_down);
free(average_wave_function_up);
free(average_wave_function_down);
return;
}

/**************************************************************************************************************************************************************/

void compute_density_matrix_chemistry_bp(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *bp_wf_up,MKL_Complex16 *bp_wf_down,MKL_Complex16 *density_matrix_up,MKL_Complex16 *density_matrix_down,int_st ist) {

  /*Compute the Back Propagated Density Matrix*/ 
  int i, j; 
  MKL_Complex16 *overlap_inverse_up, *overlap_inverse_down;      
  MKL_Complex16 *stored_product1_up, *stored_product1_down; 
  MKL_Complex16 *One,*Zero; 

  One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
  Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
  One[0].real = 1.0; One[0].imag = 0.0;
  Zero[0].real = 0.0; Zero[0].imag = 0.0;

 /*Declare Overlap Inverse Matrices*/
  overlap_inverse_up=(MKL_Complex16 *)calloc(ist.n_up_sq,sizeof(MKL_Complex16));
  overlap_inverse_down=(MKL_Complex16 *)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

  stored_product1_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
  stored_product1_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

  /*Zero Density Matrices*/
  czero_vec(density_matrix_up,ist.n_spatial_orbitals_sq); 
  czero_vec(density_matrix_down,ist.n_spatial_orbitals_sq); 

  /*Zero All Total Vectors*/
  czero_vec(stored_product1_up,ist.n_spatial_orbitals_n_up);
  czero_vec(stored_product1_down,ist.n_spatial_orbitals_n_down);

  /*Get The Overlap Inverse Matrices*/
  compute_overlap_inverse(bp_wf_up,wf_up,overlap_inverse_up,ist,0); 
  compute_overlap_inverse(bp_wf_down,wf_down,overlap_inverse_down,ist,1);  
     
  /*Get Total Stored Product From multiplying Into Trial Wf*/
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_up,ist.n_spatial_orbitals,ist.n_up,One,overlap_inverse_up,ist.n_up,bp_wf_up,ist.n_up,Zero,stored_product1_up,ist.n_spatial_orbitals);
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_down,ist.n_spatial_orbitals,ist.n_down,One,overlap_inverse_down,ist.n_down,bp_wf_down,ist.n_down,Zero,stored_product1_down,ist.n_spatial_orbitals);

  /*Now Get Total Density Matrix*/
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up,One,wf_up,ist.n_up,stored_product1_up,ist.n_spatial_orbitals,Zero,density_matrix_up,ist.n_spatial_orbitals);
  cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down,One,wf_down,ist.n_down,stored_product1_down,ist.n_spatial_orbitals,Zero,density_matrix_down,ist.n_spatial_orbitals);

free(overlap_inverse_up); free(overlap_inverse_down); 
free(stored_product1_up); free(stored_product1_down); 
return; 
}

/***********************************************************************************************************************************************************/

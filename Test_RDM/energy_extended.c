#include "afqmc.h"

/***********************************************************************************************/

void compute_energy_extended(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int *neighbors,int *number_neighbors,MKL_Complex16 *potential_energy_store, MKL_Complex16 *kinetic_energy_store, MKL_Complex16 *total_energy_store,MKL_Complex16 *total_weights,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *ke,MKL_Complex16 *te,MKL_Complex16 *tw,MKL_Complex16 *normalized_wave_function_up,MKL_Complex16 *normalized_wave_function_down) {

   /*Computes the Energy of the System*/ 
   int walker, i, j; 
   double weight_factor = 1.0; 
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;  
   MKL_Complex16 *average_wave_function_up, *average_wave_function_down; 
   MKL_Complex16 *density_matrix_up, *density_matrix_down; 
   MKL_Complex16 *overall_density_matrix_up, *overall_density_matrix_down; 
   MKL_Complex16 local_ke, local_pe_total, local_pe_Hubbard, local_pe_extended; 
   MKL_Complex16 density_prod; 
   MKL_Complex16 energy = Zero, kinetic_energy = Zero, potential_energy = Zero;
   MKL_Complex16 current_weight; 
   MKL_Complex16 weight_denominator = Zero; 

   average_wave_function_up=(MKL_Complex16 *)calloc(ist.n_sites*ist.n_up,sizeof(MKL_Complex16)); 
   average_wave_function_down=(MKL_Complex16 *)calloc(ist.n_sites*ist.n_down,sizeof(MKL_Complex16)); 

   density_matrix_up=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 
   density_matrix_down=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 
 
   overall_density_matrix_up=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));
   overall_density_matrix_down=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));

   czero_vec(average_wave_function_up,ist.n_sites*ist.n_up); 
   czero_vec(average_wave_function_down,ist.n_sites*ist.n_down);

   czero_vec(normalized_wave_function_up,ist.n_sites*ist.n_up); 
   czero_vec(normalized_wave_function_down,ist.n_sites*ist.n_down); 
 
   czero_vec(density_matrix_up,ist.n_sites_sq); 
   czero_vec(density_matrix_down,ist.n_sites_sq); 

   czero_vec(overall_density_matrix_up,ist.n_sites_sq); 
   czero_vec(overall_density_matrix_down,ist.n_sites_sq); 

   /*Obtain Population Control Rescaling Factor*/
   for (i = 0; i < 10 ; i++) {
    weight_factor *= weight_rescaling[i]; 
   } 

   for (walker = 0; walker < ist.n_walkers ; walker++ ) {

     /*Check That Walker Weights Are Non-Zero*/    
     if ( Cabs(weights[walker]) != 0 ) { 

        /*Accumulate Weight in Denominator*/
        current_weight = RCmul(weight_factor, weights[walker]); 
        weight_denominator = Cadd(weight_denominator, current_weight); 

        /*First Compute the Green's Function*/
        compute_density_matrix_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],density_matrix_up,density_matrix_down,ist);  

        /*Obtain Kinetic and Potential Energies*/ 
        local_ke = local_pe_total = local_pe_extended = local_pe_Hubbard = Zero;  
        for ( i=0; i<ist.n_sites; i++) {
         density_prod = Cmul( density_matrix_up[i*ist.n_sites+i], density_matrix_down[i*ist.n_sites+i]); 
         local_pe_Hubbard = Cadd(local_pe_Hubbard, density_prod); 
         if ( ist.n_sites_one > 2 || ist.n_sites_two > 2 ) { 
           for (j=0; j<number_neighbors[i]; j++) {
              local_ke = Csub(local_ke, Cadd(density_matrix_up[i*ist.n_sites+neighbors[4*i+j]],density_matrix_down[i*ist.n_sites+neighbors[4*i+j]])); 

              if ( j<(int)(number_neighbors[i]/2.0) ) { 
                density_prod = Cmul( Cadd(density_matrix_up[i*ist.n_sites+i], density_matrix_down[i*ist.n_sites+i]), Cadd(density_matrix_up[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]], density_matrix_down[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]]));  
                local_pe_extended = Cadd(local_pe_extended, density_prod); 
              }
 
           }
         }
         else { 
            for (j=0; j<number_neighbors[i]; j++) {
              local_ke =Csub(local_ke, RCmul(.5, Cadd(density_matrix_up[i*ist.n_sites+neighbors[4*i+j]],density_matrix_down[i*ist.n_sites+neighbors[4*i+j]])));

              if ( j<(int)(number_neighbors[i]/2.0) ) { 
                density_prod = Cmul( Cadd(density_matrix_up[i*ist.n_sites+i], density_matrix_down[i*ist.n_sites+i]), Cadd(density_matrix_up[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]], density_matrix_down[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]]));
                local_pe_extended = Cadd(local_pe_extended, density_prod); 
               }
           }
         }

          /*Construct Density Matrices*/
          for ( j=0; j<ist.n_sites; j++) {
            overall_density_matrix_up[i*ist.n_sites+j] = Cadd( overall_density_matrix_up[i*ist.n_sites+j], Cmul(density_matrix_up[i*ist.n_sites+j], current_weight)); 
            overall_density_matrix_down[i*ist.n_sites+j] = Cadd( overall_density_matrix_down[i*ist.n_sites+j], Cmul(density_matrix_down[i*ist.n_sites+j], current_weight)); 
          }
        }

        local_ke = RCmul(cns.t_fermions, Cmul(local_ke, current_weight));
        local_pe_Hubbard = RCmul(cns.U_fermions, Cmul(local_pe_Hubbard, current_weight)); 
        local_pe_extended = RCmul(cns.V_fermions, Cmul(local_pe_extended, current_weight));  
        local_pe_total = Cadd(local_pe_Hubbard, local_pe_extended); 

        /*Add to Total*/
        kinetic_energy = Cadd(local_ke, kinetic_energy); 
        potential_energy = Cadd(local_pe_total, potential_energy);          
        energy = Cadd(energy, Cadd(local_ke, local_pe_total));  

        /*Obtain Average Wavefunction As Well*/
        for (i=0; i<ist.n_sites; i++) {
         for (j=0; j<ist.n_up; j++) {
           average_wave_function_up[i*ist.n_up+j] = Cadd(average_wave_function_up[i*ist.n_up+j], Cmul(current_weight, wf_up[walker*ist.n_sites_n_up+i*ist.n_up+j]));
         }
         for (j=0; j<ist.n_down; j++) {   
           average_wave_function_down[i*ist.n_down+j] = Cadd(average_wave_function_down[i*ist.n_down+j], Cmul(current_weight, wf_down[walker*ist.n_sites_n_down+i*ist.n_down+j]));
         }
       }

      } /*Non-Zero Walker Weights*/
   }  /*Walkers*/

   /*Divide By Total Weight*/
   kinetic_energy = Cdiv(kinetic_energy, weight_denominator); 
   potential_energy = Cdiv(potential_energy, weight_denominator); 
   energy = Cdiv(energy, weight_denominator);  

   *ke = kinetic_energy; 
   *pe = potential_energy; 
   *te = energy; 
   *tw = weight_denominator; 

  /*Store the Energies for Further Averaging*/
  kinetic_energy_store[step] = kinetic_energy; 
  potential_energy_store[step] = potential_energy; 
  total_energy_store[step] = energy; 
  total_weights[step] = weight_denominator;   

  /*Obtain and Print Average WFs*/
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_up; j++) {
      average_wave_function_up[i*ist.n_up+j] = Cdiv(average_wave_function_up[i*ist.n_up+j], weight_denominator); 
    }
    for (j=0; j<ist.n_down; j++) {
      average_wave_function_down[i*ist.n_down+j] = Cdiv(average_wave_function_down[i*ist.n_down+j], weight_denominator); 
    }
  } 
  get_normalized_wf_up(normalized_wave_function_up,average_wave_function_up,ist);
  get_normalized_wf_down(normalized_wave_function_down,average_wave_function_down,ist);

free(density_matrix_up); 
free(density_matrix_down); 
free(overall_density_matrix_up); 
free(overall_density_matrix_down); 
free(average_wave_function_up); 
free(average_wave_function_down); 
return; 
}

/***********************************************************************************************************************************************************/

void compute_energy_extended_up(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int *neighbors,int *number_neighbors,MKL_Complex16 *potential_energy_store, MKL_Complex16 *kinetic_energy_store, MKL_Complex16 *total_energy_store,MKL_Complex16 *total_weights,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *ke,MKL_Complex16 *te,MKL_Complex16 *tw,MKL_Complex16 *normalized_wave_function_up) {

   /*Computes the Energy of the System*/
   int walker, i, j;
   double weight_factor = 1.0;
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;
   MKL_Complex16 *average_wave_function_up; 
   MKL_Complex16 *density_matrix_up; 
   MKL_Complex16 *overall_density_matrix_up; 
   MKL_Complex16 local_ke, local_pe; 
   MKL_Complex16 density_prod;
   MKL_Complex16 energy = Zero, kinetic_energy = Zero, potential_energy = Zero;
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator = Zero;
   FILE *pf = fopen("errors.dat", "a+");  

   average_wave_function_up=(MKL_Complex16 *)calloc(ist.n_sites*ist.n_up,sizeof(MKL_Complex16));

   density_matrix_up=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));
   overall_density_matrix_up=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));

   czero_vec(average_wave_function_up,ist.n_sites*ist.n_up);
   czero_vec(normalized_wave_function_up,ist.n_sites*ist.n_up);

   czero_vec(density_matrix_up,ist.n_sites_sq);
   czero_vec(overall_density_matrix_up,ist.n_sites_sq);

   /*Obtain Population Control Rescaling Factor*/
   for (i = 0; i < 10 ; i++) {
    weight_factor *= weight_rescaling[i];
   }

   for (walker = 0; walker < ist.n_walkers ; walker++ ) {

     /*Check That Walker Weights Are Non-Zero*/
     if ( Cabs(weights[walker]) != 0 ) {

        /*Accumulate Weight in Denominator*/
        current_weight = RCmul(weight_factor, weights[walker]);
        weight_denominator = Cadd(weight_denominator, current_weight);

        /*First Compute the Green's Function*/
        compute_density_matrix_fermions_up(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],density_matrix_up,ist);

        /*Obtain Kinetic and Potential Energies*/
        local_ke = local_pe = Zero; 
        for ( i=0; i<ist.n_sites; i++) {
         //local_pe = Cadd(local_pe, RCmul(.5, density_matrix_up[i*ist.n_sites+i])); 
         if ( ist.n_sites_one > 2 || ist.n_sites_two > 2 ) {
           for (j=0; j<number_neighbors[i]; j++) {
              local_ke = Csub(local_ke, density_matrix_up[i*ist.n_sites+neighbors[4*i+j]]); 

              if ( j<(int)(number_neighbors[i]/2.0) ) {
                density_prod = Cmul( density_matrix_up[i*ist.n_sites+i], density_matrix_up[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]]); 
                local_pe = Cadd(local_pe, density_prod);
              }

           }
         }
         else {
            for (j=0; j<number_neighbors[i]; j++) {
              local_ke =Csub(local_ke, RCmul(.5, density_matrix_up[i*ist.n_sites+neighbors[4*i+j]])); 

              if ( j<(int)(number_neighbors[i]/2.0) ) {
                density_prod = Cmul( density_matrix_up[i*ist.n_sites+i], density_matrix_up[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]]); 
                local_pe = Cadd(local_pe, density_prod);
               }
           }
         }

          /*Construct Density Matrices*/
          for ( j=0; j<ist.n_sites; j++) {
            overall_density_matrix_up[i*ist.n_sites+j] = Cadd( overall_density_matrix_up[i*ist.n_sites+j], Cmul(density_matrix_up[i*ist.n_sites+j], current_weight));
          }
        }

        local_ke = RCmul(cns.t_fermions, Cmul(local_ke, current_weight));
        local_pe = RCmul(cns.U_fermions, Cmul(local_pe, current_weight));

        /*Add to Total*/
        kinetic_energy = Cadd(local_ke, kinetic_energy);
        potential_energy = Cadd(local_pe, potential_energy);
        energy = Cadd(energy, Cadd(local_ke, local_pe));

        /*Obtain Average Wavefunction As Well*/
        for (i=0; i<ist.n_sites; i++) {
         for (j=0; j<ist.n_up; j++) {
           average_wave_function_up[i*ist.n_up+j] = Cadd(average_wave_function_up[i*ist.n_up+j], Cmul(current_weight, wf_up[walker*ist.n_sites_n_up+i*ist.n_up+j]));
         }
       }

      } /*Non-Zero Walker Weights*/
   }  /*Walkers*/

   /*Divide By Total Weight*/
   kinetic_energy = Cdiv(kinetic_energy, weight_denominator);
   potential_energy = Cdiv(potential_energy, weight_denominator);
   energy = Cdiv(energy, weight_denominator);

   *ke = kinetic_energy;
   *pe = potential_energy;
   *te = energy;
   *tw = weight_denominator;

  /*Store the Energies for Further Averaging*/
  kinetic_energy_store[step] = kinetic_energy;
  potential_energy_store[step] = potential_energy;
  total_energy_store[step] = energy;
  total_weights[step] = weight_denominator;

  /*Obtain and Print Average WFs*/
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_up; j++) {
      average_wave_function_up[i*ist.n_up+j] = Cdiv(average_wave_function_up[i*ist.n_up+j], weight_denominator);
    }
  }
  get_normalized_wf_up(normalized_wave_function_up,average_wave_function_up,ist);

free(density_matrix_up);
free(overall_density_matrix_up);
free(average_wave_function_up);
return;
}

/***********************************************************************************************/

void compute_energy_extended_down(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,int *neighbors,int *number_neighbors,MKL_Complex16 *potential_energy_store, MKL_Complex16 *kinetic_energy_store, MKL_Complex16 *total_energy_store,MKL_Complex16 *total_weights,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *ke,MKL_Complex16 *te,MKL_Complex16 *tw,MKL_Complex16 *normalized_wave_function_down) {

   /*Computes the Energy of the System*/
   int walker, i, j;
   double weight_factor = 1.0;
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;
   MKL_Complex16 *average_wave_function_down;
   MKL_Complex16 *density_matrix_down;
   MKL_Complex16 *overall_density_matrix_down;
   MKL_Complex16 local_ke, local_pe;
   MKL_Complex16 density_prod;
   MKL_Complex16 energy = Zero, kinetic_energy = Zero, potential_energy = Zero;
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator = Zero;

   average_wave_function_down=(MKL_Complex16 *)calloc(ist.n_sites*ist.n_down,sizeof(MKL_Complex16));

   density_matrix_down=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));
   overall_density_matrix_down=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));

   czero_vec(average_wave_function_down,ist.n_sites_n_down);
   czero_vec(normalized_wave_function_down,ist.n_sites_n_down);

   czero_vec(density_matrix_down,ist.n_sites_sq);
   czero_vec(overall_density_matrix_down,ist.n_sites_sq);

   /*Obtain Population Control Rescaling Factor*/
   for (i = 0; i < 10 ; i++) {
    weight_factor *= weight_rescaling[i];
   }

   for (walker = 0; walker < ist.n_walkers ; walker++ ) {

     /*Check That Walker Weights Are Non-Zero*/
     if ( Cabs(weights[walker]) != 0 ) {

        /*Accumulate Weight in Denominator*/
        current_weight = RCmul(weight_factor, weights[walker]);
        weight_denominator = Cadd(weight_denominator, current_weight);

        /*First Compute the Green's Function*/
        compute_density_matrix_fermions_down(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],density_matrix_down,ist);

        /*Obtain Kinetic and Potential Energies*/
        local_ke = local_pe = Zero;
        for ( i=0; i<ist.n_sites; i++) {
         if ( ist.n_sites_one > 2 || ist.n_sites_two > 2 ) {
           for (j=0; j<number_neighbors[i]; j++) {
              local_ke = Csub(local_ke, density_matrix_down[i*ist.n_sites+neighbors[4*i+j]]);

              if ( j<(int)(number_neighbors[i]/2.0) ) {
                density_prod = Cmul( density_matrix_down[i*ist.n_sites+i], density_matrix_down[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]]);
                local_pe = Cadd(local_pe, density_prod);
              }

           }
         }
         else {
            for (j=0; j<number_neighbors[i]; j++) {
              local_ke =Csub(local_ke, RCmul(.5, density_matrix_down[i*ist.n_sites+neighbors[4*i+j]]));

              if ( j<(int)(number_neighbors[i]/2.0) ) {
                density_prod = Cmul( density_matrix_down[i*ist.n_sites+i], density_matrix_down[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]]);
                local_pe = Cadd(local_pe, density_prod);
               }
           }
         }

          /*Construct Density Matrices*/
          for ( j=0; j<ist.n_sites; j++) {
            overall_density_matrix_down[i*ist.n_sites+j] = Cadd( overall_density_matrix_down[i*ist.n_sites+j], Cmul(density_matrix_down[i*ist.n_sites+j], current_weight));
          }
        }

        local_ke = RCmul(cns.t_fermions, Cmul(local_ke, current_weight));
        local_pe = RCmul(cns.V_fermions, Cmul(local_pe, current_weight));

        /*Add to Total*/
        kinetic_energy = Cadd(local_ke, kinetic_energy);
        potential_energy = Cadd(local_pe, potential_energy);
        energy = Cadd(energy, Cadd(local_ke, local_pe));

        /*Obtain Average Wavefunction As Well*/
        for (i=0; i<ist.n_sites; i++) {
         for (j=0; j<ist.n_down; j++) {
           average_wave_function_down[i*ist.n_down+j] = Cadd(average_wave_function_down[i*ist.n_down+j], Cmul(current_weight, wf_down[walker*ist.n_sites_n_down+i*ist.n_down+j]));
          }
       }

      } /*Non-Zero Walker Weights*/
   }  /*Walkers*/

   /*Divide By Total Weight*/
   kinetic_energy = Cdiv(kinetic_energy, weight_denominator);
   potential_energy = Cdiv(potential_energy, weight_denominator);
   energy = Cdiv(energy, weight_denominator);

   *ke = kinetic_energy;
   *pe = potential_energy;
   *te = energy;
   *tw = weight_denominator;

  /*Store the Energies for Further Averaging*/
  kinetic_energy_store[step] = kinetic_energy;
  potential_energy_store[step] = potential_energy;
  total_energy_store[step] = energy;
  total_weights[step] = weight_denominator;

  /*Obtain and Print Average WFs*/
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_down; j++) {
      average_wave_function_down[i*ist.n_down+j] = Cdiv(average_wave_function_down[i*ist.n_down+j], weight_denominator);
    }
  }
  get_normalized_wf_down(normalized_wave_function_down,average_wave_function_down,ist);

free(density_matrix_down);
free(overall_density_matrix_down);
free(average_wave_function_down);
return;
}

/***********************************************************************************************/

MKL_Complex16 compute_local_energy_extended(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,int *neighbors,int *number_neighbors,int_st ist,cns_st cns) {

   int i, j; 
   MKL_Complex16 local_ke, local_pe_total, local_pe_extended, local_pe_Hubbard;
   MKL_Complex16 density_prod;  
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;  
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0; 
   MKL_Complex16 local_energy; 
   MKL_Complex16 *density_matrix_up, *density_matrix_down; 

   density_matrix_up=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 
   density_matrix_down=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 

   /*First Compute the Green's Function*/
   compute_density_matrix_fermions(wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_inverse_up,overlap_inverse_down,density_matrix_up,density_matrix_down,ist); 

   local_ke = local_pe_total = local_pe_extended = local_pe_Hubbard = Zero;
   for ( i=0; i<ist.n_sites; i++) {
     density_prod = Cmul( density_matrix_up[i*ist.n_sites+i], density_matrix_down[i*ist.n_sites+i]);
     local_pe_Hubbard = Cadd(local_pe_Hubbard, density_prod);
     if ( ist.n_sites_one > 2 || ist.n_sites_two > 2 ) {
       for (j=0; j<number_neighbors[i]; j++) {
          local_ke = Csub(local_ke, Cadd(density_matrix_up[i*ist.n_sites+neighbors[4*i+j]],density_matrix_down[i*ist.n_sites+neighbors[4*i+j]]));

          if ( j < (int)(number_neighbors[i]/2.0) ) { 
            density_prod = Cmul(Cadd(density_matrix_up[i*ist.n_sites+i],density_matrix_down[i*ist.n_sites+i]), Cadd(density_matrix_up[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]], density_matrix_down[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]])); 
            local_pe_extended = Cadd(local_pe_extended, density_prod);  
          }
        }
      }
      else {
        for (j=0; j<number_neighbors[i]; j++) {
          local_ke =Csub(local_ke, RCmul(.5, Cadd(density_matrix_up[i*ist.n_sites+neighbors[4*i+j]],density_matrix_down[i*ist.n_sites+neighbors[4*i+j]])));
          
          if ( j<(int)(number_neighbors[i]/2.0) ) {
            density_prod = Cmul(Cadd(density_matrix_up[i*ist.n_sites+i],density_matrix_down[i*ist.n_sites+i]), Cadd(density_matrix_up[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]], density_matrix_down[neighbors[4*i+j]*ist.n_sites+neighbors[4*i+j]]));
            local_pe_extended = Cadd(local_pe_extended, density_prod); 
          }
         }
       }
    }

    local_ke = RCmul(cns.t_fermions, local_ke); 
    local_pe_Hubbard = RCmul(cns.U_fermions, local_pe_Hubbard); 
    local_pe_extended = RCmul(cns.V_fermions, local_pe_extended); 
    local_pe_total = Cadd(local_pe_Hubbard, local_pe_extended); 

    local_energy = Cadd(local_ke, local_pe_total); 

free(density_matrix_up); 
free(density_matrix_down); 
return(local_energy);  
}

/***********************************************************************************************/

MKL_Complex16 compute_shift_extended(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,int_st ist,cns_st cns) {

   int i; 
   MKL_Complex16 density_sum;
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;
   MKL_Complex16 *greens_function_up, *greens_function_down;

   greens_function_up=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));
   greens_function_down=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));

   /*First Compute the Green's Function*/
   compute_density_matrix_fermions(wf_up,wf_down,trial_wf_up,trial_wf_down,overlap_inverse_up,overlap_inverse_down,greens_function_up,greens_function_down,ist);

   density_sum = Zero;  
   for ( i=0; i<ist.n_sites; i++) {
     density_sum = Cadd( greens_function_up[i*ist.n_sites+i], greens_function_down[i*ist.n_sites+i]);
   }

   density_sum = Cmul(density_sum, cns.propagation_exponent_fermions); 
   density_sum = RCmul(-1.0, density_sum); 

free(greens_function_up);
free(greens_function_down);
return(density_sum);
}

/***********************************************************************************************/

void compute_averaged_energy_extended(MKL_Complex16 *potential_energy,MKL_Complex16 *kinetic_energy,MKL_Complex16 *total_energy,MKL_Complex16 *weights,int_st ist,char *str) {
     /*Compute the Averages of the Averages*/

     int i;
     double normalization_2 = 1.0/((double) sqrt(ist.n_steps_production));
     MKL_Complex16 square;
     MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
     MKL_Complex16 normalization;
     MKL_Complex16 average_potential_energy, average_kinetic_energy, average_total_energy;
     MKL_Complex16 std_dev_potential, std_dev_kinetic, std_dev_total;
     FILE *pf = fopen(str, "a+");

     /*Obtain Averages*/
     normalization = average_potential_energy = average_kinetic_energy = average_total_energy = Zero;
     for (i=ist.n_steps_equilibration; i<ist.n_steps; i+=ist.n_steps_energy) {
         average_potential_energy = Cadd(average_potential_energy, Cmul(potential_energy[i], weights[i]));
         average_kinetic_energy = Cadd(average_kinetic_energy, Cmul(kinetic_energy[i], weights[i]));
         average_total_energy = Cadd(average_total_energy, Cmul(total_energy[i], weights[i]));
         normalization = Cadd(normalization, weights[i]);
     }
     average_potential_energy = Cdiv(average_potential_energy, normalization);
     average_kinetic_energy = Cdiv(average_kinetic_energy, normalization);
     average_total_energy = Cdiv(average_total_energy, normalization);
    
     /*Then Obtain Standard Deviations*/
     std_dev_potential = std_dev_kinetic = std_dev_total = Zero;
     for (i=ist.n_steps_equilibration; i<ist.n_steps; i+=ist.n_steps_energy) {
        square = Csub(average_potential_energy, potential_energy[i]);
        std_dev_potential = Cadd(std_dev_potential, Cmul(weights[i], Cmul(square, square)));
        square = Cmul(square, square);

        square = Csub(average_kinetic_energy, kinetic_energy[i]);
        std_dev_kinetic = Cadd(std_dev_kinetic, Cmul(weights[i], Cmul(square, square)));

        square = Csub(average_total_energy, total_energy[i]);
        std_dev_total = Cadd(std_dev_total, Cmul(weights[i], Cmul(square, square)));
    }
    std_dev_potential = Cdiv(std_dev_potential, normalization);
    std_dev_potential = RCmul(normalization_2, std_dev_potential);

    std_dev_kinetic = Cdiv(std_dev_kinetic, normalization);
    std_dev_kinetic = RCmul(normalization_2, std_dev_kinetic);

    std_dev_total = Cdiv(std_dev_total, normalization);
    std_dev_total = RCmul(normalization_2, std_dev_total);
    

    /*File Printing*/
    fprintf(pf, "%f\t %f\t\t %f\t %f\t\t %f\t %f\n", average_potential_energy.real, std_dev_potential.real, average_kinetic_energy.real, std_dev_kinetic.real, average_total_energy.real, std_dev_total.real); fflush(pf);

return;
}

/*******************************************************************************************************/

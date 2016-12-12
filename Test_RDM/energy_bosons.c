#include "afqmc.h"

/***********************************************************************************************/

void compute_energy_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *weights,double *kinetic_matrix_bosons,double *potential_matrix_bosons,int *neighbors,int *number_neighbors,MKL_Complex16 *potential_energy_store, MKL_Complex16 *kinetic_energy_store, MKL_Complex16 *total_energy_store,MKL_Complex16 *total_weights,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *ke,MKL_Complex16 *te,MKL_Complex16 *tw) {

   /*Computes the Energy of the System*/
   int walker, i, j;
   double weight_factor = 1.0;
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;
   MKL_Complex16 *average_wave_function_bosons; 
   MKL_Complex16 local_ke, local_pe, local_singlebody_pe, local_twobody_pe;
   MKL_Complex16 normalization_one_body, normalization_two_body;  
   MKL_Complex16 energy = Zero, kinetic_energy = Zero, potential_energy = Zero;
   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator = Zero;
   /*FILE *pf = fopen(str, "a+");
   FILE *pf3 = fopen("wave_function.dat", "a+");
   FILE *pf2 = fopen("errors.dat", "a+");
   FILE *pf4 = fopen("checkingnpotentials.dat", "a+"); 
   */

   average_wave_function_bosons=(MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
   czero_vec(average_wave_function_bosons,ist.n_sites);

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

        /*Compute Permanent*/
        normalization_one_body = normalization_one_body_bosons(trial_wf_bosons,&wf_bosons[walker*ist.n_sites],ist); 
        normalization_two_body = normalization_two_body_bosons(trial_wf_bosons,&wf_bosons[walker*ist.n_sites],ist); 

        /*Obtain Kinetic and Potential Energies*/
        local_ke = compute_ke_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,kinetic_matrix_bosons,ist); 
        local_ke = Cmul(normalization_one_body, local_ke); 

        local_singlebody_pe = compute_pe_singlebody_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,potential_matrix_bosons,ist); 
        local_singlebody_pe = Cmul(normalization_one_body, local_singlebody_pe); 

        local_twobody_pe = Zero; 
        for (j=0; j<ist.n_sites; j++) {
          local_twobody_pe = Cadd(local_twobody_pe, compute_pe_twobody_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,j,ist)); 
        }
        local_twobody_pe = Cmul(normalization_two_body, local_twobody_pe);
 
        /*Normalize By Wavefunction*/
        local_pe = Cadd(local_singlebody_pe, local_twobody_pe); 
        local_ke = RCmul(cns.t_bosons, Cmul(local_ke, current_weight));
        local_pe = RCmul(.5 * cns.U_bosons, Cmul(local_pe, current_weight));

        /*Add to Total*/
        kinetic_energy = Cadd(local_ke, kinetic_energy);
        potential_energy = Cadd(local_pe, potential_energy);
        energy = Cadd(energy, Cadd(local_ke, local_pe));

        /*Obtain Average Wavefunction As Well*/
        for (i=0; i<ist.n_sites; i++) {
         average_wave_function_bosons[i] = Cadd(average_wave_function_bosons[i], Cmul(current_weight, wf_bosons[walker*ist.n_sites+i]));
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
  /*fprintf(pf3, "%d\t\t", step);
  for (i=0; i<ist.n_sites; i++) {
    average_wave_function_bosons[i] = Cdiv(average_wave_function_bosons[i], weight_denominator);
    fprintf(pf3, "%f+%fi\t", average_wave_function_bosons[i].real,average_wave_function_bosons[i].imag);
 } fprintf(pf3, "\n\n");
  */

//fclose(pf);
//fclose(pf2);
//fclose(pf3);
//fclose(pf4); 
return;
}

/********************************************************************************************/

void compute_averaged_energy_bosons(MKL_Complex16 *potential_energy,MKL_Complex16 *kinetic_energy,MKL_Complex16 *total_energy,MKL_Complex16 *weights,int_st ist,char *str) {
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

/***********************************************************************************************/

MKL_Complex16 compute_shift_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,double *potential_matrix_bosons,int_st ist,cns_st cns,int site) {

   /*Compute the Optimal Shift for Bosons to Use in the Propagator*/

   int i;
   MKL_Complex16 sqrt_U; 
   MKL_Complex16 shift; 
   MKL_Complex16 normalization_one_body; 
   MKL_Complex16 *store_product, *store_product_2;

   store_product = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
   store_product_2 = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));

   /*Get Normalization*/
   normalization_one_body = normalization_one_body_bosons(trial_wf_bosons,wf_bosons,ist);    

    /*Conjugate Trial WF*/
    for (i=0; i<ist.n_sites; i++) {
      store_product_2[i] = conjugate(trial_wf_bosons[i]);
    }

    /*Find A Psi*/
    dmat_cmat(potential_matrix_bosons,wf_bosons,store_product,ist.n_sites,ist.n_sites,1);

    /*Zero Elements That Are Not On Particular Site of Interest*/
    for (i=0; i<ist.n_sites; i++) {
      if ( i!= site ) {
        store_product[i].real = store_product[i].imag = 0.0; 
      }
    } 
    
    /*Find Full Expectation Value*/
    cmat_cmat(store_product_2,store_product,&shift,1,ist.n_sites,1);

    /*Multiply By U Which Is In Shift*/ 
    if ( cns.U_bosons > 0 ) {
        sqrt_U.imag = sqrt(cns.U_bosons); 
        sqrt_U.real = 0.0; 
    }
    else {
        sqrt_U.real = sqrt(fabs(cns.U_bosons)); 
        sqrt_U.imag = 0.0; 
    } 

    shift = Cmul(shift, sqrt_U); 

free(store_product);
free(store_product_2);
return(shift); 
}

/**********************************************************************************************/

MKL_Complex16 normalization_one_body_bosons(MKL_Complex16 *trial_wf_bosons, MKL_Complex16 *wf_bosons,int_st ist) {

   /*Compute the Normalization for Calculating Single Particle Averages*/
   MKL_Complex16 N; N.real = ist.n_bosons; N.imag = 0.0; 
   MKL_Complex16 normalization; 
   MKL_Complex16 permanent; 

   permanent = perm_bosons(trial_wf_bosons,wf_bosons,ist); 

   normalization = Cdiv(N, permanent); 

return(normalization); 
}

/**********************************************************************************************/

MKL_Complex16 normalization_two_body_bosons(MKL_Complex16 *trial_wf_bosons, MKL_Complex16 *wf_bosons,int_st ist) {

   /*Compute the Normalization for Calculating Single Particle Averages*/
   MKL_Complex16 N; N.real = ist.n_bosons *  (ist.n_bosons-1); /*N.r = 1.0; */  N.imag = 0.0;
   MKL_Complex16 normalization;
   MKL_Complex16 permanent;

   if ( ist.n_bosons > 1 ) {
     permanent = perm_bosons(trial_wf_bosons,wf_bosons,ist);
     permanent = Cmul(permanent, permanent); 

     normalization = Cdiv(N, permanent);
   }
   else {
     normalization.real = 1.0; 
     normalization.imag = 0.0; 
   }  

return(normalization);
}

/********************************************************************************************/

MKL_Complex16 compute_ke_bosons(MKL_Complex16 *wf_bosons, MKL_Complex16 *trial_wf_bosons,double *kinetic_energy,int_st ist) {

    /*Compute the Kinetic Energy of Bosons*/ 

    int i; 
    MKL_Complex16 ke;  
    MKL_Complex16 *store_product, *store_product_2; 

    store_product = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16)); 
    store_product_2 = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16)); 

    /*Conjugate Trial WF*/ 
    for (i=0; i<ist.n_sites; i++) {
      store_product_2[i] = conjugate(trial_wf_bosons[i]); 
    }

    dmat_cmat(kinetic_energy,wf_bosons,store_product,ist.n_sites,ist.n_sites,1); 
    cmat_cmat(store_product_2,store_product,&ke,1,ist.n_sites,1);  

free(store_product); 
free(store_product_2); 
return(ke); 
}

/******************************************************************************************/

MKL_Complex16 compute_pe_singlebody_bosons(MKL_Complex16 *wf_bosons, MKL_Complex16 *trial_wf_bosons,double *potential_energy,int_st ist) {

    /*Compute PE for Bosons Using Identical Orbirtal Representation Evaluation for Two-Body Operators*/

    int i;  
    MKL_Complex16 pe_single; 
    MKL_Complex16 *store_product, *store_product_2; 

    store_product = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
    store_product_2 = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16)); 

    /*Get Conjugate of Trial Wf*/
    for (i=0; i<ist.n_sites; i++) {
      store_product_2[i] = conjugate(trial_wf_bosons[i]);  
    }

    /*Obtain Single Body contribution*/
    dmat_cmat(potential_energy,wf_bosons,store_product,ist.n_sites,ist.n_sites,1);
    cmat_cmat(store_product_2,store_product,&pe_single,1,ist.n_sites,1);

free(store_product);
free(store_product_2); 
return(pe_single);
}

/********************************************************************************************/

MKL_Complex16 compute_pe_twobody_bosons(MKL_Complex16 *wf_bosons, MKL_Complex16 *trial_wf_bosons,int site,int_st ist) {

    /*Compute PE for Bosons Using Identical Orbirtal Representation Evaluation for Two-Body Operators*/
    MKL_Complex16 pe;
    
    /*Obtain Two-body Contribution*/
    pe = Cmul(conjugate(trial_wf_bosons[site]),conjugate(trial_wf_bosons[site]));
    pe = Cmul(pe, Cmul(wf_bosons[site], wf_bosons[site]));

return(pe);
}

/*********************************************************************************************************/

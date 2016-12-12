#include "afqmc.h"

/***********************************************************************************************/

void compute_energy_bose_fermi(MKL_Complex16 *wf_bosons,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *weights,double *kinetic_matrix_bosons,double *potential_matrix_bosons,int *neighbors,int *number_neighbors,MKL_Complex16 *potential_energy_store, MKL_Complex16 *kinetic_energy_store, MKL_Complex16 *coupling_energy_store, MKL_Complex16 *total_energy_store,MKL_Complex16 *total_weights,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *ke,MKL_Complex16 *ce,MKL_Complex16 *te,MKL_Complex16 *tw) {

   /*Computes the Energy of the System*/
   int walker, i, j;
   int n_sites_n_up = ist.n_sites * ist.n_up;
   int n_sites_n_down = ist.n_sites * ist.n_down;
   int n_up_sq = ist.n_up * ist.n_up;
   int n_down_sq = ist.n_down * ist.n_down;
   double weight_factor = 1.0;

   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;

   MKL_Complex16 *average_wave_function_bosons, *average_wave_function_up, *average_wave_function_down; 
   MKL_Complex16 *greens_function_up, *greens_function_down; 

   MKL_Complex16 local_ke_total, local_pe_total; 
   MKL_Complex16 local_ke_bosons, local_pe_bosons, local_singlebody_pe_bosons, local_twobody_pe_bosons;
   MKL_Complex16 local_ke_fermions, local_pe_fermions, average_local_pe_fermions = Zero;
   MKL_Complex16 local_coupling; 

   MKL_Complex16 density_prod;  
   MKL_Complex16 normalization_one_body, normalization_two_body;  

   MKL_Complex16 energy = Zero, kinetic_energy = Zero, potential_energy = Zero, coupling_energy = Zero;
   MKL_Complex16 kinetic_energy_fermions = Zero, kinetic_energy_bosons = Zero; 
   MKL_Complex16 potential_energy_fermions = Zero, potential_energy_bosons = Zero; 

   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator = Zero;

   /*FILE *pf = fopen(str, "a+");
   FILE *pf2 = fopen("errors.dat", "a+");
   FILE *pf3 = fopen("wave_function_bosons.dat", "a+");
   FILE *pf4 = fopen("wave_function_fermions_up.dat", "a+"); 
   FILE *pf5 = fopen("wave_function_fermions_down.dat", "a+"); 
   */

   average_wave_function_bosons=(MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
   average_wave_function_up=(MKL_Complex16 *)calloc(n_sites_n_up,sizeof(MKL_Complex16)); 
   average_wave_function_down=(MKL_Complex16 *)calloc(n_sites_n_down,sizeof(MKL_Complex16)); 

   greens_function_up=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));
   greens_function_down=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));

   czero_vec(average_wave_function_bosons,ist.n_sites);
   czero_vec(average_wave_function_up,n_sites_n_up); 
   czero_vec(average_wave_function_down,n_sites_n_down); 

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

        /*First Get Fermion Observables*/
        compute_density_matrix_fermions(&wf_up[walker*n_sites_n_up],&wf_down[walker*n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*n_up_sq],&overlap_inverse_down[walker*n_down_sq],greens_function_up,greens_function_down,ist);

        /*Obtain Kinetic and Potential Energies*/
        local_ke_fermions = local_pe_fermions = Zero;
        for ( i=0; i<ist.n_sites; i++) {
          density_prod = Cmul( Csub(One, greens_function_up[i*ist.n_sites+i]), Csub(One, greens_function_down[i*ist.n_sites+i]));
          local_pe_fermions = Cadd(local_pe_fermions, density_prod);
          if ( ist.n_sites_one > 2 || ist.n_sites_two > 2 ) {
            for (j=0; j<number_neighbors[i]; j++) {
               local_ke_fermions = Cadd(local_ke_fermions, Cadd(greens_function_up[i*ist.n_sites+neighbors[4*i+j]],greens_function_down[i*ist.n_sites+neighbors[4*i+j]]));
            }
          }
          else {
             for (j=0; j<number_neighbors[i]; j++) {
              local_ke_fermions =Cadd(local_ke_fermions, RCmul(.5, Cadd(greens_function_up[i*ist.n_sites+neighbors[4*i+j]],greens_function_down[i*ist.n_sites+neighbors[4*i+j]])));
            }
          }
        }

        average_local_pe_fermions = Cadd(average_local_pe_fermions, Cmul(local_pe_fermions, current_weight)); 

        local_ke_fermions = RCmul(cns.t_fermions, Cmul(local_ke_fermions, current_weight));
        local_pe_fermions = RCmul(cns.U_fermions, Cmul(local_pe_fermions, current_weight));

        /*Compute Permanent*/
        normalization_one_body = normalization_one_body_bosons(trial_wf_bosons,&wf_bosons[walker*ist.n_sites],ist); 
        normalization_two_body = normalization_two_body_bosons(trial_wf_bosons,&wf_bosons[walker*ist.n_sites],ist); 

        /*Obtain Kinetic and Potential Energies*/
        local_ke_bosons = compute_ke_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,kinetic_matrix_bosons,ist); 
        local_ke_bosons = Cmul(normalization_one_body, local_ke_bosons); 

        local_singlebody_pe_bosons = compute_pe_singlebody_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,potential_matrix_bosons,ist); 
        local_singlebody_pe_bosons = Cmul(normalization_one_body, local_singlebody_pe_bosons); 

        local_twobody_pe_bosons = Zero; 
        for (j=0; j<ist.n_sites; j++) {
          local_twobody_pe_bosons = Cadd(local_twobody_pe_bosons, compute_pe_twobody_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,j,ist)); 
        }
        local_twobody_pe_bosons = Cmul(normalization_two_body, local_twobody_pe_bosons);
 
        /*Find Local Coupling Energy*/
        local_coupling = compute_coupling_bose_fermi(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,potential_matrix_bosons,greens_function_up,greens_function_down,ist); 
        local_coupling = Cmul(normalization_one_body, local_coupling); 
        local_coupling = RCmul(cns.U_bose_fermi, Cmul(local_coupling, current_weight)); 
 
        /*Normalize By Wavefunction*/
        local_pe_bosons = Cadd(local_singlebody_pe_bosons, local_twobody_pe_bosons); 
        local_ke_bosons = RCmul(cns.t_bosons, Cmul(local_ke_bosons, current_weight));
        local_pe_bosons = RCmul(.5 * cns.U_bosons, Cmul(local_pe_bosons, current_weight));

        /*Add to Total*/
        local_pe_total = Cadd(local_pe_fermions, local_pe_bosons); 
        local_ke_total = Cadd(local_ke_fermions, local_ke_bosons); 
       
        /*Get Broken Up KEs and PEs*/ 
        kinetic_energy_bosons = Cadd(kinetic_energy_bosons, local_ke_bosons); 
        kinetic_energy_fermions = Cadd(kinetic_energy_fermions, local_ke_fermions); 

        potential_energy_bosons = Cadd(potential_energy_bosons, local_pe_bosons); 
        potential_energy_fermions = Cadd(potential_energy_fermions, local_pe_fermions); 

        coupling_energy = Cadd( local_coupling, coupling_energy); 
        kinetic_energy = Cadd( local_ke_total, kinetic_energy);
        potential_energy = Cadd( local_pe_total, potential_energy);

        energy = Cadd(energy, local_ke_total); 
        energy = Cadd(energy, local_pe_total); 
        energy = Cadd(energy, local_coupling); 

        /*Obtain Average Wavefunction As Well*/
        for (i=0; i<ist.n_sites; i++) {
          average_wave_function_bosons[i] = Cadd(average_wave_function_bosons[i], Cmul(current_weight, wf_bosons[walker*ist.n_sites+i]));
          for (j=0; j<ist.n_up; j++) {
           average_wave_function_up[i*ist.n_up+j] = Cadd(average_wave_function_up[i*ist.n_up+j], Cmul(current_weight, wf_up[walker*n_sites_n_up+i*ist.n_up+j]));
          }
          for (j=0; j<ist.n_down; j++) {
           average_wave_function_down[i*ist.n_down+j] = Cadd(average_wave_function_down[i*ist.n_down+j], Cmul(current_weight, wf_down[walker*n_sites_n_down+i*ist.n_down+j]));
          }
        }


      } /*Non-Zero Walker Weights*/
   }  /*Walkers*/

   /*Divide By Total Weight*/
   kinetic_energy_fermions = Cdiv(kinetic_energy_fermions, weight_denominator); 
   kinetic_energy_bosons = Cdiv(kinetic_energy_bosons, weight_denominator); 
   
   potential_energy_fermions = Cdiv(potential_energy_fermions, weight_denominator); 
   potential_energy_bosons = Cdiv(potential_energy_bosons, weight_denominator); 

   kinetic_energy = Cdiv(kinetic_energy, weight_denominator);
   potential_energy = Cdiv(potential_energy, weight_denominator);
   coupling_energy = Cdiv(coupling_energy, weight_denominator); 
   energy = Cdiv(energy, weight_denominator);

   //Print Desired Information*******************************************
   *ke = kinetic_energy;
   *pe = potential_energy;
   *te = energy;
   *ce = coupling_energy; 
   *tw = weight_denominator;

   average_local_pe_fermions = Cdiv(average_local_pe_fermions, weight_denominator);  

  /*Store the Energies for Further Averaging*/
  kinetic_energy_store[step] = kinetic_energy;
  potential_energy_store[step] = potential_energy;
  coupling_energy_store[step] = coupling_energy;  
  total_energy_store[step] = energy;
  total_weights[step] = weight_denominator;

  /*Obtain and Print Average WFs*/
  /*fprintf(pf3, "%d\t\t", step);
  for (i=0; i<ist.n_sites; i++) {
    average_wave_function_bosons[i] = Cdiv(average_wave_function_bosons[i], weight_denominator);
    fprintf(pf3, "%f+%fi\t", average_wave_function_bosons[i].real,average_wave_function_bosons[i].imag);
  } fprintf(pf3, "\n\n");
 
  //Average Wave Function Fermions Up  
  fprintf(pf4, "%d\t\t", step);
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_up; j++) {
      average_wave_function_up[i*ist.n_up+j] = Cdiv(average_wave_function_up[i*ist.n_up+j], weight_denominator);
      fprintf(pf4, "%f+%fi\t", average_wave_function_up[i*ist.n_up+j].real,average_wave_function_up[i*ist.n_up+j].imag);
   }
  } fprintf(pf4, "\n\n");

  //Average Wavef Function Down 
  fprintf(pf5, "%d\t\t", step); 
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_down; j++) {
      average_wave_function_down[i*ist.n_down+j] = Cdiv(average_wave_function_down[i*ist.n_down+j], weight_denominator);
      fprintf(pf5, "%f+%fi\t", average_wave_function_down[i*ist.n_down+j].real, average_wave_function_down[i*ist.n_down+j].imag);
    }
  }
  fprintf(pf5, "\n"); fflush(pf5);  
  */

free(average_wave_function_bosons); 
free(average_wave_function_up); 
free(average_wave_function_down); 

free(greens_function_up); 
free(greens_function_down); 

//fclose(pf);
//fclose(pf2);
//fclose(pf3);
//fclose(pf4); 
//fclose(pf5); 
return;
}

/********************************************************************************************/

MKL_Complex16 compute_coupling_bose_fermi(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,double *potential_matrix_bosons,MKL_Complex16 *greens_function_up,MKL_Complex16 *greens_function_down,int_st ist) {

   /*Calculate the Coupling Between the Bosons and Fermions*/
   int i;
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0; 
   MKL_Complex16 coupling;  
   MKL_Complex16 *fermion_density_potential_matrix; 
   MKL_Complex16 *store_product, *store_product_2;
   MKL_Complex16 *density_matrix_up, *density_matrix_down; 
   MKL_Complex16 *One_2, *Zero;

   One_2 = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   One_2[0].real = 1.0; One_2[0].imag = 0.0;
   Zero[0].real = 0.0; Zero[0].imag = 0.0;

   store_product = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
   store_product_2 = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 
   density_matrix_down = (MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 
   fermion_density_potential_matrix = (MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 

   /*Get Conjugate of Trial Wf*/
   for (i=0; i<ist.n_sites; i++) {
     store_product_2[i] = conjugate(trial_wf_bosons[i]);
     density_matrix_up[i*ist.n_sites+i] = Csub(One, greens_function_up[i*ist.n_sites+i]); 
     density_matrix_down[i*ist.n_sites+i] = Csub(One, greens_function_down[i*ist.n_sites+i]);  
     fermion_density_potential_matrix[i*ist.n_sites+i] = RCmul( potential_matrix_bosons[i*ist.n_sites+i], Cadd(density_matrix_up[i*ist.n_sites+i], density_matrix_down[i*ist.n_sites+i])); 
   }

   /*Obtain Single Body contribution*/
   cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_sites,ist.n_sites,1,One_2,fermion_density_potential_matrix,1,wf_bosons,ist.n_sites,Zero,store_product,ist.n_sites);
   //cmat_cmat(fermion_density_potential_matrix,wf_bosons,store_product,ist.n_sites,ist.n_sites,1);
   cmat_cmat(store_product_2,store_product,&coupling,1,ist.n_sites,1);

free(store_product);
free(store_product_2);
free(density_matrix_up); 
free(density_matrix_down);
free(fermion_density_potential_matrix);
free(One_2); 
free(Zero);   
return(coupling);
}

/**************************************************************************************************/

void compute_averaged_energy_bose_fermi(MKL_Complex16 *potential_energy,MKL_Complex16 *kinetic_energy,MKL_Complex16 *coupling_energy,MKL_Complex16 *total_energy,MKL_Complex16 *weights,int_st ist,char *str) {
     /*Compute the Averages of the Averages*/

     int i;
     double normalization_2 = 1.0/((double) sqrt(ist.n_steps_production));
     MKL_Complex16 square;
     MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
     MKL_Complex16 normalization;
     MKL_Complex16 average_potential_energy, average_kinetic_energy, average_total_energy, average_coupling_energy;
     MKL_Complex16 std_dev_potential, std_dev_kinetic, std_dev_total, std_dev_coupling;
     FILE *pf = fopen(str, "a+");

     /*Obtain Averages*/
     normalization = average_potential_energy = average_kinetic_energy = average_total_energy = average_coupling_energy = Zero;
     for (i=ist.n_steps_equilibration; i<ist.n_steps; i+=ist.n_steps_energy) {
         average_potential_energy = Cadd(average_potential_energy, Cmul(potential_energy[i], weights[i]));
         average_kinetic_energy = Cadd(average_kinetic_energy, Cmul(kinetic_energy[i], weights[i]));
         average_coupling_energy = Cadd(average_coupling_energy, Cmul(coupling_energy[i], weights[i]));  
         average_total_energy = Cadd(average_total_energy, Cmul(total_energy[i], weights[i]));
         normalization = Cadd(normalization, weights[i]);
     }
     average_potential_energy = Cdiv(average_potential_energy, normalization);
     average_kinetic_energy = Cdiv(average_kinetic_energy, normalization);
     average_coupling_energy = Cdiv(average_coupling_energy, normalization); 
     average_total_energy = Cdiv(average_total_energy, normalization);
     

     /*Then Obtain Standard Deviations*/
     std_dev_potential = std_dev_kinetic = std_dev_coupling = std_dev_total = Zero;
     for (i=ist.n_steps_equilibration; i<ist.n_steps; i+=ist.n_steps_energy) {
        square = Csub(average_potential_energy, potential_energy[i]);
        std_dev_potential = Cadd(std_dev_potential, Cmul(weights[i], Cmul(square, square)));
        square = Cmul(square, square);

        square = Csub(average_kinetic_energy, kinetic_energy[i]);
        std_dev_kinetic = Cadd(std_dev_kinetic, Cmul(weights[i], Cmul(square, square)));
        
        square = Csub(average_coupling_energy, coupling_energy[i]); 
        std_dev_coupling = Cadd(std_dev_coupling, Cmul(weights[i], Cmul(square, square)));  

        square = Csub(average_total_energy, total_energy[i]);
        std_dev_total = Cadd(std_dev_total, Cmul(weights[i], Cmul(square, square)));
    }
    std_dev_potential = Cdiv(std_dev_potential, normalization);
    std_dev_potential = RCmul(normalization_2, std_dev_potential);

    std_dev_kinetic = Cdiv(std_dev_kinetic, normalization);
    std_dev_kinetic = RCmul(normalization_2, std_dev_kinetic);

    std_dev_coupling = Cdiv(std_dev_coupling, normalization); 
    std_dev_coupling = RCmul(normalization_2, std_dev_coupling); 

    std_dev_total = Cdiv(std_dev_total, normalization);
    std_dev_total = RCmul(normalization_2, std_dev_total);

    /*File Printing*/
    fprintf(pf, "%f\t %f\t %f\t %f\t\t %f\t %f\t\t %f\t %f\n", average_potential_energy.real, std_dev_potential.real, average_kinetic_energy.real, std_dev_kinetic.real, average_coupling_energy.real, std_dev_coupling.real, average_total_energy.real, std_dev_total.real); fflush(pf);

return;
}

/***********************************************************************************************************/


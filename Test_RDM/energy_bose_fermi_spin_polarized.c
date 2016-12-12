#include "afqmc.h"

/***********************************************************************************************/

void compute_energy_bose_fermi_spin_polarized(MKL_Complex16 *wf_bosons,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,double *kinetic_matrix_bosons,double *potential_matrix_bosons,int *neighbors,int *number_neighbors,MKL_Complex16 *potential_energy_store, MKL_Complex16 *kinetic_energy_store, MKL_Complex16 *coupling_energy_store, MKL_Complex16 *total_energy_store,MKL_Complex16 *total_weights,double *weight_rescaling,int_st ist,cns_st cns,int step,MKL_Complex16 *pe,MKL_Complex16 *ke,MKL_Complex16 *ce,MKL_Complex16 *te,MKL_Complex16 *tw) {

   /*Computes the Energy of the System*/
   int walker, i, j;
   int n_sites_n_up = ist.n_sites * ist.n_up;
   int n_up_sq = ist.n_up * ist.n_up;
   double weight_factor = 1.0;

   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;

   MKL_Complex16 *average_wave_function_bosons, *average_wave_function_up; 
   MKL_Complex16 *greens_function_up; 
   MKL_Complex16 *average_density_matrix_up; 

   MKL_Complex16 local_ke_total, local_pe_total; 
   MKL_Complex16 local_ke_bosons, local_pe_bosons, local_singlebody_pe_bosons, local_twobody_pe_bosons;
   MKL_Complex16 local_ke_fermions; 
   MKL_Complex16 local_coupling; 

   MKL_Complex16 normalization_one_body, normalization_two_body;  

   MKL_Complex16 energy = Zero, kinetic_energy = Zero, potential_energy = Zero, coupling_energy = Zero;
   MKL_Complex16 kinetic_energy_fermions = Zero, kinetic_energy_bosons = Zero; 
   MKL_Complex16 potential_energy_bosons = Zero; 

   MKL_Complex16 current_weight;
   MKL_Complex16 weight_denominator = Zero;

   /*FILE *pf = fopen(str, "a+");
   FILE *pf2 = fopen("errors.dat", "a+");
   FILE *pf3 = fopen("wave_function_bosons.dat", "a+");
   FILE *pf4 = fopen("wave_function_fermions_up.dat", "a+"); 
   FILE *pf5 = fopen("density_matrix_fermions.dat", "a+");  
  */ 

   average_wave_function_bosons=(MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
   average_wave_function_up=(MKL_Complex16 *)calloc(n_sites_n_up,sizeof(MKL_Complex16)); 

   average_density_matrix_up=(MKL_Complex16 *)calloc(ist.n_sites*ist.n_sites,sizeof(MKL_Complex16)); 

   greens_function_up=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));

   czero_vec(average_wave_function_bosons,ist.n_sites);
   czero_vec(average_wave_function_up,n_sites_n_up); 

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
        compute_greens_function_fermions_spin_polarized(&wf_up[walker*n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*n_up_sq],greens_function_up,ist);

       for (i=0; i<ist.n_sites; i++) {
         for (j=0; j<ist.n_sites; j++) {
            average_density_matrix_up[i*ist.n_sites+j] = Cadd(average_density_matrix_up[i*ist.n_sites+j], Cmul(current_weight,greens_function_up[i*ist.n_sites+j])); 
          }
        }     


        /*Obtain Kinetic and Potential Energies*/
        local_ke_fermions = Zero;
        for ( i=0; i<ist.n_sites; i++) {
          if ( ist.n_sites_one > 2 || ist.n_sites_two > 2 ) {
            for (j=0; j<number_neighbors[i]; j++) {
               local_ke_fermions = Cadd(local_ke_fermions, greens_function_up[i*ist.n_sites+neighbors[4*i+j]]); 
            }
          }
          else {
             for (j=0; j<number_neighbors[i]; j++) {
              local_ke_fermions =Cadd(local_ke_fermions, RCmul(.5, greens_function_up[i*ist.n_sites+neighbors[4*i+j]])); 
            }
          }
        }

        local_ke_fermions = RCmul(cns.t_fermions, Cmul(local_ke_fermions, current_weight));

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
        local_coupling = compute_coupling_bose_fermi_spin_polarized(walker,&wf_bosons[walker*ist.n_sites],trial_wf_bosons,potential_matrix_bosons,greens_function_up,ist); 
        local_coupling = Cmul(normalization_one_body, local_coupling);
        local_coupling = RCmul(cns.U_bose_fermi, Cmul(local_coupling, current_weight)); 

        /*Normalize By Wavefunction*/
        local_pe_bosons = Cadd(local_singlebody_pe_bosons, local_twobody_pe_bosons); 
        local_ke_bosons = RCmul(cns.t_bosons, Cmul(local_ke_bosons, current_weight));
        local_pe_bosons = RCmul(.5 * cns.U_bosons, Cmul(local_pe_bosons, current_weight));

        /*Add to Total*/
        local_pe_total = local_pe_bosons; 
        local_ke_total = Cadd(local_ke_fermions, local_ke_bosons); 
       
        /*Get Broken Up KEs and PEs*/ 
        kinetic_energy_bosons = Cadd(kinetic_energy_bosons, local_ke_bosons); 
        kinetic_energy_fermions = Cadd(kinetic_energy_fermions, local_ke_fermions); 
        potential_energy_bosons = Cadd(potential_energy_bosons, local_pe_bosons); 

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
        }


      } /*Non-Zero Walker Weights*/
   }  /*Walkers*/

   /*Divide By Total Weight*/
   kinetic_energy_fermions = Cdiv(kinetic_energy_fermions, weight_denominator); 
   kinetic_energy_bosons = Cdiv(kinetic_energy_bosons, weight_denominator); 
   
   potential_energy_bosons = Cdiv(potential_energy_bosons, weight_denominator); 

   kinetic_energy = Cdiv(kinetic_energy, weight_denominator);
   potential_energy = Cdiv(potential_energy, weight_denominator);
   coupling_energy = Cdiv(coupling_energy, weight_denominator); 
   energy = Cdiv(energy, weight_denominator);

   /*Send Values Out for Printing*/
   *ke = kinetic_energy;
   *pe = potential_energy;
   *te = energy;
   *ce = coupling_energy;
   *tw = weight_denominator;

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
  */ 

  /*Average Wave Function Fermions Up*/  
  /*fprintf(pf4, "%d\t\t", step);
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_up; j++) {
      average_wave_function_up[i*ist.n_up+j] = Cdiv(average_wave_function_up[i*ist.n_up+j], weight_denominator);
      fprintf(pf4, "%f+%fi\t", average_wave_function_up[i*ist.n_up+j].real,average_wave_function_up[i*ist.n_up+j].imag);
   }
  } fprintf(pf4, "\n\n");
  */

  /*Print Out Average Fermion Density Matrix*/
 /* fprintf(pf5, "%d\t\t", step); 
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_sites; j++) {
      average_density_matrix_up[i*ist.n_sites+j] = Cdiv(average_density_matrix_up[i*ist.n_sites+j], weight_denominator);
      fprintf(pf5, "%f+%fi\t", average_density_matrix_up[i*ist.n_sites+j].real,average_density_matrix_up[i*ist.n_sites+j].imag);
   }
  } fprintf(pf5, "\n\n");
 */ 

free(average_wave_function_bosons); 
free(average_wave_function_up); 
free(average_density_matrix_up);  

free(greens_function_up); 

//fclose(pf);
//fclose(pf2);
//fclose(pf3);
//fclose(pf4);
//fclose(pf5);  
return;
}

/********************************************************************************************/

MKL_Complex16 compute_coupling_bose_fermi_spin_polarized(int walker, MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,double *potential_matrix_bosons,MKL_Complex16 *greens_function_up,int_st ist) {

   /*Calculate the Coupling Between the Bosons and Fermions*/
   int i;
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0; 
   MKL_Complex16 coupling;  
   MKL_Complex16 *fermion_density_potential_matrix; 
   MKL_Complex16 *store_product, *store_product_2;
   MKL_Complex16 *density_matrix_up; 
   MKL_Complex16 *One_2, *Zero;

   One_2 = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   One_2[0].real = 1.0; One_2[0].imag = 0.0;
   Zero[0].real = 0.0; Zero[0].imag = 0.0;

   store_product = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
   store_product_2 = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
   density_matrix_up = (MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 
   fermion_density_potential_matrix = (MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 

   czero_vec(fermion_density_potential_matrix, ist.n_sites_sq); 

   /*Get Conjugate of Trial Wf*/
   for (i=0; i<ist.n_sites; i++) {
     store_product_2[i] = conjugate(trial_wf_bosons[i]);
     density_matrix_up[i*ist.n_sites+i] = Csub(One, greens_function_up[i*ist.n_sites+i]); 
     fermion_density_potential_matrix[i*ist.n_sites+i] = RCmul( potential_matrix_bosons[i*ist.n_sites+i], density_matrix_up[i*ist.n_sites+i]); 
   }

   /*Obtain Single Body contribution*/
   cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_sites,ist.n_sites,1,One_2,fermion_density_potential_matrix,1,wf_bosons,ist.n_sites,Zero,store_product,ist.n_sites);
   //cmat_cmat(fermion_density_potential_matrix,wf_bosons,store_product,ist.n_sites,ist.n_sites,1);
   cmat_cmat(store_product_2,store_product,&coupling,1,ist.n_sites,1);

free(One_2); 
free(Zero); 
free(store_product);
free(store_product_2);
free(density_matrix_up); 
free(fermion_density_potential_matrix);  
return(coupling);
}

/**************************************************************************************************/

void compute_greens_function_fermions_spin_polarized(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *greens_function_up,int_st ist) {

   /*Computes the Spin Up and Down Greens Functions*/
   int i, j;
   MKL_Complex16 *stored_product1_up; 
   MKL_Complex16 *stored_product2_up; 

   stored_product1_up=(MKL_Complex16 *)calloc(ist.n_sites*ist.n_up,sizeof(MKL_Complex16));
   stored_product2_up=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));

   /*Get Up Green's Function First*/
   cmat_transpose_cmat(overlap_inverse_up,trial_wf_up,stored_product1_up,ist.n_up,ist.n_sites,ist.n_up);
   cmat_cmat(wf_up,stored_product1_up,stored_product2_up,ist.n_sites,ist.n_up,ist.n_sites);

   /*Now Get Greens Function By Subtracting From Unity*/
   for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_sites; j++) {
      greens_function_up[i*ist.n_sites+j] = RCmul(-1, stored_product2_up[i*ist.n_sites+j]);
    }
    greens_function_up[i*ist.n_sites+i].real += 1.0;
   }


free(stored_product1_up);
free(stored_product2_up);
}

/***********************************************************************************************************/

MKL_Complex16 compute_shift_bose_fermi_spin_polarized(MKL_Complex16 *wf_bosons,MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_inverse_up,double *potential_matrix_bosons,int_st ist,cns_st cns,int site) {

   /*Compute the Optimal Shift for Bosons to Use in the Propagator*/

   int i; 
   MKL_Complex16 sqrt_U;
   MKL_Complex16 density_bosons, density_fermions; 
   MKL_Complex16 shift;
   MKL_Complex16 normalization_one_body;
   MKL_Complex16 *store_product, *store_product_2;
   MKL_Complex16 *greens_function_up; 

   store_product = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
   store_product_2 = (MKL_Complex16 *)calloc(ist.n_sites,sizeof(MKL_Complex16));
   greens_function_up = (MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 

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
    cmat_cmat(store_product_2,store_product,&density_bosons,1,ist.n_sites,1);

    /*Find Greens Function To Get Fermion Density*/
    compute_greens_function_fermions_spin_polarized(wf_up,trial_wf_up,overlap_inverse_up,greens_function_up,ist); 

    density_fermions.real = density_fermions.imag = 0.0; 
    for ( i=0; i<ist.n_sites; i++) {
       density_fermions = Cadd(greens_function_up[i*ist.n_sites+i], density_fermions); 
    }

    /*Get Full Shift*/
    shift = Cadd(density_fermions, density_bosons); 

    /*Multiply By U Which Is In Shift*/
    if ( cns.U_bose_fermi > 0 ) {
        sqrt_U.real = sqrt(cns.U_bose_fermi) * sqrt(cns.dtau);
        sqrt_U.imag = 0.0;
    }
    else {
        sqrt_U.imag = sqrt(fabs(cns.U_bose_fermi)) * sqrt(cns.dtau);
        sqrt_U.real = 0.0;
    }

    shift = Cmul(shift, sqrt_U);

free(greens_function_up); 
free(store_product);
free(store_product_2);
return(shift);
}

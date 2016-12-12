#include "afqmc.h"

void init(int_st *ist,cns_st *cns) { 

  FILE *pf = fopen("afqmc.par", "r"); 
  FILE *po = fopen("afqmc-parameters.dat", "w+"); 

  /*Scan In Relevant Parameters*/

  fscanf (pf,"%d",&ist->n_sites_one); 
  fscanf (pf,"%d",&ist->n_sites_two); 

  fscanf (pf,"%d",&ist->n_spatial_orbitals);
  fscanf (pf,"%d",&ist->n_determinants_max); 
  fscanf (pf,"%lf",&ist->det_cutoff_energy); 
  fscanf (pf,"%lf",&ist->det_cutoff_phaseless);  

  fscanf (pf,"%d",&ist->n_up); 
  fscanf (pf,"%d",&ist->n_down); 
  fscanf (pf,"%d",&ist->n_bosons);   

  fscanf (pf,"%d",&ist->n_walkers);
  fscanf (pf,"%d",&ist->n_steps_energy); 
  fscanf (pf,"%d",&ist->n_steps_orthogonalize); 
  fscanf (pf,"%d",&ist->n_steps_population); 
  fscanf (pf,"%d",&ist->n_steps_equilibration); 
  fscanf (pf,"%d",&ist->n_steps_start); 
  fscanf (pf,"%d",&ist->n_steps_back_propagation); 

  fscanf (pf,"%d",&ist->flag_trial); 
  fscanf (pf,"%d",&ist->flag_pbc);   
  fscanf (pf,"%d",&ist->flag_cp); 
  fscanf (pf,"%d",&ist->flag_discrete); 
  fscanf (pf,"%d",&ist->flag_phaseless);
  fscanf (pf,"%d",&ist->flag_local_energy);  
  fscanf (pf,"%d",&ist->flag_meanfield); 
  fscanf (pf,"%d",&ist->flag_spin_polarized); 

  fscanf (pf,"%d",&ist->flag_simulation_type); 
  fscanf (pf,"%d",&ist->flag_read_in_trial); 
  fscanf (pf,"%d",&ist->flag_restart); 
  fscanf (pf,"%d",&ist->flag_compute_2RDM); 
  fscanf (pf,"%d",&ist->flag_back_propagation); 

  fscanf (pf,"%lg",&cns->beta);
  fscanf (pf,"%lg",&cns->dtau); 

  fscanf (pf,"%lg",&cns->U_fermions); 
  fscanf (pf,"%lg",&cns->t_fermions);
  fscanf (pf,"%lg",&cns->U_bosons); 
  fscanf (pf,"%lg",&cns->t_bosons); 
  fscanf (pf,"%lg",&cns->U_bose_fermi);
  fscanf (pf,"%lf",&cns->V_fermions);  

  fscanf (pf,"%lg",&cns->trial_density_up);
  fscanf (pf,"%lg",&cns->trial_density_down);   
  fscanf (pf,"%lg",&cns->trial_density_bosons);  

  fscanf (pf,"%lg",&cns->min_walker_weight); 
  fscanf (pf,"%lg",&cns->max_walker_weight); 
  fscanf (pf,"%lg",&cns->population_control_factor);  

  fscanf (pf,"%lg",&cns->energy_cap_constant); 

  /*Steps*/
  ist->n_sites = ist->n_sites_one * ist->n_sites_two; 

  ist->n_spatial_orbitals_sq = ist->n_spatial_orbitals * ist->n_spatial_orbitals; 
  ist->n_spatial_orbitals_third = ist->n_spatial_orbitals_sq * ist->n_spatial_orbitals; 
  ist->n_spatial_orbitals_fourth = ist->n_spatial_orbitals_sq * ist->n_spatial_orbitals_sq; 
  ist->n_spatial_orbitals_n_up = ist->n_spatial_orbitals * ist->n_up; 
  ist->n_spatial_orbitals_n_down = ist->n_spatial_orbitals * ist->n_down;

  ist->n_spin_orbitals = 2 * ist->n_spatial_orbitals; 
  ist->n_spin_orbitals_sq = ist->n_spin_orbitals * ist->n_spin_orbitals; 
  ist->n_spin_orbitals_third = ist->n_spin_orbitals_sq * ist->n_spin_orbitals; 
  ist->n_spin_orbitals_fourth = ist->n_spin_orbitals_third * ist->n_spin_orbitals;

  /*Determine The Number of Energy Determinants*/
  ist->n_determinants_trial_energy = determine_number_dets(ist->det_cutoff_energy,ist->n_determinants_max,ist->n_up+ist->n_down);  
  ist->n_determinants_trial_phaseless = determine_number_dets(ist->det_cutoff_phaseless,ist->n_determinants_max,ist->n_up+ist->n_down); 

  ist->n_determinants_energy_n_spatial_orbitals_n_up = ist->n_spatial_orbitals_n_up * ist->n_determinants_trial_energy; 
  ist->n_determinants_energy_n_spatial_orbitals_n_down = ist->n_spatial_orbitals_n_down * ist->n_determinants_trial_energy;
  ist->n_determinants_energy_n_up_sq = ist->n_determinants_trial_energy * ist->n_up * ist->n_up; 
  ist->n_determinants_energy_n_down_sq = ist->n_determinants_trial_energy * ist->n_down * ist->n_down; 

  ist->n_determinants_phaseless_n_spatial_orbitals_n_up = ist->n_spatial_orbitals_n_up * ist->n_determinants_trial_phaseless;
  ist->n_determinants_phaseless_n_spatial_orbitals_n_down = ist->n_spatial_orbitals_n_down * ist->n_determinants_trial_phaseless;
  ist->n_determinants_phaseless_n_up_sq = ist->n_determinants_trial_phaseless * ist->n_up * ist->n_up;
  ist->n_determinants_phaseless_n_down_sq = ist->n_determinants_trial_phaseless * ist->n_down * ist->n_down;

  ist->n_steps = (int)(cns->beta/cns->dtau); 

  ist->n_sites_sq = ist->n_sites * ist->n_sites;  
  ist->n_sites_third = ist->n_sites_sq * ist->n_sites; 
  ist->n_sites_fourth = ist->n_sites_sq * ist->n_sites_sq;

  ist->n_electrons = ist->n_up + ist->n_down;  
  ist->n_up_sq = ist->n_up * ist->n_up; 
  ist->n_down_sq = ist->n_down * ist->n_down; 
  ist->n_sites_n_up = ist->n_sites * ist->n_up; 
  ist->n_sites_n_down = ist->n_sites * ist->n_down; 
  ist->n_states = 10000; 

  /*Assume the Eigenvectors are Real Unless You Know Otherwise*/
  ist->flag_real_eigenvectors = 1; 
 
  /*Find Steps for Production Runs*/
  ist->n_steps_production = ist->n_steps - ist->n_steps_equilibration; 

  fprintf (po,"Sites One = %d\n", ist->n_sites_one); 
  fprintf (po,"Sites Two = %d\n", ist->n_sites_two); 
  fprintf (po,"Sites = %d\n", ist->n_sites); 
  fprintf (po, "\n"); 

  fprintf (po,"Spatial Orbitals = %d\n", ist->n_spatial_orbitals); 
  fprintf (po,"Spin Orbitals = %d\n", ist->n_spin_orbitals);
  fprintf (po,"N Determinants Max = %d\n", ist->n_determinants_max); 
  fprintf (po,"Det Cutoff Energy = %f\n", ist->det_cutoff_energy); 
  fprintf (po,"N Trial Determinants for Energy = %d\n", ist->n_determinants_trial_energy);  
  fprintf (po,"Det Cutoff Phaseless = %f\n", ist->det_cutoff_phaseless); 
  fprintf (po,"N Trial Determinants for Phaseless = %d\n", ist->n_determinants_trial_phaseless); 
  fprintf (po,"Number Electrons Up = %d\n", ist->n_up); 
  fprintf (po,"Number Electrons Down = %d\n", ist->n_down); 
  fprintf (po,"Number Bosons = %d\n", ist->n_bosons); 

  fprintf (po, "\n");  
  fprintf (po,"Number of Walkers = %d\n", ist->n_walkers); 
  fprintf (po,"Nsteps = %d\n", ist->n_steps);
  fprintf (po,"Nsteps Energy = %d\n", ist->n_steps_energy); 
  fprintf (po,"Nsteps Orthogonalize = %d\n", ist->n_steps_orthogonalize); 
  fprintf (po,"Nsteps Population Control = %d\n", ist->n_steps_population); 
  fprintf (po,"Nsteps Equilibration = %d\n", ist->n_steps_equilibration);
  fprintf (po,"Nsteps Start = %d\n", ist->n_steps_start);  
  fprintf (po,"Nsteps Production = %d\n", ist->n_steps_production);  
  fprintf (po,"Nsteps Back Propagation = %d\n", ist->n_steps_back_propagation); 

  fprintf (po, "\n"); 
  fprintf (po,"Flag CP = %d\n", ist->flag_cp);  
  fprintf (po,"Trial Wavefunction Type = %d\n", ist->flag_trial); 
  fprintf (po,"Periodic Boundary Conditions = %d\n", ist->flag_pbc); 
  fprintf (po,"Discrete Transform = %d\n", ist->flag_discrete); 
  fprintf (po,"Phaseless = %d\n", ist->flag_phaseless);
  fprintf (po,"Local Energy = %d\n", ist->flag_local_energy); 
  fprintf (po,"Mean Field = %d\n", ist->flag_meanfield); 
  fprintf (po,"Spin Polarized = %d\n", ist->flag_spin_polarized);  
  fprintf (po,"Simulation Type = %d\n", ist->flag_simulation_type);  
  fprintf (po,"Flag Real Eigenvectors = %d\n", ist->flag_real_eigenvectors);  
  fprintf (po,"Flag Read in Trial = %d\n", ist->flag_read_in_trial);  
  fprintf (po,"Flag Restart = %d\n", ist->flag_restart); 
  fprintf (po,"Flag Compute 2RDM = %d\n", ist->flag_compute_2RDM); 
  fprintf (po,"Flag Back Propagation = %d\n", ist->flag_back_propagation); 

  /*Reset Values if Spin-Polarized System*/
  if ( ist->flag_spin_polarized == 1 ) {
    ist->n_down = 0;
    cns->U_fermions = 0.0;
    if ( ist->flag_simulation_type == 2 ) {
      ist->flag_simulation_type = 3; 
    } 
  }

  /*Also Ensure Appropriate Type of Transformation If For Bosons*/
  if ( ist->flag_simulation_type != 0 ) {
    ist->flag_discrete = 0; 
  }

  /*If Fermimon, Boson, or Bose-Fermi Simulation Propagation Constants*********************************/
  if ( ist->flag_simulation_type == 0 ) {
    /*Find Constant U Exponent*/
    if ( ist->flag_discrete == 1 ) {
      cns->U_exponent_fermions = exp(-1 * cns->U_fermions * cns->dtau/2.0); 

      /*Find Gamma Given U and Delta Tau*/
      cns->gamma = find_gamma(cns->U_fermions,cns->dtau);

      cns->factor_spin_up_field_down = exp(cns->gamma);
      cns->factor_spin_up_field_up = exp(-1*cns->gamma);  
    }
    else {
      cns->U_exponent_fermions = exp(cns->U_fermions * cns->dtau/2.0); 
      cns->total_trial_density = cns->trial_density_up - cns->trial_density_down;
      cns->total_trial_density *= 1.0/((double) ist->n_sites);  
      cns->mean_field_exponential_fermions = exp(-cns->U_fermions * cns->dtau * (cns->total_trial_density)); 
      if ( cns->U_fermions > 0 ) {
        cns->propagation_exponent_fermions.real = sqrt(cns->U_fermions * cns->dtau);  
        cns->propagation_exponent_fermions.imag = 0; 
      }
      else {
        cns->propagation_exponent_fermions.imag = sqrt(fabs(cns->U_fermions)*cns->dtau); 
        cns->propagation_exponent_fermions.real = 0.0; 
      }
    } 
  }
  else if ( ist->flag_simulation_type == 1 ) {
      /*If for Bosons*/

      cns->U_exponent_bosons = 1.0; 
      cns->total_trial_density = cns->trial_density_bosons;
      cns->total_trial_density *= 1.0/((double) ist->n_sites);  
      cns->mean_field_exponential_bosons = exp(-cns->U_bosons * cns->dtau * cns->total_trial_density); 
      if ( cns->U_bosons > 0 ) {
        cns->propagation_exponent_bosons.real = 0.0;
        cns->propagation_exponent_bosons.imag = sqrt(fabs(cns->U_bosons)*cns->dtau);
      }
      else {
        cns->propagation_exponent_bosons.real = sqrt(cns->U_bosons*cns->dtau);
        cns->propagation_exponent_bosons.imag = 0.0; 
      }
      
  }
  else if ( ist->flag_simulation_type == 2 )  { 
    /*If For Mixtures*/

     /*Bosons*/
     cns->U_exponent_bosons = 1.0;
     cns->total_trial_density = cns->trial_density_bosons;
     cns->total_trial_density *= 1.0/((double) ist->n_sites);
     cns->mean_field_exponential_bosons = exp(-cns->U_bosons * cns->dtau * cns->total_trial_density);
  
     cns->total_trial_density = cns->trial_density_up - cns->trial_density_down;
     cns->total_trial_density *= 1.0/((double) ist->n_sites);
     cns->mean_field_exponential_fermions = exp(-cns->U_fermions * cns->dtau * (cns->total_trial_density));
 
     /*Mixture Decoupling*/
     if ( cns->U_bose_fermi > 0 ) {
       cns->propagation_exponent_bose_fermi.real = sqrt(cns->U_bose_fermi * cns->dtau);
       cns->propagation_exponent_bose_fermi.imag = 0.0; 
     }
     else {
       cns->propagation_exponent_bose_fermi.real = 0.0; 
       cns->propagation_exponent_bose_fermi.imag = sqrt(fabs(cns->U_bose_fermi) * cns->dtau);
     }    

     /*Fermion Decoupling*/ 
     cns->U_exponent_fermions = exp( -1.0 * ( 2.0 * cns->U_bose_fermi + cns->U_fermions) * cns->dtau/2.0);
     if ( cns->U_bose_fermi + cns->U_fermions > 0 ) {
       cns->propagation_exponent_fermions.real = sqrt( (cns->U_bose_fermi + cns->U_fermions) * cns->dtau);
       cns->propagation_exponent_fermions.imag = 0;
     }
     else {
       cns->propagation_exponent_fermions.imag = sqrt( fabs(cns->U_bose_fermi + cns->U_fermions) * cns->dtau);
       cns->propagation_exponent_fermions.real = 0.0;
     }

     /*Boson Decoupling*/
     if ( cns->U_bosons + cns->U_bose_fermi  > 0 ) {
        cns->propagation_exponent_bosons.real = 0.0; 
        cns->propagation_exponent_bosons.imag = sqrt( (cns->U_bosons + cns->U_bose_fermi) * cns->dtau ); 
     }
     else {
        cns->propagation_exponent_bosons.imag = 0.0; 
        cns->propagation_exponent_bosons.real = sqrt( fabs(cns->U_bosons + cns->U_bose_fermi) * cns->dtau ); 
     }

  }
  else if ( ist->flag_simulation_type == 3 ) {
  
    /*If For Mixtures*/

     /*Bosons*/
      cns->total_trial_density = cns->trial_density_bosons;
      cns->total_trial_density *= 1.0/((double) ist->n_sites);
      cns->mean_field_exponential_bosons = exp(-cns->U_bose_fermi * cns->dtau * cns->total_trial_density);
     if ( cns->U_bose_fermi + cns->U_bosons  > 0 ) {
       cns->propagation_exponent_bosons.real = 0.0;
       cns->propagation_exponent_bosons.imag = sqrt((cns->U_bose_fermi+cns->U_bosons)*cns->dtau);
     }
     else {
       cns->propagation_exponent_bosons.real = sqrt(fabs(cns->U_bose_fermi+cns->U_bosons)*cns->dtau);
       cns->propagation_exponent_bosons.imag = 0.0;
     }

     /*For Bose Fermi Term*/
     if ( cns->U_bose_fermi > 0 ) {
        cns->propagation_exponent_bose_fermi.real = sqrt(cns->U_bose_fermi*cns->dtau); 
        cns->propagation_exponent_bose_fermi.imag = 0.0; 
     }
     else {
        cns->propagation_exponent_bose_fermi.real = 0.0; 
        cns->propagation_exponent_bose_fermi.imag = sqrt(fabs(cns->U_bose_fermi)*cns->dtau); 
     }  

     /*Fermions*/
     cns->U_exponent_fermions = exp(-1.0 * cns->U_bose_fermi * cns->dtau/2.0);
  }
  else if ( ist->flag_simulation_type == 7 ) {
     /*For the Extended Hubbard Model*/
     cns->U_exponent_fermions = exp(cns->U_fermions * cns->dtau/2.0);
 
     /*Get Propagation Exponents*/
     if ( cns->U_fermions > 0 ) {
       cns->propagation_exponent_fermions.real = sqrt(cns->dtau * cns->U_fermions); 
       cns->propagation_exponent_fermions.imag = 0.0;    
     }
     else {
       cns->propagation_exponent_fermions.real = 0.0; 
       cns->propagation_exponent_fermions.imag = sqrt(cns->dtau * fabs(cns->U_fermions));   
     }    

     if ( cns->V_fermions > 0 ) {
       cns->propagation_exponent_extended.real = sqrt(cns->dtau * cns->V_fermions); 
       cns->propagation_exponent_extended.imag = 0.0;   
     }
     else {
       cns->propagation_exponent_extended.real = 0.0; 
       cns->propagation_exponent_extended.imag = sqrt(cns->dtau * fabs(cns->V_fermions));  
     }

  }

  /*Min and Max Number of Walkers*/
  cns->min_number_walkers = (int)(0.5 * ist->n_walkers);
  cns->max_number_walkers = (int)(8.0 * ist->n_walkers);

  fprintf (po,"Beta = %g\n", cns->beta);
  fprintf (po,"Dtau = %g\n", cns->dtau); 
  fprintf (po, "\n\n"); 

  fprintf (po,"U Fermions = %g\n", cns->U_fermions); 
  fprintf (po,"t Fermions = %g\n", cns->t_fermions);
  fprintf (po,"U Bosons = %g\n", cns->U_bosons); 
  fprintf (po,"t Bosons = %g\n", cns->t_bosons); 
  fprintf (po,"U Bose Fermi = %g\n", cns->U_bose_fermi); 
  fprintf (po,"V Fermions = %g\n", cns->V_fermions); 
  fprintf (po, "\n\n"); 

  if ( ist->flag_simulation_type!=5 && ist->flag_simulation_type!=8 && ist->flag_simulation_type!=9 ) {
    fprintf (po,"trial density up = %g\n", cns->trial_density_up);
    fprintf (po,"trial density down = %g\n", cns->trial_density_down);  
    fprintf (po,"trial density bosons = %g\n", cns->trial_density_bosons);   
    fprintf (po,"trial total density = %g\n", cns->total_trial_density); 
    fprintf (po, "\n\n"); 

    fprintf (po,"gamma = %g\n", cns->gamma);
    fprintf (po,"U exponent fermions = %g\n", cns->U_exponent_fermions);  
    fprintf (po,"U exponent bosons = %g\n", cns->U_exponent_bosons);  
    fprintf (po,"factor_spin_up_field_up = %g\n", cns->factor_spin_up_field_up); 
    fprintf (po,"factor_spin_up_field_down = %g\n", cns->factor_spin_up_field_down);  
    fprintf (po, "\n\n"); 
   
    fprintf (po,"prop exponent fermions = %f+%fi\n", cns->propagation_exponent_fermions.real, cns->propagation_exponent_fermions.imag);
    fprintf (po,"prop exponent bosons = %f+%fi\n", cns->propagation_exponent_bosons.real, cns->propagation_exponent_bosons.imag);
    fprintf (po,"prop exponent bose fermi = %f+%fi\n", cns->propagation_exponent_bose_fermi.real, cns->propagation_exponent_bose_fermi.imag);

   fprintf (po,"mean field exponential fermions = %f\n", cns->mean_field_exponential_fermions);
   fprintf (po,"mean field exponential bosons = %f\n", cns->mean_field_exponential_bosons);
   fprintf (po,"mean field exponential bose fermi = %f\n", cns->mean_field_exponential_bose_fermi);

  }


  fprintf (po,"min walker weight = %g\n", cns->min_walker_weight);
  fprintf (po,"max walker weight = %g\n", cns->max_walker_weight); 
  fprintf (po,"population control factor = %g\n", cns->population_control_factor); 
  fprintf (po,"min number walkers = %d\n", cns->min_number_walkers);   
  fprintf (po,"max number walkers = %d\n", cns->max_number_walkers); 
  fprintf (po,"energy cap constant = %f\n", cns->energy_cap_constant); 
  fprintf (po, "\n\n"); 

  fclose(pf); 
  fclose(po);  

return; 
}

/**************************************************************/

int determine_number_dets(double energy_cutoff,int number_max_dets,int number_electrons) {

   /*Determines the Number of Determinants Below a Certain Threshhold*/
    FILE *pf = fopen("trialwavefunction_readin.par", "r");
    int i, j, number_dets = 0; 
    int *electrons; 
    double det_size; 

    electrons=(int*)calloc(number_electrons,sizeof(int)); 

    for (i=0; i<number_max_dets; i++) {
       fscanf(pf, "%lf", &det_size); 

       for (j=0; j<number_electrons; j++) { 
        fscanf(pf, "%d", &electrons[j]); 
       }
   
       if (fabs(det_size) > energy_cutoff) {
         number_dets++; 
       }
    }

fclose(pf); 
free(electrons); 
return(number_dets); 
}

/*************************************************************/

void init_neighbors(int *neighbors,int *number_neighbors,int_st ist){

  int i;  

  /*Finds the Neighbors for Each Sites*/
  if (ist.flag_pbc == 1) {
   for (i=0; i<ist.n_sites; i++) {
      neighborsperiodicboundary(i,neighbors,number_neighbors,ist);
   }
  }
  else {
   for (i=0; i<ist.n_sites; i++) {
      neighborsopenboundary(i,neighbors,number_neighbors,ist);
   }
  }
 
return; 
}

/*************************************************************/

void init_walkers_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *weights,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_up, MKL_Complex16 *overlap_inverse_down, double *weight_rescaling,int_st ist,cns_st cns) {

  int iw, is, ip, iwp1, iwp2, isp1, isp2;
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;
  FILE *pf = fopen("errors.dat", "a+");
  int i, j;

  /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<cns.max_number_walkers; iw++) {

     weights[iw] = One;
     overlap_up[iw] = One;
     overlap_down[iw] = One;

     iwp1 = iw * ist.n_sites * ist.n_up;
     iwp2 = iw * ist.n_sites * ist.n_down;

     /*Run Through Sites*/
     for ( is = 0; is<ist.n_sites ; is++ ) {

       isp1 = is * ist.n_up;
       isp2 = is * ist.n_down;

       /*Run Through Electrons*/
       for ( ip = 0; ip<ist.n_up ; ip++ ) {
         wf_up[ iwp1 + isp1 + ip] = trial_wf_up[ isp1 + ip];
       }
       for ( ip = 0 ; ip<ist.n_down ; ip++ ) {
         wf_down[ iwp2 + isp2 + ip] = trial_wf_down[ isp2 + ip];
       }

     } /*is*/

     iwp1 = iw * ist.n_up_sq;
     iwp2 = iw * ist.n_down_sq;

     for ( ip = 0; ip < ist.n_up; ip++) {
      overlap_inverse_up[iwp1 + ip * ist.n_up + ip] = One;
     }
    for ( ip = 0; ip < ist.n_down; ip++) {
      overlap_inverse_down[iwp2 + ip * ist.n_down + ip] = One;
     }

  } /*iw*/

  fprintf(pf, "trial wf\n");
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_up; j++) {
     fprintf(pf, "%f\t", trial_wf_up[i*ist.n_up+j].real);
    }
    fprintf(pf, "\n");
  }
  fprintf(pf, "\n\n"); fflush(pf);

  fprintf(pf, "trial wf\n");
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_down; j++) {
     fprintf(pf, "%f\t", trial_wf_down[i*ist.n_down+j].real);
    }
    fprintf(pf, "\n");
  }
  fprintf(pf, "\n\n"); fflush(pf);

  /*Also Initialize All WEight Rescaling Factors*/
  for ( iw=0 ; iw<10; iw++) {
     weight_rescaling[iw] = 1.0;
  }

return;
}

/****************************************************************/

void init_walkers_fermions_up(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *weights,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_inverse_up,double *weight_rescaling,int_st ist,cns_st cns) {

  int iw, is, ip, iwp1, isp1; 
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0; 
  FILE *pf = fopen("errors.dat", "a+"); 
  int i, j; 

  /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<cns.max_number_walkers; iw++) {

     weights[iw] = One; 
     overlap_up[iw] = One; 

     iwp1 = iw * ist.n_sites * ist.n_up;      

     /*Run Through Sites*/
     for ( is = 0; is<ist.n_sites ; is++ ) {
              
       isp1 = is * ist.n_up;

       /*Run Through Electrons*/
       for ( ip = 0; ip<ist.n_up ; ip++ ) { 
         wf_up[ iwp1 + isp1 + ip] = trial_wf_up[ isp1 + ip]; 
       }

     } /*is*/

     iwp1 = iw * ist.n_up_sq; 

     for ( ip = 0; ip < ist.n_up; ip++) {
      overlap_inverse_up[iwp1 + ip * ist.n_up + ip] = One; 
     }
     
  } /*iw*/ 

  fprintf(pf, "trial wf\n"); 
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_up; j++) {
     fprintf(pf, "%f\t", trial_wf_up[i*ist.n_up+j].real); 
    }
    fprintf(pf, "\n"); 
  }
  fprintf(pf, "\n\n"); fflush(pf); 

  /*Also Initialize All WEight Rescaling Factors*/
  for ( iw=0 ; iw<10; iw++) {
     weight_rescaling[iw] = 1.0; 
  }   

return; 
}

/*************************************************************/

void init_walkers_fermions_down(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down,MKL_Complex16 *weights,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_down,double *weight_rescaling,int_st ist,cns_st cns) {

  int iw, is, ip, iwp1, isp1; 
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;
  FILE *pf = fopen("errors.dat", "a+");
  int i, j;

  /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<cns.max_number_walkers; iw++) {

     weights[iw] = One;
     overlap_down[iw] = One;

     iwp1 = iw * ist.n_sites * ist.n_down;

     /*Run Through Sites*/
     for ( is = 0; is<ist.n_sites ; is++ ) {

       isp1 = is * ist.n_down;

       /*Run Through Electrons*/
       for ( ip = 0; ip<ist.n_down ; ip++ ) {
         wf_down[ iwp1 + isp1 + ip] = trial_wf_down[ isp1 + ip];
       }

     } /*is*/

     iwp1 = iw * ist.n_down_sq;

     for ( ip = 0; ip < ist.n_down; ip++) {
      overlap_inverse_down[iwp1 + ip * ist.n_down + ip] = One;
     }

  } /*iw*/

  fprintf(pf, "trial wf\n");
  for (i=0; i<ist.n_sites; i++) {
    for (j=0; j<ist.n_down; j++) {
     fprintf(pf, "%f\t", trial_wf_down[i*ist.n_down+j].real);
    }
    fprintf(pf, "\n");
  }
  fprintf(pf, "\n\n"); fflush(pf);

  /*Also Initialize All WEight Rescaling Factors*/
  for ( iw=0 ; iw<10; iw++) {
     weight_rescaling[iw] = 1.0;
  }

return;
}

/**************************************************************/

void init_walkers_chemistry_no_restart_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *weights,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_up, MKL_Complex16 *overlap_inverse_down, double *weight_rescaling,int_st ist,cns_st cns) {

  int iw, is, ip, iwp1, iwp2, isp1, isp2;
  int i, j; 
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;
  FILE *pf = fopen("errors.dat", "a+"); 

   /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<ist.n_walkers; iw++) {

     weights[iw] = One;
     overlap_up[iw] = One;
     overlap_down[iw] = One;

     iwp1 = iw * ist.n_spatial_orbitals_n_up;
     iwp2 = iw * ist.n_spatial_orbitals_n_down;

     /*Run Through Sites*/
     for ( is = 0; is<ist.n_spatial_orbitals ; is++ ) {

       isp1 = is * ist.n_up;
       isp2 = is * ist.n_down;

       /*Run Through Electrons*/
       /*Assumes First Determinant Is the One To Start On with WF*/
       for ( ip = 0; ip<ist.n_up ; ip++ ) {
         wf_up[ iwp1 + isp1 + ip] = trial_wf_up[ isp1 + ip];
       }
       for ( ip = 0 ; ip<ist.n_down ; ip++ ) {
         wf_down[ iwp2 + isp2 + ip] = trial_wf_down[ isp2 + ip];
       }

     } /*is*/

     iwp1 = iw * ist.n_up_sq;
     iwp2 = iw * ist.n_down_sq; 

     for ( ip = 0; ip < ist.n_up; ip++) {
      overlap_inverse_up[iwp1 + ip * ist.n_up + ip] = One;
     }
     for ( ip = 0; ip < ist.n_down; ip++) {
      overlap_inverse_down[iwp2 + ip * ist.n_down + ip] = One;
     }

  } /*iw*/

  /*Also Initialize All WEight Rescaling Factors*/
  for ( iw=0 ; iw<10; iw++) {
     weight_rescaling[iw] = 1.0;
  }

return;
}

/***************************************************************/

void init_walkers_chemistry_no_restart_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *weights,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up, MKL_Complex16 *overlap_inverse_down, double *weight_rescaling,int_st ist,cns_st cns) {

  int iw, is, ip, iwp1, iwp2, isp1, isp2; 
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;
  FILE *pf = fopen("errors.dat", "a+"); 

   fprintf(pf, "n determinants %d %d\n", ist.n_determinants_trial_phaseless, ist.n_determinants_trial_energy); fflush(pf); 

  /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<ist.n_walkers; iw++) {

     weights[iw] = One;

     iwp1 = iw * ist.n_spatial_orbitals * ist.n_up;
     iwp2 = iw * ist.n_spatial_orbitals * ist.n_down;

     /*Run Through Sites*/
     for ( is = 0; is<ist.n_spatial_orbitals ; is++ ) {

       isp1 = is * ist.n_up;
       isp2 = is * ist.n_down;

       /*Run Through Electrons*/
       /*Assumes First Determinant Is the One To Start On with WF*/
       for ( ip = 0; ip<ist.n_up ; ip++ ) {
         wf_up[ iwp1 + isp1 + ip] = trial_wf_up_phaseless[ isp1 + ip];
       }
       for ( ip = 0 ; ip<ist.n_down ; ip++ ) {
         wf_down[ iwp2 + isp2 + ip] = trial_wf_down_phaseless[ isp2 + ip];
       }

     } /*is*/

     /*Get Overlaps and Inverses*/
     iwp1 = iw * ist.n_determinants_phaseless_n_up_sq;
     iwp2 = iw * ist.n_determinants_phaseless_n_down_sq;

     isp1 = iw * ist.n_spatial_orbitals_n_up; 
     isp2 = iw * ist.n_spatial_orbitals_n_down; 

     overlap_total[iw] = compute_overlap_inverse_multi_phaseless(trial_wf_up_phaseless,trial_wf_down_phaseless,&wf_up[isp1],&wf_down[isp2],trial_determinant_coefficients_phaseless,&overlap_inverse_up[iwp1],&overlap_inverse_down[iwp2],&det_overlap_up[iw*ist.n_determinants_trial_phaseless],&det_overlap_down[iw*ist.n_determinants_trial_phaseless],ist); 

  } /*iw*/

  /*Also Initialize All WEight Rescaling Factors*/
  for ( iw=0 ; iw<10; iw++) {
     weight_rescaling[iw] = 1.0;
  }

return;
}

/**************************************************************************/

void init_walkers_chemistry_restart_multi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *weights,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *overlap_inverse_up, MKL_Complex16 *overlap_inverse_down, double *weight_rescaling,int_st ist,cns_st cns) {

  int iw, is, ip, iwp1, iwp2, isp1, isp2;
  FILE *pf = fopen("restart.par", "r"); 

  /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<ist.n_walkers; iw++) {

     iwp1 = iw * ist.n_spatial_orbitals_n_up;
     iwp2 = iw * ist.n_spatial_orbitals_n_down;

     /*Scan in Weights and Wave functions*/
     fscanf(pf, "%lf %lf\t", &weights[iw].real, &weights[iw].imag); 
  
     for (is=0; is<ist.n_spatial_orbitals; is++) {
       isp1 = is * ist.n_up;

       for (ip = 0; ip<ist.n_up; ip++) {
         fscanf(pf, "%lf %lf\t", &wf_up[ iwp1 + isp1 + ip].real, &wf_up[ iwp1 + isp1 + ip].imag); 
       }
     } 
     for (is=0; is<ist.n_spatial_orbitals; is++) {
       isp2 = is * ist.n_down; 
 
       for (ip = 0; ip<ist.n_down; ip++) {
         fscanf(pf, "%lf %lf\t", &wf_down[ iwp1 + isp2 + ip].real, &wf_down[ iwp1 + isp2 + ip].imag);
       }
     }     

     isp1 = iw * ist.n_determinants_phaseless_n_up_sq; 
     isp2 = iw * ist.n_determinants_phaseless_n_down_sq; 

     overlap_total[iw] = compute_overlap_inverse_multi_phaseless(trial_wf_up_phaseless,trial_wf_down_phaseless,&wf_up[iwp1],&wf_down[iwp2],trial_determinant_coefficients_phaseless,&overlap_inverse_up[isp1],&overlap_inverse_down[isp2],&det_overlap_up[iw*ist.n_determinants_trial_phaseless],&det_overlap_down[iw*ist.n_determinants_trial_phaseless],ist);
 
  } /*iw*/

  /*Also Initialize All WEight Rescaling Factors*/
  for ( iw=0 ; iw<10; iw++) {
     weight_rescaling[iw] = 1.0;
  }

fclose(pf);
return;
}

/**************************************************************/

void init_walkers_chemistry_restart_restricted(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *weights,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_inverse_up, MKL_Complex16 *overlap_inverse_down, double *weight_rescaling,int_st ist,cns_st cns) {

  int iw, is, ip, iwp1, iwp2, isp1, isp2;
  FILE *pf = fopen("restart.par", "r");

  /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<ist.n_walkers; iw++) {

     iwp1 = iw * ist.n_spatial_orbitals_n_up;
     iwp2 = iw * ist.n_spatial_orbitals_n_down;

     /*Scan in Weights and Wave functions*/
     fscanf(pf, "%lf %lf\t", &weights[iw].real, &weights[iw].imag);

     for (is=0; is<ist.n_spatial_orbitals; is++) {
       isp1 = is * ist.n_up;

       for (ip = 0; ip<ist.n_up; ip++) {
         fscanf(pf, "%lf %lf\t", &wf_up[ iwp1 + isp1 + ip].real, &wf_up[ iwp1 + isp1 + ip].imag);
       }
     }
     for (is=0; is<ist.n_spatial_orbitals; is++) {
       isp2 = is * ist.n_down;

       for (ip = 0; ip<ist.n_down; ip++) {
         fscanf(pf, "%lf %lf\t", &wf_down[ iwp1 + isp2 + ip].real, &wf_down[ iwp1 + isp2 + ip].imag);
       }
     }

     /*Now Determine New Overlaps*/
     overlap_up[iw] = compute_overlap_inverse(trial_wf_up,&wf_up[iw*ist.n_spatial_orbitals_n_up],&overlap_inverse_up[iw*ist.n_up_sq],ist,0);   
     overlap_down[iw] = compute_overlap_inverse(trial_wf_down,&wf_down[iw*ist.n_spatial_orbitals_n_down],&overlap_inverse_down[iw*ist.n_down_sq],ist,1);     

  } /*iw*/

  /*Also Initialize All WEight Rescaling Factors*/
  for ( iw=0 ; iw<10; iw++) {
     weight_rescaling[iw] = 1.0;
  }

fclose(pf);
return;
}

/*************************************************************/

void init_walkers_bosons(MKL_Complex16 *wf,MKL_Complex16 *trial_wf,MKL_Complex16 *weights,MKL_Complex16 *overlap,double *weight_rescaling,int_st ist,cns_st cns) {

  int iw, is, iwp; 
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;

  /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<cns.max_number_walkers; iw++) {
     weights[iw] = One;
     overlap[iw] = One;

     iwp = iw * ist.n_sites;

     /*Run Through Sites*/
     for ( is = 0; is<ist.n_sites ; is++ ) {

       wf[ iwp + is ] = trial_wf[ is ];

     } /*is*/

  } /*iw*/

  /*Also Initialize All WEight Rescaling Factors*/
  for ( iw=0 ; iw<10; iw++) {
     weight_rescaling[iw] = 1.0;
  }

return;
}

/***********************************************************************/

void init_walkers_bose_fermi(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *weights,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_bosons,double *weight_rescaling,int_st ist,cns_st cns) {

  int iw, is, ip, iwp1, iwp2, isp1, isp2;
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;

  /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<cns.max_number_walkers; iw++) {
     weights[iw] = One;
     overlap_up[iw] = One;
     overlap_down[iw] = One;
     overlap_bosons[iw] = One; 

     iwp1 = iw * ist.n_sites * ist.n_up;
     iwp2 = iw * ist.n_sites * ist.n_down;

     /*Run Through Sites*/
     for ( is = 0; is<ist.n_sites ; is++ ) {

       isp1 = is * ist.n_up;
       isp2 = is * ist.n_down;

       /*Run Through Electrons*/
       for ( ip = 0; ip<ist.n_up ; ip++ ) {
         wf_up[ iwp1 + isp1 + ip] = trial_wf_up[ isp1 + ip];
       }
       for ( ip = 0 ; ip<ist.n_down ; ip++ ) {
         wf_down[ iwp2 + isp2 + ip] = trial_wf_down[ isp2 + ip];
       }
       wf_bosons[ iw * ist.n_sites + is] = trial_wf_bosons[ is ]; 

     } /*is*/

  } /*iw*/

  /*Also Initialize All WEight Rescaling Factors*/
  for ( iw=0 ; iw<10; iw++) {
     weight_rescaling[iw] = 1.0;
  }

return;
}

/***************************************************************/

void init_walkers_bose_fermi_spin_polarized(MKL_Complex16 *wf_up,MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *weights,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_bosons,double *weight_rescaling,int_st ist,cns_st cns) {

  int iw, is, ip, iwp1, isp1;
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;

  /*Initialize All Walker Weights and Overlaps*/
  for ( iw=0 ; iw<cns.max_number_walkers; iw++) {
     weights[iw] = One;
     overlap_up[iw] = One;
     overlap_bosons[iw] = One;

     iwp1 = iw * ist.n_sites * ist.n_up;

     /*Run Through Sites*/
     for ( is = 0; is<ist.n_sites ; is++ ) {

       isp1 = is * ist.n_up;

       /*Run Through Electrons*/
       for ( ip = 0; ip<ist.n_up ; ip++ ) {
         wf_up[ iwp1 + isp1 + ip] = trial_wf_up[ isp1 + ip];
       }
       wf_bosons[ iw * ist.n_sites + is] = trial_wf_bosons[ is ];

     } /*is*/

  } /*iw*/

  /*Also Initialize All WEight Rescaling Factors*/
  for ( iw=0 ; iw<10; iw++) {
    weight_rescaling[iw] = 1.0;
  }

return;
}

/**************************************************************/

void init_wf_fermions(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,double *kinetic_eigs_fermions,double *kinetic_eigvecs_fermions,int_st ist,cns_st cns) {

  /*Determine Which Type of Trial WF To Form*/
  if ( ist.flag_read_in_trial == 0 ) {

    if ( ist.flag_trial == 0 ) {
     /*Just an Identity Matrix*/
     /*Be Careful, Because May Not Reach Correct Steady State Because Could Be Orthogonal!!!!*/ 
     trial_identity_fermions(trial_wf_up,trial_wf_down,ist); 
    }
    else if ( ist.flag_trial == 1 ) {
      /*The Free Case*/
      trial_free_fermions(trial_wf_up,trial_wf_down,kinetic_eigs_fermions,kinetic_eigvecs_fermions,ist,cns);  
    }
    else if ( ist.flag_trial == 2 ) {
      /*The Random Case*/
     trial_random_fermions(trial_wf_up,trial_wf_down,ist);
    }

  }
  else {
    
     /*Read in Wavefunction*/
     trial_fermions_readin(trial_wf_up,trial_wf_down,ist); 
 
  }

return; 
}

/**************************************************************/

void init_wf_fermions_up(MKL_Complex16 *trial_wf_up,double *kinetic_eigs_fermions,double *kinetic_eigvecs_fermions,int_st ist,cns_st cns) {

  /*Determine Which Type of Trial WF To Form*/
  if ( ist.flag_read_in_trial == 0 ) {

    if ( ist.flag_trial == 0 ) {
     /*Just an Identity Matrix*/
     /*Be Careful, Because May Not Reach Correct Steady State Because Could Be Orthogonal!!!!*/
     trial_identity_fermions_up(trial_wf_up,ist);
    }
    else if ( ist.flag_trial == 1 ) {
      /*The Free Case*/
      trial_free_fermions_up(trial_wf_up,kinetic_eigs_fermions,kinetic_eigvecs_fermions,ist,cns);
    }
    else if ( ist.flag_trial == 2 ) {
      /*The Random Case*/
     trial_random_fermions_up(trial_wf_up,ist);
    }

  }
  else {

     /*Read in Wavefunction*/
     trial_fermions_readin_up(trial_wf_up,ist);

  }

return;
}

/**************************************************************/

void init_wf_fermions_down(MKL_Complex16 *trial_wf_down,double *kinetic_eigs_fermions,double *kinetic_eigvecs_fermions,int_st ist,cns_st cns) {

  /*Determine Which Type of Trial WF To Form*/
  if ( ist.flag_read_in_trial == 0 ) {
     
    if ( ist.flag_trial == 0 ) { 
     /*Just an Identity Matrix*/
     /*Be Careful, Because May Not Reach Correct Steady State Because Could Be Orthogonal!!!!*/
     trial_identity_fermions_down(trial_wf_down,ist);
    }
    else if ( ist.flag_trial == 1 ) {
      /*The Free Case*/
      trial_free_fermions_down(trial_wf_down,kinetic_eigs_fermions,kinetic_eigvecs_fermions,ist,cns);
    }  
    else if ( ist.flag_trial == 2 ) {
      /*The Random Case*/
     trial_random_fermions_down(trial_wf_down,ist);
    }

  }
  else {

     /*Read in Wavefunction*/
     trial_fermions_readin_down(trial_wf_down,ist);

  }

return;
}

/***************************************************************/

int init_wf_chemistry_restricted(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *kinetic_matrix_original_fermions_up,int *list_occupied_orbitals,int_st ist,cns_st cns) {

     int n_max_orbitals; 

     /*If Read-In*/
     if ( ist.flag_read_in_trial == 0 ) {

       /*Start As HF Solution*/
       n_max_orbitals = trial_chemistry_restricted_noreadin(trial_wf_up,trial_wf_down,kinetic_matrix_original_fermions_up,list_occupied_orbitals,ist);

     }
     else {
  
       /*Read In Wavefunction*/
       n_max_orbitals = trial_chemistry_restricted_readin(trial_wf_up,trial_wf_down,list_occupied_orbitals,ist); 

     } 

return(n_max_orbitals); 
}

/***************************************************************/

int init_wf_chemistry_unrestricted(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *kinetic_matrix_original_fermions_up,MKL_Complex16 *kinetic_matrix_original_fermions_down,int *list_occupied_orbitals,int_st ist,cns_st cns) {

     int n_max_orbitals_trial_energy; 

     if ( ist.flag_read_in_trial == 0 ) {

       /*Start As HF Solution*/
       n_max_orbitals_trial_energy = trial_chemistry_unrestricted_noreadin(trial_wf_up,trial_wf_down,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,list_occupied_orbitals,ist);

    }
    else {
  
       /*Read In WF*/
       n_max_orbitals_trial_energy = trial_chemistry_restricted_readin(trial_wf_up,trial_wf_down,list_occupied_orbitals,ist); 
  
    } 

return(n_max_orbitals_trial_energy); 
}

/***************************************************************/

int init_wf_chemistry_multi(MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinants_coefficients_phaseless,int *list_occupied_orbitals,int_st ist,cns_st cns) {

     int n_max_orbital_trial_energy;  

     /*Read In Wavefunction*/
     n_max_orbital_trial_energy = trial_chemistry_multi_readin(trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinants_coefficients_phaseless,list_occupied_orbitals,ist);

return(n_max_orbital_trial_energy); 
}

/***************************************************************/

void init_wf_bosons(MKL_Complex16 *trial_wf,double *kinetic_eigs_bosons,double *kinetic_eigvecs_bosons,int_st ist,cns_st cns) {

  /*Determine Which Type of Trial WF To Form*/
  if ( ist.flag_trial == 0 ) {
    /*Just an Identity Matrix*/

    /*Be Careful, Because May Not Reach Correct Steady State Because Could Be Orthogonal!!!!*/
    trial_identity_bosons(trial_wf,ist);

  }
  else if ( ist.flag_trial == 1 ) {
    /*The Free Case*/

    trial_free_bosons(trial_wf,kinetic_eigs_bosons,kinetic_eigvecs_bosons,ist,cns);

  }
  else if ( ist.flag_trial == 2 ) {
     /*The Restricted Hartree Fock Case*/

    trial_random_bosons(trial_wf,ist); 

  }
  else if ( ist.flag_trial == 3 ) {
      /*Just a Random WF*/

     trial_random_bosons(trial_wf,ist);

  }
  else {
     /*Put In Close to Exact Answer*/

     trial_wf[0].real = .666215;
     trial_wf[0].imag = 0.0;

     trial_wf[1].real = .529221;
     trial_wf[1].imag = 0.0;

     trial_wf[2].real = .427738;
     trial_wf[2].imag = 0.0;

  }

return;
}

/**************************************************************************************/

void init_kinetic_fermions(double *kinetic_full_fermions,double *kinetic_backwards_fermions,double *kinetic_forwards_half_fermions,double *kinetic_eigs_fermions,double *kinetic_eigvecs_fermions,int *neighbors,int *number_neighbors,int_st ist,cns_st cns) {

  /*Creates Kinetic Matrices********/

  int i, j, k;
  int info; 

  /*Zero Eigenvectors*/ 
  dzero_vec(kinetic_eigvecs_fermions,ist.n_sites_sq); 

  /*If Non-Zero T***************/
  if (cns.t_fermions!=0) {

    /*Creating Original Kinetic Matrix and Then Diagonalizing It*/
    for (i=0; i<ist.n_sites;  i++) {
      kinetic_eigvecs_fermions[i*ist.n_sites+i]+= cns.V_fermions;  
      //kinetic_full_fermions[i*ist.n_sites+i]+= cns.U_fermions/2.0; 
      for (j=0; j<number_neighbors[i]; j++) {
        kinetic_eigvecs_fermions[i*ist.n_sites+neighbors[i*4+j]]=-cns.t_fermions;
        kinetic_eigvecs_fermions[neighbors[i*4+j]*ist.n_sites+i]=kinetic_eigvecs_fermions[i*ist.n_sites+neighbors[i*4+j]];

        /*Add One Body Terms from V*/
        if ( j < (int)(number_neighbors[i]/2.0)) {
          kinetic_eigvecs_fermions[neighbors[i*4+j]+ist.n_sites*neighbors[i*4+j]]+= cns.V_fermions; 
        }
      }   
     }

    /*Now Compute Eigenvalues and Eigenvectors Using NRC Routines*/
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',ist.n_sites,kinetic_eigvecs_fermions,ist.n_sites,kinetic_eigs_fermions);

    dzero_vec(kinetic_full_fermions,ist.n_sites_sq); 
    dzero_vec(kinetic_backwards_fermions,ist.n_sites_sq); 
    dzero_vec(kinetic_forwards_half_fermions,ist.n_sites_sq); 

    /*Find Various Kinetic Terms*/
    for (k=0; k<ist.n_sites; k++) {
     for (i=0; i<ist.n_sites; i++) {
      for (j=0; j<ist.n_sites; j++) {
        kinetic_full_fermions[i*ist.n_sites+j]+=exp(-1*kinetic_eigs_fermions[k]*cns.dtau)*kinetic_eigvecs_fermions[k*ist.n_sites+i]*kinetic_eigvecs_fermions[k*ist.n_sites+j];
        kinetic_backwards_fermions[i*ist.n_sites+j]+=exp(kinetic_eigs_fermions[k]*cns.dtau/2.0)*kinetic_eigvecs_fermions[k*ist.n_sites+i]*kinetic_eigvecs_fermions[k*ist.n_sites+j];
        kinetic_forwards_half_fermions[i*ist.n_sites+j]+=exp(-1*kinetic_eigs_fermions[k]*cns.dtau/2.0)*kinetic_eigvecs_fermions[k*ist.n_sites+i]*kinetic_eigvecs_fermions[k*ist.n_sites+j];
      }
     }
    }

  }
  else { /*If Zero T********************************************/

    /*Make Fermions Matrices Diagonal*/
    for (i=0; i<ist.n_sites; i++) {
     for (j=0; j<ist.n_sites; j++) {
       kinetic_full_fermions[i*ist.n_sites+j]=kinetic_backwards_fermions[i*ist.n_sites+j]=kinetic_forwards_half_fermions[i*ist.n_sites+j]=0; 
       kinetic_eigvecs_fermions[i*ist.n_sites+j]=0;
     }
     kinetic_eigs_fermions[i]=0;
     kinetic_full_fermions[i*ist.n_sites+i]=kinetic_backwards_fermions[i*ist.n_sites+i]=kinetic_forwards_half_fermions[i*ist.n_sites+i]=1; 
     kinetic_eigvecs_fermions[i*ist.n_sites+i]=1;   
   }

  } /**********************************************************/

return; 
}

/******************************************************************************/

void init_kinetic_chemistry_real_restricted(MKL_Complex16 *kinetic_matrix_original_fermions, MKL_Complex16 *potential_supermatrix_fermions, double *kinetic_full_fermions,double *kinetic_backwards_fermions,double *kinetic_forwards_half_fermions,double *kinetic_eigs_fermions,double *kinetic_eigvecs_fermions,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_eigs_fermions,MKL_Complex16 *mean_field_density_first_transform,MKL_Complex16 *mean_field_density_second_transform,double *potential_one_body_matrix,int_st ist,cns_st cns) {

  /*Creates Kinetic Matrices********/

  int i, j, k, l;  
  int gamma;
  int info;  
  int i_orbital, l_orbital;  
  int i_spin, j_spin, k_spin, l_spin; 
  double subtracted_element; 
  MKL_Complex16 gamma_constant;
  MKL_Complex16 product;  

   /*Find One-Body Matrix That Derives From Vijkl Matrix*/
   dzero_vec(kinetic_eigvecs_fermions, ist.n_spatial_orbitals_sq); 
   for (i=0; i<ist.n_spatial_orbitals; i++) {
      for (k=0; k<ist.n_spatial_orbitals; k++) {
         kinetic_eigvecs_fermions[i*ist.n_spatial_orbitals+k] = kinetic_matrix_original_fermions[i*ist.n_spatial_orbitals+k].real + potential_one_body_matrix[i*ist.n_spatial_orbitals+k];  
       }
    }

     /*Subtract Out Mean Field Portion, If Necessary************************************************************************/
     if ( ist.flag_meanfield == 1 ) {

      for (gamma = 0; gamma<ist.n_spin_orbitals_sq; gamma++) {

        if ( Cabs(mean_field_density_first_transform[gamma]) > .00001 || Cabs(mean_field_density_second_transform[gamma]) > .000001 ) {

          gamma_constant = RCmul(0.25, potential_eigs_fermions[gamma]); 

          for (i=0; i<ist.n_spin_orbitals; i++) {
             i_spin = i%2;
             i_orbital = (int)(i/2.0);

             for (l=0; l<ist.n_spin_orbitals; l++) {
              l_spin = l%2;
              l_orbital = (int)(l/2.0);

              /*Only Works If Two Spins are the Same*/
              if ( i_spin == l_spin ) {

                 /*If Spin Up*/
                 if ( i_spin == 0 ) {
                    /*I AM ASSUMING UP AND DOWN MATRICES ARE SYMMETRIC HERE...CHANGE FOR MORE GENERAL CASE*/

                    /*From First Transform*/
                    product = Cmul(gamma_constant, Cmul(potential_eigvecs_fermions[i*ist.n_spin_orbitals_third+l*ist.n_spin_orbitals_sq+gamma], mean_field_density_first_transform[gamma])); 
                    
                    kinetic_eigvecs_fermions[i_orbital*ist.n_spatial_orbitals+l_orbital] += product.real;  

                    product = Cmul(gamma_constant, Cmul(conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third+i*ist.n_spin_orbitals_sq+gamma]), mean_field_density_first_transform[gamma])); 
                    kinetic_eigvecs_fermions[l_orbital*ist.n_spatial_orbitals+i_orbital] += product.real; 


                    /*From Second Transform*/
                    product = Cmul(gamma_constant, Cmul(potential_eigvecs_fermions[i*ist.n_spin_orbitals_third+l*ist.n_spin_orbitals_sq+gamma], mean_field_density_second_transform[gamma])); 
                    kinetic_eigvecs_fermions[i_orbital*ist.n_spatial_orbitals+l_orbital] += product.real;  

                    product = Cmul(gamma_constant, Cmul(conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third+i*ist.n_spin_orbitals_sq+gamma]), mean_field_density_second_transform[gamma])); 
                    kinetic_eigvecs_fermions[l_orbital*ist.n_spatial_orbitals+i_orbital] -= product.real; 

                }
              }

             } /*l*/
            } /*i*/

         } /*gamma*/
        } /*run through gamma*/
     } /*If Mean Field*********************************************************************************************************/

    /*Now Compute Eigenvalues and Eigenvectors Using NRC Routines*/
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',ist.n_spatial_orbitals,kinetic_eigvecs_fermions,ist.n_spatial_orbitals,kinetic_eigs_fermions);

    dzero_vec(kinetic_full_fermions,ist.n_spatial_orbitals_sq);
    dzero_vec(kinetic_backwards_fermions,ist.n_spatial_orbitals_sq);
    dzero_vec(kinetic_forwards_half_fermions,ist.n_spatial_orbitals_sq);

    /*Find Various Kinetic Terms*/
    for (k=0; k<ist.n_spatial_orbitals; k++) {
     for (i=0; i<ist.n_spatial_orbitals; i++) {
      for (j=0; j<ist.n_spatial_orbitals; j++) {
        kinetic_full_fermions[i*ist.n_spatial_orbitals+j]+=exp(-1*kinetic_eigs_fermions[k]*cns.dtau)*kinetic_eigvecs_fermions[k*ist.n_spatial_orbitals+i]*kinetic_eigvecs_fermions[k*ist.n_spatial_orbitals+j];
        kinetic_backwards_fermions[i*ist.n_spatial_orbitals+j]+=exp(kinetic_eigs_fermions[k]*cns.dtau/2.0)*kinetic_eigvecs_fermions[k*ist.n_spatial_orbitals+i]*kinetic_eigvecs_fermions[k*ist.n_spatial_orbitals+j];
        kinetic_forwards_half_fermions[i*ist.n_spatial_orbitals+j]+=exp(-1*kinetic_eigs_fermions[k]*cns.dtau/2.0)*kinetic_eigvecs_fermions[k*ist.n_spatial_orbitals+i]*kinetic_eigvecs_fermions[k*ist.n_spatial_orbitals+j];

      }
     }
    }

return; 
}

/*******************************************************************************/

void init_kinetic_chemistry_real_unrestricted(MKL_Complex16 *kinetic_matrix_original_fermions_up, MKL_Complex16 *kinetic_matrix_original_fermions_down, MKL_Complex16 *potential_supermatrix_fermions, double *kinetic_full_fermions_up,double *kinetic_backwards_fermions_up,double *kinetic_forwards_half_fermions_up,double *kinetic_eigs_fermions_up,double *kinetic_eigvecs_fermions_up,double *kinetic_full_fermions_down,double *kinetic_backwards_fermions_down,double *kinetic_forwards_half_fermions_down,double *kinetic_eigs_fermions_down,double *kinetic_eigvecs_fermions_down,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_eigs_fermions,MKL_Complex16 *mean_field_density_first_transform,MKL_Complex16 *mean_field_density_second_transform,double *potential_one_body_matrix_up,double *potential_one_body_matrix_down,int_st ist,cns_st cns) {

  /*Creates Kinetic Matrices********/

  int i, j, k, l;
  int info; 
  int gamma;
  int i_orbital, l_orbital;
  int i_spin, j_spin, k_spin, l_spin;
  double subtracted_element;
  MKL_Complex16 gamma_constant;
  MKL_Complex16 product;
  FILE *pf = fopen("errors.dat", "a+"); 

   /*Find One-Body Matrix That Derives From Vijkl Matrix*/
   dzero_vec(kinetic_eigvecs_fermions_up, ist.n_spatial_orbitals_sq);
   dzero_vec(kinetic_eigvecs_fermions_down, ist.n_spatial_orbitals_sq); 
   for (i=0; i<ist.n_spatial_orbitals; i++) {
      for (k=0; k<ist.n_spatial_orbitals; k++) {
         kinetic_eigvecs_fermions_up[i*ist.n_spatial_orbitals+k] = kinetic_matrix_original_fermions_up[i*ist.n_spatial_orbitals+k].real + potential_one_body_matrix_up[i*ist.n_spatial_orbitals+k];
         kinetic_eigvecs_fermions_down[i*ist.n_spatial_orbitals+k] = kinetic_matrix_original_fermions_down[i*ist.n_spatial_orbitals+k].real + potential_one_body_matrix_down[i*ist.n_spatial_orbitals+k]; 
       }
    }

     /*Subtract Out Mean Field Portion, If Necessary************************************************************************/
     if ( ist.flag_meanfield == 1 ) {

      for (gamma = 0; gamma<ist.n_spin_orbitals_sq; gamma++) {

        if ( Cabs(mean_field_density_first_transform[gamma]) > .00001 || Cabs(mean_field_density_second_transform[gamma]) > .000001 ) {

          gamma_constant = RCmul(0.25, potential_eigs_fermions[gamma]);

          for (i=0; i<ist.n_spin_orbitals; i++) {
             i_spin = i%2;
             i_orbital = (int)(i/2.0);

             for (l=0; l<ist.n_spin_orbitals; l++) {
              l_spin = l%2;
              l_orbital = (int)(l/2.0);
 
              /*Only Works If Two Spins are the Same*/
              if ( i_spin == l_spin ) {

                 /*If Spin Up*/
                 if ( i_spin == 0 ) {

                    /*From First Transform*/
                    product = Cmul(gamma_constant, Cmul(potential_eigvecs_fermions[i*ist.n_spin_orbitals_third+l*ist.n_spin_orbitals_sq+gamma], mean_field_density_first_transform[gamma]));
                    kinetic_eigvecs_fermions_up[i_orbital*ist.n_spatial_orbitals+l_orbital] += product.real;

                    product = Cmul(gamma_constant, Cmul(conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third+i*ist.n_spin_orbitals_sq+gamma]), mean_field_density_first_transform[gamma]));
                    kinetic_eigvecs_fermions_up[l_orbital*ist.n_spatial_orbitals+i_orbital] += product.real;

                    /*From Second Transform*/
                    product = Cmul(gamma_constant, Cmul(potential_eigvecs_fermions[i*ist.n_spin_orbitals_third+l*ist.n_spin_orbitals_sq+gamma], mean_field_density_second_transform[gamma]));
                    kinetic_eigvecs_fermions_up[i_orbital*ist.n_spatial_orbitals+l_orbital] += product.real;

                    product = Cmul(gamma_constant, Cmul(conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third+i*ist.n_spin_orbitals_sq+gamma]), mean_field_density_second_transform[gamma]));
                    kinetic_eigvecs_fermions_up[l_orbital*ist.n_spatial_orbitals+i_orbital] -= product.real;

                } /*If Spin Up*/
                else { 
                    /*From First Transform*/
                    product = Cmul(gamma_constant, Cmul(potential_eigvecs_fermions[i*ist.n_spin_orbitals_third+l*ist.n_spin_orbitals_sq+gamma], mean_field_density_first_transform[gamma]));
                    kinetic_eigvecs_fermions_down[i_orbital*ist.n_spatial_orbitals+l_orbital] += product.real;

                    product = Cmul(gamma_constant, Cmul(conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third+i*ist.n_spin_orbitals_sq+gamma]), mean_field_density_first_transform[gamma]));
                    kinetic_eigvecs_fermions_down[l_orbital*ist.n_spatial_orbitals+i_orbital] += product.real;

                    /*From Second Transform*/
                    product = Cmul(gamma_constant, Cmul(potential_eigvecs_fermions[i*ist.n_spin_orbitals_third+l*ist.n_spin_orbitals_sq+gamma], mean_field_density_second_transform[gamma]));
                    kinetic_eigvecs_fermions_down[i_orbital*ist.n_spatial_orbitals+l_orbital] += product.real;

                    product = Cmul(gamma_constant, Cmul(conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third+i*ist.n_spin_orbitals_sq+gamma]), mean_field_density_second_transform[gamma]));
                    kinetic_eigvecs_fermions_down[l_orbital*ist.n_spatial_orbitals+i_orbital] -= product.real;

                } /*If Spin Down*/

              }

             } /*l*/
            } /*i*/

         } /*gamma*/
        } /*run through gamma*/
     } /*If Mean Field*********************************************************************************************************/

    /*Now Compute Eigenvalues and Eigenvectors Using NRC Routines*/
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',ist.n_spatial_orbitals,kinetic_eigvecs_fermions_up,ist.n_spatial_orbitals,kinetic_eigs_fermions_up);
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',ist.n_spatial_orbitals,kinetic_eigvecs_fermions_down,ist.n_spatial_orbitals,kinetic_eigs_fermions_down);

    dzero_vec(kinetic_full_fermions_up,ist.n_spatial_orbitals_sq);
    dzero_vec(kinetic_backwards_fermions_up,ist.n_spatial_orbitals_sq);
    dzero_vec(kinetic_forwards_half_fermions_up,ist.n_spatial_orbitals_sq);

    dzero_vec(kinetic_full_fermions_down,ist.n_spatial_orbitals_sq);
    dzero_vec(kinetic_backwards_fermions_down,ist.n_spatial_orbitals_sq);
    dzero_vec(kinetic_forwards_half_fermions_down,ist.n_spatial_orbitals_sq);

    /*Find Various Kinetic Terms*/
    for (k=0; k<ist.n_spatial_orbitals; k++) {
     for (i=0; i<ist.n_spatial_orbitals; i++) {
      for (j=0; j<ist.n_spatial_orbitals; j++) {
        kinetic_full_fermions_up[i*ist.n_spatial_orbitals+j]+=exp(-1*kinetic_eigs_fermions_up[k]*cns.dtau)*kinetic_eigvecs_fermions_up[k*ist.n_spatial_orbitals+i]*kinetic_eigvecs_fermions_up[k*ist.n_spatial_orbitals+j];

        kinetic_backwards_fermions_up[i*ist.n_spatial_orbitals+j]+=exp(kinetic_eigs_fermions_up[k]*cns.dtau/2.0)*kinetic_eigvecs_fermions_up[k*ist.n_spatial_orbitals+i]*kinetic_eigvecs_fermions_up[k*ist.n_spatial_orbitals+j];

        kinetic_forwards_half_fermions_up[i*ist.n_spatial_orbitals+j]+=exp(-1*kinetic_eigs_fermions_up[k]*cns.dtau/2.0)*kinetic_eigvecs_fermions_up[k*ist.n_spatial_orbitals+i]*kinetic_eigvecs_fermions_up[k*ist.n_spatial_orbitals+j];

        kinetic_full_fermions_down[i*ist.n_spatial_orbitals+j]+=exp(-1*kinetic_eigs_fermions_down[k]*cns.dtau)*kinetic_eigvecs_fermions_down[k*ist.n_spatial_orbitals+i]*kinetic_eigvecs_fermions_down[k*ist.n_spatial_orbitals+j];

        kinetic_backwards_fermions_down[i*ist.n_spatial_orbitals+j]+=exp(kinetic_eigs_fermions_down[k]*cns.dtau/2.0)*kinetic_eigvecs_fermions_down[k*ist.n_spatial_orbitals+i]*kinetic_eigvecs_fermions_down[k*ist.n_spatial_orbitals+j];

        kinetic_forwards_half_fermions_down[i*ist.n_spatial_orbitals+j]+=exp(-1*kinetic_eigs_fermions_down[k]*cns.dtau/2.0)*kinetic_eigvecs_fermions_down[k*ist.n_spatial_orbitals+i]*kinetic_eigvecs_fermions_down[k*ist.n_spatial_orbitals+j];
      }
     }
    }

return;
}

/*******************************************************************************/

void init_kinetic_bosons(double *kinetic_full_bosons,double *kinetic_backwards_bosons,double *kinetic_forwards_half_bosons,double *kinetic_eigs_bosons,double *kinetic_eigvecs_bosons,double *kinetic_matrix_bosons,double *potential_matrix_bosons,int *neighbors,int *number_neighbors,int_st ist,cns_st cns) {

  /*Creates Kinetic Matrices********/

  int i, j, k;
  int info; 

  dzero_vec(kinetic_eigvecs_bosons,ist.n_sites_sq); 
  /*If Non-Zero T***************/
  if (cns.t_bosons!=0) {

    /*Creating Original Kinetic Matrix and Then Diagonalizing It*/
    for (i=0; i<ist.n_sites;  i++) {
      potential_matrix_bosons[i*ist.n_sites+i] = 1.0; 
      for (j=0; j<number_neighbors[i]; j++) {
        kinetic_eigvecs_bosons[i*ist.n_sites+neighbors[i*4+j]]=-cns.t_bosons;
        kinetic_eigvecs_bosons[neighbors[i*4+j]*ist.n_sites+i]=kinetic_eigvecs_bosons[i*ist.n_sites+neighbors[i*4+j]];
        
        kinetic_matrix_bosons[i*ist.n_sites+neighbors[i*4+j]] = -1.0; 
        kinetic_matrix_bosons[neighbors[i*4+j]*ist.n_sites+i] = kinetic_matrix_bosons[i*ist.n_sites+neighbors[i*4+j]]; 
      }
     }

    /*Now Compute Eigenvalues and Eigenvectors Using NRC Routines*/
    info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',ist.n_sites,kinetic_eigvecs_bosons,ist.n_sites,kinetic_eigs_bosons);

    dzero_vec(kinetic_full_bosons,ist.n_sites_sq);
    dzero_vec(kinetic_backwards_bosons,ist.n_sites_sq);
    dzero_vec(kinetic_forwards_half_bosons,ist.n_sites_sq);

    /*Find Various Kinetic Terms*/
    for (k=0; k<ist.n_sites; k++) {
     for (i=0; i<ist.n_sites; i++) {
      for (j=0; j<ist.n_sites; j++) {
        kinetic_full_bosons[i*ist.n_sites+j]+=exp(-1*kinetic_eigs_bosons[k]*cns.dtau)*kinetic_eigvecs_bosons[k*ist.n_sites+i]*kinetic_eigvecs_bosons[k*ist.n_sites+j];
        kinetic_backwards_bosons[i*ist.n_sites+j]+=exp(kinetic_eigs_bosons[k]*cns.dtau/2.0)*kinetic_eigvecs_bosons[k*ist.n_sites+i]*kinetic_eigvecs_bosons[k*ist.n_sites+j];
        kinetic_forwards_half_bosons[i*ist.n_sites+j]+=exp(-1*kinetic_eigs_bosons[k]*cns.dtau/2.0)*kinetic_eigvecs_bosons[k*ist.n_sites+i]*kinetic_eigvecs_bosons[k*ist.n_sites+j];
      }
     }
    }

  }
  else { /*If Zero T********************************************/

    /*Make Fermions Matrices Diagonal*/
    for (i=0; i<ist.n_sites; i++) {
     for (j=0; j<ist.n_sites; j++) {
       kinetic_full_bosons[i*ist.n_sites+j]=kinetic_backwards_bosons[i*ist.n_sites+j]=kinetic_forwards_half_bosons[i*ist.n_sites+j]=0;
       kinetic_eigvecs_bosons[i*ist.n_sites+j]=0;
     }
     kinetic_eigs_bosons[i]=0;
     kinetic_full_bosons[i*ist.n_sites+i]=kinetic_backwards_bosons[i*ist.n_sites+i]=kinetic_forwards_half_bosons[i*ist.n_sites+i]=1;
     kinetic_eigvecs_bosons[i*ist.n_sites+i]=1;
   }

  } /**********************************************************/

return;
}

/***********************************************************/

void init_kinetic_many(double *kinetic_many_full,double *kinetic_eigs,double *kinetic_eigvecs,int *neighbors,int *number_neighbors,int_st ist,cns_st cns) {

    /*Create A Kinetic Matrix with Several Large Time Steps To Create Free Trial Wavefunction*/
   
    int i, j, k; 
    double full_time_step = cns.dtau * ist.n_steps_orthogonalize; 
 
    /*If Non-Zero T*/
    if (cns.t_fermions != 0 ) {  
  
       dzero_vec(kinetic_many_full,ist.n_sites_sq);
 
       /*Find Various Kinetic Terms*/
       for (k=0; k<ist.n_sites; k++) {
         for (i=0; i<ist.n_sites; i++) {
           for (j=0; j<ist.n_sites; j++) {
             kinetic_many_full[i*ist.n_sites+j]+=exp(-1*kinetic_eigs[k]*full_time_step)*kinetic_eigvecs[i*ist.n_sites+k]*kinetic_eigvecs[j*ist.n_sites+k];
           }
         }
       }  

      print_dmat(kinetic_many_full,ist.n_sites,ist.n_sites,"errors.dat"); 

    }
    else {  /*If T=0***************************************************************/
   
       /*Make Fermions Matrices Diagonal*/
       for (i=0; i<ist.n_sites; i++) {
        for (j=0; j<ist.n_sites; j++) {
         kinetic_many_full[i*ist.n_sites+j]=0.0; 
        }
       kinetic_many_full[i*ist.n_sites+i]=1.0; 
     } 

   }

return; 
}

/*****************************************************************************************/

double get_trial_energy_density_restricted(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *trial_density_up_total,MKL_Complex16 *trial_density_down_total,MKL_Complex16 *kinetic_original_matrix,MKL_Complex16 *potential_original_matrix,MKL_Complex16 *energy_shifts,int_st ist) {

   /*Obtains the Trial Energy Based Upon the Trial Wavefunction for Use in Estimating Walker Weights*/
   MKL_Complex16 *One;
   MKL_Complex16 trial_ke, trial_pe;
   MKL_Complex16 *trial_density_down, *trial_density_up;
   double trial_energy;

   One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
   One[0].real = 1.0; One[0].imag = 0.0;

   /*Density Matrices Down and Up*/
   trial_density_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
   trial_density_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

   /*Calculate Density Matrix*/
   cmat_transpose_cmat(trial_wf_up,trial_wf_up,trial_density_up,ist.n_spatial_orbitals,ist.n_up,ist.n_spatial_orbitals);
   cmat_transpose_cmat(trial_wf_down,trial_wf_down,trial_density_down,ist.n_spatial_orbitals,ist.n_down,ist.n_spatial_orbitals);

   /*Calculate Energy Based upon Density Matrix*/
   trial_ke = compute_kinetic_energy_density_matrix_restricted(trial_density_up,trial_density_down,kinetic_original_matrix,ist);
   trial_pe = compute_potential_energy_density_matrix_restricted(trial_density_up,trial_density_down,potential_original_matrix,ist);
   trial_energy = trial_ke.real + trial_pe.real; 

free(trial_density_up);
free(trial_density_down);
free(One);
return(trial_energy);
}

/***************************************************************************************/

double get_trial_energy_density_fermions(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *trial_density_up,MKL_Complex16 *trial_density_down,int *neighbors,int *number_neighbors,int_st ist,cns_st cns) {

   /*Obtains the Trial Energy Based Upon the Trial Wavefunction for Use in Estimating Walker Weights*/
   int i, j; 
   MKL_Complex16 trial_ke, trial_pe;
   MKL_Complex16 density_prod; 
   double trial_energy;

   /*Calculate Density Matrix*/
   cmat_transpose_cmat(trial_wf_up,trial_wf_up,trial_density_up,ist.n_sites,ist.n_up,ist.n_sites);
   cmat_transpose_cmat(trial_wf_down,trial_wf_down,trial_density_down,ist.n_sites,ist.n_down,ist.n_sites);

   /*Calculate Energy Based upon Density Matrix*/
   trial_ke.real = trial_ke.imag = trial_pe.real = trial_pe.imag = 0.0;
   for ( i=0; i<ist.n_sites; i++) {
      density_prod = Cmul( trial_density_up[i*ist.n_sites+i], trial_density_down[i*ist.n_sites+i]);
      trial_pe = Cadd(trial_pe, density_prod);
      if ( ist.n_sites_one > 2 || ist.n_sites_two > 2 ) {
        for (j=0; j<number_neighbors[i]; j++) {
           trial_ke = Csub(trial_ke, Cadd(trial_density_up[i*ist.n_sites+neighbors[4*i+j]],trial_density_down[i*ist.n_sites+neighbors[4*i+j]]));
         }
      }
      else {
        for (j=0; j<number_neighbors[i]; j++) {
           trial_ke =Csub(trial_ke, RCmul(.5, Cadd(trial_density_up[i*ist.n_sites+neighbors[4*i+j]],trial_density_down[i*ist.n_sites+neighbors[4*i+j]])));
        }
      }
    }

    trial_ke = RCmul(cns.t_fermions, trial_ke); 
    trial_pe = RCmul(cns.U_fermions, trial_pe); 
    trial_energy = trial_ke.real + trial_pe.real;  

free(trial_density_down); 
free(trial_density_up);  
return(trial_energy);
}

/**************************************************************************************/

double get_trial_energy_density_fermions_up(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_density_up,int *neighbors,int *number_neighbors,int_st ist,cns_st cns) {

   /*Obtains the Trial Energy Based Upon the Trial Wavefunction for Use in Estimating Walker Weights*/
   int i, j;
   MKL_Complex16 trial_ke; 
   double trial_energy;

   /*Calculate Density Matrix*/
   cmat_transpose_cmat(trial_wf_up,trial_wf_up,trial_density_up,ist.n_sites,ist.n_up,ist.n_sites);

   /*Calculate Energy Based upon Density Matrix*/
   trial_ke.real = trial_ke.imag = 0.0;
   for ( i=0; i<ist.n_sites; i++) {
      if ( ist.n_sites_one > 2 || ist.n_sites_two > 2 ) {
        for (j=0; j<number_neighbors[i]; j++) {
           trial_ke = Csub(trial_ke, trial_density_up[i*ist.n_sites+neighbors[4*i+j]]); 
         }
      }
      else {
        for (j=0; j<number_neighbors[i]; j++) {
           trial_ke =Csub(trial_ke, RCmul(.5, trial_density_up[i*ist.n_sites+neighbors[4*i+j]])); 
        }
      }
    }

    trial_ke = RCmul(cns.t_fermions, trial_ke);
    trial_energy = trial_ke.real; 

free(trial_density_up);  
return(trial_energy);
}

/**************************************************************************************/

double get_trial_energy_density_fermions_down(MKL_Complex16 *trial_wf_down,MKL_Complex16 *trial_density_down,int *neighbors,int *number_neighbors,int_st ist,cns_st cns) {

   /*Obtains the Trial Energy Based Upon the Trial Wavefunction for Use in Estimating Walker Weights*/
   int i, j;
   MKL_Complex16 trial_ke; 
   double trial_energy;

   /*Calculate Density Matrix*/
   cmat_transpose_cmat(trial_wf_down,trial_wf_down,trial_density_down,ist.n_sites,ist.n_down,ist.n_sites);

   /*Calculate Energy Based upon Density Matrix*/
   trial_ke.real = trial_ke.imag = 0.0;
   for ( i=0; i<ist.n_sites; i++) {
      if ( ist.n_sites_one > 2 || ist.n_sites_two > 2 ) {
        for (j=0; j<number_neighbors[i]; j++) {
           trial_ke = Csub(trial_ke, trial_density_down[i*ist.n_sites+neighbors[4*i+j]]);
         }
      }
      else {
        for (j=0; j<number_neighbors[i]; j++) {
           trial_ke =Csub(trial_ke, RCmul(.5, trial_density_down[i*ist.n_sites+neighbors[4*i+j]]));
        }
      }
    }

    trial_ke = RCmul(cns.t_fermions, trial_ke);
    trial_energy = trial_ke.real; 

free(trial_density_down);
return(trial_energy);
}

/***************************************************************************************/
double get_trial_energy_density_unrestricted(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *trial_density_up_total,MKL_Complex16 *trial_density_down_total,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *kinetic_original_matrix_down, MKL_Complex16 *potential_original_matrix_up, MKL_Complex16 *potential_original_matrix_down, MKL_Complex16 *potential_original_matrix_updown, MKL_Complex16 *energy_shifts, int_st ist) {

   /*Obtains the Trial Energy Based Upon the Trial Wavefunction for Use in Estimating Walker Weights*/
   int i, j; 
   MKL_Complex16 *One;
   MKL_Complex16 *trial_density_up, *trial_density_down; 
   MKL_Complex16 trial_ke, trial_pe;
   double trial_energy;

   trial_density_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16)); 
   trial_density_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16)); 

   czero_vec(trial_density_up_total,ist.n_spatial_orbitals_sq); 
   czero_vec(trial_density_down_total,ist.n_spatial_orbitals_sq); 
   for (i=0; i<ist.n_determinants_trial_energy; i++) {
     for (j=0; j<ist.n_determinants_trial_energy; j++) { 

      /*Calculate Density Matrix*/
      cmat_transpose_cmat(&trial_wf_up[i*ist.n_spatial_orbitals_n_up],&trial_wf_up[j*ist.n_spatial_orbitals_n_up],trial_density_up,ist.n_spatial_orbitals,ist.n_up,ist.n_spatial_orbitals);
      cmat_transpose_cmat(&trial_wf_down[i*ist.n_spatial_orbitals_n_down],&trial_wf_down[j*ist.n_spatial_orbitals_n_down],trial_density_down,ist.n_spatial_orbitals,ist.n_down,ist.n_spatial_orbitals);

      /*Add To Total Density Matrix*/
      cblas_zaxpy(ist.n_spatial_orbitals_sq,One,trial_density_up,1,trial_density_up_total,1);
      cblas_zaxpy(ist.n_spatial_orbitals_sq,One,trial_density_down,1,trial_density_down_total,1);

     }
   }

   /*Calculate Energy Based upon Density Matrix*/
   trial_ke = compute_kinetic_energy_density_matrix_unrestricted(trial_density_up_total,trial_density_down_total,kinetic_original_matrix_up,kinetic_original_matrix_down,ist);
   trial_pe = compute_potential_energy_density_matrix_unrestricted(trial_density_up_total,trial_density_down_total,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,ist);
   trial_energy = trial_ke.real + trial_pe.real ; //+ energy_shifts[0].r + energy_shifts[1].r + energy_shifts[2].r + energy_shifts[3].r; 

free(One); 
free(trial_density_up); 
free(trial_density_down); 
return(trial_energy);
}

/**************************************************************************************/

void get_mean_field_densities(MKL_Complex16 *trial_density_up,MKL_Complex16 *trial_density_down,MKL_Complex16 *mean_field_density_first_transform,MKL_Complex16 *mean_field_density_second_transform,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_eigs_first_transform,MKL_Complex16 *potential_eigs_second_transform,MKL_Complex16 *mean_field_first_transform_constant, MKL_Complex16 *mean_field_second_transform_constant,int_st ist){

    /*Forms the Mean Field Operators Used in the HS Transforms If Mean Field Portion Is Subtracted Out*/
    int gamma, i, l, j;   
    int i_spin, l_spin; 
    int i_orbital, l_orbital;  
    long int tidum; 

    Randomize(); tidum = -random();

    /*Zero the Mean Field Density Matrices*/
    czero_vec(mean_field_density_first_transform,ist.n_spin_orbitals_sq); 
    czero_vec(mean_field_density_second_transform,ist.n_spin_orbitals_sq); 

    /*Run Through Sites*/
    for (gamma=0; gamma<ist.n_spin_orbitals_sq; gamma++) {

           /*Zero These Totals*/
           mean_field_density_first_transform[gamma].real = mean_field_density_first_transform[gamma].imag = 0.0; 
           mean_field_density_second_transform[gamma].real = mean_field_density_second_transform[gamma].imag = 0.0; 

           for (i=0; i<ist.n_spin_orbitals; i++) {
            i_spin = i%2;
            i_orbital = (int)(i/2.0);

            for (l=0; l<ist.n_spin_orbitals; l++) {
             l_spin = l%2;
             l_orbital = (int)(l/2.0);

             /*Only Works If Two Spins are the Same*/
             if ( i_spin == l_spin ) {

                 /*If Spin Up*/
                 if ( i_spin == 0 ) {

                   /*Get From First Transform For Up*/
                   mean_field_density_first_transform[gamma] = Cadd(mean_field_density_first_transform[gamma], Cmul(potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma], trial_density_up[i_orbital*ist.n_spatial_orbitals+l_orbital])); 

                   mean_field_density_first_transform[gamma] = Cadd(mean_field_density_first_transform[gamma], Cmul(conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]), trial_density_up[l_orbital*ist.n_spatial_orbitals+i_orbital])); 

                   /*Get Second Transform for Up*/
                   mean_field_density_second_transform[gamma] = Cadd(mean_field_density_second_transform[gamma], Cmul(potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma], trial_density_up[i_orbital*ist.n_spatial_orbitals+l_orbital]));
                   mean_field_density_second_transform[gamma] = Csub(mean_field_density_second_transform[gamma], Cmul(conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]), trial_density_up[l_orbital*ist.n_spatial_orbitals+i_orbital]));

               } /*Spin Up*/
               else { /*If Spin Down*/

                   /*Get From First Transform For Down*/
                   mean_field_density_first_transform[gamma] = Cadd(mean_field_density_first_transform[gamma], Cmul(potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma], trial_density_down[i_orbital*ist.n_spatial_orbitals+l_orbital]));
                   mean_field_density_first_transform[gamma] = Cadd(mean_field_density_first_transform[gamma], Cmul(conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]), trial_density_down[l_orbital*ist.n_spatial_orbitals+i_orbital]));

                   /*Get Second Transform for Down*/                   
                   mean_field_density_second_transform[gamma] = Cadd(mean_field_density_second_transform[gamma], Cmul(potential_eigvecs_fermions[i*ist.n_spin_orbitals_third + l * ist.n_spin_orbitals_sq + gamma], trial_density_down[i_orbital*ist.n_spatial_orbitals+l_orbital]));
                   mean_field_density_second_transform[gamma] = Csub(mean_field_density_second_transform[gamma], Cmul(conjugate(potential_eigvecs_fermions[l*ist.n_spin_orbitals_third + i * ist.n_spin_orbitals_sq + gamma]), trial_density_down[l_orbital*ist.n_spatial_orbitals+i_orbital]));

               } /*Spin Down*/

           } /*Equal Spins*/
 

         } /*l*/
        } /*i*/



        /*Get Mean Field Constants*/
        if ( Cabs(mean_field_density_first_transform[gamma]) > .00001 ) {
          mean_field_first_transform_constant[gamma] = RCmul(-1,Cmul(potential_eigs_first_transform[gamma], mean_field_density_first_transform[gamma]));
        }
        else {
          mean_field_first_transform_constant[gamma].real = mean_field_first_transform_constant[gamma].imag = 0.0; 
        }

        if ( Cabs(mean_field_density_second_transform[gamma]) > .00001 ) {  
          mean_field_second_transform_constant[gamma] = RCmul(-1,Cmul(potential_eigs_second_transform[gamma], mean_field_density_second_transform[gamma])); 
        }
        else {
          mean_field_second_transform_constant[gamma].real = mean_field_second_transform_constant[gamma].imag = 0.0; 
        } 

       } /*Gamma*/

return; 
}

/*****************************************************************************************/

void init_wf_fermions_spin_polarized(MKL_Complex16 *trial_wf_up,double *kinetic_eigs_fermions,double *kinetic_eigvecs_fermions,int_st ist,cns_st cns) {

  /*Determine Which Type of Trial WF To Form*/
  if ( ist.flag_trial == 0 ) {
    /*Just an Identity Matrix*/

    /*Be Careful, Because May Not Reach Correct Steady State Because Could Be Orthogonal!!!!*/
    trial_identity_fermions_spin_polarized(trial_wf_up,ist);

  }
  else if ( ist.flag_trial == 1 ) {
    /*The Free Case*/

    trial_free_fermions_spin_polarized(trial_wf_up,kinetic_eigs_fermions,kinetic_eigvecs_fermions,ist,cns);

  }
  else if ( ist.flag_trial == 2 ) {
     /*The Restricted Hartree Fock Case*/

   trial_rhf_fermions_spin_polarized(trial_wf_up,ist);  

  }
  else if ( ist.flag_trial == 3 ) {
      /*Just a Random WF*/

     trial_random_fermions_spin_polarized(trial_wf_up,ist);

  }
 else {
     /*Put In Close to Exact Answer*/

     trial_wf_up[0].real = .666215;
     trial_wf_up[0].imag = 0.0;

     trial_wf_up[1].real = .529221;
     trial_wf_up[1].imag = 0.0;

     trial_wf_up[2].real = .427738;
     trial_wf_up[2].imag = 0.0;

  }

return;
}

/**************************************************************************************************************/

MKL_Complex16 compute_overlap_inverse_multi_energy(MKL_Complex16 *trial_wf_up_energy, MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,int_st ist) {

     MKL_Complex16 total_overlap; 

     if ( ist.n_up > 0 ) {
       if ( ist.n_down > 0 ) {
           total_overlap = compute_overlap_inverse_multi_energy_both(trial_wf_up_energy,trial_wf_down_energy,wf_up,wf_down,trial_determinant_coefficients_energy,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,ist); 
       }
       else {
            total_overlap = compute_overlap_inverse_multi_energy_up(trial_wf_up_energy,wf_up,trial_determinant_coefficients_energy,overlap_inverse_up,det_overlap_up,ist);  
       }
    }
    else {
      total_overlap = compute_overlap_inverse_multi_energy_down(trial_wf_down_energy,wf_down,trial_determinant_coefficients_energy,overlap_inverse_down,det_overlap_down,ist);
   }

return total_overlap; 
}

/************************************************************************************************************/

MKL_Complex16 compute_overlap_inverse_multi_phaseless(MKL_Complex16 *trial_wf_up_phaseless, MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,int_st ist) {

     MKL_Complex16 total_overlap;
     FILE *pf = fopen("errors.dat", "a+"); 

     if ( ist.n_up > 0 ) {
       if ( ist.n_down > 0 ) {
   //fprintf(pf, "in first up and down\n"); fflush(pf); 
           total_overlap = compute_overlap_inverse_multi_phaseless_both(trial_wf_up_phaseless,trial_wf_down_phaseless,wf_up,wf_down,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,ist);
   //fprintf(pf, "after in first\n"); fflush(pf); 
       }
       else {
            total_overlap = compute_overlap_inverse_multi_phaseless_up(trial_wf_up_phaseless,wf_up,trial_determinant_coefficients_phaseless,overlap_inverse_up,det_overlap_up,ist);
       }
    }
    else {
      total_overlap = compute_overlap_inverse_multi_phaseless_down(trial_wf_down_phaseless,wf_down,trial_determinant_coefficients_phaseless,overlap_inverse_down,det_overlap_down,ist);
   }

return total_overlap;
}

/************************************************************************************************************/

MKL_Complex16 compute_overlap_inverse_multi_energy_both(MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,int_st ist) {

    int i, j, k;
    MKL_Complex16 total_overlap;
    MKL_Complex16 *Zero, *One;
    MKL_Complex16 *stored_product_up, *stored_product_down;

    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero[0].real = Zero[0].imag = 0.0; 
    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;

    /*If Up Overlaps*******************************************************/
    stored_product_up = (MKL_Complex16 *)calloc(ist.n_up_sq,sizeof(MKL_Complex16));
    stored_product_down = (MKL_Complex16 *)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

    /*Set Total Overlap To Zero*/
    total_overlap = Zero[0];

    /*Build Up Overlap Matrix for All Determinants*/
    for (i=0; i<ist.n_determinants_trial_energy; i++) {

        cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_spatial_orbitals,One,&trial_wf_up_energy[i*ist.n_spatial_orbitals_n_up],ist.n_up,wf_up,ist.n_up,Zero,stored_product_up,ist.n_up);
        det_overlap_up[i]=complex_inverse_det_fermions(stored_product_up,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up);

        cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_spatial_orbitals,One,&trial_wf_down_energy[i*ist.n_spatial_orbitals_n_down],ist.n_down,wf_down,ist.n_down,Zero,stored_product_down,ist.n_down);
        det_overlap_down[i]=complex_inverse_det_fermions(stored_product_down,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down);  

        total_overlap = Cadd(total_overlap, Cmul(conjugate(trial_determinant_coefficients_energy[i]), Cmul(det_overlap_up[i],det_overlap_down[i])));
     }


//free(stored_product_up);
//free(stored_product_down); 
return(total_overlap);
}

/************************************************************************************************************/

MKL_Complex16 compute_overlap_inverse_multi_phaseless_both(MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,int_st ist) {

    int i, j, k;
    MKL_Complex16 total_overlap;
    MKL_Complex16 *Zero, *One;
    MKL_Complex16 *stored_product_up, *stored_product_down;
    FILE *pf = fopen("errors.dat", "a+"); 

    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero[0].real = Zero[0].imag = 0.0; 
    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;

    /*If Up Overlaps*******************************************************/
    stored_product_up = (MKL_Complex16 *)calloc(ist.n_up_sq,sizeof(MKL_Complex16));
    stored_product_down = (MKL_Complex16 *)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

    /*Set Total Overlap To Zero*/
    total_overlap = Zero[0];

    /*Build Up Overlap Matrix for All Determinants*/
    for (i=0; i<ist.n_determinants_trial_phaseless; i++) {

        cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_spatial_orbitals,One,&trial_wf_up_phaseless[i*ist.n_spatial_orbitals_n_up],ist.n_up,wf_up,ist.n_up,Zero,stored_product_up,ist.n_up);
        det_overlap_up[i]=complex_inverse_det_fermions(stored_product_up,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up);

        cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_spatial_orbitals,One,&trial_wf_down_phaseless[i*ist.n_spatial_orbitals_n_down],ist.n_down,wf_down,ist.n_down,Zero,stored_product_down,ist.n_down);
        det_overlap_down[i]=complex_inverse_det_fermions(stored_product_down,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down);

        total_overlap = Cadd(total_overlap, Cmul(conjugate(trial_determinant_coefficients_phaseless[i]), Cmul(det_overlap_up[i],det_overlap_down[i])));
     }


//free(stored_product_up);
////free(stored_product_down); 
return(total_overlap);
}

/***********************************************************************************************************/

MKL_Complex16 compute_overlap_inverse_multi_phaseless_up(MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *wf_up,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *det_overlap_up,int_st ist) {

    int i, j, k;
    MKL_Complex16 total_overlap;
    MKL_Complex16 *Zero, *One;
    MKL_Complex16 *stored_product_up; 

    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero[0].real = Zero[0].imag = 0.0;
    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;

    /*If Up Overlaps*******************************************************/
    stored_product_up = (MKL_Complex16 *)calloc(ist.n_up_sq,sizeof(MKL_Complex16));

    /*Set Total Overlap To Zero*/
    total_overlap = Zero[0];

    /*Build Up Overlap Matrix for All Determinants*/
    for (i=0; i<ist.n_determinants_trial_phaseless; i++) {
        cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_spatial_orbitals,One,&trial_wf_up_phaseless[i*ist.n_spatial_orbitals_n_up],ist.n_up,wf_up,ist.n_up,Zero,stored_product_up,ist.n_up);
        det_overlap_up[i]=complex_inverse_det_fermions(stored_product_up,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up);

        total_overlap = Cadd(total_overlap, Cmul(conjugate(trial_determinant_coefficients_phaseless[i]), det_overlap_up[i])); 
     }

free(stored_product_up);
return(total_overlap);
}

/***********************************************************************************************************/

MKL_Complex16 compute_overlap_inverse_multi_energy_up(MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *wf_up,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *det_overlap_up,int_st ist) {

    int i, j, k;
    MKL_Complex16 total_overlap;
    MKL_Complex16 *Zero, *One;
    MKL_Complex16 *stored_product_up;

    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero[0].real = Zero[0].imag = 0.0;
    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;

    /*If Up Overlaps*******************************************************/
    stored_product_up = (MKL_Complex16 *)calloc(ist.n_up_sq,sizeof(MKL_Complex16));

    /*Set Total Overlap To Zero*/
    total_overlap = Zero[0];

    /*Build Up Overlap Matrix for All Determinants*/
    for (i=0; i<ist.n_determinants_trial_energy; i++) {
        cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_spatial_orbitals,One,&trial_wf_up_energy[i*ist.n_spatial_orbitals_n_up],ist.n_up,wf_up,ist.n_up,Zero,stored_product_up,ist.n_up);
        det_overlap_up[i]=complex_inverse_det_fermions(stored_product_up,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up);

        total_overlap = Cadd(total_overlap, Cmul(conjugate(trial_determinant_coefficients_energy[i]), det_overlap_up[i]));
     }

free(stored_product_up);
return(total_overlap);
}

/******************************************************************************************************************************/

MKL_Complex16 compute_overlap_inverse_multi_energy_down(MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *wf_down,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_down,int_st ist) {

    int i, j, k;
    MKL_Complex16 total_overlap;
    MKL_Complex16 *Zero, *One;
    MKL_Complex16 *stored_product_down;

    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero[0].real = Zero[0].imag = 0.0;
    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;

    /*If Up Overlaps*******************************************************/
    stored_product_down = (MKL_Complex16 *)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

    /*Set Total Overlap To Zero*/
    total_overlap = Zero[0];

    /*Build Up Overlap Matrix for All Determinants*/
    for (i=0; i<ist.n_determinants_trial_energy; i++) {
        cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_spatial_orbitals,One,&trial_wf_down_energy[i*ist.n_spatial_orbitals_n_down],ist.n_down,wf_down,ist.n_down,Zero,stored_product_down,ist.n_down);
        det_overlap_down[i]=complex_inverse_det_fermions(stored_product_down,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down);

        total_overlap = Cadd(total_overlap, Cmul(conjugate(trial_determinant_coefficients_energy[i]), det_overlap_down[i]));
     }

free(stored_product_down);
return(total_overlap);
}

/******************************************************************************************************************************/

MKL_Complex16 compute_overlap_inverse_multi_phaseless_down(MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *wf_down,MKL_Complex16 *trial_determinant_coefficients_phaseless,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *det_overlap_down,int_st ist) {

    int i, j, k;
    MKL_Complex16 total_overlap;
    MKL_Complex16 *Zero, *One;
    MKL_Complex16 *stored_product_down;

    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero[0].real = Zero[0].imag = 0.0;
    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;

    /*If Up Overlaps*******************************************************/
    stored_product_down = (MKL_Complex16 *)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

    /*Set Total Overlap To Zero*/
    total_overlap = Zero[0];

    /*Build Up Overlap Matrix for All Determinants*/
    for (i=0; i<ist.n_determinants_trial_phaseless; i++) {
        cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_spatial_orbitals,One,&trial_wf_down_phaseless[i*ist.n_spatial_orbitals_n_down],ist.n_down,wf_down,ist.n_down,Zero,stored_product_down,ist.n_down);
        det_overlap_down[i]=complex_inverse_det_fermions(stored_product_down,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down);

        total_overlap = Cadd(total_overlap, Cmul(conjugate(trial_determinant_coefficients_phaseless[i]), det_overlap_down[i]));
     }

free(stored_product_down);
return(total_overlap);
}

/************************************************************************************************/

MKL_Complex16 compute_overlap_inverse(MKL_Complex16 *trial_wf, MKL_Complex16 *wf,MKL_Complex16 *overlap_inverse,int_st ist,int flag_up_down) {
      
    int i;
    MKL_Complex16 overlap; 
    MKL_Complex16 *Zero, *One;
    MKL_Complex16 *stored_product_up, *stored_product_down;
    MKL_Complex16 *stored_product_up_total, *stored_product_down_total;
    
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero[0].real = 0.0; Zero[0].imag = 0.0;
      
    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
      
    /*If Up Overlaps*******************************************************/
    if ( flag_up_down == 0 ) {
      
      stored_product_up = (MKL_Complex16 *)calloc(ist.n_up_sq,sizeof(MKL_Complex16));
      stored_product_up_total = (MKL_Complex16 *)calloc(ist.n_up_sq,sizeof(MKL_Complex16));
        
      /*Zero Total Overlap Product*/
      czero_vec(stored_product_up_total,ist.n_up_sq);
      
      /*Build Up Overlap Matrix for All Determinants*/
      for (i=0; i<ist.n_determinants_trial_energy; i++) {
        cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_spatial_orbitals,One,&trial_wf[i*ist.n_spatial_orbitals_n_up],ist.n_up,wf,ist.n_up,Zero,stored_product_up,ist.n_up);
        cblas_zaxpy(ist.n_up_sq,One,stored_product_up,1,stored_product_up_total,1);
      }

      overlap=complex_inverse_det_fermions(stored_product_up_total,overlap_inverse,ist.n_up);

      free(stored_product_up);
      free(stored_product_up_total);  
    }
    else {

      stored_product_down = (MKL_Complex16 *)calloc(ist.n_down_sq,sizeof(MKL_Complex16));
      stored_product_down_total = (MKL_Complex16 *)calloc(ist.n_down_sq,sizeof(MKL_Complex16));   

      /*Zero Total Overlap Product*/
      czero_vec(stored_product_down_total,ist.n_down_sq);

      /*Build Up Overlap Matrix for All Determinants*/
      for (i=0; i<ist.n_determinants_trial_energy; i++) {
        cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_spatial_orbitals,One,&trial_wf[i*ist.n_spatial_orbitals_n_down],ist.n_down,wf,ist.n_down,Zero,stored_product_down,ist.n_down);
        cblas_zaxpy(ist.n_down_sq,One,stored_product_down,1,stored_product_down_total,1);  
      }

      overlap=complex_inverse_det_fermions(stored_product_down_total,overlap_inverse,ist.n_down);

      free(stored_product_down); 
      free(stored_product_down_total); 
    } 

free(One); 
free(Zero); 
return(overlap); 
}

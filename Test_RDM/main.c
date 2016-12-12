/*This Ig an Example  Ground-State AFQMC Code Meant To Be Used on the Hubbard Model. Many Coding Aspects of It Can Be Improved; It is Meant to Exemplify the Physics. Started on September 18, 2013*/  

#include "afqmc.h"

int main() {

  FILE *pout; 
  FILE *pf = fopen("energy.dat", "a+"); 
  FILE *pf9 = fopen("energy_only.dat", "a+");  
  FILE *pf2 = fopen("average_density_matrix_up.dat", "a+"); 
  FILE *pf3 = fopen("average_density_matrix_down.dat", "a+"); 
  FILE *pf4 = fopen("average_two_body_density_matrix.dat", "a+");  
  FILE *pf5 = fopen("average_wave_function_up.dat", "a+"); 
  FILE *pf6 = fopen("average_wave_function_down.dat", "a+");  
  FILE *pf7 = fopen("overlaps.dat", "a+"); 
  FILE *pf25 = fopen("checkoverlap.dat", "a+");
  FILE *pf26 = fopen("checkweights.dat", "a+"); 
  FILE *pf27 = fopen("checkdensitymatrix_up.dat", "a+"); 
  FILE *pf28 = fopen("checkdensitymatrix_down.dat", "a+");    
  FILE *pf29 = fopen("checkwavefunction_up.dat", "a+"); 
  FILE *pf30 = fopen("checkwavefunction_down.dat", "a+");  
  FILE *pf31 = fopen("checkkineticphases.dat", "a+");  
  FILE *pf32 = fopen("checkpotentialphases.dat", "a+");  
  FILE *pf33 = fopen("check_flags.dat", "a+"); 
  int steps; 
  int i, j, k, l; 
  int *neighbors, *number_neighbors; 
  int *walkers_big, *walkers_small, *walkers_copied; 
  int *kinetic_ij_sparse_up, *kinetic_ij_sparse_down; 
  int *potential_ij_sparse_up, *potential_ij_sparse_down, *potential_ij_sparse_updown; 
  int number_kinetic_sparse_up, number_kinetic_sparse_down, number_potential_sparse_up, number_potential_sparse_down, number_potential_sparse_updown; 
  int *list_orbitals_trial_energy;  
  MKL_Complex16 max_overlap, min_overlap, max_weight, min_weight; 
  double max_overlap_real, min_overlap_real, max_weight_real, min_weight_real; 
  double *kinetic_matrix_sparse_up, *kinetic_matrix_sparse_down;
  double *potential_matrix_sparse_up, *potential_matrix_sparse_down, *potential_matrix_sparse_updown; 
  MKL_Complex16 accumulated_energy, accumulated_weight; 
  MKL_Complex16 ke_mixed, ke_mixed_error, pe_mixed, pe_mixed_error, te_mixed, te_mixed_error, ce, tw_mixed;
  MKL_Complex16 ke_exact, ke_exact_error, pe_exact, pe_exact_error, te_exact, te_exact_error, tw_exact;
  MKL_Complex16 ke_bp, ke_bp_error, pe_bp, pe_bp_error, te_bp, te_bp_error, tw_bp; 
  int *local_energy_flag, *total_energy_flag, *field_flag;  
  MKL_Complex16 *weights, *phases, *phases2, *energies, *energies2, *energies3, *overlaps, *density_matrices_up, *density_matrices_down,*overlap_up, *overlap_down, *overlap_bosons, *overlap_total;
  MKL_Complex16 *det_overlap_up, *det_overlap_down;
  MKL_Complex16 *overlap_inverse_up, *overlap_inverse_down; 
  MKL_Complex16 *trial_determinant_coefficients_energy, *trial_determinant_coefficients_phaseless;  
  MKL_Complex16 *wf_up, *wf_down, *wf_bosons;
  MKL_Complex16 *bp_wf_up, *bp_wf_down;  
  MKL_Complex16 *prev_wf_up, *prev_wf_down;  
  MKL_Complex16 *average_wf_up, *average_wf_down; 
  MKL_Complex16 *average_density_matrix_up, *average_density_matrix_down; 
  MKL_Complex16 *average_two_body_density_matrix;
  MKL_Complex16 *exact_two_body_density_matrix; 
  MKL_Complex16 *average_wave_function_up, *average_wave_function_down;  
  MKL_Complex16 *trial_wf_up_energy, *trial_wf_down_energy, *trial_wf_bosons;
  MKL_Complex16 *trial_wf_up_phaseless, *trial_wf_down_phaseless;  
  MKL_Complex16 *trial_density_up, *trial_density_down; 
  MKL_Complex16 *new_walker_weights, *new_overlap_up, *new_overlap_down, *new_overlap_bosons, *new_overlap_total;
  MKL_Complex16 *new_det_overlap_up, *new_det_overlap_down; 
  MKL_Complex16 *new_overlap_inverse_up, *new_overlap_inverse_down; 
  MKL_Complex16 *new_wf_up, *new_wf_down, *new_wf_bosons;
  MKL_Complex16 *potential_supermatrix_fermions; 
  MKL_Complex16 *potential_matrices_back_propagation_up, *potential_matrices_back_propagation_down; 
  MKL_Complex16 *energy_shifts;
  MKL_Complex16 *kinetic_matrix_original_fermions_up, *kinetic_matrix_original_fermions_down; 
  MKL_Complex16 *potential_matrix_original_fermions_up, *potential_matrix_original_fermions_down, *potential_matrix_original_fermions_updown; 
  MKL_Complex16 *potential_matrix_original_fermions_spin; 
  MKL_Complex16 *potential_eigs_fermions, *potential_eigvecs_fermions; 
  MKL_Complex16 *potential_eigs_fermions_first_transform, *potential_eigs_fermions_second_transform;   
  MKL_Complex16 *mean_field_density_first_transform, *mean_field_density_second_transform; 
  MKL_Complex16 *mean_field_first_transform_constant, *mean_field_second_transform_constant; 
  double *kinetic_full_fermions_real_up, *kinetic_backwards_half_fermions_real_up, *kinetic_forwards_half_fermions_real_up; 
  double *kinetic_eigs_fermions_real_up, *kinetic_eigvecs_fermions_real_up;
  double *kinetic_full_fermions_real_down, *kinetic_backwards_half_fermions_real_down, *kinetic_forwards_half_fermions_real_down;
  double *kinetic_eigs_fermions_real_down, *kinetic_eigvecs_fermions_real_down;
  double *kinetic_full_bosons, *kinetic_backwards_half_bosons, *kinetic_forwards_half_bosons; 
  double *kinetic_eigs_bosons, *kinetic_eigvecs_bosons; 
  double *kinetic_matrix_bosons, *potential_matrix_bosons;
  double *weight_rescaling;  
  double *potential_one_body_matrix_up, *potential_one_body_matrix_down; 
  MKL_Complex16 *R_up, *R_down; 
  MKL_Complex16 *potential_energy, *kinetic_energy, *coupling_energy, *total_energy, *total_weights;  
  MKL_Complex16 *total_potential_energy, *total_kinetic_energy, *total_energy_2, *total_weight; 
  long idum; 
  cns_st cns; int_st ist;

  /*******************************************************************************************/

  pout = fopen("errors.dat", "a+"); 

  /********************************************************************************************/

  /*Obtain Parameters*/
  init(&ist,&cns); 
  Randomize(); idum = -random();

  /*******************************************************************************************/

  /*If a Pure Fermion Simulation**************************************************************/

  if ( ist.flag_simulation_type == 0 ) { /**************************************************************************************/
    neighbors=(int*)calloc(4*ist.n_sites,sizeof(int));
    number_neighbors=(int*)calloc(ist.n_sites,sizeof(int));

    weights=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    /*For Population Control*/
    new_walker_weights = (MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    weight_rescaling=(double*)calloc(10,sizeof(double));
    walkers_big = (int*)calloc(cns.max_number_walkers,sizeof(int));
    walkers_small = (int*)calloc(cns.max_number_walkers,sizeof(int));
    walkers_copied = (int*)calloc(cns.max_number_walkers,sizeof(int));

    potential_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    kinetic_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    total_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    total_weights=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));

    overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    overlap_down=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));
    overlap_inverse_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));

    wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_up,sizeof(MKL_Complex16));
    wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_down,sizeof(MKL_Complex16));

    trial_wf_up_energy=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));
    trial_wf_down_energy=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));
    trial_density_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
    trial_density_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
 
    kinetic_full_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_backwards_half_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_forwards_half_fermions_real_up=(double *)calloc(ist.n_sites_sq,sizeof(double)); 
 
    kinetic_eigs_fermions_real_up=(double*)calloc(ist.n_sites,sizeof(double));
    kinetic_eigvecs_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double)); 
 
    /*For Population Control*/
    new_wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_up,sizeof(MKL_Complex16));
    new_wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_down,sizeof(MKL_Complex16));
    new_overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));
    new_overlap_inverse_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));
    new_overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    new_overlap_down=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
 
    /*For Orthogonalization*/
    R_up=(MKL_Complex16 *)calloc(ist.n_up*ist.n_up,sizeof(MKL_Complex16));
    R_down=(MKL_Complex16 *)calloc(ist.n_down*ist.n_down,sizeof(MKL_Complex16));
 
    /******************************************************************************************/
 
    /*Obtain Matrices, Wavefunctions, and Initialize Walkers*/
    init_neighbors(neighbors,number_neighbors,ist); 
    init_kinetic_fermions(kinetic_full_fermions_real_up,kinetic_backwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_up,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,neighbors,number_neighbors,ist,cns);
    init_wf_fermions(trial_wf_up_energy,trial_wf_down_energy,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,ist,cns); 
    init_walkers_fermions(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,weights,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weight_rescaling,ist,cns); 

    cns.trial_energy = get_trial_energy_density_fermions(trial_wf_up_energy,trial_wf_down_energy,trial_density_up,trial_density_down,neighbors,number_neighbors,ist,cns); 
    cns.exp_trial_energy = exp(cns.dtau * cns.trial_energy);

    /*******************************************************************************************/

    /*First Regress the Walkers By Half a Kinetic Propagator*/
    propagate_half_backwards_kinetic_fermions(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);  

    /*Production Phase*************************************************************************/
    for (steps = 0; steps < ist.n_steps; steps++) {
      if ( steps%20== 0 ) {
        fprintf(pout, "steps: %d\n", steps); fflush(pout); 
      }

      /*First Propagate a Full Step Forwards*/
      propagate_forwards_kinetic_fermions(kinetic_full_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist);

      /*Then Propagate Forwards a Full U Step*/
      if ( cns.U_fermions != 0 ) {
        /*Use Discrete or Continuous Transform*/
        if ( ist.flag_discrete == 1 ) {
           //propagate_forwards_potential_discrete_fermions(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist,cns,&idum);
        }
        else {
          if ( ist.flag_phaseless == 0 ) {
            //propagate_forwards_potential_continuous_fermions(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist,cns,&idum); 
            propagate_forwards_potential_continuous_allshifts_fermions(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist,cns,&idum);
          }
          else {
             //propagate_forwards_potential_continuous_phaseless_allshifts_fermions(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist,cns,&idum);
          }
        } 
      }

      /*Normalize By Trial Energy*/
      normalize(weights,ist,cns); 

      /*Now Collect Energies and Other Observables*/
      if ( steps%ist.n_steps_energy==0 ) { 
        /*First Propagate One Half Step Forward*/
        propagate_half_forwards_kinetic_fermions(kinetic_forwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

        /*Compute Energ*/
        compute_energy_fermions(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,neighbors,number_neighbors,potential_energy,kinetic_energy,total_energy,total_weights,weight_rescaling,ist,cns,steps,&pe_mixed,&ke_mixed,&te_mixed,&tw_mixed);
        cns.exp_trial_energy = exp(cns.dtau * (pe_mixed.real + ke_mixed.real)); 
        
        /*Print Energies*/
        fprintf(pf, "%d\t\t  %d\t\t  %f\t\t   %f\t\t  %f\t\t  %f\t\t %f\n", steps, ist.n_walkers, steps*cns.dtau, ke_mixed.real, pe_mixed.real, te_mixed.real, tw_mixed.real); fflush(pf); 
        fprintf(pf9, "%f\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9); 

        /*Propagate Back Backwards*/
        propagate_half_backwards_kinetic_fermions(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist); 
      }   

      /*Orthogonalize Q and R if Necessary*/
      if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
        orthogonalize_fermions(wf_up,wf_down,overlap_up,overlap_down,weights,R_up,R_down,ist);   
      }

      /*Apply Population Control If Necessary*/
      if ( cns.U_fermions != 0 ) { 
        if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
          //ist.n_walkers = population_control_fermions(weights,wf_up,wf_down,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,new_walker_weights,new_wf_up,new_wf_down,new_overlap_inverse_up,new_overlap_inverse_down,new_overlap_up,new_overlap_down,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns); 
        } 
      }

    } /*Steps*/

    /*Now Compute Average Overall Energy*/
    compute_averaged_energy_fermions(potential_energy,kinetic_energy,total_energy,total_weights,ist,"final_average_energy.dat"); 

    /******************************************************************************************/ 

    free(weights);
    free(new_walker_weights);

    free(potential_energy); free(kinetic_energy); free(total_energy); free(total_weights);

    free(overlap_up); free(overlap_down); 
    free(overlap_inverse_up); free(overlap_inverse_down); 
    free(wf_up); free(wf_down); 
    free(trial_wf_up_energy); free(trial_wf_down_energy); 
    free(kinetic_full_fermions_real_up); free(kinetic_backwards_half_fermions_real_up); free(kinetic_forwards_half_fermions_real_up);  
    free(kinetic_eigs_fermions_real_up); free(kinetic_eigvecs_fermions_real_up);  

    free(new_overlap_up); free(new_overlap_down); 
    free(new_overlap_inverse_up); free(new_overlap_inverse_down); 
    free(new_wf_up); free(new_wf_down); 

    free(R_up); free(R_down); 

  } /*Fermion Simulation*****************************************************************************************************************************/
  else if ( ist.flag_simulation_type == 1 ) { /*If Bosons*********************************************************************************************/

    neighbors=(int*)calloc(4*ist.n_sites,sizeof(int));
    number_neighbors=(int*)calloc(ist.n_sites,sizeof(int));

    weights=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    /*For Population Control*/
    new_walker_weights = (MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    weight_rescaling=(double*)calloc(10,sizeof(double));
    walkers_big = (int*)calloc(cns.max_number_walkers,sizeof(int));
    walkers_small = (int*)calloc(cns.max_number_walkers,sizeof(int));
    walkers_copied = (int*)calloc(cns.max_number_walkers,sizeof(int));

    potential_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    kinetic_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    total_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    total_weights=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));

    overlap_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    wf_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites,sizeof(MKL_Complex16));
    trial_wf_bosons=(MKL_Complex16*)calloc(ist.n_sites,sizeof(MKL_Complex16));

    kinetic_matrix_bosons=(double *)calloc(ist.n_sites_sq,sizeof(double)); 
    potential_matrix_bosons=(double *)calloc(ist.n_sites_sq,sizeof(double)); 

    kinetic_full_bosons=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_backwards_half_bosons=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_forwards_half_bosons=(double *)calloc(ist.n_sites_sq,sizeof(double));

    kinetic_eigs_bosons=(double*)calloc(ist.n_sites,sizeof(double));
    kinetic_eigvecs_bosons=(double*)calloc(ist.n_sites_sq,sizeof(double));

    /*For Population Control*/
    new_wf_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites,sizeof(MKL_Complex16));
    new_overlap_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    /******************************************************************************************/

    /*Obtain Matrices, Wavefunctions, and Initialize Walkers*/
    init_neighbors(neighbors,number_neighbors,ist);

    init_kinetic_bosons(kinetic_full_bosons,kinetic_backwards_half_bosons,kinetic_forwards_half_bosons,kinetic_eigs_bosons,kinetic_eigvecs_bosons,kinetic_matrix_bosons,potential_matrix_bosons,neighbors,number_neighbors,ist,cns);
    init_wf_bosons(trial_wf_bosons,kinetic_eigs_bosons,kinetic_eigvecs_bosons,ist,cns);
    init_walkers_bosons(wf_bosons,trial_wf_bosons,weights,overlap_bosons,weight_rescaling,ist,cns);

    /*******************************************************************************************/

    /*First Regress the Walkers By Half a Kinetic Propagator*/
    propagate_half_backwards_kinetic_bosons(kinetic_backwards_half_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);

    /*Production Phase*************************************************************************/
    for (steps = 0; steps < ist.n_steps; steps++) {
      if ( steps%20== 0 ) {
        fprintf(pout, "steps: %d\n", steps); fflush(pout);
      }

      /*First Propagate a Full Step Forwards*/
      propagate_forwards_kinetic_bosons(kinetic_full_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);


      /*Then Propagate Forwards a Full U Step*/
      if ( cns.U_bosons != 0 ) {
        /*Use Discrete or Continuous Transform*/
        if ( ist.flag_phaseless == 0 ) {
            //propagate_forwards_potential_continuous_bosons(wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist,cns,&idum); 
        }
        else {
             //propagate_forwards_potential_continuous_phaseless_bosons(wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist,cns,&idum);
        }
      }

      /*Normalize By Trial Energy*/
      normalize(weights,ist,cns);

      /*Now Collect Energies and Other Observables*/
      if ( steps%ist.n_steps_energy==0 ) {
        /*First Propagate One Half Step Forward*/
        propagate_half_forwards_kinetic_bosons(kinetic_forwards_half_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);

        /*Compute Energ*/
        compute_energy_bosons(wf_bosons,trial_wf_bosons,weights,kinetic_matrix_bosons,potential_matrix_bosons,neighbors,number_neighbors,potential_energy,kinetic_energy,total_energy,total_weights,weight_rescaling,ist,cns,steps,&pe_mixed,&ke_mixed,&te_mixed,&tw_mixed);

        /*Print Desired Information*******************************************/
        fprintf(pf, "%d\t\t  %d\t\t  %f\t\t   %f\t\t  %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, ke_mixed.real, pe_mixed.real, te_mixed.real, tw_mixed.real); fflush(pf); 
        fprintf(pf9, "%f\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9);   

        /*Propagate Back Backwards*/
        propagate_half_backwards_kinetic_bosons(kinetic_backwards_half_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);
      }

      /*Normalize Wave Function If Necessary*/
      if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
        // normalize_bosons(wf_bosons,overlap_bosons,weights,ist); 
      } 

      /*Apply Population Control If Necessary*/
      if ( cns.U_bosons != 0 ) {
        if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
          //ist.n_walkers = population_control_bosons(weights,wf_bosons,overlap_bosons,new_walker_weights,new_wf_bosons,new_overlap_bosons,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns);
        }
      }

    } /*Steps*/

    /*Now Compute Average Overall Energy*/
    compute_averaged_energy_bosons(potential_energy,kinetic_energy,total_energy,total_weights,ist,"final_average_energy.dat");

    /******************************************************************************************/

    free(weights);
    free(new_walker_weights);

    free(potential_energy); free(kinetic_energy); free(total_energy); free(total_weights);
    free(overlap_bosons);
    free(wf_bosons);
    free(trial_wf_bosons);
    free(kinetic_full_bosons); free(kinetic_backwards_half_bosons); free(kinetic_forwards_half_bosons);
    free(kinetic_eigs_bosons); free(kinetic_eigvecs_bosons);
    free(kinetic_matrix_bosons); free(potential_matrix_bosons);  

    free(new_overlap_bosons);
    free(new_wf_bosons);

  } /*Boson Simulation******************************************************************************************************************************/
  else if ( ist.flag_simulation_type == 2 ) { /**Bose-Fermi Mixtures********************************************************************************/

    /*Declarations**************************************************/
    neighbors=(int*)calloc(4*ist.n_sites,sizeof(int));
    number_neighbors=(int*)calloc(ist.n_sites,sizeof(int));

    weights=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    /*For Population Control*/
    new_walker_weights = (MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    weight_rescaling=(double*)calloc(10,sizeof(double));
    walkers_big = (int*)calloc(cns.max_number_walkers,sizeof(int));
    walkers_small = (int*)calloc(cns.max_number_walkers,sizeof(int));
    walkers_copied = (int*)calloc(cns.max_number_walkers,sizeof(int));

    potential_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    kinetic_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    coupling_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    total_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    total_weights=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));

    kinetic_matrix_bosons=(double *)calloc(ist.n_sites_sq,sizeof(double));
    potential_matrix_bosons=(double *)calloc(ist.n_sites_sq,sizeof(double));

    overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    overlap_down=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    overlap_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16)); 

    overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));
    overlap_inverse_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));

    wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_up,sizeof(MKL_Complex16));
    wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_down,sizeof(MKL_Complex16));
    wf_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites,sizeof(MKL_Complex16));

    trial_wf_up_energy=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));
    trial_wf_down_energy=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));
    trial_wf_bosons=(MKL_Complex16*)calloc(ist.n_sites,sizeof(MKL_Complex16));

    kinetic_full_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_backwards_half_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_forwards_half_fermions_real_up=(double *)calloc(ist.n_sites_sq,sizeof(double));

    kinetic_eigs_fermions_real_up=(double*)calloc(ist.n_sites,sizeof(double));
    kinetic_eigvecs_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));

    kinetic_full_bosons=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_backwards_half_bosons=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_forwards_half_bosons=(double *)calloc(ist.n_sites_sq,sizeof(double));

    kinetic_eigs_bosons=(double*)calloc(ist.n_sites,sizeof(double));
    kinetic_eigvecs_bosons=(double*)calloc(ist.n_sites_sq,sizeof(double));

    /*For Population Control*/
    new_wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_up,sizeof(MKL_Complex16));
    new_wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_down,sizeof(MKL_Complex16));
    new_wf_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites,sizeof(MKL_Complex16));

    new_overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));
    new_overlap_inverse_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));

    new_overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    new_overlap_down=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    new_overlap_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16)); 
   
    /*For Orthogonalization*/
    R_up=(MKL_Complex16 *)calloc(ist.n_up*ist.n_up,sizeof(MKL_Complex16));
    R_down=(MKL_Complex16 *)calloc(ist.n_down*ist.n_down,sizeof(MKL_Complex16)); 

    /**************************************************************************************************************************/

    /*Initialize Walkers and Matrices*/
    init_neighbors(neighbors,number_neighbors,ist);

    init_kinetic_fermions(kinetic_full_fermions_real_up,kinetic_backwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_up,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,neighbors,number_neighbors,ist,cns);
    init_kinetic_bosons(kinetic_full_bosons,kinetic_backwards_half_bosons,kinetic_forwards_half_bosons,kinetic_eigs_bosons,kinetic_eigvecs_bosons,kinetic_matrix_bosons,potential_matrix_bosons,neighbors,number_neighbors,ist,cns); 

    init_wf_fermions(trial_wf_up_energy,trial_wf_down_energy,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,ist,cns);
    init_wf_bosons(trial_wf_bosons,kinetic_eigs_bosons,kinetic_eigvecs_bosons,ist,cns);  

    init_walkers_bose_fermi(wf_up,wf_down,wf_bosons,trial_wf_up_energy,trial_wf_down_energy,trial_wf_bosons,weights,overlap_up,overlap_down,overlap_bosons,weight_rescaling,ist,cns);

    /**********************************/
    
    /*First Regress the Walkers By Half a Kinetic Propagator*/
    propagate_half_backwards_kinetic_bosons(kinetic_backwards_half_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);
    propagate_half_backwards_kinetic_fermions(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist); 

    /*Production Phase*************************************************************************/
    for (steps = 0; steps < ist.n_steps; steps++) {
      if ( steps%20== 0 ) {
        fprintf(pout, "steps: %d\n", steps); fflush(pout);
      }

      /*First Propagate a Full Step Forwards*/
      propagate_forwards_kinetic_bosons(kinetic_full_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);
      propagate_forwards_kinetic_fermions(kinetic_full_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist); 

      /*Then Propagate Forwards a Full U Step*/
      if ( cns.U_bosons != 0 || cns.U_fermions != 0 || cns.U_bose_fermi != 0 ) {
        /*Use Discrete or Continuous Transform*/
        if ( ist.flag_phaseless == 0 ) {
          propagate_forwards_potential_continuous_bose_fermi(wf_up,wf_down,wf_bosons,trial_wf_up_energy,trial_wf_down_energy,trial_wf_bosons,overlap_up,overlap_down,overlap_bosons,overlap_inverse_up,overlap_inverse_down,weights,ist,cns,&idum); 
        }
        else {
          propagate_forwards_potential_continuous_phaseless_allshifts_bose_fermi(wf_up,wf_down,wf_bosons,trial_wf_up_energy,trial_wf_down_energy,trial_wf_bosons,overlap_up,overlap_down,overlap_bosons,overlap_inverse_up,overlap_inverse_down,weights,ist,cns,&idum);
        }
      }

      /*Normalize By Trial Energy*/
      normalize(weights,ist,cns);
      
      /*Now Collect Energies and Other Observables*/
      if ( steps%ist.n_steps_energy==0 ) {

        /*First Propagate One Half Step Forward*/
        propagate_half_forwards_kinetic_bosons(kinetic_forwards_half_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);
        propagate_half_forwards_kinetic_fermions(kinetic_forwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);      
 
        /*Compute Energ*/
        compute_energy_bose_fermi(wf_bosons,wf_up,wf_down,trial_wf_bosons,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_bosons,potential_matrix_bosons,neighbors,number_neighbors,potential_energy,kinetic_energy,coupling_energy,total_energy,total_weights,weight_rescaling,ist,cns,steps,&pe_mixed,&ke_mixed,&ce,&te_mixed,&tw_mixed);      

         /*Print Desired Information*******************************************/
         fprintf(pf, "%d\t %d %f\t\t\t %f\t %f\t %f\t %f\t %f\n", steps, ist.n_walkers, steps*cns.dtau, ke_mixed.real, pe_mixed.real, ce.real, te_mixed.real, tw_mixed.real); fflush(pf); 
         fprintf(pf9, "%f\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9); 


        /*Propagate Back Backwards*/
        propagate_half_backwards_kinetic_bosons(kinetic_backwards_half_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);
        propagate_half_backwards_kinetic_fermions(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist); 
      }
      
      /*Normalize Wave Function If Necessary*/
      if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
         normalize_bosons(wf_bosons,overlap_bosons,weights,ist);
         orthogonalize_fermions(wf_up,wf_down,overlap_up,overlap_down,weights,R_up,R_down,ist);
      }
    
         /*Apply Population Control If Necessary*/
      if ( cns.U_bosons != 0 || cns.U_fermions != 0 || cns.U_bose_fermi != 0) {
        if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
          ist.n_walkers = population_control_bose_fermi(weights,wf_up,wf_down,wf_bosons,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,overlap_bosons,new_walker_weights,new_wf_up,new_wf_down,new_wf_bosons,new_overlap_inverse_up,new_overlap_inverse_down,new_overlap_up,new_overlap_down,new_overlap_bosons,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns); 
        }
      }

    } /*Steps*/

    /*Now Compute Average Overall Energy*/
    compute_averaged_energy_bose_fermi(potential_energy,kinetic_energy,coupling_energy,total_energy,total_weights,ist,"final_average_energy.dat");
  
    /*Free the Matrices*/
    free(neighbors); free(number_neighbors); 

    free(weights); 
    free(new_walker_weights); free(weight_rescaling); 
    free(walkers_big); free(walkers_small); free(walkers_copied);  
 
    free(total_weights); free(total_energy); free(kinetic_energy); free(potential_energy); free(coupling_energy);  

    free(overlap_up); free(overlap_down); free(overlap_bosons);  
    free(overlap_inverse_up); free(overlap_inverse_down);

    free(wf_up); free(wf_down); free(wf_bosons); 
    free(trial_wf_up_energy); free(trial_wf_down_energy); free(trial_wf_bosons); 

    free(kinetic_full_fermions_real_up); free(kinetic_backwards_half_fermions_real_up); free(kinetic_forwards_half_fermions_real_up); 
    free(kinetic_eigs_fermions_real_up); free(kinetic_eigvecs_fermions_real_up); 

    free(kinetic_full_bosons); free(kinetic_backwards_half_bosons); free(kinetic_forwards_half_bosons); 
    free(kinetic_eigs_bosons); free(kinetic_eigvecs_bosons); 

    free(new_wf_up); free(new_wf_down); free(new_wf_bosons); 
    free(new_overlap_inverse_up); free(new_overlap_inverse_down); 
    free(new_overlap_up); free(new_overlap_down); free(new_overlap_bosons); 

    free(R_up); free(R_down); 


  } /*Bose Fermi simulation*************************************************************************************************************************/ 
  else if ( ist.flag_simulation_type == 3 ) { /*Bose Fermi Simulation Spin-Polarized*************************/

    /*Declarations**************************************************/
    neighbors=(int*)calloc(4*ist.n_sites,sizeof(int));
    number_neighbors=(int*)calloc(ist.n_sites,sizeof(int));

    weights=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    /*For Population Control*/
    new_walker_weights = (MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    weight_rescaling=(double*)calloc(10,sizeof(double));
    walkers_big = (int*)calloc(cns.max_number_walkers,sizeof(int));
    walkers_small = (int*)calloc(cns.max_number_walkers,sizeof(int));
    walkers_copied = (int*)calloc(cns.max_number_walkers,sizeof(int));

    potential_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    kinetic_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    coupling_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    total_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
    total_weights=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));

    kinetic_matrix_bosons=(double *)calloc(ist.n_sites_sq,sizeof(double));
    potential_matrix_bosons=(double *)calloc(ist.n_sites_sq,sizeof(double));

    overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    overlap_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));

    wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_up,sizeof(MKL_Complex16));
    wf_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites,sizeof(MKL_Complex16));

    trial_wf_up_energy=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));
    trial_wf_bosons=(MKL_Complex16*)calloc(ist.n_sites,sizeof(MKL_Complex16));

    kinetic_full_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_backwards_half_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_forwards_half_fermions_real_up=(double *)calloc(ist.n_sites_sq,sizeof(double));

    kinetic_eigs_fermions_real_up=(double*)calloc(ist.n_sites,sizeof(double));
    kinetic_eigvecs_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));

    kinetic_full_bosons=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_backwards_half_bosons=(double*)calloc(ist.n_sites_sq,sizeof(double));
    kinetic_forwards_half_bosons=(double *)calloc(ist.n_sites_sq,sizeof(double));

    kinetic_eigs_bosons=(double*)calloc(ist.n_sites,sizeof(double));
    kinetic_eigvecs_bosons=(double*)calloc(ist.n_sites_sq,sizeof(double));

    /*For Population Control*/
    new_wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_up,sizeof(MKL_Complex16));
    new_wf_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites,sizeof(MKL_Complex16));

    new_overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));

    new_overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    new_overlap_bosons=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    /*For Orthogonalization*/
    R_up=(MKL_Complex16 *)calloc(ist.n_up_sq,sizeof(MKL_Complex16));
 
    /*Initialize Walkers and Matrices*/
    init_neighbors(neighbors,number_neighbors,ist);

    init_kinetic_fermions(kinetic_full_fermions_real_up,kinetic_backwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_up,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,neighbors,number_neighbors,ist,cns);
    init_kinetic_bosons(kinetic_full_bosons,kinetic_backwards_half_bosons,kinetic_forwards_half_bosons,kinetic_eigs_bosons,kinetic_eigvecs_bosons,kinetic_matrix_bosons,potential_matrix_bosons,neighbors,number_neighbors,ist,cns);

    init_wf_fermions_spin_polarized(trial_wf_up_energy,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,ist,cns);
    init_wf_bosons(trial_wf_bosons,kinetic_eigs_bosons,kinetic_eigvecs_bosons,ist,cns);

    init_walkers_bose_fermi_spin_polarized(wf_up,wf_bosons,trial_wf_up_energy,trial_wf_bosons,weights,overlap_up,overlap_bosons,weight_rescaling,ist,cns);

    /**********************************/

    /*First Regress the Walkers By Half a Kinetic Propagator*/
    propagate_half_backwards_kinetic_bosons(kinetic_backwards_half_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);
    propagate_half_backwards_kinetic_fermions_spin_polarized(kinetic_backwards_half_fermions_real_up,wf_up,trial_wf_up_energy,overlap_up,overlap_inverse_up,weights,ist);

    /*Production Phase*************************************************************************/
    for (steps = 0; steps < ist.n_steps; steps++) {
      if ( steps%20== 0 ) {
        fprintf(pout, "steps: %d\n", steps); fflush(pout);
      }

      /*First Propagate a Full Step Forwards*/
      propagate_forwards_kinetic_bosons(kinetic_full_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);
      propagate_forwards_kinetic_fermions_spin_polarized(kinetic_full_fermions_real_up,wf_up,trial_wf_up_energy,overlap_up,overlap_inverse_up,weights,ist);

      /*Then Propagate Forwards a Full U Step*/
      if ( cns.U_bosons != 0 || cns.U_fermions != 0 || cns.U_bose_fermi != 0 ) {
        /*Use Discrete or Continuous Transform*/
        if ( ist.flag_phaseless == 0 ) {
            propagate_forwards_potential_continuous_meanfield_bose_fermi_spin_polarized(wf_up,wf_bosons,trial_wf_up_energy,trial_wf_bosons,overlap_up,overlap_bosons,overlap_inverse_up,weights,ist,cns,&idum); 
        }
        else {
            propagate_forwards_potential_continuous_phaseless_allshifts_bose_fermi_spin_polarized(wf_up,wf_bosons,trial_wf_up_energy,trial_wf_bosons,overlap_up,overlap_bosons,overlap_inverse_up,weights,potential_matrix_bosons,ist,cns,&idum); 
        }
      } 

      /*Normalize By Trial Energy*/
      normalize(weights,ist,cns);

      /*Now Collect Energies and Other Observables*/
      if ( steps%ist.n_steps_energy==0 ) {

        /*First Propagate One Half Step Forward*/
        propagate_half_forwards_kinetic_bosons(kinetic_forwards_half_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);
        propagate_half_forwards_kinetic_fermions_spin_polarized(kinetic_forwards_half_fermions_real_up,wf_up,trial_wf_up_energy,overlap_up,overlap_inverse_up,weights,ist);

        /*Normalize the Bosons Before Measuring*/
        normalize_bosons(wf_bosons,overlap_bosons,weights,ist);

        /*Compute Energ*/
        compute_energy_bose_fermi_spin_polarized(wf_bosons,wf_up,trial_wf_bosons,trial_wf_up_energy,overlap_inverse_up,weights,kinetic_matrix_bosons,potential_matrix_bosons,neighbors,number_neighbors,potential_energy,kinetic_energy,coupling_energy,total_energy,total_weights,weight_rescaling,ist,cns,steps,&pe_mixed,&ke_mixed,&ce,&te_mixed,&tw_mixed);

         /*Print Desired Information*******************************************/
         fprintf(pf, "%d\t %d %f\t\t\t %f\t %f\t %f\t %f\t %f\n", steps, ist.n_walkers, steps*cns.dtau, ke_mixed.real, pe_mixed.real, ce.real,te_mixed.real,tw_mixed.real); fflush(pf); 
         fprintf(pf9, "%f\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9);  
 

        /*Propagate Back Backwards*/
        propagate_half_backwards_kinetic_bosons(kinetic_backwards_half_bosons,wf_bosons,trial_wf_bosons,overlap_bosons,weights,ist);
        propagate_half_backwards_kinetic_fermions_spin_polarized(kinetic_backwards_half_fermions_real_up,wf_up,trial_wf_up_energy,overlap_up,overlap_inverse_up,weights,ist);

       }

      /*Normalize Wave Function If Necessary*/
      if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
         normalize_bosons(wf_bosons,overlap_bosons,weights,ist);
         orthogonalize_fermions_spin_polarized(wf_up,overlap_up,weights,R_up,ist);
      }

         /*Apply Population Control If Necessary*/
      if ( cns.U_bosons != 0 || cns.U_bose_fermi != 0 || cns.U_fermions != 0  ) {
        if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
           ist.n_walkers = population_control_bose_fermi_spin_polarized(weights,wf_up,wf_bosons,overlap_inverse_up,overlap_up,overlap_bosons,new_walker_weights,new_wf_up,new_wf_bosons,new_overlap_inverse_up,new_overlap_up,new_overlap_bosons,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns,steps);
        }
      }


    } /*Steps*/

    /*Now Compute Average Overall Energy*/
    compute_averaged_energy_bose_fermi(potential_energy,kinetic_energy,coupling_energy,total_energy,total_weights,ist,"final_average_energy.dat");

    /*Free the Matrices*/
    free(neighbors); free(number_neighbors);

    free(weights);
    free(new_walker_weights); free(weight_rescaling);
    free(walkers_big); free(walkers_small); free(walkers_copied);

    free(total_weights); free(total_energy); free(kinetic_energy); free(potential_energy); free(coupling_energy);

    free(overlap_up); free(overlap_bosons);
    free(overlap_inverse_up);

    free(wf_up); free(wf_bosons);
    free(trial_wf_up_energy); free(trial_wf_bosons);

    free(kinetic_full_fermions_real_up); free(kinetic_backwards_half_fermions_real_up); free(kinetic_forwards_half_fermions_real_up);
    free(kinetic_eigs_fermions_real_up); free(kinetic_eigvecs_fermions_real_up);

    free(kinetic_full_bosons); free(kinetic_backwards_half_bosons); free(kinetic_forwards_half_bosons);
    free(kinetic_eigs_bosons); free(kinetic_eigvecs_bosons);

    free(new_wf_up); free(new_wf_bosons);
    free(new_overlap_inverse_up);
    free(new_overlap_up); free(new_overlap_bosons);

    free(R_up); 

   } /*If Bose Fermi Spin Polarized******************************************************************************/
   else if ( ist.flag_simulation_type == 5 ) { 

    weights=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));

    /*For Population Control*/
    new_walker_weights = (MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));
    weight_rescaling=(double*)malloc(10*sizeof(double));
    walkers_big = (int*)malloc(cns.max_number_walkers*sizeof(int));
    walkers_small = (int*)malloc(cns.max_number_walkers*sizeof(int));
    walkers_copied = (int*)malloc(cns.max_number_walkers*sizeof(int));

    total_energy_2=(MKL_Complex16*)malloc(3*sizeof(MKL_Complex16)); 
    total_potential_energy=(MKL_Complex16*)malloc(3*sizeof(MKL_Complex16)); 
    total_kinetic_energy=(MKL_Complex16*)malloc(3*sizeof(MKL_Complex16));  
    total_weight=(MKL_Complex16*)malloc(1*sizeof(MKL_Complex16));  

    overlap_up=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));
    overlap_down=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));

    overlap_inverse_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_up_sq*sizeof(MKL_Complex16));
    overlap_inverse_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_down_sq*sizeof(MKL_Complex16));

    average_density_matrix_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16)); 
    average_density_matrix_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
    average_two_body_density_matrix=(MKL_Complex16*)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));  

    average_wave_function_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16)); 
    average_wave_function_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16)); 

    wf_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    wf_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

    trial_wf_up_energy=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    trial_wf_down_energy=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));
    trial_density_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16)); 
    trial_density_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16)); 

    kinetic_full_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
    kinetic_backwards_half_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
    kinetic_forwards_half_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));

    kinetic_eigs_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals*sizeof(double));
    kinetic_eigvecs_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));

    energy_shifts=(MKL_Complex16 *)malloc(4*sizeof(MKL_Complex16)); 

    kinetic_matrix_original_fermions_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

    potential_matrix_original_fermions_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals_fourth*sizeof(MKL_Complex16));
    potential_matrix_original_fermions_spin=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16)); 

    potential_supermatrix_fermions=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));
    potential_eigs_fermions=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
    potential_eigvecs_fermions=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));

    potential_eigs_fermions_first_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
    potential_eigs_fermions_second_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));

    potential_one_body_matrix_up = (double *)calloc(ist.n_spatial_orbitals_sq,sizeof(double)); 

    mean_field_density_first_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16)); 
    mean_field_density_second_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16)); 
   
    mean_field_first_transform_constant=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16)); 
    mean_field_second_transform_constant=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16)); 

   list_orbitals_trial_energy = (int *)malloc(ist.n_spatial_orbitals*sizeof(int));

    /*For Population Control*/
    new_wf_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
    new_wf_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));
    new_overlap_inverse_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_up_sq*sizeof(MKL_Complex16));
    new_overlap_inverse_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_down_sq*sizeof(MKL_Complex16));
    new_overlap_up=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));
    new_overlap_down=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));

    /*For Orthogonalization*/
    R_up=(MKL_Complex16 *)malloc(ist.n_up*ist.n_up*sizeof(MKL_Complex16));
    R_down=(MKL_Complex16 *)malloc(ist.n_down*ist.n_down*sizeof(MKL_Complex16));

    /*If There Is Back Propagation*/
    if ( ist.flag_back_propagation == 1 ) {
      potential_matrices_back_propagation_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_sq*ist.n_steps_back_propagation*sizeof(MKL_Complex16)); 
      potential_matrices_back_propagation_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_sq*ist.n_steps_back_propagation*sizeof(MKL_Complex16));

      bp_wf_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
      bp_wf_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

      prev_wf_up=(MKL_Complex16*)malloc(ist.n_steps_back_propagation*cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
      prev_wf_down=(MKL_Complex16*)malloc(ist.n_steps_back_propagation*cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));
    }

    /******************************************************************************************/

    /*Zero the Total Energy Vectors*/
    czero_vec(total_potential_energy,3); 
    czero_vec(total_kinetic_energy,3); 
    czero_vec(total_energy_2,3); 
    czero_vec(total_weight,1);

    /*Obtain Trial Wave Function, Trial Energy, and Trial Density Matrices from WF*/
    ist.n_max_orbital_trial_energy = init_wf_chemistry_restricted(trial_wf_up_energy,trial_wf_down_energy,kinetic_matrix_original_fermions_up,list_orbitals_trial_energy,ist,cns);

    /*Input the Vijkl Elements and Resort Them to Find Eigenvalues*/
    get_potential_matrix_elements_molpro_restricted(kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,potential_matrix_original_fermions_spin,energy_shifts,list_orbitals_trial_energy,ist,cns,"potential_matrix_elements.par");

    if ( ist.flag_restart == 0 ) {
      init_walkers_chemistry_no_restart_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,weights,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weight_rescaling,ist,cns);
    }
    else {
      init_walkers_chemistry_restart_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,weights,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weight_rescaling,ist,cns);
    } 

    /*Obtain Initial Trial Energy*/
    cns.trial_energy = get_trial_energy_density_restricted(trial_wf_up_energy,trial_wf_down_energy,trial_density_up,trial_density_down,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,ist);
    fprintf(pout, "trial energy %f\n", cns.trial_energy+energy_shifts[0].real); fflush(pout); 
    cns.init_trial_energy = cns.trial_energy; 
    cns.exp_trial_energy = exp(cns.dtau * cns.trial_energy);

    /*Form the Supermatrix and Diagonalize It to Use in Transforms*/
    init_potential_chemistry_restricted(potential_matrix_original_fermions_up,potential_matrix_original_fermions_spin,potential_supermatrix_fermions,potential_eigs_fermions,potential_eigvecs_fermions,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_one_body_matrix_up,ist,cns);


    /*If A Mean Field Is Subtracted from the Transform, Then Need Mean Fields to Subtract*/
    if ( ist.flag_meanfield == 1 ) {
       get_mean_field_densities(trial_density_up,trial_density_down,mean_field_density_first_transform,mean_field_density_second_transform,potential_eigvecs_fermions,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,ist); 
    }

    /*Obtain Matrices, Wavefunctions, and Initialize Walkers*/
    init_kinetic_chemistry_real_restricted(kinetic_matrix_original_fermions_up,potential_supermatrix_fermions,kinetic_full_fermions_real_up,kinetic_backwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_up,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,potential_eigvecs_fermions,potential_eigs_fermions,mean_field_density_first_transform,mean_field_density_second_transform,potential_one_body_matrix_up,ist,cns);

    /*First Regress the Walkers By Half a Kinetic Propagator*/
    propagate_half_backwards_kinetic_chemistry_real_restricted(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

    /*Production Phase*************************************************************************/
    for (steps = 0; steps < ist.n_steps; steps++) {

      if ( steps%20== 0 ) {
        fprintf(pout, "steps: %d\n", steps); fflush(pout);
      }

      /*If Back Propagation Possible*/
      if ( ist.flag_back_propagation == 0 ) { /**************************************************************************************/ 

        /*First Propagate a Full Step Forwards*/
        propagate_forwards_kinetic_chemistry_real_restricted(kinetic_full_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist);

        /*Then Propagate Forwards a Full U Step******************************************************************************************/
        /*Use Discrete or Continuous Transform*/         
        if ( ist.flag_real_eigenvectors == 0 ) { /*Determine If the Eigenctors Are Real or Not*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {
              propagate_forwards_potential_continuous_mixed_chemistry_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum);
            } /*Mean Field Flag*/
            else { /*If Mean Field*/
              propagate_forwards_potential_continuous_mixed_chemistry_meanfield_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,ist,cns,&idum); 
           }
          } /*If Phaseless*/
          else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
                propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,ist,cns,&idum);
            }
            else {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,ist,cns,&idum);
            }
          }   /*If Mean Not Field*/
          else { /*If Mean Field*/
           if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_meanfield_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum); 
           }
           else {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum); 
           }
          } /*If Mean Field*/
         } /*Phaseless Else*/ 
        }
        else {   /*If Eigenvectors Are Real*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {
                propagate_forwards_potential_continuous_mixed_chemistry_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum); 
             } /*If Not Mean Field*/
             else { /*If Mean Field*/
                 propagate_forwards_potential_continuous_mixed_chemistry_meanfield_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,ist,cns,&idum); 
             } 
           }
           else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,ist,cns,&idum);
             }
            else {
              propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,ist,cns,&idum);
            }
           }
           else { /*If Mean Field*/
             if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_meanfield_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum); 
             }
             else {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum); 
             }
           }
          } /*Phaseless Else*/
        } /*End Eigenvectors Real********************************************************************************************************/

         /*Normalize the Walker Weights by the Trial Energy*/
         normalize(weights,ist,cns); 

        /*Now Collect Energies and Other Observables*/
        if ( steps%ist.n_steps_energy==0 ) {
          /*First Propagate One Half Step Forward*/

          propagate_half_forwards_kinetic_chemistry_real_restricted(kinetic_forwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

          /*If You Want to Calculate 2 RDMS*/ 
          if ( ist.flag_compute_2RDM == 0 ) {

              /*Compute Energ*/
              compute_energy_chemistry_mixed_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);

          }
          else {   /*2RDMS*******************************************************************************************************************************************************/ 

/*in her with mixed*/

            /*Compute Energy and Density Matrices*/
            compute_energy_chemistry_mixed_two_body_density_matrix_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down);

           /*Print RDMs******************************/
           fprintf(pf2, "%d\n", steps); 
           fprintf(pf3, "%d\n", steps); 
           for (i=0; i<ist.n_spatial_orbitals; i++) {
            for (j=0; j<ist.n_spatial_orbitals; j++) {
              //if ( Cabs(average_density_matrix_up[i*ist.n_spatial_orbitals+j]) > .0001 || Cabs(average_density_matrix_down[i*ist.n_spatial_orbitals+j]) > .0001 ) {
                fprintf(pf2, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_up[i*ist.n_spatial_orbitals+j].real, average_density_matrix_up[i*ist.n_spatial_orbitals+j].imag); 
                fprintf(pf3, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_down[i*ist.n_spatial_orbitals+j].real, average_density_matrix_down[i*ist.n_spatial_orbitals+j].imag);

              //}
             }
            } 
            fprintf(pf2, "\n\n"); fflush(pf2); 
            fprintf(pf3, "\n\n"); fflush(pf3);  

            /*Print 2-Body RDM*/
            fprintf(pf4, "%d\n", steps); 
            for (i=0; i<ist.n_spin_orbitals; i++) {
             for (j=0; j<ist.n_spin_orbitals; j++) {
              for (k=0; k<ist.n_spin_orbitals; k++) {
               for (l=0; l<ist.n_spin_orbitals; l++) {
                if (  Cabs(average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]) > .0001 ) { 
                 fprintf(pf4, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].imag); 
                }
               }
              }
             }
            }
            fprintf(pf4, "\n\n"); fflush(pf4);            

          } 
          cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));

          /*Print Energies*/
          fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, energy_shifts[0].real, te_mixed.real, te_mixed_error.real, tw_mixed.real, tw_mixed.imag); fflush(pf);  
          fprintf(pf9, "%f\t\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9);  

          /*Print Wave Functions*/
          fprintf(pf5, "%d\n", steps); 
          fprintf(pf6, "%d\n", steps); 
          for (i=0; i<ist.n_spatial_orbitals; i++) {
            for (j=0; j<ist.n_up; j++) {
             fprintf(pf5, "%f+%fi\t", average_wave_function_up[i*ist.n_up+j].real, average_wave_function_up[i*ist.n_up+j].imag); 
           }
           for (j=0; j<ist.n_down; j++) {
             fprintf(pf6, "%f+%fi\t", average_wave_function_down[i*ist.n_down+j].real, average_wave_function_down[i*ist.n_down+j].imag); 
           }
           fprintf(pf5, "\n"); 
           fprintf(pf6, "\n"); 
          }
          fprintf(pf5, "\n\n"); fflush(pf5); 
          fprintf(pf6, "\n\n"); fflush(pf6);  

          /*Propagate Back Backwards*/
           propagate_half_backwards_kinetic_chemistry_real_restricted(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);
 
        }/*Energy Calculations***************************************************************************************************/
  
        /*Orthogonalize Q and R if Necessary*/
        if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
          orthogonalize_chemistry(wf_up,wf_down,overlap_up,overlap_down,weights,R_up,R_down,ist);
        }
        /*Apply Population Control If Necessary*/
        if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
            ist.n_walkers = population_control_chemistry(weights,wf_up,wf_down,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,new_walker_weights,new_wf_up,new_wf_down,new_overlap_inverse_up,new_overlap_inverse_down,new_overlap_up,new_overlap_down,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns);
        }


     } /*If No Back Propagation*/
     else { /*If Back Propagation*/

        /*First Propagate a Full Step Forwards*/
        propagate_forwards_kinetic_chemistry_real_restricted(kinetic_full_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist);

        /*Then Propagate Forwards a Full U Step******************************************************************************************/
        /*Use Discrete or Continuous Transform*/
        if ( ist.flag_real_eigenvectors == 0 ) { /*Determine If the Eigenctors Are Real or Not*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {  
              propagate_forwards_potential_continuous_bp_chemistry_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            } /*Mean Field Flag*/
            else { /*If Mean Field*/
              propagate_forwards_potential_continuous_bp_chemistry_meanfield_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
           }
          } /*If Phaseless*/
          else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
                 propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            }
            else {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            }
          }   /*If Mean Not Field*/
          else { /*If Mean Field*/
           if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_meanfield_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
           }
           else {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_meanfield_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
           }
          } /*If Mean Field*/
         } /*Phaseless Else*/
        }
        else {   /*If Eigenvectors Are Real*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {  
                propagate_forwards_potential_continuous_bp_chemistry_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             } /*If Not Mean Field*/
             else { /*If Mean Field*/
                 propagate_forwards_potential_continuous_bp_chemistry_meanfield_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
           }
           else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
                propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
            else {
              propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            }
           }
           else { /*If Mean Field*/
             if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_meanfield_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
             else {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_meanfield_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
           }
          } /*Phaseless Else*/
        } /*End Eigenvectors Real********************************************************************************************************/

         /*Store WF, If Needed for Back Prop*/
         store_prev_wf(wf_up,wf_down,prev_wf_up,prev_wf_down,kinetic_forwards_half_fermions_real_up,ist);

         /*Normalize the Walker Weights by the Trial Energy*/
         normalize(weights,ist,cns);

        /*Now Collect Energies and Other Observables*/
        if ( steps%ist.n_steps_energy==0 ) {
          /*First Propagate One Half Step Forward*/

          propagate_half_forwards_kinetic_chemistry_real_restricted(kinetic_forwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);


          /*If You Want to Calculate 2 RDMS*/
          if ( ist.flag_compute_2RDM == 0 ) {

            /*If You Want to Use Back Propagated Wave Function*/
            if (  ist.flag_back_propagation == 0 || steps < ist.n_steps_back_propagation ) {

  fprintf(pout, "in not back propagation\n"); fflush(pout); 

              /*Compute Energ*/
              compute_energy_chemistry_mixed_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);

             //compute_energy_chemistry_restricted_exact(wf_up,wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_exact,&ke_exact_error,&pe_exact,&pe_exact_error,&te_exact,&te_exact_error,&tw_exact,average_wave_function_up,average_wave_function_down);

             cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));

             /*Print Energies*/
             fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, energy_shifts[0].real, te_mixed.real, te_mixed_error.real, tw_mixed.real, tw_mixed.imag); fflush(pf);
             fprintf(pf9, "%f\t\t %f\t\t %f\n", steps*cns.dtau, te_mixed.real, te_exact.real); fflush(pf9);

            }
            else {
 fprintf(pout, "in back propagation\n"); fflush(pout); 

              /*Back Propagate First*/
              back_propagation_chemistry_restricted(trial_wf_up_energy,trial_wf_down_energy,bp_wf_up,bp_wf_down,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,kinetic_backwards_half_fermions_real_up,kinetic_full_fermions_real_up,kinetic_forwards_half_fermions_real_up,ist,cns);

              /*Compute Energy Second*/
              compute_energy_chemistry_bp_restricted(wf_up,wf_down,bp_wf_up,bp_wf_down,prev_wf_up,prev_wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_bp,&ke_bp_error,&pe_bp,&pe_bp_error,&te_bp,&te_bp_error,&tw_bp,average_wave_function_up,average_wave_function_down);
 
              fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_bp.real, pe_bp_error.real, ke_bp.real, ke_bp_error.real, energy_shifts[0].real, te_bp.real, te_bp_error.real, tw_bp.real, tw_bp.imag); fflush(pf);
              fprintf(pf9, "%f\t\t %f\t\t %f\n", steps*cns.dtau, te_bp.real, te_exact.real); fflush(pf9);
              cns.exp_trial_energy = exp(cns.dtau * (ke_bp.real + pe_bp.real)); 


              //compute_energy_chemistry_restricted_exact(wf_up,wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_exact,&ke_exact_error,&pe_exact,&pe_exact_error,&te_exact,&te_exact_error,&tw_exact,average_wave_function_up,average_wave_function_down);

              //compute_energy_chemistry_mixed_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);

              /*Print Energies*/
             //fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, energy_shifts[0].real, te_mixed.real, te_mixed_error.real, tw_mixed.real, tw_mixed.imag); fflush(pf);
             //fprintf(pf9, "%f\t\t %f\t\t %f\n", steps*cns.dtau, te_mixed.real, te_exact.real); fflush(pf9);
             // cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));
 

           }
          }
          else {   /*2RDMS*******************************************************************************************************************************************************/

           /*If You Want to Back Propagate Or Not*/
           if ( ist.flag_back_propagation == 0 || steps < ist.n_steps_back_propagation ) {

 fprintf(pout, "in 2RDM back prop\n"); fflush(pout); 

              /*Compute Energy and Density Matrices*/
              compute_energy_chemistry_mixed_two_body_density_matrix_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down);

             /*Print RDMs******************************/
             fprintf(pf2, "%d\n", steps);
             fprintf(pf3, "%d\n", steps); 
             for (i=0; i<ist.n_spatial_orbitals; i++) {
               for (j=0; j<ist.n_spatial_orbitals; j++) {
                if ( Cabs(average_density_matrix_up[i*ist.n_spatial_orbitals+j]) > .0001 || Cabs(average_density_matrix_down[i*ist.n_spatial_orbitals+j]) > .0001 ) {
                   fprintf(pf2, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_up[i*ist.n_spatial_orbitals+j].real, average_density_matrix_up[i*ist.n_spatial_orbitals+j].imag);
                   fprintf(pf3, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_down[i*ist.n_spatial_orbitals+j].real, average_density_matrix_down[i*ist.n_spatial_orbitals+j].imag);
                }
               }
              }
              fprintf(pf2, "\n\n"); fflush(pf2);
              fprintf(pf3, "\n\n"); fflush(pf3);   

              /*Print 2-Body RDM*/
              fprintf(pf4, "%d\n", steps);
              for (i=0; i<ist.n_spin_orbitals; i++) {
               for (j=0; j<ist.n_spin_orbitals; j++) {
                for (k=0; k<ist.n_spin_orbitals; k++) {
                 for (l=0; l<ist.n_spin_orbitals; l++) {
                  if ( Cabs(average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]) > .0001 ) {  
                   fprintf(pf4, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].imag);
                  }
                 }
                }
               }
              }
              fprintf(pf4, "\n\n"); fflush(pf4);

              cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));

              /*Print Energies*/
              fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, energy_shifts[0].real, te_mixed.real, te_mixed_error.real, tw_mixed.real, tw_mixed.imag); fflush(pf);
              fprintf(pf9, "%f\t\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9);

          }
          else {

   fprintf(pout, "in back prop + 2RDM\n"); fflush(pout); 

              /*Back Propagate First*/
              back_propagation_chemistry_restricted(trial_wf_up_energy,trial_wf_down_energy,bp_wf_up,bp_wf_down,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,kinetic_backwards_half_fermions_real_up,kinetic_full_fermions_real_up,kinetic_forwards_half_fermions_real_up,ist,cns);

 compute_energy_chemistry_mixed_two_body_density_matrix_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_bp,&ke_bp_error,&pe_bp,&pe_bp_error,&te_bp,&te_bp_error,&tw_bp,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down);

              /*Compute Energy and Density Matrices*/
              //compute_energy_chemistry_bp_two_body_density_matrix_restricted(wf_up,wf_down,bp_wf_up,bp_wf_down,prev_wf_up,prev_wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_bp,&ke_bp_error,&pe_bp,&pe_bp_error,&te_bp,&te_bp_error,&tw_bp,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down);

             /*Print RDMs******************************/
             fprintf(pf2, "%d\n", steps);
             fprintf(pf3, "%d\n", steps); 
             for (i=0; i<ist.n_spatial_orbitals; i++) {
              for (j=0; j<ist.n_spatial_orbitals; j++) {
               if ( Cabs(average_density_matrix_up[i*ist.n_spatial_orbitals+j]) > .0001 || Cabs(average_density_matrix_down[i*ist.n_spatial_orbitals+j]) >.0001 ) {

                   fprintf(pf2, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_up[i*ist.n_spatial_orbitals+j].real, average_density_matrix_up[i*ist.n_spatial_orbitals+j].imag);
                   fprintf(pf3, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_down[i*ist.n_spatial_orbitals+j].real, average_density_matrix_down[i*ist.n_spatial_orbitals+j].imag);
                }
               }
              }
              fprintf(pf2, "\n\n"); fflush(pf2);
              fprintf(pf3, "\n\n"); fflush(pf3); 


              /*Print 2-Body RDM*/
              fprintf(pf4, "%d\n", steps);
              for (i=0; i<ist.n_spin_orbitals; i++) {
               for (j=0; j<ist.n_spin_orbitals; j++) {
                for (k=0; k<ist.n_spin_orbitals; k++) {
                 for (l=0; l<ist.n_spin_orbitals; l++) {
                   if ( Cabs(average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]) > .0001 ) {
                      fprintf(pf4, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].imag);
                   }
                 }
                }
               }
              }
              fprintf(pf4, "\n\n"); fflush(pf4);

             cns.exp_trial_energy = exp(cns.dtau * (ke_bp.real + pe_bp.real));

            /*Print Energies*/
            fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_bp.real, pe_bp_error.real, ke_bp.real, ke_bp_error.real, energy_shifts[0].real, te_bp.real, te_bp_error.real, tw_bp.real, tw_bp.imag); fflush(pf);
            fprintf(pf9, "%f\t\t %f\n", steps*cns.dtau, te_bp.real); fflush(pf9);

           } /*Else BP and two RDMS*/
          }

          /*Print Wave Functions*/
          fprintf(pf5, "%d\n", steps);
          fprintf(pf6, "%d\n", steps);
          for (i=0; i<ist.n_spatial_orbitals; i++) {
            for (j=0; j<ist.n_up; j++) {
             fprintf(pf5, "%f+%fi\t", average_wave_function_up[i*ist.n_up+j].real, average_wave_function_up[i*ist.n_up+j].imag);
           }
           for (j=0; j<ist.n_down; j++) {
             fprintf(pf6, "%f+%fi\t", average_wave_function_down[i*ist.n_down+j].real, average_wave_function_down[i*ist.n_down+j].imag);
           }
           fprintf(pf5, "\n");
           fprintf(pf6, "\n");
          }
          fprintf(pf5, "\n\n"); fflush(pf5);
          fprintf(pf6, "\n\n"); fflush(pf6);

          /*Propagate Back Backwards*/
           propagate_half_backwards_kinetic_chemistry_real_restricted(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

        }/*Energy Calculations***************************************************************************************************/

        /*Copy Wave Functions Back If Using Back Propagation*/
        copy_back_propagation_matrices_forward(potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns);
        copy_prev_wf_forward(prev_wf_up,prev_wf_down,ist,cns);

        /*Orthogonalize Q and R if Necessary*/
        if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
          orthogonalize_chemistry(wf_up,wf_down,overlap_up,overlap_down,weights,R_up,R_down,ist);
        }

        /*Apply Population Control If Necessary*/
	 if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
           //ist.n_walkers = population_control_chemistry(weights,wf_up,wf_down,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,new_walker_weights,new_wf_up,new_wf_down,new_overlap_inverse_up,new_overlap_inverse_down,new_overlap_up,new_overlap_down,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns);
         }

     } /*If Back Propagation */

    } /*Steps*/

    /******************************************************************************************/

    free(weights);
    free(new_walker_weights);

    free(total_potential_energy); free(total_kinetic_energy); free(total_energy_2); free(total_weight);  

    free(energy_shifts); 
    free(kinetic_matrix_original_fermions_up);
    free(potential_matrix_original_fermions_up); 
    free(potential_supermatrix_fermions);
    free(potential_eigs_fermions); free(potential_eigvecs_fermions);
    free(potential_eigs_fermions_first_transform); free(potential_eigs_fermions_second_transform);

    free(mean_field_density_first_transform);  free(mean_field_density_second_transform);
    free(mean_field_first_transform_constant); free(mean_field_second_transform_constant); 

    free(overlap_up); free(overlap_down);
    free(overlap_inverse_up); free(overlap_inverse_down);
    free(wf_up); free(wf_down);
    free(trial_wf_up_energy); free(trial_wf_down_energy);
    free(trial_density_up); free(trial_density_down); 

    free(kinetic_full_fermions_real_up); free(kinetic_backwards_half_fermions_real_up); free(kinetic_forwards_half_fermions_real_up);
    free(kinetic_eigs_fermions_real_up); free(kinetic_eigvecs_fermions_real_up);

    free(new_overlap_up); free(new_overlap_down);
    free(new_overlap_inverse_up); free(new_overlap_inverse_down);
    free(new_wf_up); free(new_wf_down);

    free(average_density_matrix_up); free(average_density_matrix_down); 
    free(average_two_body_density_matrix); 

    free(average_wave_function_up); free(average_wave_function_down); 

    free(R_up); free(R_down);

    if ( ist.flag_back_propagation == 1 ) {
       free(bp_wf_up); free(bp_wf_down); 
       free(prev_wf_up); free(prev_wf_down); 

       free(potential_matrices_back_propagation_up);
       free(potential_matrices_back_propagation_down);  
    }

    fclose(pf); fclose(pf9); 
   }
   else if ( ist.flag_simulation_type == 6 ) { /*Unrestricted HF*/

    weights=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    /*For Population Control*/
    new_walker_weights = (MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    weight_rescaling=(double*)calloc(10,sizeof(double));
    walkers_big = (int*)calloc(cns.max_number_walkers,sizeof(int));
    walkers_small = (int*)calloc(cns.max_number_walkers,sizeof(int));
    walkers_copied = (int*)calloc(cns.max_number_walkers,sizeof(int));

    total_energy_2=(MKL_Complex16*)calloc(3,sizeof(MKL_Complex16));
    total_potential_energy=(MKL_Complex16*)calloc(3,sizeof(MKL_Complex16));
    total_kinetic_energy=(MKL_Complex16*)calloc(3,sizeof(MKL_Complex16));
    total_weight=(MKL_Complex16*)calloc(1,sizeof(MKL_Complex16));

    overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    overlap_down=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_up_sq,sizeof(MKL_Complex16));
    overlap_inverse_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_down_sq,sizeof(MKL_Complex16));

    wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));

    average_density_matrix_up=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16)); 
    average_density_matrix_down=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16)); 
    
    average_two_body_density_matrix=(MKL_Complex16*)calloc(ist.n_spin_orbitals_fourth,sizeof(MKL_Complex16)); 

    average_wave_function_up=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16)); 
    average_wave_function_down=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16)); 

    trial_wf_up_energy=(MKL_Complex16*)calloc(ist.n_determinants_energy_n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    trial_wf_down_energy=(MKL_Complex16*)calloc(ist.n_determinants_energy_n_spatial_orbitals_n_down,sizeof(MKL_Complex16));
    trial_density_up=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
    trial_density_down=(MKL_Complex16*)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

    kinetic_full_fermions_real_up=(double*)calloc(ist.n_spatial_orbitals_sq,sizeof(double));
    kinetic_backwards_half_fermions_real_up=(double*)calloc(ist.n_spatial_orbitals_sq,sizeof(double));
    kinetic_forwards_half_fermions_real_up=(double*)calloc(ist.n_spatial_orbitals_sq,sizeof(double));

    kinetic_eigs_fermions_real_up=(double*)calloc(ist.n_spatial_orbitals,sizeof(double));
    kinetic_eigvecs_fermions_real_up=(double*)calloc(ist.n_spatial_orbitals_sq,sizeof(double));

    kinetic_full_fermions_real_down=(double*)calloc(ist.n_spatial_orbitals_sq,sizeof(double));
    kinetic_backwards_half_fermions_real_down=(double*)calloc(ist.n_spatial_orbitals_sq,sizeof(double));
    kinetic_forwards_half_fermions_real_down=(double*)calloc(ist.n_spatial_orbitals_sq,sizeof(double));

    kinetic_eigs_fermions_real_down=(double*)calloc(ist.n_spatial_orbitals,sizeof(double));
    kinetic_eigvecs_fermions_real_down=(double*)calloc(ist.n_spatial_orbitals_sq,sizeof(double));

    energy_shifts=(MKL_Complex16 *)calloc(3,sizeof(MKL_Complex16));

    kinetic_matrix_original_fermions_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));
    kinetic_matrix_original_fermions_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_sq,sizeof(MKL_Complex16));

    kinetic_matrix_sparse_up=(double *)calloc(ist.n_kinetic_sparse_up,sizeof(double)); 
    kinetic_matrix_sparse_down=(double *)calloc(ist.n_kinetic_sparse_down,sizeof(double));

    kinetic_ij_sparse_up=(int*)calloc(ist.n_kinetic_sparse_up,sizeof(int)); 
    kinetic_ij_sparse_down=(int*)calloc(ist.n_kinetic_sparse_down,sizeof(int));  

    potential_one_body_matrix_up=(double *)calloc(ist.n_spatial_orbitals_sq,sizeof(double)); 
    potential_one_body_matrix_down=(double *)calloc(ist.n_spatial_orbitals_sq,sizeof(double));  

    potential_matrix_sparse_up=(double*)calloc(ist.n_potential_sparse_up,sizeof(double)); 
    potential_matrix_sparse_down=(double*)calloc(ist.n_potential_sparse_down,sizeof(double)); 
    potential_matrix_sparse_updown=(double*)calloc(2*ist.n_potential_sparse_updown,sizeof(double));

    potential_ij_sparse_up=(int*)calloc(ist.n_potential_sparse_up*4,sizeof(int)); 
    potential_ij_sparse_down=(int*)calloc(ist.n_potential_sparse_down*4,sizeof(int)); 
    potential_ij_sparse_updown=(int*)calloc(ist.n_potential_sparse_updown*2,sizeof(int));  

    potential_matrix_original_fermions_up=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_fourth,sizeof(MKL_Complex16));
    potential_matrix_original_fermions_down=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_fourth,sizeof(MKL_Complex16));
    potential_matrix_original_fermions_updown=(MKL_Complex16 *)calloc(ist.n_spatial_orbitals_fourth,sizeof(MKL_Complex16));
    potential_matrix_original_fermions_spin=(MKL_Complex16 *)calloc(ist.n_spin_orbitals_fourth,sizeof(MKL_Complex16)); 

    potential_supermatrix_fermions=(MKL_Complex16 *)calloc(ist.n_spin_orbitals_fourth,sizeof(MKL_Complex16));
    potential_eigs_fermions=(MKL_Complex16 *)calloc(ist.n_spin_orbitals_sq,sizeof(MKL_Complex16));
    potential_eigvecs_fermions=(MKL_Complex16*)calloc(ist.n_spin_orbitals_fourth,sizeof(MKL_Complex16));

    potential_eigs_fermions_first_transform=(MKL_Complex16 *)calloc(ist.n_spin_orbitals_sq,sizeof(MKL_Complex16));
    potential_eigs_fermions_second_transform=(MKL_Complex16 *)calloc(ist.n_spin_orbitals_sq,sizeof(MKL_Complex16));

    mean_field_density_first_transform=(MKL_Complex16 *)calloc(ist.n_spin_orbitals_sq,sizeof(MKL_Complex16));
    mean_field_density_second_transform=(MKL_Complex16 *)calloc(ist.n_spin_orbitals_sq,sizeof(MKL_Complex16));

    mean_field_first_transform_constant=(MKL_Complex16 *)calloc(ist.n_spin_orbitals_sq,sizeof(MKL_Complex16));
    mean_field_second_transform_constant=(MKL_Complex16 *)calloc(ist.n_spin_orbitals_sq,sizeof(MKL_Complex16));

    list_orbitals_trial_energy = (int *)malloc(ist.n_spatial_orbitals*sizeof(int));

    /*For Population Control*/
    new_wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
    new_wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));
    new_overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_up_sq,sizeof(MKL_Complex16));
    new_overlap_inverse_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_down_sq,sizeof(MKL_Complex16));
    new_overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
    new_overlap_down=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

    /*For Orthogonalization*/
    R_up=(MKL_Complex16 *)calloc(ist.n_up*ist.n_up,sizeof(MKL_Complex16));
    R_down=(MKL_Complex16 *)calloc(ist.n_down*ist.n_down,sizeof(MKL_Complex16));

    /*If There Is Back Propagation*/
    if ( ist.flag_back_propagation == 1 ) {
       bp_wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16)); 
       bp_wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16)); 
    } 

    /******************************************************************************************/

    /*Zero the Accumulated Energies and Weights*/
    czero_vec(total_potential_energy,3); 
    czero_vec(total_kinetic_energy,3); 
    czero_vec(total_energy_2,3); 
    czero_vec(total_weight,1); 

    /*Obtain Trial Wavefunction*/ 
    ist.n_max_orbital_trial_energy = init_wf_chemistry_unrestricted(trial_wf_up_energy,trial_wf_down_energy,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,list_orbitals_trial_energy,ist,cns);

    /*Input the Vijkl Elements and Resort Them to Find Eigenvalues*/
    get_potential_matrix_elements_molpro_unrestricted(kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,potential_matrix_original_fermions_spin,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,list_orbitals_trial_energy,energy_shifts,&number_kinetic_sparse_up,&number_kinetic_sparse_down,&number_potential_sparse_up,&number_potential_sparse_down,&number_potential_sparse_updown,ist,cns,"potential_matrix_elements.par");
    ist.n_kinetic_sparse_up = number_kinetic_sparse_up;
    ist.n_kinetic_sparse_down = number_kinetic_sparse_down;
    ist.n_potential_sparse_up = number_potential_sparse_up;
    ist.n_potential_sparse_down = number_potential_sparse_down;
    ist.n_potential_sparse_updown = number_potential_sparse_updown;

    /*Obtain Trial Wave Function, Trial Energy, and Trial Density Matrices from WF*/
    if ( ist.flag_restart == 0 ) {
      init_walkers_chemistry_no_restart_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,weights,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weight_rescaling,ist,cns);
    }
    else {
      init_walkers_chemistry_restart_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,weights,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weight_rescaling,ist,cns);
    }

    /*Get Trial Energies*/
    cns.trial_energy = get_trial_energy_density_unrestricted(trial_wf_up_energy,trial_wf_down_energy,trial_density_up,trial_density_down,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,energy_shifts,ist);
    cns.exp_trial_energy = exp(cns.dtau * cns.trial_energy);

/*Form the Supermatrix and Diagonalize It to Use in Transforms*/
    init_potential_chemistry_unrestricted(potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,potential_supermatrix_fermions,potential_eigs_fermions,potential_eigvecs_fermions,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_one_body_matrix_up,potential_one_body_matrix_down,ist,cns);

    /*If A Mean Field Is Subtracted from the Transform, Then Need Mean Fields to Subtract*/
    if ( ist.flag_meanfield == 1 ) {
       get_mean_field_densities(trial_density_up,trial_density_down,mean_field_density_first_transform,mean_field_density_second_transform,potential_eigvecs_fermions,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,ist);
    }

    /*Obtain Matrices, Wavefunctions, and Initialize Walkers*/
    init_kinetic_chemistry_real_unrestricted(kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_supermatrix_fermions,kinetic_full_fermions_real_up,kinetic_backwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_up,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,kinetic_full_fermions_real_down,kinetic_backwards_half_fermions_real_down,kinetic_forwards_half_fermions_real_down,kinetic_eigs_fermions_real_down,kinetic_eigvecs_fermions_real_down,potential_eigvecs_fermions,potential_eigs_fermions,mean_field_density_first_transform,mean_field_density_second_transform,potential_one_body_matrix_up,potential_one_body_matrix_down,ist,cns);
 
    /*First Regress the Walkers By Half a Kinetic Propagator*/
    propagate_half_backwards_kinetic_chemistry_real_unrestricted(kinetic_backwards_half_fermions_real_up,kinetic_backwards_half_fermions_real_down,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

    /*Production Phase*************************************************************************/
    for (steps = 0; steps < ist.n_steps; steps++) {

      if ( steps%20== 0 ) {
        fprintf(pout, "steps: %d\n", steps); fflush(pout);
      }

      /*First Propagate a Full Step Forwards*/
      propagate_forwards_kinetic_chemistry_real_unrestricted(kinetic_full_fermions_real_up,kinetic_full_fermions_real_down,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist);

      /*Then Propagate Forwards a Full U Step**********************************************************************************/
      /*Use Discrete or Continuous Transform*/
      if ( ist.flag_real_eigenvectors == 0 ) { /*Determine If the Eigenctors Are Real or Not*/
        if ( ist.flag_phaseless == 0 ) {
         if ( ist.flag_meanfield == 0 ) {
          propagate_forwards_potential_continuous_mixed_chemistry_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum);
         } /*Mean Field Flag*/
         else { /*If Mean Field*/
           propagate_forwards_potential_continuous_mixed_chemistry_meanfield_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,ist,cns,&idum);
         }
         } /*If Phaseless*/
        else {
          if ( ist.flag_meanfield == 0 ) {
           if ( ist.flag_local_energy == 0 ) {
            propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,ist,cns,&idum);
          }
          else {
            propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_unrestricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,ist,cns,&idum);
          }
        } /*If Mean Not Field*/
        else { /*If Mean Field*/
         if ( ist.flag_local_energy == 0 ) {
           propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_meanfield_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
         }
         else {
           propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_unrestricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
         }
        } /*If Mean Field*/
       } /*Phaseless Else*/
      }
      else {   /*If Eigenvectors Are Real*/
         if ( ist.flag_phaseless == 0 ) {
          if ( ist.flag_meanfield == 0 ) {
              propagate_forwards_potential_continuous_mixed_chemistry_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum);
           } /*If Not Mean Field*/
           else { /*If Mean Field*/
               propagate_forwards_potential_continuous_mixed_chemistry_meanfield_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,ist,cns,&idum);
           }
         }
         else {
          if ( ist.flag_meanfield == 0 ) {
            if ( ist.flag_local_energy == 0 ) {
             propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,ist,cns,&idum);
           }
          else {
            propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_real_unrestricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,ist,cns,&idum);
          }
         }
         else { /*If Mean Field*/
           if ( ist.flag_local_energy == 0 ) {
             propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_meanfield_real_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
           }
           else {
             propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_unrestricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
}
         }
        } /*Phaseless Else*/
      } /*End Eigenvectors Real********************************************************************************************************/

     /*Normalize by Trial Energy*/
     normalize(weights,ist,cns);

      /*Now Collect Energies and Other Observables*/
      if ( steps%ist.n_steps_energy==0 ) {
        /*First Propagate One Half Step Forward*/
        propagate_half_forwards_kinetic_chemistry_real_unrestricted(kinetic_forwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_down,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

        /*Compute Energ*/
        if ( ist.flag_compute_2RDM == 0 ) {
          compute_energy_chemistry_mixed_unrestricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&pe_mixed,&pe_mixed_error,&ke_mixed,&ke_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);
        }
        else {  
           compute_energy_chemistry_mixed_two_body_density_matrix_unrestricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&pe_mixed,&pe_mixed_error,&ke_mixed,&ke_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down); 
       
          /*Print RDMs******************************/
          fprintf(pf2, "%d\n", steps);
          fprintf(pf3, "%d\n", steps); 
          for (i=0; i<ist.n_spatial_orbitals; i++) {
            for (j=0; j<ist.n_spatial_orbitals; j++) {
             if ( Cabs(average_density_matrix_up[i*ist.n_spatial_orbitals+j]) > .0001 || Cabs(average_density_matrix_down[i*ist.n_spatial_orbitals+j]) > .0001 ) {
               fprintf(pf2, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_up[i*ist.n_spatial_orbitals+j].real, average_density_matrix_up[i*ist.n_spatial_orbitals+j].imag);
               fprintf(pf3, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_down[i*ist.n_spatial_orbitals+j].real, average_density_matrix_down[i*ist.n_spatial_orbitals+j].imag);
             }
            }
           }
           fprintf(pf2, "\n\n"); fflush(pf2);
           fprintf(pf3, "\n\n"); fflush(pf3); 

           /*Print 2-Body RDM*/
           fprintf(pf4, "%d\n", steps);
           for (i=0; i<ist.n_spin_orbitals; i++) {
            for (j=0; j<ist.n_spin_orbitals; j++) {
             for (k=0; k<ist.n_spin_orbitals; k++) {
              for (l=0; l<ist.n_spin_orbitals; l++) {
                if ( Cabs(average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]) > .0001 ) {
                  fprintf(pf4, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].imag);
                }
              }
             }
            }
           }          
           fprintf(pf4, "\n\n"); fflush(pf4);
          
        }
        cns.exp_trial_energy = exp(cns.dtau * (pe_mixed.real + ke_mixed.real));

        /*Print Out Energies*/
        fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, energy_shifts[0].real, te_mixed.real, te_mixed_error.real, tw_mixed.real, tw_mixed.imag); fflush(pf);
        fprintf(pf9, "%f\t\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9);

        /*Print Out Wave Functions*/
        fprintf(pf5, "step %d\n", steps); 
        fprintf(pf6, "step %d\n", steps); 
        for (i=0; i<ist.n_spatial_orbitals; i++) {
         for (j=0; j<ist.n_up; j++) {
            fprintf(pf5, "%f+%fi\t", average_wave_function_up[i*ist.n_up+j].real, average_wave_function_up[i*ist.n_up+j].imag); 
         }
         for (j=0; j<ist.n_down; j++) {
            fprintf(pf6, "%f+%fi\t", average_wave_function_down[i*ist.n_down+j].real, average_wave_function_down[i*ist.n_down+j].imag); 
          }
         fprintf(pf5, "\n"); 
         fprintf(pf6, "\n"); 
        }
        fprintf(pf5, "\n\n"); fflush(pf5); 
        fprintf(pf6, "\n\n"); fflush(pf6); 

        /*Propagate Back Backwards*/
         propagate_half_backwards_kinetic_chemistry_real_unrestricted(kinetic_backwards_half_fermions_real_up,kinetic_backwards_half_fermions_real_down,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);
      }

      /*Orthogonalize Q and R if Necessary*/
      if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
        orthogonalize_chemistry(wf_up,wf_down,overlap_up,overlap_down,weights,R_up,R_down,ist);
      }

      /*Apply Population Control If Necessary*/
      if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
          //ist.n_walkers = population_control_chemistry(weights,wf_up,wf_down,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,new_walker_weights,new_wf_up,new_wf_down,new_overlap_inverse_up,new_overlap_inverse_down,new_overlap_up,new_overlap_down,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns);
      }

    } /*Steps*/

    /******************************************************************************************/

    free(weights);
    free(new_walker_weights);

    free(energy_shifts);
    free(kinetic_matrix_original_fermions_up); free(kinetic_matrix_original_fermions_down);

    free(potential_matrix_original_fermions_up); free(potential_matrix_original_fermions_down); free(potential_matrix_original_fermions_updown);
    free(potential_supermatrix_fermions);
    free(potential_eigs_fermions); free(potential_eigvecs_fermions);
    free(potential_eigs_fermions_first_transform); free(potential_eigs_fermions_second_transform);

    free(potential_one_body_matrix_up); free(potential_one_body_matrix_down); 

    free(mean_field_density_first_transform);  free(mean_field_density_second_transform);
    free(mean_field_first_transform_constant); free(mean_field_second_transform_constant);

    free(overlap_up); free(overlap_down);
    free(overlap_inverse_up); free(overlap_inverse_down);
    free(wf_up); free(wf_down);

    free(trial_wf_up_energy); free(trial_wf_down_energy);
    free(trial_density_up); free(trial_density_down);

    free(average_density_matrix_up); free(average_density_matrix_down); 
    free(average_two_body_density_matrix); 

    free(average_wave_function_up); free(average_wave_function_down); 

    free(kinetic_full_fermions_real_up); free(kinetic_backwards_half_fermions_real_up); free(kinetic_forwards_half_fermions_real_up);
    free(kinetic_eigs_fermions_real_up); free(kinetic_eigvecs_fermions_real_up);
 
    free(kinetic_full_fermions_real_down); free(kinetic_backwards_half_fermions_real_down); 
    free(kinetic_forwards_half_fermions_real_down);
    free(kinetic_eigs_fermions_real_down); free(kinetic_eigvecs_fermions_real_down);

    free(new_overlap_up); free(new_overlap_down);
    free(new_overlap_inverse_up); free(new_overlap_inverse_down);
    free(new_wf_up); free(new_wf_down);
    free(R_up); free(R_down);

    /*Free Back Propagation*/
    if ( ist.flag_back_propagation == 1 ) {
      free(bp_wf_up); free(bp_wf_down); 
    }

   } /*Else If Unrestricted HF Simulation*/
   else if ( ist.flag_simulation_type == 7 ) { /*Extended Hubbard Model************************************/ 

    /*If Just One Up Electron */
    if ( ist.n_up >=1 && ist.n_down >=1 ) {

      neighbors=(int*)calloc(4*ist.n_sites,sizeof(int));
      number_neighbors=(int*)calloc(ist.n_sites,sizeof(int));

      weights=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

      /*For Population Control*/
      new_walker_weights = (MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
      weight_rescaling=(double*)calloc(10,sizeof(double));
      walkers_big = (int*)calloc(cns.max_number_walkers,sizeof(int));
      walkers_small = (int*)calloc(cns.max_number_walkers,sizeof(int));
      walkers_copied = (int*)calloc(cns.max_number_walkers,sizeof(int));

      potential_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
      kinetic_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
      total_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
      total_weights=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));

      overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
      overlap_down=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

      overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));
      overlap_inverse_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));

      wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_up,sizeof(MKL_Complex16));
      wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_down,sizeof(MKL_Complex16));

      average_wf_up=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16)); 
      average_wf_down=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16)); 

      trial_wf_up_energy=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));
      trial_wf_down_energy=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));

      trial_density_up=(MKL_Complex16*)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 
      trial_density_down=(MKL_Complex16*)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 
 
      kinetic_full_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));
      kinetic_backwards_half_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));
      kinetic_forwards_half_fermions_real_up=(double *)calloc(ist.n_sites_sq,sizeof(double));

      kinetic_eigs_fermions_real_up=(double*)calloc(ist.n_sites,sizeof(double));
      kinetic_eigvecs_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));

      /*For Population Control*/
      new_wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_up,sizeof(MKL_Complex16));
      new_wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_down,sizeof(MKL_Complex16));
      new_overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));
      new_overlap_inverse_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));
      new_overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
      new_overlap_down=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

      /*For Orthogonalization*/
      R_up=(MKL_Complex16 *)calloc(ist.n_up*ist.n_up,sizeof(MKL_Complex16));
      R_down=(MKL_Complex16 *)calloc(ist.n_down*ist.n_down,sizeof(MKL_Complex16));

      /******************************************************************************************/

      /*Obtain Matrices, Wavefunctions, and Initialize Walkers*/
      init_neighbors(neighbors,number_neighbors,ist);

      init_kinetic_fermions(kinetic_full_fermions_real_up,kinetic_backwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_up,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,neighbors,number_neighbors,ist,cns);
      init_wf_fermions(trial_wf_up_energy,trial_wf_down_energy,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,ist,cns);
      init_walkers_fermions(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,weights,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weight_rescaling,ist,cns);

      cns.trial_energy = get_trial_energy_density_fermions(trial_wf_up_energy,trial_wf_down_energy,trial_density_up,trial_density_down,neighbors,number_neighbors,ist,cns);
      cns.exp_trial_energy = exp(cns.dtau * cns.trial_energy);

      /*******************************************************************************************/

      /*First Regress the Walkers By Half a Kinetic Propagator*/
      propagate_half_backwards_kinetic_fermions(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

      /*Production Phase*************************************************************************/
      for (steps = 0; steps < ist.n_steps; steps++) {
        if ( steps%20== 0 ) {
          fprintf(pout, "steps: %d\n", steps); fflush(pout);
        }

        /*First Propagate a Full Step Forwards*/
        propagate_forwards_kinetic_fermions(kinetic_full_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,weights,ist);

        /*Then Propagate Forwards a Full U Step*/
        if ( cns.U_fermions != 0 || cns.V_fermions != 0 ) {
          /*Use Discrete or Continuous Transform*/
              //propagate_forwards_potential_continuous_extended(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,neighbors,number_neighbors,weights,ist,cns,&idum); 
        }

        /*Now Collect Energies and Other Observables*/
        if ( steps%ist.n_steps_energy==0 ) {
          /*First Propagate One Half Step Forward*/
          propagate_half_forwards_kinetic_fermions(kinetic_forwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

          /*Compute Energ*/
          compute_energy_extended(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_inverse_up,overlap_inverse_down,weights,neighbors,number_neighbors,potential_energy,kinetic_energy,total_energy,total_weights,weight_rescaling,ist,cns,steps,&pe_mixed,&ke_mixed,&te_mixed,&tw_mixed,average_wf_up,average_wf_down);
          cns.exp_trial_energy = exp(cns.dtau * (pe_mixed.real + ke_mixed.real));
 
          /*Print Energies*/
          fprintf(pf, "%d\t\t  %d\t\t  %f\t\t   %f\t\t  %f\t\t  %f\t\t %f\n", steps, ist.n_walkers, steps*cns.dtau, ke_mixed.real, pe_mixed.real, te_mixed.real, tw_mixed.real); fflush(pf);
          fprintf(pf9, "%f\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9);
  
          /*Print Wavefunctions*/
          print_cmat(average_wf_up,ist.n_sites,ist.n_up, "average_wave_function_up.dat"); 
          print_cmat(average_wf_down,ist.n_sites,ist.n_down,"average_wave_function_down.dat"); 
  
          /*Propagate Back Backwards*/
          propagate_half_backwards_kinetic_fermions(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,overlap_up,overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);
        }
  
        /*Orthogonalize Q and R if Necessary*/
        if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
          orthogonalize_fermions(wf_up,wf_down,overlap_up,overlap_down,weights,R_up,R_down,ist);
        }
  
        /*Normalize By Trial Energy*/
        normalize(weights,ist,cns);
 
        /*Apply Population Control If Necessary*/
        if ( cns.U_fermions != 0 || cns.V_fermions != 0 ) {
          if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
            //ist.n_walkers = population_control_fermions(weights,wf_up,wf_down,overlap_inverse_up,overlap_inverse_down,overlap_up,overlap_down,new_walker_weights,new_wf_up,new_wf_down,new_overlap_inverse_up,new_overlap_inverse_down,new_overlap_up,new_overlap_down,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns); 
          }
        }
 
      } /*Steps*/

      /*Now Compute Average Overall Energy*/
      compute_averaged_energy_fermions(potential_energy,kinetic_energy,total_energy,total_weights,ist,"final_average_energy.dat");

      /******************************************************************************************/

      free(weights);
      free(new_walker_weights);

      free(potential_energy); free(kinetic_energy); free(total_energy); free(total_weights);
  
      free(overlap_up); free(overlap_down);
      free(overlap_inverse_up); free(overlap_inverse_down);
      free(wf_up); free(wf_down);
      free(average_wf_up); free(average_wf_down); 
      free(trial_wf_up_energy); free(trial_wf_down_energy);
      free(trial_density_up); free(trial_density_down); 
      free(kinetic_full_fermions_real_up); free(kinetic_backwards_half_fermions_real_up); free(kinetic_forwards_half_fermions_real_up);
      free(kinetic_eigs_fermions_real_up); free(kinetic_eigvecs_fermions_real_up);

      free(new_overlap_up); free(new_overlap_down);
      free(new_overlap_inverse_up); free(new_overlap_inverse_down);
      free(new_wf_up); free(new_wf_down);
 
      free(R_up); free(R_down);
     } /*If Both Electrons*****************************************************/
     else if ( ist.n_up >=1 && ist.n_down == 0 ) { /*If Just Up Electrons*******************************/

      neighbors=(int*)calloc(4*ist.n_sites,sizeof(int));
      number_neighbors=(int*)calloc(ist.n_sites,sizeof(int));

      weights=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

      /*For Population Control*/
      new_walker_weights = (MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
      weight_rescaling=(double*)calloc(10,sizeof(double));
      walkers_big = (int*)calloc(cns.max_number_walkers,sizeof(int));
      walkers_small = (int*)calloc(cns.max_number_walkers,sizeof(int));
      walkers_copied = (int*)calloc(cns.max_number_walkers,sizeof(int));

      potential_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
      kinetic_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
      total_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
      total_weights=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));

      overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
      overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));

      wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_up,sizeof(MKL_Complex16));
      average_wf_up=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));

      trial_wf_up_energy=(MKL_Complex16*)calloc(ist.n_sites_n_up,sizeof(MKL_Complex16));
      trial_density_up=(MKL_Complex16*)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));

      kinetic_full_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));
      kinetic_backwards_half_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));
      kinetic_forwards_half_fermions_real_up=(double *)calloc(ist.n_sites_sq,sizeof(double));

      kinetic_eigs_fermions_real_up=(double*)calloc(ist.n_sites,sizeof(double));
      kinetic_eigvecs_fermions_real_up=(double*)calloc(ist.n_sites_sq,sizeof(double));

      /*For Population Control*/
      new_wf_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_up,sizeof(MKL_Complex16));
      new_overlap_inverse_up=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));
      new_overlap_up=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

      /*For Orthogonalization*/
      R_up=(MKL_Complex16 *)calloc(ist.n_up*ist.n_up,sizeof(MKL_Complex16));

      /******************************************************************************************/

      /*Obtain Matrices, Wavefunctions, and Initialize Walkers*/
      init_neighbors(neighbors,number_neighbors,ist);

      init_kinetic_fermions(kinetic_full_fermions_real_up,kinetic_backwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_up,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,neighbors,number_neighbors,ist,cns);
      init_wf_fermions_up(trial_wf_up_energy,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,ist,cns);
      init_walkers_fermions_up(wf_up,trial_wf_up_energy,weights,overlap_up,overlap_inverse_up,weight_rescaling,ist,cns);

      cns.trial_energy = get_trial_energy_density_fermions_up(trial_wf_up_energy,trial_density_up,neighbors,number_neighbors,ist,cns);
      cns.exp_trial_energy = exp(cns.dtau * cns.trial_energy);

      /*******************************************************************************************/

      /*First Regress the Walkers By Half a Kinetic Propagator*/
      propagate_half_backwards_kinetic_fermions_up(kinetic_backwards_half_fermions_real_up,wf_up,trial_wf_up_energy,overlap_up,overlap_inverse_up,weights,ist);

      /*Production Phase*************************************************************************/
      for (steps = 0; steps < ist.n_steps; steps++) {
        if ( steps%20== 0 ) {
          fprintf(pout, "steps: %d\n", steps); fflush(pout);
        }

        /*First Propagate a Full Step Forwards*/
        propagate_forwards_kinetic_fermions_up(kinetic_full_fermions_real_up,wf_up,trial_wf_up_energy,overlap_inverse_up,overlap_up,weights,ist);

         /*Then Propagate Forwards a Full U Step*/
        if ( cns.U_fermions != 0 || cns.V_fermions != 0 ) {
          /*Use Discrete or Continuous Transform*/
            propagate_forwards_potential_continuous_extended_up(wf_up,trial_wf_up_energy,overlap_inverse_up,overlap_up,neighbors,number_neighbors,weights,ist,cns,&idum); 
        }

        /*Now Collect Energies and Other Observables*/
        if ( steps%ist.n_steps_energy==0 ) {
          /*First Propagate One Half Step Forward*/
          propagate_half_forwards_kinetic_fermions_up(kinetic_forwards_half_fermions_real_up,wf_up,trial_wf_up_energy,overlap_up,overlap_inverse_up,weights,ist);

          /*Compute Energ*/
          compute_energy_extended_up(wf_up,trial_wf_up_energy,overlap_inverse_up,weights,neighbors,number_neighbors,potential_energy,kinetic_energy,total_energy,total_weights,weight_rescaling,ist,cns,steps,&pe_mixed,&ke_mixed,&te_mixed,&tw_mixed,average_wf_up);
          cns.exp_trial_energy = exp(cns.dtau * (pe_mixed.real + ke_mixed.real));

          /*Print Energies*/
          fprintf(pf, "%d\t\t  %d\t\t  %f\t\t   %f\t\t  %f\t\t  %f\t\t %f\n", steps, ist.n_walkers, steps*cns.dtau, ke_mixed.real, pe_mixed.real, te_mixed.real, tw_mixed.real); fflush(pf);
          fprintf(pf9, "%f\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9);

          /*Print Wavefunctions*/
          print_cmat(average_wf_up,ist.n_sites,ist.n_up, "average_wave_function_up.dat");

          /*Propagate Back Backwards*/
          propagate_half_backwards_kinetic_fermions_up(kinetic_backwards_half_fermions_real_up,wf_up,trial_wf_up_energy,overlap_up,overlap_inverse_up,weights,ist);
        }

        /*Orthogonalize Q and R if Necessary*/
        if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
          orthogonalize_fermions_up(wf_up,overlap_up,weights,R_up,ist);
        }

        /*Normalize By Trial Energy*/
        normalize(weights,ist,cns);

        /*Apply Population Control If Necessary*/
        if ( cns.U_fermions != 0 || cns.V_fermions != 0 ) {
          if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
            //ist.n_walkers = population_control_fermions(weights,wf_up,overlap_inverse_up,overlap_up,new_walker_weights,new_wf_up,new_overlap_inverse_up,new_overlap_up,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns); 
          }
        }

      } /*Steps*/

      /*Now Compute Average Overall Energy*/
      compute_averaged_energy_fermions(potential_energy,kinetic_energy,total_energy,total_weights,ist,"final_average_energy.dat");

      /******************************************************************************************/

      free(weights);
      free(new_walker_weights);

      free(potential_energy); free(kinetic_energy); free(total_energy); free(total_weights);

      free(overlap_up); 
      free(overlap_inverse_up); 
      free(wf_up); 
      free(average_wf_up); 
      free(trial_wf_up_energy); 
      free(trial_density_up); 
      free(kinetic_full_fermions_real_up); free(kinetic_backwards_half_fermions_real_up); free(kinetic_forwards_half_fermions_real_up);
      free(kinetic_eigs_fermions_real_up); free(kinetic_eigvecs_fermions_real_up);

      free(new_overlap_up); 
      free(new_overlap_inverse_up); 
      free(new_wf_up); 

      free(R_up); 
     }
     else if ( ist.n_up == 0 && ist.n_down >= 1 ) { /*If All Down Electrons***************************************************************/ 
      neighbors=(int*)calloc(4*ist.n_sites,sizeof(int));
      number_neighbors=(int*)calloc(ist.n_sites,sizeof(int));

      weights=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

      /*For Population Control*/
      new_walker_weights = (MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
      weight_rescaling=(double*)calloc(10,sizeof(double));
      walkers_big = (int*)calloc(cns.max_number_walkers,sizeof(int));
      walkers_small = (int*)calloc(cns.max_number_walkers,sizeof(int));
      walkers_copied = (int*)calloc(cns.max_number_walkers,sizeof(int));

      potential_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
      kinetic_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
      total_energy=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));
      total_weights=(MKL_Complex16*)calloc(ist.n_steps,sizeof(MKL_Complex16));

      overlap_down=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));
      overlap_inverse_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));

      wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_down,sizeof(MKL_Complex16));
      average_wf_down=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));

      trial_wf_down_energy=(MKL_Complex16*)calloc(ist.n_sites_n_down,sizeof(MKL_Complex16));
      trial_density_down=(MKL_Complex16*)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));

      kinetic_full_fermions_real_down=(double*)calloc(ist.n_sites_sq,sizeof(double));
      kinetic_backwards_half_fermions_real_down=(double*)calloc(ist.n_sites_sq,sizeof(double));
      kinetic_forwards_half_fermions_real_down=(double *)calloc(ist.n_sites_sq,sizeof(double));

      kinetic_eigs_fermions_real_down=(double*)calloc(ist.n_sites,sizeof(double));
      kinetic_eigvecs_fermions_real_down=(double*)calloc(ist.n_sites_sq,sizeof(double));

      /*For Population Control*/
      new_wf_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_n_down,sizeof(MKL_Complex16));
      new_overlap_inverse_down=(MKL_Complex16*)calloc(cns.max_number_walkers*ist.n_sites_sq,sizeof(MKL_Complex16));
      new_overlap_down=(MKL_Complex16*)calloc(cns.max_number_walkers,sizeof(MKL_Complex16));

      /*For Orthogonalization*/
      R_down=(MKL_Complex16 *)calloc(ist.n_down*ist.n_down,sizeof(MKL_Complex16));

      /******************************************************************************************/

      /*Obtain Matrices, Wavefunctions, and Initialize Walkers*/
      init_neighbors(neighbors,number_neighbors,ist);

      init_kinetic_fermions(kinetic_full_fermions_real_down,kinetic_backwards_half_fermions_real_down,kinetic_forwards_half_fermions_real_down,kinetic_eigs_fermions_real_down,kinetic_eigvecs_fermions_real_down,neighbors,number_neighbors,ist,cns);
      init_wf_fermions_down(trial_wf_down_energy,kinetic_eigs_fermions_real_down,kinetic_eigvecs_fermions_real_down,ist,cns);
      init_walkers_fermions_down(wf_down,trial_wf_down_energy,weights,overlap_down,overlap_inverse_down,weight_rescaling,ist,cns);

      cns.trial_energy = get_trial_energy_density_fermions_down(trial_wf_down_energy,trial_density_down,neighbors,number_neighbors,ist,cns);
      cns.exp_trial_energy = exp(cns.dtau * cns.trial_energy);

      /*******************************************************************************************/

      /*First Regress the Walkers By Half a Kinetic Propagator*/
      propagate_half_backwards_kinetic_fermions_down(kinetic_backwards_half_fermions_real_down,wf_down,trial_wf_down_energy,overlap_down,overlap_inverse_down,weights,ist);

      /*Production Phase*************************************************************************/
      for (steps = 0; steps < ist.n_steps; steps++) {
        if ( steps%20== 0 ) {
          fprintf(pout, "steps: %d\n", steps); fflush(pout);
        }

        /*First Propagate a Full Step Forwards*/
        propagate_forwards_kinetic_fermions_down(kinetic_full_fermions_real_down,wf_down,trial_wf_down_energy,overlap_inverse_down,overlap_down,weights,ist);

         /*Then Propagate Forwards a Full U Step*/
        if ( cns.U_fermions != 0 || cns.V_fermions != 0 ) {
          /*Use Discrete or Continuous Transform*/
          propagate_forwards_potential_continuous_extended_down(wf_down,trial_wf_down_energy,overlap_inverse_down,overlap_down,neighbors,number_neighbors,weights,ist,cns,&idum);
        }

        /*Now Collect Energies and Other Observables*/
        if ( steps%ist.n_steps_energy==0 ) {
          /*First Propagate One Half Step Forward*/
          propagate_half_forwards_kinetic_fermions_down(kinetic_forwards_half_fermions_real_down,wf_down,trial_wf_down_energy,overlap_down,overlap_inverse_down,weights,ist);

          /*Compute Energ*/
          compute_energy_extended_down(wf_down,trial_wf_down_energy,overlap_inverse_down,weights,neighbors,number_neighbors,potential_energy,kinetic_energy,total_energy,total_weights,weight_rescaling,ist,cns,steps,&pe_mixed,&ke_mixed,&te_mixed,&tw_mixed,average_wf_down);
          cns.exp_trial_energy = exp(cns.dtau * (pe_mixed.real + ke_mixed.real));

          /*Print Energies*/
          fprintf(pf, "%d\t\t  %d\t\t  %f\t\t   %f\t\t  %f\t\t  %f\t\t %f\n", steps, ist.n_walkers, steps*cns.dtau, ke_mixed.real, pe_mixed.real, te_mixed.real, tw_mixed.real); fflush(pf);
          fprintf(pf9, "%f\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9);

          /*Print Wavefunctions*/
          print_cmat(average_wf_down,ist.n_sites,ist.n_down, "average_wave_function_down.dat");

          /*Propagate Back Backwards*/
          propagate_half_backwards_kinetic_fermions_down(kinetic_backwards_half_fermions_real_down,wf_down,trial_wf_down_energy,overlap_down,overlap_inverse_down,weights,ist);
        }

        /*Orthogonalize Q and R if Necessary*/
        if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
          orthogonalize_fermions_down(wf_down,overlap_down,weights,R_down,ist);
        }

        /*Normalize By Trial Energy*/
        normalize(weights,ist,cns);

        /*Apply Population Control If Necessary*/
        if ( cns.U_fermions != 0 || cns.V_fermions != 0 ) {
          if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
            //ist.n_walkers = population_control_fermions(weights,wf_down,overlap_inverse_down,overlap_down,new_walker_weights,new_wf_down,new_overlap_inverse_down,new_overlap_down,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns); 
          }
        }

      } /*Steps*/

      /*Now Compute Average Overall Energy*/
      compute_averaged_energy_fermions(potential_energy,kinetic_energy,total_energy,total_weights,ist,"final_average_energy.dat");

      /******************************************************************************************/

      free(weights);
      free(new_walker_weights);

      free(potential_energy); free(kinetic_energy); free(total_energy); free(total_weights);

      free(overlap_down);
      free(overlap_inverse_down);
      free(wf_down);
      free(average_wf_down);
      free(trial_wf_down_energy);
      free(trial_density_down);
      free(kinetic_full_fermions_real_down); free(kinetic_backwards_half_fermions_real_down); free(kinetic_forwards_half_fermions_real_down);
      free(kinetic_eigs_fermions_real_down); free(kinetic_eigvecs_fermions_real_down);

      free(new_overlap_down);
      free(new_overlap_inverse_down);
      free(new_wf_down);

      free(R_down); 
   
     } /*If No Up Electrons*/ 
    } /*End of Function 7**********************************************************************************************************************************************************************************/ 
    else if ( ist.flag_simulation_type == 8 ) { /*MultiDeterminant Trial Wave Functions RHF*************************************************************************************************************************/ 

      weights=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));

      /*For Population Control*/
      new_walker_weights = (MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));
      weight_rescaling=(double*)malloc(10*sizeof(double));
      walkers_big = (int*)malloc(cns.max_number_walkers*sizeof(int));
      walkers_small = (int*)malloc(cns.max_number_walkers*sizeof(int));
      walkers_copied = (int*)malloc(cns.max_number_walkers*sizeof(int));

      total_energy_2=(MKL_Complex16*)malloc(3*sizeof(MKL_Complex16));
      total_potential_energy=(MKL_Complex16*)malloc(3*sizeof(MKL_Complex16));
      total_kinetic_energy=(MKL_Complex16*)malloc(3*sizeof(MKL_Complex16));
      total_weight=(MKL_Complex16*)malloc(1*sizeof(MKL_Complex16));

      det_overlap_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*sizeof(MKL_Complex16)); 
      det_overlap_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*sizeof(MKL_Complex16)); 

      overlap_total=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));

      overlap_inverse_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*ist.n_up_sq*sizeof(MKL_Complex16));
      overlap_inverse_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*ist.n_down_sq*sizeof(MKL_Complex16));

      average_density_matrix_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
      average_density_matrix_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
      average_two_body_density_matrix=(MKL_Complex16*)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));

      exact_two_body_density_matrix=(MKL_Complex16*)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));

      average_wave_function_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
      average_wave_function_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

      wf_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
      wf_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

      trial_determinant_coefficients_energy=(MKL_Complex16*)malloc(ist.n_determinants_trial_energy*sizeof(MKL_Complex16));
      trial_wf_up_energy=(MKL_Complex16*)malloc(ist.n_determinants_energy_n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
      trial_wf_down_energy=(MKL_Complex16*)malloc(ist.n_determinants_energy_n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

      trial_determinant_coefficients_phaseless=(MKL_Complex16*)malloc(ist.n_determinants_trial_phaseless*sizeof(MKL_Complex16)); 
      trial_wf_up_phaseless=(MKL_Complex16*)malloc(ist.n_determinants_phaseless_n_spatial_orbitals_n_up*sizeof(MKL_Complex16)); 
      trial_wf_down_phaseless=(MKL_Complex16*)malloc(ist.n_determinants_phaseless_n_spatial_orbitals_n_down*sizeof(MKL_Complex16)); 

      trial_density_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
      trial_density_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

      kinetic_full_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
      kinetic_backwards_half_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
      kinetic_forwards_half_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));

      kinetic_eigs_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals*sizeof(double));
      kinetic_eigvecs_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));

      kinetic_full_fermions_real_down=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
      kinetic_backwards_half_fermions_real_down=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
      kinetic_forwards_half_fermions_real_down=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));

      kinetic_eigs_fermions_real_down=(double*)malloc(ist.n_spatial_orbitals*sizeof(double));
      kinetic_eigvecs_fermions_real_down=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));

      energy_shifts=(MKL_Complex16 *)malloc(4*sizeof(MKL_Complex16));
  
      kinetic_matrix_original_fermions_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

      potential_matrix_original_fermions_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals_fourth*sizeof(MKL_Complex16));
      potential_matrix_original_fermions_spin=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16)); 
      potential_matrix_original_fermions_down=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals_fourth*sizeof(MKL_Complex16)); 

      potential_supermatrix_fermions=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));

      potential_eigs_fermions=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
      potential_eigvecs_fermions=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));

      potential_eigs_fermions_first_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
      potential_eigs_fermions_second_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
  
      potential_one_body_matrix_up = (double *)calloc(ist.n_spatial_orbitals_sq,sizeof(double));

      mean_field_density_first_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
      mean_field_density_second_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));

      mean_field_first_transform_constant=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
      mean_field_second_transform_constant=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));

      /*For Population Control*/
      new_wf_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
      new_wf_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

      new_overlap_inverse_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_energy_n_up_sq*sizeof(MKL_Complex16));
      new_overlap_inverse_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_energy_n_down_sq*sizeof(MKL_Complex16));

      new_overlap_total=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));
    
      new_det_overlap_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*sizeof(MKL_Complex16)); 
      new_det_overlap_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*sizeof(MKL_Complex16));  

      /*List Orbitals*/ 
      list_orbitals_trial_energy = (int *)malloc(ist.n_spatial_orbitals*sizeof(int));

      /*For Orthogonalization*/
      R_up=(MKL_Complex16 *)malloc(ist.n_up_sq*sizeof(MKL_Complex16));
      R_down=(MKL_Complex16 *)malloc(ist.n_down_sq*sizeof(MKL_Complex16));

      /*If There Is Back Propagation*/
      if ( ist.flag_back_propagation == 1 ) {
        potential_matrices_back_propagation_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_sq*ist.n_steps_back_propagation*sizeof(MKL_Complex16));
        potential_matrices_back_propagation_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_sq*ist.n_steps_back_propagation*sizeof(MKL_Complex16));

        bp_wf_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
        bp_wf_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

        prev_wf_up=(MKL_Complex16*)malloc(ist.n_steps_back_propagation*cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
        prev_wf_down=(MKL_Complex16*)malloc(ist.n_steps_back_propagation*cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));
      }

      /******************************************************************************************/

    /*Zero the Total Energy Vectors*/
    czero_vec(total_potential_energy,3);
    czero_vec(total_kinetic_energy,3);
    czero_vec(total_energy_2,3);
    czero_vec(total_weight,1);

    /*Obtain Trial Wave Function, Trial Energy, and Trial Density Matrices from WF*/
    ist.n_max_orbital_trial_energy = init_wf_chemistry_multi(trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,list_orbitals_trial_energy,ist,cns);

    /*Input the Vijkl Elements and Resort Them to Find Eigenvalues*/
    get_potential_matrix_elements_molpro_restricted(kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,potential_matrix_original_fermions_spin,energy_shifts,list_orbitals_trial_energy,ist,cns,"potential_matrix_elements.par");

    if ( ist.flag_restart == 0 ) {
      init_walkers_chemistry_no_restart_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,weights,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weight_rescaling,ist,cns);
    }
    else {
      init_walkers_chemistry_restart_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,weights,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weight_rescaling,ist,cns);
    }

    /*Obtain Initial Trial Energy*/
    cns.trial_energy = get_trial_energy_density_multi_restricted(trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,trial_density_up,trial_density_down,exact_two_body_density_matrix,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,potential_matrix_original_fermions_spin,energy_shifts,ist);
    fprintf(pout, "trial energy %f\n", cns.trial_energy+energy_shifts[0].real); fflush(pout);
    cns.init_trial_energy = cns.trial_energy;
    cns.exp_trial_energy = exp(cns.dtau * cns.trial_energy);

    /*Form the Supermatrix and Diagonalize It to Use in Transforms*/
    init_potential_chemistry_restricted(potential_matrix_original_fermions_up,potential_matrix_original_fermions_spin,potential_supermatrix_fermions,potential_eigs_fermions,potential_eigvecs_fermions,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_one_body_matrix_up,ist,cns);

    /*If A Mean Field Is Subtracted from the Transform, Then Need Mean Fields to Subtract*/
    if ( ist.flag_meanfield == 1 ) {
      get_mean_field_densities(trial_density_up,trial_density_down,mean_field_density_first_transform,mean_field_density_second_transform,potential_eigvecs_fermions,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,ist);
    }

    /*Obtain Matrices, Wavefunctions, and Initialize Walkers*/
    init_kinetic_chemistry_real_restricted(kinetic_matrix_original_fermions_up,potential_supermatrix_fermions,kinetic_full_fermions_real_up,kinetic_backwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_up,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,potential_eigvecs_fermions,potential_eigs_fermions,mean_field_density_first_transform,mean_field_density_second_transform,potential_one_body_matrix_up,ist,cns);

    /*First Regress the Walkers By Half a Kinetic Propagator*/
    propagate_half_backwards_kinetic_chemistry_real_multi_restricted(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

    /*Production Phase*************************************************************************/
    for (steps = 0; steps < ist.n_steps; steps++) {

      if ( steps%20== 0 ) {
        fprintf(pout, "steps: %d\n", steps); fflush(pout);
      }

      /*If Back Propagation Possible*/
      if ( ist.flag_back_propagation == 0 ) { /**************************************************************************************/

        /*First Propagate a Full Step Forwards*/
        propagate_forwards_kinetic_chemistry_real_multi_restricted(kinetic_full_fermions_real_up,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,ist);

        /*Then Propagate Forwards a Full U Step******************************************************************************************/
        /*Use Discrete or Continuous Transform*/
        if ( ist.flag_real_eigenvectors == 0 ) { /*Determine If the Eigenctors Are Real or Not*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {
              propagate_forwards_potential_continuous_mixed_chemistry_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum);
            } /*Mean Field Flag*/
            else { /*If Mean Field*/
              propagate_forwards_potential_continuous_mixed_chemistry_meanfield_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,ist,cns,&idum);
           }
          } /*If Phaseless*/
          else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
                propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum);
            }
            else {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_multi_restricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,ist,cns,&idum);
            }
          }   /*If Mean Not Field*/
          else { /*If Mean Field*/
           if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_meanfield_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
           }
           else {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_multi_restricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
           }
          } /*If Mean Field*/
         } /*Phaseless Else*/
        }
        else {   /*If Eigenvectors Are Real*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {
                propagate_forwards_potential_continuous_mixed_chemistry_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum);
             } /*If Not Mean Field*/
             else { /*If Mean Field*/
                 propagate_forwards_potential_continuous_mixed_chemistry_meanfield_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,ist,cns,&idum);
             }
           }
           else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,ist,cns,&idum);
             }
            else {
              propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_real_multi_restricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,ist,cns,&idum);
            }
           }
           else { /*If Mean Field*/
             if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_meanfield_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
             }
             else {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_restricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
             }
           }
          } /*Phaseless Else*/
        } /*End Eigenvectors Real********************************************************************************************************/

         /*Normalize the Walker Weights by the Trial Energy*/
         normalize(weights,ist,cns); 

        /*Now Collect Energies and Other Observables*/
        if ( steps%ist.n_steps_energy==0 ) {
          /*First Propagate One Half Step Forward*/
          propagate_half_forwards_kinetic_chemistry_real_multi_restricted(kinetic_forwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,weights,ist);

          /*If You Want to Calculate 2 RDMS*/
          if ( ist.flag_compute_2RDM == 0 ) {
              /*Compute Energ*/

              compute_energy_chemistry_mixed_multi_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);


          }
          else {   /*2RDMS*******************************************************************************************************************************************************/

            /*Compute Energy and Density Matrices*/
            compute_energy_chemistry_mixed_two_body_density_matrix_multi_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_total,det_overlap_up,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down);

           /*Print RDMs******************************/
           fprintf(pf2, "%d\n", steps);
           fprintf(pf3, "%d\n", steps); 
           for (i=0; i<ist.n_spatial_orbitals; i++) {
            for (j=0; j<ist.n_spatial_orbitals; j++) {
              if ( Cabs(average_density_matrix_up[i*ist.n_spatial_orbitals+j]) > .0001 || Cabs(average_density_matrix_down[i*ist.n_spatial_orbitals+j]) > .0001 ) {
                 fprintf(pf2, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_up[i*ist.n_spatial_orbitals+j].real, average_density_matrix_up[i*ist.n_spatial_orbitals+j].imag);
                  fprintf(pf3, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_down[i*ist.n_spatial_orbitals+j].real, average_density_matrix_down[i*ist.n_spatial_orbitals+j].imag);
              }
             }
            }
            fprintf(pf2, "\n\n"); fflush(pf2);
            fprintf(pf3, "\n\n"); fflush(pf3); 

            /*Print 2-Body RDM*/
            fprintf(pf4, "%d\n", steps);
            for (i=0; i<ist.n_spin_orbitals; i++) {
             for (j=0; j<ist.n_spin_orbitals; j++) {
              for (k=0; k<ist.n_spin_orbitals; k++) {
               for (l=0; l<ist.n_spin_orbitals; l++) {
                if ( Cabs(average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l] ) > .0001 ) {  
                 fprintf(pf4, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].imag);
                }
               }
              }
             }
            }
            fprintf(pf4, "\n\n"); fflush(pf4);

          }
          cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));

          /*Print Energies*/
          fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, energy_shifts[0].real, te_mixed.real, te_mixed_error.real,tw_mixed.real, tw_mixed.imag); fflush(pf);
          fprintf(pf9, "%f\t\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9);

          /*Print Wave Functions*/
          fprintf(pf5, "%d\n", steps);
          fprintf(pf6, "%d\n", steps);
          for (i=0; i<ist.n_spatial_orbitals; i++) {
            for (j=0; j<ist.n_up; j++) {
             fprintf(pf5, "%f+%fi\t", average_wave_function_up[i*ist.n_up+j].real, average_wave_function_up[i*ist.n_up+j].imag);
           }
           for (j=0; j<ist.n_down; j++) {
             fprintf(pf6, "%f+%fi\t", average_wave_function_down[i*ist.n_down+j].real, average_wave_function_down[i*ist.n_down+j].imag);
           }
           fprintf(pf5, "\n");
           fprintf(pf6, "\n");
          }
          fprintf(pf5, "\n\n"); fflush(pf5);
          fprintf(pf6, "\n\n"); fflush(pf6);

          /*Propagate Back Backwards*/
          propagate_half_backwards_kinetic_chemistry_real_multi_restricted(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

        }/*Energy Calculations***************************************************************************************************/

        /*Orthogonalize Q and R if Necessary*/
        if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
          orthogonalize_chemistry_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,weights,R_up,R_down,ist);  
        }

        /*Apply Population Control If Necessary*/
        if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
            ist.n_walkers = population_control_chemistry_multi(weights,wf_up,wf_down,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,new_walker_weights,new_wf_up,new_wf_down,new_overlap_inverse_up,new_overlap_inverse_down,new_overlap_total,new_det_overlap_up,new_det_overlap_down,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns);
        }

     } /*If No Back Propagation*/
     else { /*If Back Propagation*/

        /*First Propagate a Full Step Forwards*/
        propagate_forwards_kinetic_chemistry_real_multi_restricted(kinetic_full_fermions_real_up,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,ist);

        /*Then Propagate Forwards a Full U Step******************************************************************************************/
        /*Use Discrete or Continuous Transform*/
        if ( ist.flag_real_eigenvectors == 0 ) { /*Determine If the Eigenctors Are Real or Not*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {   
              propagate_forwards_potential_continuous_bp_chemistry_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            } /*Mean Field Flag*/
            else { /*If Mean Field*/
              propagate_forwards_potential_continuous_bp_chemistry_meanfield_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
           }
          } /*If Phaseless*/
          else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
                propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            }
            else {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_multi_restricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            }
          }   /*If Mean Not Field*/
          else { /*If Mean Field*/
           if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_meanfield_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
}
           else {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_meanfield_multi_restricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
           }
          } /*If Mean Field*/
         } /*Phaseless Else*/
        }
        else {   /*If Eigenvectors Are Real*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {   
                propagate_forwards_potential_continuous_bp_chemistry_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             } /*If Not Mean Field*/
             else { /*If Mean Field*/
                 propagate_forwards_potential_continuous_bp_chemistry_meanfield_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
           }
           else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {  
                 propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
            else {
              propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_real_multi_restricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            }
           }
           else { /*If Mean Field*/
             if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_meanfield_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
             else {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_meanfield_real_multi_restricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up, potential_matrix_original_fermions_up,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
           }
          } /*Phaseless Else*/
        } /*End Eigenvectors Real********************************************************************************************************/

         /*Store WF, If Needed for Back Prop*/
         store_prev_wf(wf_up,wf_down,prev_wf_up,prev_wf_down,kinetic_forwards_half_fermions_real_up,ist);

         /*Normalize the Walker Weights by the Trial Energy*/
         normalize(weights,ist,cns);

        /*Now Collect Energies and Other Observables*/
        if ( steps%ist.n_steps_energy==0 ) {
          /*First Propagate One Half Step Forward*/

          propagate_half_forwards_kinetic_chemistry_real_multi_restricted(kinetic_forwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);


          /*If You Want to Calculate 2 RDMS*/
          if ( ist.flag_compute_2RDM == 0 ) {

            /*If You Want to Use Back Propagated Wave Function*/
            if (  steps < ist.n_steps_back_propagation ) {
              /*Compute Energ*/
              //compute_energy_chemistry_mixed_multi_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_up,overlap_down,det_overlap_total,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);

               //compute_energy_chemistry_restricted_exact(wf_up,wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_exact,&ke_exact_error,&pe_exact,&pe_exact_error,&te_exact,&te_exact_error,&tw_exact,average_wave_function_up,average_wave_function_down);

             cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));

             /*Print Energies*/
             fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, energy_shifts[0].real, te_mixed.real, te_mixed_error.real, te_exact.real, te_exact_error.real, tw_mixed.real, tw_mixed.imag); fflush(pf);
             fprintf(pf9, "%f\t\t %f\t\t %f\n", steps*cns.dtau, te_mixed.real, te_exact.real); fflush(pf9);

            }
            else {

              /*Back Propagate First*/
              //back_propagation_chemistry_restricted(trial_wf_up_phaseless,trial_wf_down_phaseless,bp_wf_up,bp_wf_down,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,kinetic_backwards_half_fermions_real_up,kinetic_full_fermions_real_up,kinetic_forwards_half_fermions_real_up,ist,cns);

              /*Compute Energy Second*/
              //compute_energy_chemistry_restricted_bp(wf_up,wf_down,bp_wf_up,bp_wf_down,prev_wf_up,prev_wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_bp,&ke_bp_error,&pe_bp,&pe_bp_error,&te_bp,&te_bp_error,&tw_bp,average_wave_function_up,average_wave_function_down);

              compute_energy_chemistry_mixed_multi_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);

              //compute_energy_chemistry_restricted_exact(wf_up,wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_exact,&ke_exact_error,&pe_exact,&pe_exact_error,&te_exact,&te_exact_error,&tw_exact,average_wave_function_up,average_wave_function_down); 
cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));

             /*Print Energies*/
             fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, energy_shifts[0].real, te_mixed.real, te_mixed_error.real, te_exact.real, te_exact.imag, tw_mixed.real, tw_mixed.imag); fflush(pf);
             fprintf(pf9, "%f\t\t %f\t\t %f\n", steps*cns.dtau, te_mixed.real, te_exact.real); fflush(pf9);

           }
          }
          else {   /*2RDMS*******************************************************************************************************************************************************/

           /*If You Want to Back Propagate Or Not*/
           if ( steps < ist.n_steps_back_propagation ) {

              /*Compute Energy and Density Matrices*/
              compute_energy_chemistry_mixed_two_body_density_matrix_multi_restricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_bp,&ke_bp_error,&pe_bp,&pe_bp_error,&te_bp,&te_bp_error,&tw_bp,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down);

             /*Print RDMs******************************/
             fprintf(pf2, "%d\n", steps);
             fprintf(pf3, "%d\n", steps);  
             for (i=0; i<ist.n_spatial_orbitals; i++) {
               for (j=0; j<ist.n_spatial_orbitals; j++) {
                if ( Cabs(average_density_matrix_up[i*ist.n_spatial_orbitals+j]) > .0001 || Cabs(average_density_matrix_down[i*ist.n_spatial_orbitals+j]) > .0001 ) { 

                  fprintf(pf2, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_up[i*ist.n_spatial_orbitals+j].real, average_density_matrix_up[i*ist.n_spatial_orbitals+j].imag);
                  fprintf(pf3, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_down[i*ist.n_spatial_orbitals+j].real, average_density_matrix_down[i*ist.n_spatial_orbitals+j].imag);
                }
               }
              }
              fprintf(pf2, "\n\n"); fflush(pf2);
              fprintf(pf3, "\n\n"); fflush(pf3); 

              /*Print 2-Body RDM*/
              fprintf(pf4, "%d\n", steps);
              for (i=0; i<ist.n_spin_orbitals; i++) {
               for (j=0; j<ist.n_spin_orbitals; j++) {
                for (k=0; k<ist.n_spin_orbitals; k++) {
                 for (l=0; l<ist.n_spin_orbitals; l++) {
                   if ( Cabs(average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]) > .0001 ) {
                      fprintf(pf4, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].imag);
                   }
                 }
                }
               }
              }
              fprintf(pf4, "\n\n"); fflush(pf4);

              cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));

              /*Print Energies*/
              fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, energy_shifts[0].real, te_mixed.real, te_mixed_error.real, tw_mixed.real, tw_mixed.imag); fflush(pf);
              fprintf(pf9, "%f\t\t %f\n", steps*cns.dtau, te_mixed.real); fflush(pf9);

          }
          else {

              /*Back Propagate First*/
              back_propagation_chemistry_multi(trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,bp_wf_up,bp_wf_down,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,kinetic_backwards_half_fermions_real_up,kinetic_full_fermions_real_up,kinetic_forwards_half_fermions_real_up,ist,cns);

              /*Compute Energy and Density Matrices*/
              compute_energy_chemistry_bp_two_body_density_matrix_restricted(wf_up,wf_down,bp_wf_up,bp_wf_down,prev_wf_up,prev_wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_bp,&ke_bp_error,&pe_bp,&pe_bp_error,&te_bp,&te_bp_error,&tw_bp,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down);

             /*Print RDMs******************************/
             fprintf(pf2, "%d\n", steps);
             fprintf(pf3, "%d\n", steps); 
             for (i=0; i<ist.n_spatial_orbitals; i++) {
              for (j=0; j<ist.n_spatial_orbitals; j++) {
               if ( Cabs(average_density_matrix_up[i*ist.n_spatial_orbitals+j]) > .0001 || Cabs(average_density_matrix_down[i*ist.n_spatial_orbitals+j]) > .0001 ) {
                     fprintf(pf2, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_up[i*ist.n_spatial_orbitals+j].real, average_density_matrix_up[i*ist.n_spatial_orbitals+j].imag);
                     fprintf(pf3, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_down[i*ist.n_spatial_orbitals+j].real, average_density_matrix_down[i*ist.n_spatial_orbitals+j].imag);
                }
               }
              }
              fprintf(pf2, "\n\n"); fflush(pf2);
              fprintf(pf3, "\n\n"); fflush(pf3); 

              /*Print 2-Body RDM*/
              fprintf(pf4, "%d\n", steps);
              for (i=0; i<ist.n_spin_orbitals; i++) {
               for (j=0; j<ist.n_spin_orbitals; j++) {
                for (k=0; k<ist.n_spin_orbitals; k++) {
                 for (l=0; l<ist.n_spin_orbitals; l++) {
                  if ( Cabs(average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]) > .0001 ) {
                   fprintf(pf4, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].imag);
                  }
                 }
                }
               }
              }
              fprintf(pf4, "\n\n"); fflush(pf4);

             cns.exp_trial_energy = exp(cns.dtau * (ke_bp.real + pe_bp.real));

            /*Print Energies*/
            fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_bp.real, pe_bp_error.real, ke_bp.real, ke_bp_error.real, energy_shifts[0].real, te_bp.real, te_bp_error.real, tw_bp.real, tw_bp.imag); fflush(pf);
            fprintf(pf9, "%f\t\t %f\n", steps*cns.dtau, te_bp.real); fflush(pf9);

           } /*Else BP and two RDMS*/
          }

          /*Print Wave Functions*/
          fprintf(pf5, "%d\n", steps);
          fprintf(pf6, "%d\n", steps);
          for (i=0; i<ist.n_spatial_orbitals; i++) {
            for (j=0; j<ist.n_up; j++) {
             fprintf(pf5, "%f+%fi\t", average_wave_function_up[i*ist.n_up+j].real, average_wave_function_up[i*ist.n_up+j].imag);
           }
           for (j=0; j<ist.n_down; j++) {
             fprintf(pf6, "%f+%fi\t", average_wave_function_down[i*ist.n_down+j].real, average_wave_function_down[i*ist.n_down+j].imag);
           }
           fprintf(pf5, "\n");
           fprintf(pf6, "\n");
          }
          fprintf(pf5, "\n\n"); fflush(pf5);
          fprintf(pf6, "\n\n"); fflush(pf6);

          /*Propagate Back Backwards*/
           propagate_half_backwards_kinetic_chemistry_real_multi_restricted(kinetic_backwards_half_fermions_real_up,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

        }/*Energy Calculations***************************************************************************************************/

        /*Copy Wave Functions Back If Using Back Propagation*/
        copy_back_propagation_matrices_forward(potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns);
        copy_prev_wf_forward(prev_wf_up,prev_wf_down,ist,cns);

        /*Orthogonalize Q and R if Necessary*/
        if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
          orthogonalize_chemistry_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,weights,R_up,R_down,ist);
        }

        /*Apply Population Control If Necessary*/
         if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
           ist.n_walkers = population_control_chemistry_multi(weights,wf_up,wf_down,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,new_walker_weights,new_wf_up,new_wf_down,new_overlap_inverse_up,new_overlap_inverse_down,new_overlap_total,new_det_overlap_up,new_det_overlap_down,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns);
         }

     } /*If Back Propagation */
    } /*Steps*/

    /******************************************************************************************/

    free(weights);
    free(new_walker_weights);

    free(total_potential_energy); free(total_kinetic_energy); free(total_energy_2); free(total_weight);

    free(energy_shifts);
    free(kinetic_matrix_original_fermions_up);
    free(kinetic_matrix_original_fermions_down); 

    free(potential_matrix_original_fermions_up);
    free(potential_matrix_original_fermions_down); 
    free(potential_matrix_original_fermions_spin);
 
    free(potential_supermatrix_fermions);
    free(potential_eigs_fermions); free(potential_eigvecs_fermions);
    free(potential_eigs_fermions_first_transform); free(potential_eigs_fermions_second_transform);
    free(potential_one_body_matrix_up); 

    free(mean_field_density_first_transform);  free(mean_field_density_second_transform);
    free(mean_field_first_transform_constant); free(mean_field_second_transform_constant);

    free(overlap_total);
    free(det_overlap_up); free(det_overlap_down); 
    free(overlap_inverse_up); free(overlap_inverse_down);
    free(wf_up); free(wf_down);

    free(trial_determinant_coefficients_energy); free(trial_determinant_coefficients_phaseless); 
 
    free(trial_wf_up_energy); free(trial_wf_down_energy);
    free(trial_wf_up_phaseless); free(trial_wf_down_phaseless); 

    free(trial_density_up); free(trial_density_down);

    free(kinetic_full_fermions_real_up); free(kinetic_backwards_half_fermions_real_up); free(kinetic_forwards_half_fermions_real_up);
    free(kinetic_eigs_fermions_real_up); free(kinetic_eigvecs_fermions_real_up);

    free(kinetic_full_fermions_real_down); free(kinetic_backwards_half_fermions_real_down); free(kinetic_forwards_half_fermions_real_down); 
    free(kinetic_eigs_fermions_real_down); free(kinetic_eigvecs_fermions_real_down); 

    free(new_overlap_total);
    free(new_det_overlap_up); free(new_det_overlap_down); 
    free(new_overlap_inverse_up); free(new_overlap_inverse_down);
    free(new_wf_up); free(new_wf_down);

    free(average_density_matrix_up); free(average_density_matrix_down);
    free(average_two_body_density_matrix);

    free(exact_two_body_density_matrix); 

    free(average_wave_function_up); free(average_wave_function_down);

    free(R_up); free(R_down);

    if ( ist.flag_back_propagation == 1 ) {
       free(bp_wf_up); free(bp_wf_down);
       free(prev_wf_up); free(prev_wf_down);

       free(potential_matrices_back_propagation_up);
       free(potential_matrices_back_propagation_down);
    }

   fclose(pf); fclose(pf9);
   }
   else if ( ist.flag_simulation_type == 9 ) { /*Unrestricted with Many Determinants******************************************************************************/

      /*Determine Number Kinetic and Potential Sparse Terms*/
      get_number_matrix_elements_molpro_unrestricted(&number_kinetic_sparse_up,&number_kinetic_sparse_down,&number_potential_sparse_up,&number_potential_sparse_down,&number_potential_sparse_updown,"potential_matrix_elements.par");  
      ist.n_kinetic_sparse_up = number_kinetic_sparse_up; 
      ist.n_kinetic_sparse_down = number_kinetic_sparse_down; 
      ist.n_potential_sparse_up = number_potential_sparse_up; 
      ist.n_potential_sparse_down = number_potential_sparse_down; 
      ist.n_potential_sparse_updown = number_potential_sparse_updown; 

      accumulated_energy.real = accumulated_energy.imag = 0.0; 
      accumulated_weight.real = accumulated_weight.imag = 0.0;  

      /*DECLARE********************************************************************************************/
        /*Weights of Walkers*/
        weights=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));
        energies=(MKL_Complex16*)malloc(3*cns.max_number_walkers*sizeof(MKL_Complex16)); 
        energies2=(MKL_Complex16*)malloc(3*cns.max_number_walkers*sizeof(MKL_Complex16)); 
        energies3=(MKL_Complex16*)malloc(3*cns.max_number_walkers*sizeof(MKL_Complex16)); 
        phases=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16)); 
        phases2=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16)); 
        overlaps=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16)); 
        local_energy_flag = (int *)malloc(cns.max_number_walkers*sizeof(int )); 
        total_energy_flag = (int *)malloc(cns.max_number_walkers*sizeof(int )); 
        field_flag = (int *)malloc(cns.max_number_walkers*sizeof(int ));  
        density_matrices_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*cns.max_number_walkers*sizeof(MKL_Complex16)); 
        density_matrices_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*cns.max_number_walkers*sizeof(MKL_Complex16)); 

        /*For Population Control*/
        new_walker_weights = (MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));
        weight_rescaling=(double*)malloc(10*sizeof(double));
        walkers_big = (int*)malloc(cns.max_number_walkers*sizeof(int));
        walkers_small = (int*)malloc(cns.max_number_walkers*sizeof(int));
        walkers_copied = (int*)malloc(cns.max_number_walkers*sizeof(int));

        total_energy_2=(MKL_Complex16*)malloc(3*sizeof(MKL_Complex16));
        total_potential_energy=(MKL_Complex16*)malloc(3*sizeof(MKL_Complex16));
        total_kinetic_energy=(MKL_Complex16*)malloc(3*sizeof(MKL_Complex16));
        total_weight=(MKL_Complex16*)malloc(1*sizeof(MKL_Complex16));

        det_overlap_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*sizeof(MKL_Complex16));
        det_overlap_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*sizeof(MKL_Complex16));

        overlap_total=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));

        overlap_inverse_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*ist.n_up_sq*sizeof(MKL_Complex16));
        overlap_inverse_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*ist.n_down_sq*sizeof(MKL_Complex16));

        average_density_matrix_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
        average_density_matrix_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
        average_two_body_density_matrix=(MKL_Complex16*)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));

        exact_two_body_density_matrix=(MKL_Complex16*)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));

        average_wave_function_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
        average_wave_function_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

        wf_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
        wf_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

        trial_determinant_coefficients_energy=(MKL_Complex16*)malloc(ist.n_determinants_trial_energy*sizeof(MKL_Complex16));
        trial_wf_up_energy=(MKL_Complex16*)malloc(ist.n_determinants_energy_n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
        trial_wf_down_energy=(MKL_Complex16*)malloc(ist.n_determinants_energy_n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

        trial_determinant_coefficients_phaseless=(MKL_Complex16*)malloc(ist.n_determinants_trial_phaseless*sizeof(MKL_Complex16));
        trial_wf_up_phaseless=(MKL_Complex16*)malloc(ist.n_determinants_phaseless_n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
        trial_wf_down_phaseless=(MKL_Complex16*)malloc(ist.n_determinants_phaseless_n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

        trial_density_up=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
        trial_density_down=(MKL_Complex16*)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));

        kinetic_full_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
        kinetic_backwards_half_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
        kinetic_forwards_half_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));

        kinetic_eigs_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals*sizeof(double));
        kinetic_eigvecs_fermions_real_up=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
 
        kinetic_full_fermions_real_down=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
        kinetic_backwards_half_fermions_real_down=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
        kinetic_forwards_half_fermions_real_down=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));
 
        kinetic_eigs_fermions_real_down=(double*)malloc(ist.n_spatial_orbitals*sizeof(double));
        kinetic_eigvecs_fermions_real_down=(double*)malloc(ist.n_spatial_orbitals_sq*sizeof(double));

        energy_shifts=(MKL_Complex16 *)malloc(4*sizeof(MKL_Complex16));

        kinetic_matrix_original_fermions_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16));
        kinetic_matrix_original_fermions_down=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals_sq*sizeof(MKL_Complex16)); 

        kinetic_matrix_sparse_up=(double *)malloc(ist.n_kinetic_sparse_up*sizeof(double)); 
        kinetic_matrix_sparse_down=(double *)malloc(ist.n_kinetic_sparse_down*sizeof(double)); 

        kinetic_ij_sparse_up=(int*)malloc(ist.n_kinetic_sparse_up*sizeof(int)); 
        kinetic_ij_sparse_down=(int*)malloc(ist.n_kinetic_sparse_down*sizeof(int)); 

        potential_matrix_original_fermions_up=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals_fourth*sizeof(MKL_Complex16));
        potential_matrix_original_fermions_down=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals_fourth*sizeof(MKL_Complex16)); 
        potential_matrix_original_fermions_updown=(MKL_Complex16 *)malloc(ist.n_spatial_orbitals_fourth*sizeof(MKL_Complex16)); 

        potential_matrix_sparse_up=(double *)malloc(ist.n_potential_sparse_up*sizeof(double)); 
        potential_matrix_sparse_down=(double *)malloc(ist.n_potential_sparse_down*sizeof(double)); 
        potential_matrix_sparse_updown=(double *)malloc(2*ist.n_potential_sparse_updown*sizeof(double)); 

        potential_ij_sparse_up=(int *)malloc(4*ist.n_potential_sparse_up*sizeof(int)); 
        potential_ij_sparse_down=(int *)malloc(4*ist.n_potential_sparse_down*sizeof(int)); 
        potential_ij_sparse_updown=(int *)malloc(2*ist.n_potential_sparse_updown*sizeof(int)); 

        potential_matrix_original_fermions_spin=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));
        potential_supermatrix_fermions=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));
        potential_eigs_fermions=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
        potential_eigvecs_fermions=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_fourth*sizeof(MKL_Complex16));

        potential_eigs_fermions_first_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
        potential_eigs_fermions_second_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));

        potential_one_body_matrix_up = (double *)calloc(ist.n_spatial_orbitals_sq,sizeof(double));
        potential_one_body_matrix_down = (double *)calloc(ist.n_spatial_orbitals_sq,sizeof(double)); 
 
        mean_field_density_first_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
        mean_field_density_second_transform=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));

        mean_field_first_transform_constant=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));
        mean_field_second_transform_constant=(MKL_Complex16 *)malloc(ist.n_spin_orbitals_sq*sizeof(MKL_Complex16));

        list_orbitals_trial_energy = (int *)malloc(ist.n_spatial_orbitals*sizeof(int)); 

        /*For Population Control*/
        new_wf_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
        new_wf_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

        new_overlap_inverse_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_energy_n_up_sq*sizeof(MKL_Complex16));
        new_overlap_inverse_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_energy_n_down_sq*sizeof(MKL_Complex16));

        new_overlap_total=(MKL_Complex16*)malloc(cns.max_number_walkers*sizeof(MKL_Complex16));

        new_det_overlap_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*sizeof(MKL_Complex16));
        new_det_overlap_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_determinants_trial_energy*sizeof(MKL_Complex16));

        /*For Orthogonalization*/
        R_up=(MKL_Complex16 *)malloc(ist.n_up_sq*sizeof(MKL_Complex16));
        R_down=(MKL_Complex16 *)malloc(ist.n_down_sq*sizeof(MKL_Complex16));

        /*If There Is Back Propagation*/
        if ( ist.flag_back_propagation == 1 ) {
          potential_matrices_back_propagation_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_sq*ist.n_steps_back_propagation*sizeof(MKL_Complex16));
          potential_matrices_back_propagation_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_sq*ist.n_steps_back_propagation*sizeof(MKL_Complex16));

          bp_wf_up=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
          bp_wf_down=(MKL_Complex16*)malloc(cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));

          prev_wf_up=(MKL_Complex16*)malloc(ist.n_steps_back_propagation*cns.max_number_walkers*ist.n_spatial_orbitals_n_up*sizeof(MKL_Complex16));
          prev_wf_down=(MKL_Complex16*)malloc(ist.n_steps_back_propagation*cns.max_number_walkers*ist.n_spatial_orbitals_n_down*sizeof(MKL_Complex16));
        }

       /******************************************************************************************/

    /*Zero the Total Energy Vectors*/
    czero_vec(total_potential_energy,3);
    czero_vec(total_kinetic_energy,3);
    czero_vec(total_energy_2,3);
    czero_vec(total_weight,1);


    /*Obtain Trial Wave Function, Trial Energy, and Trial Density Matrices from WF*/
    ist.n_max_orbital_trial_energy = init_wf_chemistry_multi(trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,list_orbitals_trial_energy,ist,cns);

    /*Input the Vijkl Elements and Resort Them to Find Eigenvalues*/
    get_potential_matrix_elements_molpro_unrestricted(kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,potential_matrix_original_fermions_spin,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,list_orbitals_trial_energy,energy_shifts,&number_kinetic_sparse_up,&number_kinetic_sparse_down,&number_potential_sparse_up,&number_potential_sparse_down,&number_potential_sparse_updown,ist,cns,"potential_matrix_elements.par");

    ist.n_kinetic_sparse_up = number_kinetic_sparse_up;
    ist.n_kinetic_sparse_down = number_kinetic_sparse_down;
    ist.n_potential_sparse_up = number_potential_sparse_up;
    ist.n_potential_sparse_down = number_potential_sparse_down;
    ist.n_potential_sparse_updown = number_potential_sparse_updown;
    fprintf(pout, "number kinetic %d %d number potential %d %d %d\n", number_kinetic_sparse_up,number_kinetic_sparse_down,number_potential_sparse_up,number_potential_sparse_down,number_potential_sparse_updown); fflush(pout);  

    if ( ist.flag_restart == 0 ) {
      init_walkers_chemistry_no_restart_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,weights,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weight_rescaling,ist,cns);
    }
    else {
      init_walkers_chemistry_restart_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,weights,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weight_rescaling,ist,cns);
    }

    /*Obtain Initial Trial Energy*/
    cns.trial_energy = get_trial_energy_density_multi_unrestricted(trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,trial_density_up,trial_density_down,exact_two_body_density_matrix,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,potential_matrix_original_fermions_spin,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,energy_shifts,ist);
    cns.total_energy_shift = energy_shifts[0].real + energy_shifts[1].real + energy_shifts[2].real + energy_shifts[3].real; 
    fprintf(pout, "trial energy %f\n", cns.trial_energy+cns.total_energy_shift); fflush(pout); 
    cns.init_trial_energy = cns.trial_energy;
    cns.exp_trial_energy = exp(cns.dtau * cns.trial_energy);

    /*Form the Supermatrix and Diagonalize It to Use in Transforms*/
    init_potential_chemistry_unrestricted(potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,potential_supermatrix_fermions,potential_eigs_fermions,potential_eigvecs_fermions,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_one_body_matrix_up,potential_one_body_matrix_down,ist,cns);

    /*If A Mean Field Is Subtracted from the Transform, Then Need Mean Fields to Subtract*/
    if ( ist.flag_meanfield == 1 ) {
      get_mean_field_densities(trial_density_up,trial_density_down,mean_field_density_first_transform,mean_field_density_second_transform,potential_eigvecs_fermions,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,ist);
    }

    /*Obtain Matrices, Wavefunctions, and Initialize Walkers*/
    init_kinetic_chemistry_real_unrestricted(kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_supermatrix_fermions,kinetic_full_fermions_real_up,kinetic_backwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_up,kinetic_eigs_fermions_real_up,kinetic_eigvecs_fermions_real_up,kinetic_full_fermions_real_down,kinetic_backwards_half_fermions_real_down,kinetic_forwards_half_fermions_real_down,kinetic_eigs_fermions_real_down,kinetic_eigvecs_fermions_real_down,potential_eigvecs_fermions,potential_eigs_fermions,mean_field_density_first_transform,mean_field_density_second_transform,potential_one_body_matrix_up,potential_one_body_matrix_down,ist,cns); 

    /*First Regress the Walkers By Half a Kinetic Propagator*/
    propagate_half_backwards_kinetic_chemistry_real_multi_unrestricted(kinetic_backwards_half_fermions_real_up,kinetic_backwards_half_fermions_real_down,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

    /*Production Phase*************************************************************************/
    for (steps = 0; steps < ist.n_steps; steps++) {

      if ( steps%20== 0 ) {
        fprintf(pout, "steps: %d\n", steps); fflush(pout);
      }

      /*If Back Propagation Possible*/
      if ( ist.flag_back_propagation == 0 ) { /**************************************************************************************/

        /*First Propagate a Full Step Forwards*/
        propagate_forwards_kinetic_chemistry_real_multi_unrestricted(kinetic_full_fermions_real_up,kinetic_full_fermions_real_down,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,phases,ist);

        /*Then Propagate Forwards a Full U Step******************************************************************************************/
        /*Use Discrete or Continuous Transform*/
        if ( ist.flag_real_eigenvectors == 0 ) { /*Determine If the Eigenctors Are Real or Not*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {
              propagate_forwards_potential_continuous_mixed_chemistry_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum);
            } /*Mean Field Flag*/
            else { /*If Mean Field*/
              propagate_forwards_potential_continuous_mixed_chemistry_meanfield_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,ist,cns,&idum);
           }
          } /*If Phaseless*/
          else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
                propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum);
            }
            else {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_multi_unrestricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,ist,cns,&idum);
            }
          }   /*If Mean Not Field*/
          else { /*If Mean Field*/
           if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_meanfield_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
           }
           else {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_multi_unrestricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
           }
          } /*If Mean Field*/
         } /*Phaseless Else*/
        }
        else {   /*If Eigenvectors Are Real*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {
                propagate_forwards_potential_continuous_mixed_chemistry_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum);
             } /*If Not Mean Field*/
             else { /*If Mean Field*/
                 propagate_forwards_potential_continuous_mixed_chemistry_meanfield_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,ist,cns,&idum);
             }
           }
           else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,ist,cns,&idum);
             }
            else {
              propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_real_multi_unrestricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,ist,cns,&idum);
}
           }
           else { /*If Mean Field*/
             if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_mixed_chemistry_phaseless_squaredshift_meanfield_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
             }
             else {
               if ( steps < ist.n_steps_equilibration + ist.n_steps_energy ) {
                 propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted_sparse_equilibration(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,phases2,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
               }
               else { 
                 propagate_forwards_potential_continuous_mixed_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted_sparse_production(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,phases2,local_energy_flag,field_flag,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,mean_field_first_transform_constant,mean_field_second_transform_constant,ist,cns,&idum);
               } 
             }
           }
          } /*Phaseless Else*/
        } /*End Eigenvectors Real********************************************************************************************************/

         /*Normalize the Walker Weights by the Trial Energy*/
         normalize(weights,ist,cns);

        /*Now Collect Energies and Other Observables*/
        if ( steps%ist.n_steps_energy==0 ) {
          /*First Propagate One Half Step Forward*/
          propagate_half_forwards_kinetic_chemistry_real_multi_unrestricted(kinetic_forwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_down,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,weights,ist);

          /*If You Want to Calculate 2 RDMS*/
          if ( ist.flag_compute_2RDM == 0 ) {

            if ( steps < ist.n_steps_equilibration) { 

    printf("before calculate energy in no RDMs\n"); fflush(NULL); 

              /*Compute Energ*/
              compute_energy_chemistry_mixed_multi_unrestricted_sparse_equilibration(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,energies,energies2,energies3,overlaps,density_matrices_up,density_matrices_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&pe_mixed,&pe_mixed_error,&ke_mixed,&ke_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down); 

             /*Determine Accumulated Energy*/
             if ( steps > ist.n_steps_start ) {
               accumulated_energy = Cadd(accumulated_energy, Cmul(te_mixed, tw_mixed)); 
               accumulated_weight = Cadd(accumulated_weight, tw_mixed); 
             }

            }
            else { 

             if ( steps == ist.n_steps_equilibration ) {
                  accumulated_energy = Cdiv(accumulated_energy, accumulated_weight);
                  cns.energy_mean = accumulated_energy.real; 
                  cns.energy_max = cns.energy_mean + sqrt(cns.energy_cap_constant/cns.dtau);
                  cns.energy_min = cns.energy_mean - sqrt(cns.energy_cap_constant/cns.dtau);

                  fprintf(pout, "energy mean %f energy max %f energy min %f\n", cns.energy_mean+cns.total_energy_shift, cns.energy_max+cns.total_energy_shift, cns.energy_min+cns.total_energy_shift); fflush(pout);

               }

               /*Compute Energy Tossing Out Any Energies That Are Spurious Based Upon Mean*/
               compute_energy_chemistry_mixed_multi_unrestricted_sparse_production(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,energies,energies2,energies3,overlaps,total_energy_flag,density_matrices_up,density_matrices_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&pe_mixed,&pe_mixed_error,&ke_mixed,&ke_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);   

            }

          }
          else {   /*2RDMS*******************************************************************************************************************************************************/

            if ( steps < ist.n_steps_equilibration) {

               printf("before calculate energy in no RDMs\n"); fflush(NULL);

               /*Compute Energ*/
               compute_energy_chemistry_mixed_multi_unrestricted_sparse_equilibration(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,energies,energies2,energies3,overlaps,density_matrices_up,density_matrices_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&pe_mixed,&pe_mixed_error,&ke_mixed,&ke_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);

              /*Determine Accumulated Energy*/
              if ( steps > ist.n_steps_start ) {
                accumulated_energy = Cadd(accumulated_energy, Cmul(te_mixed, tw_mixed));
                accumulated_weight = Cadd(accumulated_weight, tw_mixed);
              }

            }
            else {

             if ( steps == ist.n_steps_equilibration ) {
                  accumulated_energy = Cdiv(accumulated_energy, accumulated_weight);
                  cns.energy_mean = accumulated_energy.real;
                  cns.energy_max = cns.energy_mean + sqrt(cns.energy_cap_constant/cns.dtau);
                  cns.energy_min = cns.energy_mean - sqrt(cns.energy_cap_constant/cns.dtau);

                  fprintf(pout, "energy mean %f energy max %f energy min %f\n", cns.energy_mean+cns.total_energy_shift, cns.energy_max+cns.total_energy_shift, cns.energy_min+cns.total_energy_shift); fflush(pout);

               }   

              /*Compute Energy and Density Matrices*/
              /*Right Now, This Routine Does Not Exist - 12-8-16*/
              compute_energy_chemistry_mixed_two_body_density_matrix_multi_unrestricted_sparse(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down);

            printf("in compute energy 2RDMs\n"); fflush(NULL); 

            //compute_energy_chemistry_mixed_multi_unrestricted_sparse_production(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,energies,energies2,energies3,overlaps,total_energy_flag,density_matrices_up,density_matrices_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&pe_mixed,&pe_mixed_error,&ke_mixed,&ke_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);

            printf("after compute energy 2RDMs\n"); fflush(NULL); 

           /*Print RDMs******************************/
           fprintf(pf2, "%d\n", steps);
           fprintf(pf3, "%d\n", steps);
           for (i=0; i<ist.n_spatial_orbitals; i++) {
            for (j=0; j<ist.n_spatial_orbitals; j++) {
              if ( Cabs(average_density_matrix_up[i*ist.n_spatial_orbitals+j]) > .0001 || Cabs(average_density_matrix_down[i*ist.n_spatial_orbitals+j]) > .0001 ) {
                 fprintf(pf2, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_up[i*ist.n_spatial_orbitals+j].real, average_density_matrix_up[i*ist.n_spatial_orbitals+j].imag);
                  fprintf(pf3, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_down[i*ist.n_spatial_orbitals+j].real, average_density_matrix_down[i*ist.n_spatial_orbitals+j].imag);
              }
             }
            }
            fprintf(pf2, "\n\n"); fflush(pf2);
            fprintf(pf3, "\n\n"); fflush(pf3);

            /*Print 2-Body RDM*/
            fprintf(pf4, "%d\n", steps);
            for (i=0; i<ist.n_spin_orbitals; i++) {
             for (j=0; j<ist.n_spin_orbitals; j++) {
              for (k=0; k<ist.n_spin_orbitals; k++) {
               for (l=0; l<ist.n_spin_orbitals; l++) {
                if ( Cabs(average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l] ) > .0001 ) {
                 fprintf(pf4, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].imag);
                }
               }
              }
             }
            }
            fprintf(pf4, "\n\n"); fflush(pf4);

          }
          cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));

         /*Print Energies*/
          fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, cns.total_energy_shift, te_mixed.real+cns.total_energy_shift, te_mixed_error.real,tw_mixed.real, tw_mixed.imag); fflush(pf);
          fprintf(pf9, "%f\t\t %f\n", steps*cns.dtau, te_mixed.real+cns.total_energy_shift); fflush(pf9); 

          fprintf(pf7, "%f\t\t %f+%fi\t %f+%fi\t %f+%fi\n", steps*cns.dtau, overlap_total[0].real, overlap_total[0].imag, det_overlap_up[0].real, det_overlap_up[0].imag, det_overlap_down[0].real, det_overlap_down[0].imag); fflush(pf7); 


         max_weight_real = 0; min_weight_real = 100000; max_overlap_real = 0; min_overlap_real = 1000000;
         max_weight.real = max_weight.imag = 0; min_weight.real = 100000; min_weight.imag = 0.0; max_overlap.real = max_overlap.imag = 0.0; min_overlap.real = 1000000.0; min_overlap.imag = 0.0;
        for (i=0; i<ist.n_walkers; i++) {

          if ( Cabs(weights[i]) > Cabs(max_weight) ) {
             max_weight.real = weights[i].real;
             max_weight.imag = weights[i].imag;
          }
          if ( Cabs(weights[i]) < Cabs(min_weight) ) {
            min_weight.real = weights[i].real;
            min_weight.imag = weights[i].imag;
          }
          if ( Cabs(overlap_total[i]) > Cabs(max_overlap) ) {
             max_overlap.real = overlap_total[i].real;
             max_overlap.imag = overlap_total[i].imag;
          }
          if ( Cabs(overlap_total[i]) < Cabs(min_overlap) ) {
             min_overlap.real = overlap_total[i].real;
             min_overlap.imag = overlap_total[i].imag;
          }

          if ( weights[i].real > max_weight_real ) {
             max_weight_real = weights[i].real;
          }
          if ( weights[i].real < min_weight_real ) {
              min_weight_real = weights[i].real;
          }
          if ( overlap_total[i].real > max_overlap_real ) {
             max_overlap_real = overlap_total[i].real;
          }
          if ( overlap_total[i].real < min_overlap_real ) {
             min_overlap_real = overlap_total[i].real;
          }

        }
        fprintf(pf25, "%d\t %f\t %f\t %f\t %f\t %f\t %f\t Real %f\t %f\t %f\t %f\n", steps, steps*cns.dtau, te_mixed.real+cns.total_energy_shift,max_weight.real, min_weight.real, max_overlap.real, min_overlap.real, max_weight_real, min_weight_real, max_overlap_real, min_overlap_real); fflush(pf25);

        fprintf(pf26, "%d\t %f\t %f\n", steps, steps*cns.dtau, te_mixed.real+cns.total_energy_shift); 
        fprintf(pf27, "%d\t %f\t %f\n", steps, steps*cns.dtau, te_mixed.real+cns.total_energy_shift);
        fprintf(pf28, "%d\t %f\t %f\n", steps, steps*cns.dtau, te_mixed.real+cns.total_energy_shift);
        fprintf(pf29, "%d\t %f\t %f\n", steps, steps*cns.dtau, te_mixed.real+cns.total_energy_shift);
        fprintf(pf30, "%d\t %f\t %f\n", steps, steps*cns.dtau, te_mixed.real+cns.total_energy_shift);
        fprintf(pf31, "%d\t %f\t %f\n", steps, steps*cns.dtau, te_mixed.real+cns.total_energy_shift); 
        fprintf(pf32, "%d\t %f\t %f\n", steps, steps*cns.dtau, te_mixed.real+cns.total_energy_shift);  
        fprintf(pf33, "%d\t %f\t %f\n", steps, steps*cns.dtau, te_mixed.real+cns.total_energy_shift);  
        for (i=0; i<ist.n_walkers; i++) {
          fprintf(pf26, "%d %f+%fi\t\t  %f+%fi\t %f+%fi\t %f+%fi\t\t  %f+%fi\t %f+%fi\t %f+%fi\t\t  %f+%fi\n", i, weights[i].real, weights[i].imag, energies[i*3].real, energies[i*3].imag, energies[i*3+1].real, energies[i*3+1].imag, energies[i*3+2].real, energies[i*3+2].imag, energies2[i*3].real, energies2[i*3].imag, energies2[i*3+1].real, energies2[i*3+1].imag, energies2[i*3+2].real, energies2[i*3+2].imag, overlaps[i].real, overlaps[i].imag);

          fprintf(pf27, "%d\n", i); 
          for (j=0; j<ist.n_spatial_orbitals_sq; j++) {
            fprintf(pf27, "%f+%fi\t", density_matrices_up[i*ist.n_spatial_orbitals_sq+j].real, density_matrices_up[i*ist.n_spatial_orbitals_sq+j].imag); 
          }
          fprintf(pf27, "\n"); fflush(pf27); 

          fprintf(pf28, "%d\n", i); 
          for (j=0; j<ist.n_spatial_orbitals_sq; j++) {
            fprintf(pf28, "%f+%fi\t", density_matrices_down[i*ist.n_spatial_orbitals_sq+j].real, density_matrices_down[i*ist.n_spatial_orbitals_sq+j].imag); 
          }
          fprintf(pf28, "\n"); fflush(pf28); 

          fprintf(pf29, "%d\n", i); 
          for (j=0; j<ist.n_spatial_orbitals_n_up; j++) {
            fprintf(pf29, "%f+%fi\t", wf_up[i*ist.n_spatial_orbitals_n_up+j].real, wf_up[i*ist.n_spatial_orbitals_n_up+j].imag); 
          }
          fprintf(pf29, "\n"); fflush(pf29); 

          fprintf(pf30, "%d\n", i);
          for (j=0; j<ist.n_spatial_orbitals_n_down; j++) {
            fprintf(pf30, "%f+%fi\t", wf_down[i*ist.n_spatial_orbitals_n_down+j].real, wf_down[i*ist.n_spatial_orbitals_n_down+j].imag);
          }
          fprintf(pf30, "\n"); fflush(pf30);

          if ( fabs(phases[i].real) > 3.1415*.5) {  
            fprintf(pf31, "%d\t %f\n", i, phases[i].real); fflush(pf31);  
          }         

         if ( fabs(phases2[3*i].real) > 3.1415*.5) {
         fprintf(pf32, "%d\t %f\t %f+%fi\t %f+%fi\n", i, phases2[i*3].real, phases2[i*3+1].real, phases2[i*3+1].imag, phases2[i*3+2].real, phases2[i*3+2].imag); 
          }

         if ( local_energy_flag[i]!=0 || total_energy_flag[i]!=0 || field_flag[i]!=0 ) {
           fprintf(pf33, "%d\t %d\t %d\t %d\n", i, local_energy_flag[i], total_energy_flag[i], field_flag[i]); fflush(pf33); 
         } 

        }
       fprintf(pf26, "\n"); fflush(pf26);  
       fprintf(pf31, "\n"); fflush(pf31); 
       fprintf(pf32, "\n"); fflush(pf32); 
       fprintf(pf33, "\n"); fflush(pf33); 



          /*Print Wave Functions*/
          fprintf(pf5, "%d\n", steps);
          fprintf(pf6, "%d\n", steps);
          for (i=0; i<ist.n_spatial_orbitals; i++) {
            for (j=0; j<ist.n_up; j++) {
             fprintf(pf5, "%f+%fi\t", average_wave_function_up[i*ist.n_up+j].real, average_wave_function_up[i*ist.n_up+j].imag);
           }
           for (j=0; j<ist.n_down; j++) {
             fprintf(pf6, "%f+%fi\t", average_wave_function_down[i*ist.n_down+j].real, average_wave_function_down[i*ist.n_down+j].imag);
           }
           fprintf(pf5, "\n");
           fprintf(pf6, "\n");
          }
          fprintf(pf5, "\n\n"); fflush(pf5);
          fprintf(pf6, "\n\n"); fflush(pf6);

        }

          /*Propagate Back Backwards*/
          propagate_half_backwards_kinetic_chemistry_real_multi_unrestricted(kinetic_backwards_half_fermions_real_up,kinetic_backwards_half_fermions_real_down,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

        }/*Energy Calculations***************************************************************************************************/

        /*Orthogonalize Q and R if Necessary*/
        if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
          orthogonalize_chemistry_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,weights,R_up,R_down,ist);
        }

        /*Apply Population Control If Necessary*/
        if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
            ist.n_walkers = population_control_chemistry_multi(weights,wf_up,wf_down,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,new_walker_weights,new_wf_up,new_wf_down,new_overlap_inverse_up,new_overlap_inverse_down,new_overlap_total,new_det_overlap_up,new_det_overlap_down,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns);
        }

     } /*If No Back Propagation*/
     else { /*If Back Propagation*/

        /*First Propagate a Full Step Forwards*/
        propagate_forwards_kinetic_chemistry_real_multi_unrestricted(kinetic_full_fermions_real_up,kinetic_full_fermions_real_down,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,phases,ist);

        /*Then Propagate Forwards a Full U Step******************************************************************************************/
        /*Use Discrete or Continuous Transform*/
        if ( ist.flag_real_eigenvectors == 0 ) { /*Determine If the Eigenctors Are Real or Not*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {
              propagate_forwards_potential_continuous_bp_chemistry_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            } /*Mean Field Flag*/
            else { /*If Mean Field*/
              propagate_forwards_potential_continuous_bp_chemistry_meanfield_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
           }
          } /*If Phaseless*/
          else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
                propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            }
            else {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_multi_unrestricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            }
          }   /*If Mean Not Field*/
          else { /*If Mean Field*/
           if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_meanfield_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
}
           else {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_meanfield_multi_unrestricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
           }
          } /*If Mean Field*/
         } /*Phaseless Else*/
        }
        else {   /*If Eigenvectors Are Real*/
           if ( ist.flag_phaseless == 0 ) {
            if ( ist.flag_meanfield == 0 ) {
                propagate_forwards_potential_continuous_bp_chemistry_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             } /*If Not Mean Field*/
             else { /*If Mean Field*/
                 propagate_forwards_potential_continuous_bp_chemistry_meanfield_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
           }
           else {
            if ( ist.flag_meanfield == 0 ) {
             if ( ist.flag_local_energy == 0 ) {
                 propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
            else {
              propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_real_multi_unrestricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
            }
           }
           else { /*If Mean Field*/
             if ( ist.flag_local_energy == 0 ) {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_squaredshift_meanfield_real_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
             else {
               propagate_forwards_potential_continuous_bp_chemistry_phaseless_localenergy_meanfield_real_multi_unrestricted(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,weights,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,potential_eigvecs_fermions,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_updown,mean_field_first_transform_constant,mean_field_second_transform_constant,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns,&idum);
             }
           }
          } /*Phaseless Else*/
        } /*End Eigenvectors Real********************************************************************************************************/

         /*Store WF, If Needed for Back Prop*/
         store_prev_wf(wf_up,wf_down,prev_wf_up,prev_wf_down,kinetic_forwards_half_fermions_real_up,ist);

         /*Normalize the Walker Weights by the Trial Energy*/
         normalize(weights,ist,cns);

        /*Now Collect Energies and Other Observables*/
        if ( steps%ist.n_steps_energy==0 ) {
          /*First Propagate One Half Step Forward*/

          propagate_half_forwards_kinetic_chemistry_real_multi_unrestricted(kinetic_forwards_half_fermions_real_up,kinetic_forwards_half_fermions_real_down,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);


          /*If You Want to Calculate 2 RDMS*/
          if ( ist.flag_compute_2RDM == 0 ) {

            /*If You Want to Use Back Propagated Wave Function*/
            if (  steps < ist.n_steps_back_propagation ) {
              /*Compute Energ*/
              /*compute_energy_chemistry_mixed_multi_unrestricted(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_up,overlap_down,det_overlap_total,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_original_fermions_up,kinetic_matrix_original_fermions_down,potential_matrix_original_fermions_up,potential_matrix_original_fermions_down,potential_matrix_original_fermions_up_down,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);*/

              /*compute_energy_chemistry_restricted_exact(wf_up,wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_exact,&ke_exact_error,&pe_exact,&pe_exact_error,&te_exact,&te_exact_error,&tw_exact,average_wave_function_up,average_wave_function_down);*/

             cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));

             /*Print Energies*/
             fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, cns.total_energy_shift, te_mixed.real+cns.total_energy_shift, te_mixed_error.real, te_exact.real, te_exact_error.real, tw_mixed.real, tw_mixed.imag); fflush(pf);
             fprintf(pf9, "%f\t\t %f\t\t %f\n", steps*cns.dtau, te_mixed.real+cns.total_energy_shift, te_exact.real); fflush(pf9); 

            }
            else {

              /*Back Propagate First*/
              /*back_propagation_chemistry_restricted(trial_wf_up_phaseless,trial_wf_down_phaseless,bp_wf_up,bp_wf_down,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,kinetic_backwards_half_fermions_real_up,kinetic_full_fermions_real_up,kinetic_forwards_half_fermions_real_up,ist,cns);*/

              /*Compute Energy Second*/
              /*compute_energy_chemistry_restricted_bp(wf_up,wf_down,bp_wf_up,bp_wf_down,prev_wf_up,prev_wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_bp,&ke_bp_error,&pe_bp,&pe_bp_error,&te_bp,&te_bp_error,&tw_bp,average_wave_function_up,average_wave_function_down);*/

              compute_energy_chemistry_mixed_multi_unrestricted_sparse_equilibration(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_energy,trial_determinant_coefficients_energy,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,energies,energies2,energies3,overlaps,density_matrices_up,density_matrices_down,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_mixed,&ke_mixed_error,&pe_mixed,&pe_mixed_error,&te_mixed,&te_mixed_error,&tw_mixed,average_wave_function_up,average_wave_function_down);

              /*compute_energy_chemistry_restricted_exact(wf_up,wf_down,weights,kinetic_matrix_original_fermions_up,potential_matrix_original_fermions_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_exact,&ke_exact_error,&pe_exact,&pe_exact_error,&te_exact,&te_exact_error,&tw_exact,average_wave_function_up,average_wave_function_down); 
 * cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));*/

             /*Print Energies*/
             fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, cns.total_energy_shift, te_mixed.real+cns.total_energy_shift, te_mixed_error.real, te_exact.real, te_exact.imag, tw_mixed.real, tw_mixed.imag); fflush(pf);
             fprintf(pf9, "%f\t\t %f\t\t %f\n", steps*cns.dtau, te_mixed.real+cns.total_energy_shift, te_exact.real); fflush(pf9);

           }
          }
          else {   /*2RDMS*******************************************************************************************************************************************************/

           /*If You Want to Back Propagate Or Not*/
           if ( steps < ist.n_steps_back_propagation ) {

              /*Compute Energy and Density Matrices*/
              compute_energy_chemistry_mixed_two_body_density_matrix_multi_unrestricted_sparse(wf_up,wf_down,trial_wf_up_energy,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,kinetic_matrix_sparse_up,kinetic_matrix_sparse_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,kinetic_ij_sparse_up,kinetic_ij_sparse_down,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_bp,&ke_bp_error,&pe_bp,&pe_bp_error,&te_bp,&te_bp_error,&tw_bp,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down);

             /*Print RDMs******************************/
             fprintf(pf2, "%d\n", steps);
             fprintf(pf3, "%d\n", steps);
             for (i=0; i<ist.n_spatial_orbitals; i++) {
               for (j=0; j<ist.n_spatial_orbitals; j++) {
                if ( Cabs(average_density_matrix_up[i*ist.n_spatial_orbitals+j]) > .0001 || Cabs(average_density_matrix_down[i*ist.n_spatial_orbitals+j]) > .0001 ) {

                  fprintf(pf2, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_up[i*ist.n_spatial_orbitals+j].real, average_density_matrix_up[i*ist.n_spatial_orbitals+j].imag);
                  fprintf(pf3, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_down[i*ist.n_spatial_orbitals+j].real, average_density_matrix_down[i*ist.n_spatial_orbitals+j].imag);
                }
               }
              }
              fprintf(pf2, "\n\n"); fflush(pf2);
              fprintf(pf3, "\n\n"); fflush(pf3);

              /*Print 2-Body RDM*/
              fprintf(pf4, "%d\n", steps);
              for (i=0; i<ist.n_spin_orbitals; i++) {
               for (j=0; j<ist.n_spin_orbitals; j++) {
                for (k=0; k<ist.n_spin_orbitals; k++) {
                 for (l=0; l<ist.n_spin_orbitals; l++) {
                   if ( Cabs(average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]) > .0001 ) {
                      fprintf(pf4, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].imag);
                   }
                 }
                }
               }
              }
fprintf(pf4, "\n\n"); fflush(pf4);

              cns.exp_trial_energy = exp(cns.dtau * (ke_mixed.real + pe_mixed.real));

              /*Print Energies*/
              fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_mixed.real, pe_mixed_error.real, ke_mixed.real, ke_mixed_error.real, cns.total_energy_shift, te_mixed.real+cns.total_energy_shift, te_mixed_error.real, tw_mixed.real, tw_mixed.imag); fflush(pf);
              fprintf(pf9, "%f\t\t %f\n", steps*cns.dtau, te_mixed.real+cns.total_energy_shift); fflush(pf9);

          }
          else {

              /*Back Propagate First*/
              back_propagation_chemistry_multi(trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,bp_wf_up,bp_wf_down,potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,kinetic_backwards_half_fermions_real_up,kinetic_full_fermions_real_up,kinetic_forwards_half_fermions_real_up,ist,cns);

              /*Compute Energy and Density Matrices*/
              compute_energy_chemistry_bp_two_body_density_matrix_restricted_sparse(wf_up,wf_down,bp_wf_up,bp_wf_down,prev_wf_up,prev_wf_down,weights,kinetic_matrix_sparse_up,potential_matrix_sparse_up,kinetic_ij_sparse_up,potential_ij_sparse_up,energy_shifts,total_potential_energy,total_kinetic_energy,total_energy_2,total_weight,weight_rescaling,ist,cns,steps,&ke_bp,&ke_bp_error,&pe_bp,&pe_bp_error,&te_bp,&te_bp_error,&tw_bp,average_density_matrix_up,average_density_matrix_down,average_two_body_density_matrix,average_wave_function_up,average_wave_function_down);

             /*Print RDMs******************************/
             fprintf(pf2, "%d\n", steps);
             fprintf(pf3, "%d\n", steps);
             for (i=0; i<ist.n_spatial_orbitals; i++) {
              for (j=0; j<ist.n_spatial_orbitals; j++) {
               if ( Cabs(average_density_matrix_up[i*ist.n_spatial_orbitals+j]) > .0001 || Cabs(average_density_matrix_down[i*ist.n_spatial_orbitals+j]) > .0001 ) {
                     fprintf(pf2, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_up[i*ist.n_spatial_orbitals+j].real, average_density_matrix_up[i*ist.n_spatial_orbitals+j].imag);
                     fprintf(pf3, "%d\t %d\t %f+%fi\n", i, j, average_density_matrix_down[i*ist.n_spatial_orbitals+j].real, average_density_matrix_down[i*ist.n_spatial_orbitals+j].imag);
                }
               }
              }
              fprintf(pf2, "\n\n"); fflush(pf2);
              fprintf(pf3, "\n\n"); fflush(pf3);

              /*Print 2-Body RDM*/
              fprintf(pf4, "%d\n", steps);
              for (i=0; i<ist.n_spin_orbitals; i++) {
               for (j=0; j<ist.n_spin_orbitals; j++) {
                for (k=0; k<ist.n_spin_orbitals; k++) {
                 for (l=0; l<ist.n_spin_orbitals; l++) {
                  if ( Cabs(average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]) > .0001 ) {
                   fprintf(pf4, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real, average_two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].imag);
                  }
                 }
                }
               }
              }
              fprintf(pf4, "\n\n"); fflush(pf4);

             cns.exp_trial_energy = exp(cns.dtau * (ke_bp.real + pe_bp.real));

            /*Print Energies*/
            fprintf(pf, "%d\t\t  %d\t\t  %f\t\t  %f  %f\t\t  %f %f\t\t  %f\t\t %f %f\t\t %f+%fi\n", steps, ist.n_walkers, steps*cns.dtau, pe_bp.real, pe_bp_error.real, ke_bp.real, ke_bp_error.real, energy_shifts[0].real, te_bp.real, te_bp_error.real, tw_bp.real, tw_bp.imag); fflush(pf);
            fprintf(pf9, "%f\t\t %f\n", steps*cns.dtau, te_bp.real); fflush(pf9);

           } /*Else BP and two RDMS*/
          }

          /*Print Wave Functions*/
          fprintf(pf5, "%d\n", steps);
          fprintf(pf6, "%d\n", steps);
          for (i=0; i<ist.n_spatial_orbitals; i++) {
            for (j=0; j<ist.n_up; j++) {
             fprintf(pf5, "%f+%fi\t", average_wave_function_up[i*ist.n_up+j].real, average_wave_function_up[i*ist.n_up+j].imag);
           }
           for (j=0; j<ist.n_down; j++) {
             fprintf(pf6, "%f+%fi\t", average_wave_function_down[i*ist.n_down+j].real, average_wave_function_down[i*ist.n_down+j].imag);
           }
           fprintf(pf5, "\n");
           fprintf(pf6, "\n");
          }
          fprintf(pf5, "\n\n"); fflush(pf5);
          fprintf(pf6, "\n\n"); fflush(pf6);

/*Propagate Back Backwards*/
           propagate_half_backwards_kinetic_chemistry_real_multi_unrestricted(kinetic_backwards_half_fermions_real_up,kinetic_backwards_half_fermions_real_down,wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,det_overlap_up,det_overlap_down,overlap_inverse_up,overlap_inverse_down,weights,ist);

        }/*Energy Calculations***************************************************************************************************/

        /*Copy Wave Functions Back If Using Back Propagation*/
        copy_back_propagation_matrices_forward(potential_matrices_back_propagation_up,potential_matrices_back_propagation_down,ist,cns);
        copy_prev_wf_forward(prev_wf_up,prev_wf_down,ist,cns);

        /*Orthogonalize Q and R if Necessary*/
        if ( steps%ist.n_steps_orthogonalize == 0 && steps != 0 ) {
          orthogonalize_chemistry_multi(wf_up,wf_down,trial_wf_up_phaseless,trial_wf_down_phaseless,trial_determinant_coefficients_phaseless,overlap_total,overlap_inverse_up,overlap_inverse_down,det_overlap_up,det_overlap_down,weights,R_up,R_down,ist);
        }

        /*Apply Population Control If Necessary*/
         if ( steps%ist.n_steps_population == 0 && steps != 0 ) {
           //ist.n_walkers = population_control_chemistry_multi(weights,wf_up,wf_down,overlap_inverse_up,overlap_inverse_down,overlap_total,det_overlap_up,det_overlap_down,new_walker_weights,new_wf_up,new_wf_down,new_overlap_inverse_up,new_overlap_inverse_down,new_overlap_total,new_det_overlap_up,new_det_overlap_down,walkers_big,walkers_small,walkers_copied,weight_rescaling,ist,cns);
         }

     } /*If Back Propagation */
    } /*Steps*/

    /******************************************************************************************/

    free(weights);
    free(new_walker_weights);

    free(total_potential_energy); free(total_kinetic_energy); free(total_energy_2); free(total_weight);

    free(energy_shifts);
    free(kinetic_matrix_original_fermions_up);
    free(kinetic_matrix_original_fermions_down); 

    free(kinetic_matrix_sparse_up); 
    free(kinetic_matrix_sparse_down); 
    free(kinetic_ij_sparse_up); 
    free(kinetic_ij_sparse_down); 

    free(potential_matrix_original_fermions_up);
    free(potential_matrix_original_fermions_down); 
    free(potential_matrix_original_fermions_updown);  

    free(potential_matrix_sparse_up); 
    free(potential_matrix_sparse_down); 
    free(potential_matrix_sparse_updown); 
    free(potential_ij_sparse_up); 
    free(potential_ij_sparse_down); 
    free(potential_ij_sparse_updown); 

    free(potential_matrix_original_fermions_spin);
    free(potential_supermatrix_fermions);
    free(potential_eigs_fermions); free(potential_eigvecs_fermions);
    free(potential_eigs_fermions_first_transform); free(potential_eigs_fermions_second_transform);
    free(potential_one_body_matrix_up); free(potential_one_body_matrix_down); 

    free(list_orbitals_trial_energy); 

    free(mean_field_density_first_transform);  free(mean_field_density_second_transform);
    free(mean_field_first_transform_constant); free(mean_field_second_transform_constant);

    free(overlap_total);
    free(det_overlap_up); free(det_overlap_down);
    free(overlap_inverse_up); free(overlap_inverse_down);
    free(wf_up); free(wf_down);

    free(trial_determinant_coefficients_energy); free(trial_determinant_coefficients_phaseless);

    free(trial_wf_up_energy); free(trial_wf_down_energy);
    free(trial_wf_up_phaseless); free(trial_wf_down_phaseless);

    free(trial_density_up); free(trial_density_down);

    free(kinetic_full_fermions_real_up); free(kinetic_backwards_half_fermions_real_up); free(kinetic_forwards_half_fermions_real_up);
    free(kinetic_eigs_fermions_real_up); free(kinetic_eigvecs_fermions_real_up);

    free(kinetic_full_fermions_real_down); free(kinetic_backwards_half_fermions_real_down); free(kinetic_forwards_half_fermions_real_down);
    free(kinetic_eigs_fermions_real_down); free(kinetic_eigvecs_fermions_real_down);

    free(new_overlap_total);
    free(new_det_overlap_up); free(new_det_overlap_down);
    free(new_overlap_inverse_up); free(new_overlap_inverse_down);
    free(new_wf_up); free(new_wf_down);

    free(average_density_matrix_up); free(average_density_matrix_down);
    free(average_two_body_density_matrix);

    free(exact_two_body_density_matrix);

    free(average_wave_function_up); free(average_wave_function_down);

    free(R_up); free(R_down);

    if ( ist.flag_back_propagation == 1 ) {
       free(bp_wf_up); free(bp_wf_down);
       free(prev_wf_up); free(prev_wf_down);
       free(potential_matrices_back_propagation_up);
       free(potential_matrices_back_propagation_down);
    }

   fclose(pf); fclose(pf9);
   }
  fclose(pout);
  fclose(pf); 
  fclose(pf2); 
  fclose(pf3); 
  fclose(pf4); 

  /*******************************************************************************************/

return 0; 
}

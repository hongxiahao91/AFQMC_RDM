#include "afqmc.h"

/*Functions for Determining the Trial Energy of a Multi Configuration*********************************************************/
double get_trial_energy_density_multi_restricted(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *trial_determinant_coefficients,MKL_Complex16 *one_body_density_matrix_up,MKL_Complex16 *one_body_density_matrix_down,MKL_Complex16 *two_body_density_matrix,MKL_Complex16 *kinetic_original_matrix,MKL_Complex16 *potential_original_matrix,MKL_Complex16 *potential_original_matrix_spin,MKL_Complex16 *energy_shifts,int_st ist) {

   /*Obtains the Trial Energy Based Upon the Trial Wavefunction for Use in Estimating Walker Weights*/
   MKL_Complex16 coefficient;
   MKL_Complex16 trial_ke, trial_ke_total, trial_pe, trial_pe_total;
   double trial_ke_density_matrix, trial_pe_density_matrix; 
   MKL_Complex16 trial_ke_density_matrix_old, trial_pe_density_matrix_old; 
   int count, sign_reorder = 1;
   int check_sign;  
   int *config_1, *config_2, *original_config_1; 
   double sum_coefficients;
   double trial_energy_total;
   int i, j, k, l;
   FILE *pf = fopen("errors.dat", "a+"); 
   FILE *pf1 = fopen("exact_one_body_density_matrix_up.dat", "a+"); 
   FILE *pf2 = fopen("exact_one_body_density_matrix_down.dat", "a+");  
   FILE *pf3 = fopen("exact_two_body_density_matrix.dat", "a+");  

   /*Density Matrices Down and Up*/
   config_1 = (int *)calloc(ist.n_electrons,sizeof(int)); 
   config_2 = (int *)calloc(ist.n_electrons,sizeof(int)); 
   original_config_1 = (int *)calloc(ist.n_electrons,sizeof(int)); 

   /*Zero Objects*/
   czero_vec(one_body_density_matrix_up,ist.n_spatial_orbitals_sq); 
   czero_vec(one_body_density_matrix_down,ist.n_spatial_orbitals_sq);
   czero_vec(two_body_density_matrix,ist.n_spin_orbitals_fourth);  


   trial_ke_total.real = trial_ke_total.imag = 0.0;
   trial_pe_total.real = trial_pe_total.imag = 0.0;
   trial_energy_total = 0.0;
   sum_coefficients = 0.0;
   for (i=0; i<ist.n_determinants_trial_energy; i++) {

     /*Determine the Orbital Configuration*/
     determine_configuration(original_config_1,&trial_wf_up[i*ist.n_spatial_orbitals_n_up],&trial_wf_down[i*ist.n_spatial_orbitals_n_down],ist); 

 
     for (j=0; j<ist.n_determinants_trial_energy; j++) {

        /*Determine Orbital Configuration for J*/
        determine_configuration(config_2,&trial_wf_up[j*ist.n_spatial_orbitals_n_up],&trial_wf_down[j*ist.n_spatial_orbitals_n_down],ist);


        for (k=0; k<ist.n_electrons; k++) {
           config_1[k] = original_config_1[k]; 
        } 

        /*Get How Many Orbitals Are Different*/
        count = diff_orbitals(config_1,config_2,ist.n_electrons);

        /*Only Consider Energies Of Configurations With Less Than 2 Differences*/
       if ( count < 3 ) {

          sign_reorder = 1; 
          if (count != 0 ) {
            sign_reorder = reorder_configurations(config_1,config_2,ist.n_electrons);
          }

          /*Get Coefficient of Contributions from Given Determinant*/
          coefficient = Cmul(conjugate(trial_determinant_coefficients[i]), trial_determinant_coefficients[j]);
        //fprintf(pf, "trial coefficients %d %f\n", i, trial_determinant_coefficients[i].real); fflush(pf); 

          /*Only Diagonal Coefficient Contribute to Overlap*/
          if ( i== j ) {
            sum_coefficients += coefficient.real; 
          }
          coefficient = RCmul(sign_reorder, coefficient);

          /*Now Use Configuration To Determine One and Two Body Elements*/
          trial_ke.real = get_fci_matrix_element_kinetic_restricted(config_1,config_2,one_body_density_matrix_up,one_body_density_matrix_down,kinetic_original_matrix,coefficient,ist,count); 
          trial_pe.real = get_fci_matrix_element_potential_restricted(config_1,config_2,two_body_density_matrix,potential_original_matrix,coefficient,ist,count); 
         //fprintf(pf, "trial ke %f\n", trial_ke.real); fflush(pf); 

          trial_ke_total = Cadd(trial_ke_total, trial_ke);
          trial_pe_total = Cadd(trial_pe_total, trial_pe);

        } /*Checking Count Less Than Three*/

      }
    }
    sum_coefficients = 1.0/sum_coefficients;

    /*Normalize Energies*/
    trial_pe_total = RCmul(sum_coefficients, trial_pe_total);
    trial_ke_total = RCmul(sum_coefficients, trial_ke_total);
    trial_energy_total = trial_pe_total.real + trial_ke_total.real;
    
    /*Normalize One-Body DMs*/
    /*Should Use LAPACK*/
    trial_ke_density_matrix = trial_pe_density_matrix = 0.0; 
    for (i=0; i<ist.n_spatial_orbitals; i++) {
     for (j=0; j<ist.n_spatial_orbitals; j++) {
       one_body_density_matrix_up[i*ist.n_spatial_orbitals+j] = RCmul(sum_coefficients,one_body_density_matrix_up[i*ist.n_spatial_orbitals+j]); 
       one_body_density_matrix_down[i*ist.n_spatial_orbitals+j] = RCmul(sum_coefficients,one_body_density_matrix_down[i*ist.n_spatial_orbitals+j]); 

       trial_ke_density_matrix += kinetic_original_matrix[i*ist.n_spatial_orbitals+j].real * (one_body_density_matrix_up[i*ist.n_spatial_orbitals+j].real + one_body_density_matrix_down[i*ist.n_spatial_orbitals+j].real); 

       if ( one_body_density_matrix_up[i*ist.n_spatial_orbitals+j].real != 0 ) { 
          fprintf(pf1, "%d %d %f\n", i, j, one_body_density_matrix_up[i*ist.n_spatial_orbitals+j].real); 
       }
       if ( one_body_density_matrix_down[i*ist.n_spatial_orbitals+j].real != 0 ) {
         fprintf(pf2, "%d %d %f\n", i, j, one_body_density_matrix_down[i*ist.n_spatial_orbitals+j].real);  
       }

     } 
    } 

    /*Normalize Two-Body DMs*/ 
    for (i=0; i<ist.n_spin_orbitals; i++) {
     for (j=0; j<ist.n_spin_orbitals; j++) {
      for (k=0; k<ist.n_spin_orbitals; k++) {
       for (l=0; l<ist.n_spin_orbitals; l++) {

          /*REPLACE WITH LAPACKE Routine!!!*/
          two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l] = RCmul(sum_coefficients,two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]);
  
          if (  two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real *  potential_original_matrix_spin[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real != 0 ) {
            fprintf(pf3, "%d %d %d %d %f\n", i, j, k, l, two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real); fflush(pf3); 
          }

          /*Minus To Get Proper Sign for cicjclck when matrix is cicjckcl*/
          trial_pe_density_matrix -= two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real * potential_original_matrix_spin[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real;

       }
      }
     }
    }
    trial_pe_density_matrix *= .5; 
    fflush(pf3); 

    trial_ke_density_matrix_old = compute_kinetic_energy_density_matrix_restricted(one_body_density_matrix_up,one_body_density_matrix_down,kinetic_original_matrix,ist);
    trial_pe_density_matrix_old = compute_potential_energy_density_matrix_restricted(one_body_density_matrix_up,one_body_density_matrix_down,potential_original_matrix,ist);


    /*Compute The Kinetic Energy Based Off the One-Body Density Matrix Instdead*/
    fprintf(pf, "trial ke %f %f %f\n", trial_ke_total.real, trial_ke_density_matrix, trial_ke_density_matrix_old.real); fflush(pf);
    fprintf(pf, "trial pe %f %f %f\n", trial_pe_total.real, trial_pe_density_matrix, trial_pe_density_matrix_old.real); fflush(pf); 


free(config_1); 
free(config_2);    
fclose(pf1); 
fclose(pf2); 
fclose(pf3); 
return(trial_energy_total);
}

/******************************************************************************************************************************************/

double get_trial_energy_density_multi_unrestricted(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *trial_determinant_coefficients,MKL_Complex16 *one_body_density_matrix_up,MKL_Complex16 *one_body_density_matrix_down,MKL_Complex16 *two_body_density_matrix,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *kinetic_original_matrix_down,MKL_Complex16 *potential_original_matrix_up,MKL_Complex16 *potential_original_matrix_down,MKL_Complex16 *potential_original_matrix_updown, MKL_Complex16 *potential_original_matrix_spin,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *potential_ij_sparse_up,int *potential_ij_sparse_down,int *potential_ij_sparse_updown,MKL_Complex16 *energy_shifts,int_st ist) {
 
   /*Obtains the Trial Energy Based Upon the Trial Wavefunction for Use in Estimating Walker Weights*/
   MKL_Complex16 coefficient;
   MKL_Complex16 trial_ke, trial_ke_total, trial_pe, trial_pe_total;
   double trial_ke_density_matrix, trial_pe_density_matrix;
   MKL_Complex16 trial_ke_density_matrix_old, trial_pe_density_matrix_old;
   int count, sign_reorder = 1;
   int check_sign;
   int *config_1, *config_2, *original_config_1;
   double sum_coefficients;
   double trial_energy_total;
   int i, j, k, l;
   FILE *pf = fopen("errors.dat", "a+");
   FILE *pf1 = fopen("exact_one_body_density_matrix_up.dat", "a+");
   FILE *pf2 = fopen("exact_one_body_density_matrix_down.dat", "a+");
   FILE *pf3 = fopen("exact_two_body_density_matrix.dat", "a+");

   /*Density Matrices Down and Up*/
   config_1 = (int *)calloc(ist.n_electrons,sizeof(int));
   config_2 = (int *)calloc(ist.n_electrons,sizeof(int));
   original_config_1 = (int *)calloc(ist.n_electrons,sizeof(int));

   /*Zero Objects*/
   czero_vec(one_body_density_matrix_up,ist.n_spatial_orbitals_sq);
   czero_vec(one_body_density_matrix_down,ist.n_spatial_orbitals_sq);
   czero_vec(two_body_density_matrix,ist.n_spin_orbitals_fourth);

   trial_ke_total.real = trial_ke_total.imag = 0.0;
   trial_pe_total.real = trial_pe_total.imag = 0.0;
   trial_energy_total = 0.0;
   sum_coefficients = 0.0;
   for (i=0; i<ist.n_determinants_trial_energy; i++) {

     /*Determine the Orbital Configuration*/
     determine_configuration(original_config_1,&trial_wf_up[i*ist.n_spatial_orbitals_n_up],&trial_wf_down[i*ist.n_spatial_orbitals_n_down],ist);
     for (j=0; j<ist.n_determinants_trial_energy; j++) {

        /*Determine Orbital Configuration for J*/
        determine_configuration(config_2,&trial_wf_up[j*ist.n_spatial_orbitals_n_up],&trial_wf_down[j*ist.n_spatial_orbitals_n_down],ist);

        for (k=0; k<ist.n_electrons; k++) {
           config_1[k] = original_config_1[k];
        }

        /*Get How Many Orbitals Are Different*/
        count = diff_orbitals(config_1,config_2,ist.n_electrons);

        /*Only Consider Energies Of Configurations With Less Than 2 Differences*/
       if ( count < 3 ) {

          sign_reorder = 1;
          if (count != 0 ) {
            sign_reorder = reorder_configurations(config_1,config_2,ist.n_electrons);
          }

          /*Get Coefficient of Contributions from Given Determinant*/
          coefficient = Cmul(conjugate(trial_determinant_coefficients[i]), trial_determinant_coefficients[j]);
          /*fprintf(pf, "trial coefficients %d %f\n", i, trial_determinant_coefficients[i].real); fflush(pf);*/

          /*Only Diagonal Coefficient Contribute to Overlap*/
          if ( i== j ) {
            sum_coefficients += coefficient.real;
          }
          coefficient = RCmul(sign_reorder, coefficient);

          /*Now Use Configuration To Determine One and Two Body Elements*/
          trial_ke.real = get_fci_matrix_element_kinetic_unrestricted(config_1,config_2,one_body_density_matrix_up,one_body_density_matrix_down,kinetic_original_matrix_up,kinetic_original_matrix_down,coefficient,ist,count);
          trial_pe.real = get_fci_matrix_element_potential_unrestricted(config_1,config_2,two_body_density_matrix,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,coefficient,ist,count);

          trial_ke_total = Cadd(trial_ke_total, trial_ke);
          trial_pe_total = Cadd(trial_pe_total, trial_pe);

        } /*Checking Count Less Than Three*/

      }
    }
    sum_coefficients = 1.0/sum_coefficients;

    /*Normalize Energies*/
    trial_pe_total = RCmul(sum_coefficients, trial_pe_total);
    trial_ke_total = RCmul(sum_coefficients, trial_ke_total);
    trial_energy_total = trial_pe_total.real + trial_ke_total.real;

    /*Normalize One-Body DMs*/
    /*Should Use LAPACK*/
    trial_ke_density_matrix = trial_pe_density_matrix = 0.0;
    for (i=0; i<ist.n_spatial_orbitals; i++) {
     for (j=0; j<ist.n_spatial_orbitals; j++) {
       one_body_density_matrix_up[i*ist.n_spatial_orbitals+j] = RCmul(sum_coefficients,one_body_density_matrix_up[i*ist.n_spatial_orbitals+j]);
       one_body_density_matrix_down[i*ist.n_spatial_orbitals+j] = RCmul(sum_coefficients,one_body_density_matrix_down[i*ist.n_spatial_orbitals+j]);
       
       trial_ke_density_matrix += kinetic_original_matrix_up[i*ist.n_spatial_orbitals+j].real * one_body_density_matrix_up[i*ist.n_spatial_orbitals+j].real + kinetic_original_matrix_down[i*ist.n_spatial_orbitals+j].real * one_body_density_matrix_down[i*ist.n_spatial_orbitals+j].real;

       if ( one_body_density_matrix_up[i*ist.n_spatial_orbitals+j].real != 0 ) {
          fprintf(pf1, "%d %d %f\n", i, j, one_body_density_matrix_up[i*ist.n_spatial_orbitals+j].real);
       }
       if ( one_body_density_matrix_down[i*ist.n_spatial_orbitals+j].real != 0 ) {
         fprintf(pf2, "%d %d %f\n", i, j, one_body_density_matrix_down[i*ist.n_spatial_orbitals+j].real);
       }

     }
    }

    /*Normalize Two-Body DMs*/
    for (i=0; i<ist.n_spin_orbitals; i++) {
     for (j=0; j<ist.n_spin_orbitals; j++) {
      for (k=0; k<ist.n_spin_orbitals; k++) {
       for (l=0; l<ist.n_spin_orbitals; l++) {

          /*REPLACE WITH LAPACKE Routine!!!*/
          two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l] = RCmul(sum_coefficients,two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]);

          if (  two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real *  potential_original_matrix_spin[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real != 0 ) {
            fprintf(pf3, "%d %d %d %d %f\n", i, j, k, l, two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real); fflush(pf3);
          }

          /*Minus To Get Proper Sign for cicjclck when matrix is cicjckcl*/
         trial_pe_density_matrix -= two_body_density_matrix[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real * potential_original_matrix_spin[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l].real;

       }
      }
     }
    }
    trial_pe_density_matrix *= .5;
    fflush(pf3);

   
    trial_ke_density_matrix_old = compute_kinetic_energy_density_matrix_unrestricted(one_body_density_matrix_up,one_body_density_matrix_down,kinetic_original_matrix_up,kinetic_original_matrix_down,ist);
    trial_pe_density_matrix_old = compute_potential_energy_density_matrix_unrestricted(one_body_density_matrix_up,one_body_density_matrix_down,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,ist);

   fprintf(pf, "original pe %f\n", trial_pe_density_matrix_old.real); fflush(pf); 
    trial_pe_density_matrix_old = compute_potential_energy_density_matrix_unrestricted_sparse(one_body_density_matrix_up,one_body_density_matrix_down,potential_matrix_sparse_up,potential_matrix_sparse_down,potential_matrix_sparse_updown,potential_ij_sparse_up,potential_ij_sparse_down,potential_ij_sparse_updown,ist); 
   fprintf(pf, "sparse pe %f\n", trial_pe_density_matrix_old.real); fflush(pf); 


    /*Compute The Kinetic Energy Based Off the One-Body Density Matrix Instdead*/
    fprintf(pf, "trial ke %f %f %f\n", trial_ke_total.real, trial_ke_density_matrix, trial_ke_density_matrix_old.real); fflush(pf);
    fprintf(pf, "trial pe %f %f %f\n", trial_pe_total.real, trial_pe_density_matrix, trial_pe_density_matrix_old.real); fflush(pf);
    fprintf(pf, "sum %f %f %f\n", trial_ke_total.real+trial_pe_total.real+energy_shifts[0].real+energy_shifts[1].real+energy_shifts[2].real+energy_shifts[3].real, trial_ke_density_matrix+trial_pe_density_matrix+energy_shifts[0].real+energy_shifts[1].real+energy_shifts[2].real+energy_shifts[3].real, trial_ke_density_matrix_old.real+trial_pe_density_matrix_old.real+energy_shifts[0].real+energy_shifts[1].real+energy_shifts[2].real+energy_shifts[3].real); fflush(pf); 


free(config_1);
free(config_2);
fclose(pf1);
fclose(pf2);
fclose(pf3);
return(trial_energy_total);
}

/************************************************************************************************************************************************/

void determine_configuration(int *total_config,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,int_st ist) {

     /*Determine the Spatial Orbitals Given a Trial WF*/
     int i, j, count = 0;

     /*Get Spin Up Configuration*/
     for (i=0; i<ist.n_up; i++) {
      for (j=0; j<ist.n_spatial_orbitals; j++) {

         if ( trial_wf_up[j*ist.n_up+i].real == 1.0 ) {
            total_config[count] = j*2;
            count++;
         }

      }
     }


     /*Get Spin Down Configuration*/
     for (i=0; i<ist.n_down; i++) {
       for (j=0; j<ist.n_spatial_orbitals; j++) {

         if ( trial_wf_down[j*ist.n_down+i].real == 1.0 ) {
             total_config[count] = j*2+1;
             count++;
         }

       }
     }

return;
}

/******************************************************************************************************************************************/

double get_fci_matrix_element_kinetic_restricted(int *config_vector_1,int *config_vector_2,MKL_Complex16 *trial_density_up,MKL_Complex16 *trial_density_down,MKL_Complex16 *kinetic_matrix_original_fermions,MKL_Complex16 coefficient,int_st ist,int count) {

   /*Determines the Matrix Elements from the Sum of the Electron Integrals That Remain from the Overlap of the Two WFs*/
   int i;
   int i_orbital, i_spin; 
   int change_in_sign = 1;
   int first_orbital_1, first_orbital_2;
   int first_spin_1, first_spin_2; 
   double matrix_element = 0.0;  

   /*Determine Matrix Element Based Upon Orbitals Occupied*/ 
   if ( count == 1 ) { /*If Only One Orbitals Differs*/

      /*Find Orbitals*/
      first_orbital_1 = (int)(config_vector_1[ist.n_electrons-1]/2.0);
      first_orbital_2 = (int)(config_vector_2[ist.n_electrons-1]/2.0);

      /*Get KE Matrix Element*/
      matrix_element = kinetic_matrix_original_fermions[first_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real;

      /*Find Spins To Put in Trial Density Matrices*/
      first_spin_1 = config_vector_1[ist.n_electrons-1]%2; 
      first_spin_2 = config_vector_2[ist.n_electrons-1]%2;  

      /*Add One-Body Terms to Up Matrix -- Has To Have Same Spin Because Spin Is Conserved*/
      if ( first_spin_1 == 0 ) {
        trial_density_up[first_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real += coefficient.real;    
      }
      else {  /*Add To Down Matrix*/
        trial_density_down[first_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real += coefficient.real;
      } 


    }else if ( count == 0 ) {  /*If Same Determinant*/

       /*Run Through All Electrons*/
       for (i=0; i<ist.n_electrons; i++) {

         /*Find I Spin and Orbital*/
         i_orbital = (int)(config_vector_1[i]/2.0);
         i_spin = config_vector_1[i]%2; 

         /*Add One Electron Piece*/
         matrix_element += kinetic_matrix_original_fermions[i_orbital*ist.n_spatial_orbitals+i_orbital].real;

         /*Add To Trial Density Matrices*/
         if ( i_spin == 0 ) {
           trial_density_up[i_orbital*ist.n_spatial_orbitals+i_orbital].real += coefficient.real; 
         }
         else {
           trial_density_down[i_orbital*ist.n_spatial_orbitals+i_orbital].real += coefficient.real; 
         }

       }
   }

   /*Multiply Ke By Appropriate Matrix Element*/
   matrix_element *= coefficient.real; 

return(matrix_element);
}

/******************************************************************************************************************************************/

double get_fci_matrix_element_kinetic_unrestricted(int *config_vector_1,int *config_vector_2,MKL_Complex16 *trial_density_up,MKL_Complex16 *trial_density_down,MKL_Complex16 *kinetic_matrix_original_fermions_up,MKL_Complex16 *kinetic_matrix_original_fermions_down,MKL_Complex16 coefficient,int_st ist,int count) {

   /*Determines the Matrix Elements from the Sum of the Electron Integrals That Remain from the Overlap of the Two WFs*/
   int i;
   int i_orbital, i_spin;
   int change_in_sign = 1;
   int first_orbital_1, first_orbital_2;
   int first_spin_1, first_spin_2;
   double matrix_element = 0.0;

   /*Determine Matrix Element Based Upon Orbitals Occupied*/
   if ( count == 1 ) { /*If Only One Orbitals Differs*/

      /*Find Orbitals*/
      first_orbital_1 = (int)(config_vector_1[ist.n_electrons-1]/2.0);
      first_orbital_2 = (int)(config_vector_2[ist.n_electrons-1]/2.0);

      /*Find Spins To Put in Trial Density Matrices*/
      first_spin_1 = config_vector_1[ist.n_electrons-1]%2;
      first_spin_2 = config_vector_2[ist.n_electrons-1]%2;

      /*Add One-Body Terms to Up Matrix -- Has To Have Same Spin Because Spin Is Conserved*/
      if ( first_spin_1 == 0 ) {
        matrix_element = kinetic_matrix_original_fermions_up[first_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real;
        trial_density_up[first_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real += coefficient.real;
      }
      else {  /*Add To Down Matrix*/
        matrix_element = kinetic_matrix_original_fermions_down[first_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real;
        trial_density_down[first_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real += coefficient.real;
      }


    }else if ( count == 0 ) {  /*If Same Determinant*/

       /*Run Through All Electrons*/
       for (i=0; i<ist.n_electrons; i++) {

         /*Find I Spin and Orbital*/
         i_orbital = (int)(config_vector_1[i]/2.0);
         i_spin = config_vector_1[i]%2;

         /*Add To Trial Density Matrices*/
         if ( i_spin == 0 ) {
           trial_density_up[i_orbital*ist.n_spatial_orbitals+i_orbital].real += coefficient.real;
           matrix_element += kinetic_matrix_original_fermions_up[i_orbital*ist.n_spatial_orbitals+i_orbital].real;
         }
         else {
           trial_density_down[i_orbital*ist.n_spatial_orbitals+i_orbital].real += coefficient.real;
           matrix_element += kinetic_matrix_original_fermions_down[i_orbital*ist.n_spatial_orbitals+i_orbital].real;
         }

       }
   }

   /*Multiply Ke By Appropriate Matrix Element*/
   matrix_element *= coefficient.real;

return(matrix_element);
}

/**************************************************************************************************************************/

int diff_orbitals(int *config_vector_1,int *config_vector_2,int size) {
  /*Examine How Many Different Orbitals*/

  int i, j;
  int orbital_1;
  int diff_orbital_flag;
  int count = 0;

  /*Run Through Electrons*/
  for (i=0; i<size; i++) {

     /*Find Given Orbital*/
     orbital_1 = config_vector_1[i];
     diff_orbital_flag = 1;

      /*Determine If Same Orbital Appears in Config 2*/
      for (j=0; j<size; j++) {
         if ( orbital_1 == config_vector_2[j] ) {
           diff_orbital_flag = 0;
         }
       }

      if ( diff_orbital_flag == 1 ) {
        count++;
      }
   }

return(count);
}

/******************************************************************************************************************************************/

int reorder_configurations(int *config_vector_1,int *config_vector_2,int size) {

    /*Reorders the Configurations so That the Same Orbitals Are in the Same Order*/
    /*CAN BE DONE FASTER!! - Slow Stupid Way for Now*/
    int sign = 1;  
    int i, j;
    int flag_same;
    int orbital_1;
    int *temp_config_vector_1, *temp_config_vector_2;
    int count_same = 0, count_diff = 0;

    temp_config_vector_1 = (int *)calloc(size,sizeof(int));
    temp_config_vector_2 = (int *)calloc(size,sizeof(int));

    /*First Reorder Temp Config of First So That All Orbitals That Are Same Are At Beginning*/
    for (i=0; i<size; i++) {

       orbital_1 = config_vector_1[i];
       flag_same = 0;
       for (j=0; j<size; j++) {
         if ( orbital_1 == config_vector_2[j] ) {
           flag_same = 1;
         }
       }

       /*Store In New Order Depending on Flags*/
       if ( flag_same == 0 ) {
         temp_config_vector_1[size-1-count_diff] = orbital_1;
         count_diff++;
       }
       else{
         temp_config_vector_2[count_same] = temp_config_vector_1[count_same] = orbital_1;
         count_same++;
       }

    }

    /*Now Reorder Second Configuration to Reflect First*/
    /*Find the One or Two Differing Orbitals*/
    count_diff = 0;
    for (i=0; i<size; i++) {
      orbital_1 = config_vector_2[i];

      flag_same = 0;
      for (j=0; j<count_same; j++) {
         if ( orbital_1 == temp_config_vector_2[j] ) {
          flag_same = 1;
         }
       }

       if ( flag_same == 0 ) {
         temp_config_vector_2[size-1-count_diff] = orbital_1;
         count_diff++;
         count_same++;
       }
     }

    
     sign = get_sign_configurations(config_vector_1, temp_config_vector_1,size) * get_sign_configurations(config_vector_2, temp_config_vector_2,size); 

     /*Now Restore the Vectors*/
     for (i=0; i<size; i++) {
      config_vector_1[i] = temp_config_vector_1[i];
      config_vector_2[i] = temp_config_vector_2[i];
     }

free(temp_config_vector_1);
free(temp_config_vector_2);
return(sign);
}

/*************************************************************************************************************************************/

int get_sign_configurations(int *config_vector_1,int *config_vector_2,int size) {
 
    /*Determines Sign Difference Between Two Configurations - A LIttle Bit Hacky For Now*/
    /*Improve Routine in the Future*/ 
    int sign = 1; 
    int i, j, k; 
    int hold_temp, hold_final; 
    int *temp_config_vector_1; 
 
    temp_config_vector_1 = (int*)calloc(size,sizeof(int)); 

    /*First Copy Initial Vector Into Temporary One*/
    for (i=0; i<size; i++) {
      temp_config_vector_1[i] = config_vector_1[i]; 
    }

    /*Now See What Sites Differ Between Same and Different And Attempt To Alter*/
    for (i=0; i<size; i++) {

      while (temp_config_vector_1[i] != config_vector_2[i]) {

        hold_temp = temp_config_vector_1[i]; 

        /*Where Is Value in Config 2 Vector*/
        for (j=0; j<size; j++) {  
           if ( config_vector_2[j] == hold_temp ) {
               hold_final = j;
               j = size;  
           }
        }   

        /*Now Move Temp Vector Back and Replace With Correct Site*/
        for (j=0; j<hold_final; j++) {
          if ( hold_final != size ) {
           temp_config_vector_1[j] = temp_config_vector_1[j+1]; 
          }
        }
        temp_config_vector_1[hold_final] = hold_temp;   

        sign *= pow(-1, hold_final-i); 

     } /*Check If Equal or Not*/      


   } /*Running Through all Sites*/

free(temp_config_vector_1); 
return(sign); 
}

/*************************************************************************************************************************************/

double get_fci_matrix_element_potential_restricted(int *config_vector_1,int *config_vector_2,MKL_Complex16 *two_body_density_matrix,MKL_Complex16 *potential_matrix_original_fermions,MKL_Complex16 coefficient,int_st ist,int count) {

   /*Determines the Matrix Elements from the Sum of the Electron Integrals That Remain from the Overlap of the Two WFs*/
   int i, j;
   int j_orbital, j_spin, j_spin_orbital;
   int i_orbital, i_spin, i_spin_orbital; 
   int hold;
   int electrons_minus_one = ist.n_electrons-1;
   int total_spin_1, total_spin_2;
   int first_orbital_1=0, first_orbital_2=0;
   int second_orbital_1=0, second_orbital_2=0;
   int first_spin_1=0, first_spin_2=0;
   int second_spin_1=0, second_spin_2=0;
   int first_spin_orbital_1=0, first_spin_orbital_2=0; 
   int second_spin_orbital_1=0, second_spin_orbital_2=0; 
   double matrix_element = 0.0;
   FILE *pf = fopen("errors.dat", "a+"); 


   /*Run Through Different Orbital Counts*/
   if ( count == 2 ) { /*If Two Different Electrons**********************/


        first_orbital_1 = (int)(config_vector_1[electrons_minus_one]/2.0);
        first_orbital_2 = (int)(config_vector_2[electrons_minus_one]/2.0);
        first_spin_orbital_1 = config_vector_1[electrons_minus_one]; 
        first_spin_orbital_2 = config_vector_2[electrons_minus_one]; 
        first_spin_1 = config_vector_1[electrons_minus_one]%2;
        first_spin_2 = config_vector_2[electrons_minus_one]%2;

        second_orbital_1 = (int)(config_vector_1[ist.n_electrons-2]/2.0);
        second_orbital_2 = (int)(config_vector_2[ist.n_electrons-2]/2.0);
        second_spin_orbital_1 = config_vector_1[ist.n_electrons-2]; 
        second_spin_orbital_2 = config_vector_2[ist.n_electrons-2]; 
        second_spin_1 = config_vector_1[ist.n_electrons-2]%2;
        second_spin_2 = config_vector_2[ist.n_electrons-2]%2;

        /*Two-Body Density Matrices*/
        two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+second_spin_orbital_1*ist.n_spin_orbitals_sq+second_spin_orbital_2*ist.n_spin_orbitals+first_spin_orbital_2].real += coefficient.real; 
        two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+second_spin_orbital_1*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+second_spin_orbital_2].real -= coefficient.real ;
        two_body_density_matrix[second_spin_orbital_1*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+second_spin_orbital_2*ist.n_spin_orbitals+first_spin_orbital_2].real -= coefficient.real; 
        two_body_density_matrix[second_spin_orbital_1*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+second_spin_orbital_2].real += coefficient.real; 

        /*Check If Spins the Same for Exchange Case*/
       if ( first_spin_1 == first_spin_2 && second_spin_1 == second_spin_2 ) {

         /*Check Spin and Zero If Total Spin Is Not the Same*/
         matrix_element += potential_matrix_original_fermions[first_orbital_1*ist.n_spatial_orbitals_third+first_orbital_2*ist.n_spatial_orbitals_sq+second_orbital_1*ist.n_spatial_orbitals+second_orbital_2].real;

       }
       if ( first_spin_1 == second_spin_2 && first_spin_2 == second_spin_1 ) {

         /*Having Issues with This Sign in the Get Sign Function - Should Be Negative*/
         matrix_element -= potential_matrix_original_fermions[first_orbital_1*ist.n_spatial_orbitals_third+second_orbital_2*ist.n_spatial_orbitals_sq+second_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real;

       }

     }
     else if ( count == 1 ) { /*If Only One Orbitals Differs********************/

          first_orbital_1 = (int)config_vector_1[electrons_minus_one]/2.0;
          first_orbital_2 = (int)config_vector_2[electrons_minus_one]/2.0;
          first_spin_orbital_1 = config_vector_1[electrons_minus_one]; 
          first_spin_orbital_2 = config_vector_2[electrons_minus_one]; 
          first_spin_1 = config_vector_1[electrons_minus_one]%2;
          first_spin_2 = config_vector_2[electrons_minus_one]%2;

          /*Run Through All Different Orbitals*/
          for (i=0; i<electrons_minus_one; i++) {
              i_orbital = (int)config_vector_1[i]/2.0;
              i_spin = config_vector_1[i]%2;
              i_spin_orbital = config_vector_1[i]; 

              /*Get Two Body Matrix Elements*/
              two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+i_spin_orbital*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+i_spin_orbital].real -= coefficient.real; 
              two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+i_spin_orbital*ist.n_spin_orbitals_sq+i_spin_orbital*ist.n_spin_orbitals+first_spin_orbital_2].real += coefficient.real; 
              two_body_density_matrix[i_spin_orbital*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+i_spin_orbital].real += coefficient.real; 
              two_body_density_matrix[i_spin_orbital*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+i_spin_orbital*ist.n_spin_orbitals+first_spin_orbital_2].real -= coefficient.real; 

              if ( first_spin_1 == first_spin_2 ) {
                 matrix_element += potential_matrix_original_fermions[first_orbital_1*ist.n_spatial_orbitals_third+first_orbital_2*ist.n_spatial_orbitals_sq+i_orbital*ist.n_spatial_orbitals+i_orbital].real;
               }
               if ( (first_spin_1 == i_spin) && (first_spin_2 == i_spin) ) {
                 matrix_element -= potential_matrix_original_fermions[first_orbital_1*ist.n_spatial_orbitals_third+i_orbital*ist.n_spatial_orbitals_sq+i_orbital*ist.n_spatial_orbitals+first_orbital_2].real; 
               }

           }/*sites i*/

     }else if ( count == 0 ) {  /*If Same Determinant*/

         /*Run Through All One and Two Electron Integrals*/
         for (i=0; i<ist.n_electrons; i++) {

           /*Find I Spin and Orbital*/
           i_spin = config_vector_1[i]%2;
           i_orbital = (int)(config_vector_1[i]/2.0);
           i_spin_orbital = config_vector_1[i]; 

           for (j=0; j<ist.n_electrons; j++) {

               /*Find J Spin*/
               j_spin = config_vector_2[j]%2;
               j_orbital = (int)(config_vector_2[j]/2.0);
               j_spin_orbital = config_vector_2[j]; 

               /*Two Body Density Matrix*/
               two_body_density_matrix[i_spin_orbital*ist.n_spin_orbitals_third+j_spin_orbital*ist.n_spin_orbitals_sq+i_spin_orbital*ist.n_spin_orbitals+j_spin_orbital].real -= coefficient.real; 
               two_body_density_matrix[i_spin_orbital*ist.n_spin_orbitals_third+j_spin_orbital*ist.n_spin_orbitals_sq+j_spin_orbital*ist.n_spin_orbitals+i_spin_orbital].real += coefficient.real;
 
               matrix_element += potential_matrix_original_fermions[i_orbital*ist.n_spatial_orbitals_third+i_orbital*ist.n_spatial_orbitals_sq+j_orbital*ist.n_spatial_orbitals+j_orbital].real;

               if ( i_spin == j_spin ) {

                  matrix_element -= potential_matrix_original_fermions[i_orbital*ist.n_spatial_orbitals_third+j_orbital*ist.n_spatial_orbitals_sq+j_orbital*ist.n_spatial_orbitals+i_orbital].real;

              }
           }
        }

       /*Now Multiply All by 1/2*/
      matrix_element *= .5;  

     }/*else*/
    
     /*Multiply By Weighting Coefficient*/
     matrix_element *= coefficient.real; 

return(matrix_element);
}

/****************************************************************************************************************************/

double get_fci_matrix_element_potential_unrestricted(int *config_vector_1,int *config_vector_2,MKL_Complex16 *two_body_density_matrix,MKL_Complex16 *potential_matrix_original_fermions_up,MKL_Complex16 *potential_matrix_original_fermions_down,MKL_Complex16 *potential_matrix_original_fermions_updown,MKL_Complex16 coefficient,int_st ist,int count) {

   /*Determines the Matrix Elements from the Sum of the Electron Integrals That Remain from the Overlap of the Two WFs*/
   int i, j;
   int j_orbital, j_spin, j_spin_orbital;
   int i_orbital, i_spin, i_spin_orbital;
   int hold;
   int electrons_minus_one = ist.n_electrons-1;
   int total_spin_1, total_spin_2;
   int first_orbital_1=0, first_orbital_2=0;
   int second_orbital_1=0, second_orbital_2=0;
   int first_spin_1=0, first_spin_2=0;
   int second_spin_1=0, second_spin_2=0;
   int first_spin_orbital_1=0, first_spin_orbital_2=0;
   int second_spin_orbital_1=0, second_spin_orbital_2=0;
   double matrix_element = 0.0;

   /*Run Through Different Orbital Counts*/
   if ( count == 2 ) { /*If Two Different Electrons**********************/

        first_orbital_1 = (int)(config_vector_1[electrons_minus_one]/2.0);
        first_orbital_2 = (int)(config_vector_2[electrons_minus_one]/2.0);
        first_spin_orbital_1 = config_vector_1[electrons_minus_one];
        first_spin_orbital_2 = config_vector_2[electrons_minus_one];
        first_spin_1 = config_vector_1[electrons_minus_one]%2;
        first_spin_2 = config_vector_2[electrons_minus_one]%2;

        second_orbital_1 = (int)(config_vector_1[ist.n_electrons-2]/2.0);
        second_orbital_2 = (int)(config_vector_2[ist.n_electrons-2]/2.0);
        second_spin_orbital_1 = config_vector_1[ist.n_electrons-2];
        second_spin_orbital_2 = config_vector_2[ist.n_electrons-2];
        second_spin_1 = config_vector_1[ist.n_electrons-2]%2;
        second_spin_2 = config_vector_2[ist.n_electrons-2]%2;

        /*Two-Body Density Matrices*/
        if ( first_spin_1 == 1 && second_spin_1 == 0 ) { 
          two_body_density_matrix[second_spin_orbital_1*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+second_spin_orbital_2].real += coefficient.real;
        }
        else {
          two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+second_spin_orbital_1*ist.n_spin_orbitals_sq+second_spin_orbital_2*ist.n_spin_orbitals+first_spin_orbital_2].real += coefficient.real;
        }

        if ( first_spin_1 == 1 && second_spin_1 == 0 ) {
         two_body_density_matrix[second_spin_orbital_1*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+second_spin_orbital_2*ist.n_spin_orbitals+first_spin_orbital_2].real -= coefficient.real ;  
        }
        else {
          two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+second_spin_orbital_1*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+second_spin_orbital_2].real -= coefficient.real ;
        }

        if ( second_spin_1 == 1 && first_spin_1 == 0 ) {
          two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+second_spin_orbital_1*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+second_spin_orbital_2].real -= coefficient.real;
        }
        else {
          two_body_density_matrix[second_spin_orbital_1*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+second_spin_orbital_2*ist.n_spin_orbitals+first_spin_orbital_2].real -= coefficient.real; 
        } 

        if ( second_spin_1 ==1 && first_spin_1 == 0 ) {
           two_body_density_matrix[first_spin_orbital_2*ist.n_spin_orbitals_third+second_spin_orbital_1*ist.n_spin_orbitals_sq+second_spin_orbital_2*ist.n_spin_orbitals+first_spin_orbital_1].real += coefficient.real;
        }
        else {
          two_body_density_matrix[second_spin_orbital_1*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+second_spin_orbital_2].real += coefficient.real;
        }        



       /*Coulomb Case**************************************************************************************************/
       if ( first_spin_1 == first_spin_2 && second_spin_1 == second_spin_2 ) {

          /*Both Spin Up*/
          if ( first_spin_1 == 0 && second_spin_1 == 0 ) {
             matrix_element += potential_matrix_original_fermions_up[first_orbital_1*ist.n_spatial_orbitals_third+first_orbital_2*ist.n_spatial_orbitals_sq+second_orbital_1*ist.n_spatial_orbitals+second_orbital_2].real;
          }
          else if ( first_spin_1 == 1 && second_spin_1 == 1 ) {
             matrix_element += potential_matrix_original_fermions_down[first_orbital_1*ist.n_spatial_orbitals_third+first_orbital_2*ist.n_spatial_orbitals_sq+second_orbital_1*ist.n_spatial_orbitals+second_orbital_2].real;
          }
          else { 
               if ( first_spin_1 == 0 ) {  /*Because Elements Are Arranged As Up-Up/Down-Down, Get Order Rights So That Up-Up Is First*/
                 matrix_element += potential_matrix_original_fermions_updown[first_orbital_1*ist.n_spatial_orbitals_third+first_orbital_2*ist.n_spatial_orbitals_sq+second_orbital_1*ist.n_spatial_orbitals+second_orbital_2].real;
               }
               else {
                 matrix_element += potential_matrix_original_fermions_updown[second_orbital_1*ist.n_spatial_orbitals_third+second_orbital_2*ist.n_spatial_orbitals_sq+first_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real; 
                }
          }

       }/**************************************************************************************************************/

       /*Exchange Case************************************************************************************************/  
       if ( first_spin_1 == second_spin_2 && first_spin_2 == second_spin_1 ) {

         /*Having Issues with This Sign in the Get Sign Function - Should Be Negative*/

         /*If Spin Up*/ 
         if ( first_spin_1 == 0 && first_spin_2 == 0 ) {
               matrix_element -= potential_matrix_original_fermions_up[first_orbital_1*ist.n_spatial_orbitals_third+second_orbital_2*ist.n_spatial_orbitals_sq+second_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real;
         }
         else if ( first_spin_1 == 1 && first_spin_2 == 1 ) {
               matrix_element -= potential_matrix_original_fermions_down[first_orbital_1*ist.n_spatial_orbitals_third+second_orbital_2*ist.n_spatial_orbitals_sq+second_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real;
         }
         else {
              if ( first_spin_1 == 0 ) {
                 matrix_element -= potential_matrix_original_fermions_updown[first_orbital_1*ist.n_spatial_orbitals_third+second_orbital_2*ist.n_spatial_orbitals_sq+second_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real;
              }
              else {
                 matrix_element -= potential_matrix_original_fermions_updown[second_orbital_1*ist.n_spatial_orbitals_third+first_orbital_2*ist.n_spatial_orbitals_sq+first_orbital_1*ist.n_spatial_orbitals+second_orbital_2].real;  
              }
         }      


       }

     }
     else if ( count == 1 ) { /*If Only One Orbitals Differs********************/

          first_orbital_1 = (int)config_vector_1[electrons_minus_one]/2.0;
          first_orbital_2 = (int)config_vector_2[electrons_minus_one]/2.0;
          first_spin_orbital_1 = config_vector_1[electrons_minus_one];
          first_spin_orbital_2 = config_vector_2[electrons_minus_one];
          first_spin_1 = config_vector_1[electrons_minus_one]%2;
          first_spin_2 = config_vector_2[electrons_minus_one]%2;

          /*Run Through All Different Orbitals*/
          for (i=0; i<electrons_minus_one; i++) {
              i_orbital = (int)config_vector_1[i]/2.0;
              i_spin = config_vector_1[i]%2;
              i_spin_orbital = config_vector_1[i];

              /*Get Two Body Matrix Elements*/
              if ( first_spin_1 == 1 && i_spin == 0 ) {  
               two_body_density_matrix[i_spin_orbital*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+i_spin_orbital*ist.n_spin_orbitals+first_spin_orbital_2].real -= coefficient.real;
               two_body_density_matrix[i_spin_orbital*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+i_spin_orbital].real += coefficient.real;
              }
              else {
                two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+i_spin_orbital*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+i_spin_orbital].real -= coefficient.real;
                two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+i_spin_orbital*ist.n_spin_orbitals_sq+i_spin_orbital*ist.n_spin_orbitals+first_spin_orbital_2].real += coefficient.real;
              }

              if ( i_spin == 1 && first_spin_1 == 0 ) {
                 two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+i_spin_orbital*ist.n_spin_orbitals_sq+i_spin_orbital*ist.n_spin_orbitals+first_spin_orbital_2].real += coefficient.real; 
                 two_body_density_matrix[first_spin_orbital_1*ist.n_spin_orbitals_third+i_spin_orbital*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+i_spin_orbital].real -= coefficient.real;
              }
              else {
                 two_body_density_matrix[i_spin_orbital*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+first_spin_orbital_2*ist.n_spin_orbitals+i_spin_orbital].real += coefficient.real;
                 two_body_density_matrix[i_spin_orbital*ist.n_spin_orbitals_third+first_spin_orbital_1*ist.n_spin_orbitals_sq+i_spin_orbital*ist.n_spin_orbitals+first_spin_orbital_2].real -= coefficient.real;
              } 

              /*Coulomb Term*****************************************************************************************************/
              if ( first_spin_1 == first_spin_2 ) {

                /*Up Spin*/
                if ( first_spin_1 == 0 && i_spin == 0 ) {   
                    matrix_element += potential_matrix_original_fermions_up[first_orbital_1*ist.n_spatial_orbitals_third+first_orbital_2*ist.n_spatial_orbitals_sq+i_orbital*ist.n_spatial_orbitals+i_orbital].real;

                } /*Down Spin*/
                else if ( first_spin_1 == 1 && i_spin == 1 ) {
                    matrix_element += potential_matrix_original_fermions_down[first_orbital_1*ist.n_spatial_orbitals_third+first_orbital_2*ist.n_spatial_orbitals_sq+i_orbital*ist.n_spatial_orbitals+i_orbital].real;
                }
                else {
                    if ( first_spin_1 == 0 ) {
                      matrix_element += potential_matrix_original_fermions_updown[first_orbital_1*ist.n_spatial_orbitals_third+first_orbital_2*ist.n_spatial_orbitals_sq+i_orbital*ist.n_spatial_orbitals+i_orbital].real;
                    }
                    else {
                      matrix_element += potential_matrix_original_fermions_updown[i_orbital*ist.n_spatial_orbitals_third+i_orbital*ist.n_spatial_orbitals_sq+first_orbital_1*ist.n_spatial_orbitals+first_orbital_2].real; 
                    } 
                }

              }
             
              /*Exchange Term**********************************************************************************************/
              if ( (first_spin_1 == i_spin) && (first_spin_2 == i_spin) ) {

                /*If Up */
                if ( first_spin_1 == 0 ) {
                     matrix_element -= potential_matrix_original_fermions_up[first_orbital_1*ist.n_spatial_orbitals_third+i_orbital*ist.n_spatial_orbitals_sq+i_orbital*ist.n_spatial_orbitals+first_orbital_2].real;

                }  
                else {
                     matrix_element -= potential_matrix_original_fermions_down[first_orbital_1*ist.n_spatial_orbitals_third+i_orbital*ist.n_spatial_orbitals_sq+i_orbital*ist.n_spatial_orbitals+first_orbital_2].real;
                } 

              }/************************************************************************************************************/

           }/*sites i*/

     }else if ( count == 0 ) {  /*If Same Determinant*/

         /*Run Through All One and Two Electron Integrals*/
         for (i=0; i<ist.n_electrons; i++) {

           /*Find I Spin and Orbital*/
           i_spin = config_vector_1[i]%2;
           i_orbital = (int)(config_vector_1[i]/2.0);
           i_spin_orbital = config_vector_1[i];

           for (j=0; j<ist.n_electrons; j++) {

               /*Find J Spin*/
               j_spin = config_vector_2[j]%2;
               j_orbital = (int)(config_vector_2[j]/2.0);
               j_spin_orbital = config_vector_2[j];

               /*Two Body Density Matrix*/
               
               if ( i_spin == 1 && j_spin == 0 )  {
                 two_body_density_matrix[j_spin_orbital*ist.n_spin_orbitals_third+i_spin_orbital*ist.n_spin_orbitals_sq+j_spin_orbital*ist.n_spin_orbitals+i_spin_orbital].real -= coefficient.real;
               }
               else {
                 two_body_density_matrix[i_spin_orbital*ist.n_spin_orbitals_third+j_spin_orbital*ist.n_spin_orbitals_sq+i_spin_orbital*ist.n_spin_orbitals+j_spin_orbital].real -= coefficient.real;
               }
               two_body_density_matrix[i_spin_orbital*ist.n_spin_orbitals_third+j_spin_orbital*ist.n_spin_orbitals_sq+j_spin_orbital*ist.n_spin_orbitals+i_spin_orbital].real += coefficient.real;

               /*Coulomb Terms************************************************************************************/
               /*If Up Spin*/
               if ( i_spin == 0 && j_spin == 0 ) {
                  matrix_element += potential_matrix_original_fermions_up[i_orbital*ist.n_spatial_orbitals_third+i_orbital*ist.n_spatial_orbitals_sq+j_orbital*ist.n_spatial_orbitals+j_orbital].real;
               } /*If Spin Down*/
               else if ( i_spin == 1 && j_spin == 1 ) {
                  matrix_element += potential_matrix_original_fermions_down[i_orbital*ist.n_spatial_orbitals_third+i_orbital*ist.n_spatial_orbitals_sq+j_orbital*ist.n_spatial_orbitals+j_orbital].real;
               }
               else {

                  if ( i_spin == 0 ) {
                    matrix_element += potential_matrix_original_fermions_updown[i_orbital*ist.n_spatial_orbitals_third+i_orbital*ist.n_spatial_orbitals_sq+j_orbital*ist.n_spatial_orbitals+j_orbital].real;
                  }
                  else {
                   matrix_element += potential_matrix_original_fermions_updown[j_orbital*ist.n_spatial_orbitals_third+j_orbital*ist.n_spatial_orbitals_sq+i_orbital*ist.n_spatial_orbitals+i_orbital].real;
                 } 
               }  

               /*Exchange Terms**********************************************************************************/ 
               if ( i_spin == j_spin ) {

                  /*If Spin Up*/
                  if ( i_spin == 0 ) {
                     matrix_element -= potential_matrix_original_fermions_up[i_orbital*ist.n_spatial_orbitals_third+j_orbital*ist.n_spatial_orbitals_sq+j_orbital*ist.n_spatial_orbitals+i_orbital].real;
                  }
                  else {
                     matrix_element -= potential_matrix_original_fermions_down[i_orbital*ist.n_spatial_orbitals_third+j_orbital*ist.n_spatial_orbitals_sq+j_orbital*ist.n_spatial_orbitals+i_orbital].real;
                  }
              }/****************************************************************************************************/

           }
        }

       /*Now Multiply All by 1/2*/
      matrix_element *= .5;

     }/*else*/

     /*Multiply By Weighting Coefficient*/
     matrix_element *= coefficient.real;

return(matrix_element);
}

/****************************************************************************************************************************************/

void compute_walker_energy_chemistry_mixed_multi_energy_unrestricted_trial(MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *kinetic_original_matrix_up,MKL_Complex16 *kinetic_original_matrix_down,MKL_Complex16 *potential_original_matrix_up,MKL_Complex16 *potential_original_matrix_down,MKL_Complex16 *potential_original_matrix_updown,MKL_Complex16 *total_ke,MKL_Complex16 *total_pe,int_st ist) {

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

    /*Run Through Determinants*/
    ke.real = ke.imag = 0.0;
    pe.real = pe.imag = 0.0;

    /*Determine Full Overlap*/
    for (i=0; i<ist.n_determinants_trial_energy; i++) {

       /*Zero All Total Vectors*/
       //cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_up,ist.n_spatial_orbitals,ist.n_up,One,&overlap_inverse_up[i*ist.n_up_sq],ist.n_up,&trial_wf_up_energy[i*ist.n_spatial_orbitals_n_up],ist.n_up,Zero,stored_product1_up,ist.n_spatial_orbitals);
       //cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasTrans,ist.n_down,ist.n_spatial_orbitals,ist.n_down,One,&overlap_inverse_down[i*ist.n_down_sq],ist.n_down,&trial_wf_down_energy[i*ist.n_spatial_orbitals_n_down],ist.n_down,Zero,stored_product1_down,ist.n_spatial_orbitals);

       /*Now Get Total Density Matrix*/
       //cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_up,One,wf_up,ist.n_up,stored_product1_up,ist.n_spatial_orbitals,Zero,density_matrix_up,ist.n_spatial_orbitals);
       //cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_spatial_orbitals,ist.n_spatial_orbitals,ist.n_down,One,wf_down,ist.n_down,stored_product1_down,ist.n_spatial_orbitals,Zero,density_matrix_down,ist.n_spatial_orbitals);

       /*Get Ke and PE*/
       local_ke = compute_kinetic_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,kinetic_original_matrix_up,kinetic_original_matrix_down,ist);
       local_pe = compute_potential_energy_density_matrix_unrestricted(density_matrix_up,density_matrix_down,potential_original_matrix_up,potential_original_matrix_down,potential_original_matrix_updown,ist);

       /*Get Product of dets*/
       product_dets = Cmul(conjugate(trial_determinant_coefficients_energy[i]), Cmul(det_overlap_up[i], det_overlap_down[i]));

       /*Scale By Overlap and Coefficient*/
       local_ke = Cmul(local_ke, product_dets);
       local_pe = Cmul(local_pe, product_dets);

       /*Add Up Local Pe and Ke*/
       ke = Cadd(ke, local_ke);
       pe = Cadd(pe, local_pe);


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


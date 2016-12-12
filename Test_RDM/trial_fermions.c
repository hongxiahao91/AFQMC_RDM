#include "afqmc.h"

/***************************************************************************************************************/

void trial_random_fermions(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,int_st ist) {

  /*Sets the Trial Wavefunction to the Identity Matrix*/
  int ip, ip2;
  double value; 
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;  
  long int tidum;

  Randomize(); tidum = -random(); 

  czero_vec(trial_wf_up,ist.n_sites*ist.n_up); 
  czero_vec(trial_wf_down,ist.n_sites*ist.n_down); 

  for ( ip = 0; ip < ist.n_up ; ip++ ) {
     trial_wf_up[ip * ist.n_up + ip] = One; 
  }
  for ( ip = 0; ip < ist.n_down ; ip++ ) {
     trial_wf_down[ip * ist.n_down + ip] = One; 
  }   

  /*Try Random Instead*/
   for (ip2 = 0; ip2 < ist.n_up; ip2++ ) {

     value = 0.0;
     for (ip = 0; ip<ist.n_sites; ip++) {
      trial_wf_up[ip * ist.n_up + ip2].real = ran1(&tidum);
       value += trial_wf_up[ip*ist.n_up+ip2].real * trial_wf_up[ip*ist.n_up+ip2].real;
     }
     value = 1.0/sqrt(value);

     for (ip = 0; ip<ist.n_sites; ip++ ) {
      trial_wf_up[ip * ist.n_up + ip2].real *= value;
     }

   }

   for ( ip2 = 0; ip2 < ist.n_down; ip2++ ) {

     value = 0.0;
     for (ip = 0; ip<ist.n_sites; ip++) {
      trial_wf_down[ip * ist.n_down + ip2].real = ran1(&tidum);
       value += trial_wf_down[ip*ist.n_down+ip2].real * trial_wf_down[ip*ist.n_down+ip2].real;
     }
     value = 1.0/sqrt(value);

     for (ip = 0; ip<ist.n_sites; ip++ ) {
      trial_wf_down[ip * ist.n_down + ip2].real *= value;
     }
   }

return; 
}

/***************************************************************************************************************/

void trial_random_fermions_up(MKL_Complex16 *trial_wf_up,int_st ist) {

  /*Sets the Trial Wavefunction to the Identity Matrix*/
  int ip, ip2;
  double value;
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;
  long int tidum;

  Randomize(); tidum = -random();

  czero_vec(trial_wf_up,ist.n_sites*ist.n_up);

  for ( ip = 0; ip < ist.n_up ; ip++ ) {
     trial_wf_up[ip * ist.n_up + ip] = One;
  }

  /*Try Random Instead*/
   for (ip2 = 0; ip2 < ist.n_up; ip2++ ) {

     value = 0.0;
     for (ip = 0; ip<ist.n_sites; ip++) {
      trial_wf_up[ip * ist.n_up + ip2].real = ran1(&tidum);
       value += trial_wf_up[ip*ist.n_up+ip2].real * trial_wf_up[ip*ist.n_up+ip2].real;
     }
     value = 1.0/sqrt(value);

     for (ip = 0; ip<ist.n_sites; ip++ ) {
      trial_wf_up[ip * ist.n_up + ip2].real *= value;
     }

   }

return;
}

/***************************************************************************************************************/

void trial_random_fermions_down(MKL_Complex16 *trial_wf_down,int_st ist) {

  /*Sets the Trial Wavefunction to the Identity Matrix*/
  int ip, ip2;
  double value;
  MKL_Complex16 One; One.real = 1.0; One.imag = 0.0;
  long int tidum;
   
  Randomize(); tidum = -random();
   
  czero_vec(trial_wf_down,ist.n_sites_n_down);
   
  for ( ip = 0; ip < ist.n_down ; ip++ ) {
     trial_wf_down[ip * ist.n_down + ip] = One;
  }
      
  /*Try Random Instead*/
   for (ip2 = 0; ip2 < ist.n_down; ip2++ ) {

     value = 0.0;
     for (ip = 0; ip<ist.n_sites; ip++) {
      trial_wf_down[ip * ist.n_down + ip2].real = ran1(&tidum);
       value += trial_wf_down[ip*ist.n_down+ip2].real * trial_wf_down[ip*ist.n_down+ip2].real;
     }
     value = 1.0/sqrt(value);
   
     for (ip = 0; ip<ist.n_sites; ip++ ) {
      trial_wf_down[ip * ist.n_down + ip2].real *= value;
     }

   }

return;
}

/***********************************************************************************/

void trial_identity_fermions(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,int_st ist) {

   /*Sets the Trial Wavefunction to the Restricted Hartree-Fock Wavefunction Using Coefficients from MolPro*/ 
   int ip; 

   czero_vec(trial_wf_up,ist.n_sites*ist.n_up);
   czero_vec(trial_wf_down,ist.n_sites*ist.n_down); 

   /*Get Up Components First*/
   for (ip = 0; ip < ist.n_up; ip++) {
      trial_wf_up[ip * ist.n_up + ip].real = 1.0; 
   }

   /*Get Down*/
   for (ip = 0; ip < ist.n_down; ip++) {
      trial_wf_down[ip*ist.n_down+ip].real = 1.0; 
   } 

return; 
}

/***********************************************************************************/

void trial_identity_fermions_up(MKL_Complex16 *trial_wf_up,int_st ist) {

   /*Sets the Trial Wavefunction to the Restricted Hartree-Fock Wavefunction Using Coefficients from MolPro*/
   int ip;
   
  czero_vec(trial_wf_up,ist.n_sites*ist.n_up);

   /*Get Up Components First*/
   for (ip = 0; ip < ist.n_up; ip++) {
      trial_wf_up[ip * ist.n_up + ip].real = 1.0;
   }

return;
}

/***********************************************************************************/

void trial_identity_fermions_down(MKL_Complex16 *trial_wf_down,int_st ist) {

   /*Sets the Trial Wavefunction to the Restricted Hartree-Fock Wavefunction Using Coefficients from MolPro*/
   int ip;
   
   czero_vec(trial_wf_down,ist.n_sites_n_down);

   /*Get Up Components First*/
   for (ip = 0; ip < ist.n_down; ip++) {
      trial_wf_down[ip * ist.n_down + ip].real = 1.0;
   }

return;
}

/**************************************************************************************/

int trial_chemistry_multi_readin(MKL_Complex16 *trial_wf_up_energy,MKL_Complex16 *trial_wf_down_energy,MKL_Complex16 *trial_determinant_coefficients_energy,MKL_Complex16 *trial_wf_up_phaseless,MKL_Complex16 *trial_wf_down_phaseless,MKL_Complex16 *trial_determinant_coefficients_phaseless,int *list_occupied_orbitals,int_st ist) {

   int l, i, j, k;
   int orbital;
   MKL_Complex16 p;  
   int *trial_wf_electrons_temp; 
   MKL_Complex16 *trial_determinant_coefficients_energy_temp;  
   int *temp_electron_config_storage; 
   int n_max_trial_orbital;
   int store_orbital, check_orbital;  
   FILE *input = fopen("trialwavefunction_readin.par", "r");
   FILE *output = fopen("trialwavefunction.dat", "a+");

   /*First Zero Trial Wavefunctions*/
   trial_wf_electrons_temp = (int*)calloc((ist.n_determinants_max+1)*ist.n_electrons,sizeof(int)); 
   trial_determinant_coefficients_energy_temp = (MKL_Complex16*)calloc(ist.n_determinants_max+1,sizeof(MKL_Complex16)); 
   temp_electron_config_storage = (int*)calloc(ist.n_electrons,sizeof(int)); 

   izero_vec(list_occupied_orbitals,ist.n_spatial_orbitals); 

   czero_vec(trial_wf_up_energy,ist.n_determinants_energy_n_spatial_orbitals_n_up); 
   czero_vec(trial_wf_down_energy,ist.n_determinants_energy_n_spatial_orbitals_n_down); 

   czero_vec(trial_wf_up_phaseless,ist.n_determinants_phaseless_n_spatial_orbitals_n_up);
   czero_vec(trial_wf_down_phaseless,ist.n_determinants_phaseless_n_spatial_orbitals_n_down);

   /*Read In Inputs*/
   for (i=0; i<ist.n_determinants_max; i++) {
      fscanf(input, "%lf", &trial_determinant_coefficients_energy_temp[i+1].real);
      
      for (j=0; j<ist.n_electrons; j++) {
        fscanf(input, "%d", &trial_wf_electrons_temp[(i+1)*ist.n_electrons+j]); 
      }
   }    

   /*Now Reorganize Determinants In Terms of Coefficient Magnitude*/
    for (i=1;i<ist.n_determinants_max;i++) {
       k=i;
       p=trial_determinant_coefficients_energy_temp[k=i];
       for (l=0; l<ist.n_electrons; l++) {
        temp_electron_config_storage[l] = trial_wf_electrons_temp[i*ist.n_electrons+l];  
       } 

       for (j=i+1;j<=ist.n_determinants_max;j++) {
         if (Cabs(trial_determinant_coefficients_energy_temp[j]) >= Cabs(p)) {
           p=trial_determinant_coefficients_energy_temp[k=j]; 
           for (l=0; l<ist.n_electrons; l++) {
             temp_electron_config_storage[l] = trial_wf_electrons_temp[j*ist.n_electrons+l]; 
           }
          } 
         }
  
          if (k != i) {
            trial_determinant_coefficients_energy_temp[k] = trial_determinant_coefficients_energy_temp[i]; 
            for (l=0; l<ist.n_electrons; l++) {
               trial_wf_electrons_temp[k*ist.n_electrons+l] = trial_wf_electrons_temp[i*ist.n_electrons+l]; 
             }
             trial_determinant_coefficients_energy_temp[i] = p; 
             for (l=0; l<ist.n_electrons; l++) {
                trial_wf_electrons_temp[i*ist.n_electrons+l] = temp_electron_config_storage[l]; 
            }
          }
          
     }

    /*Now Copy Over To Appropriate Wavefunctions*/
    n_max_trial_orbital = 0; 
    for (i=0; i<ist.n_determinants_trial_energy; i++) {

     if ( i >= ist.n_determinants_trial_phaseless) {
      trial_determinant_coefficients_energy[i] = trial_determinant_coefficients_energy_temp[i+1]; 
      fprintf(output, "%f\t %f\n", trial_determinant_coefficients_energy[i].real, trial_determinant_coefficients_energy[i].imag);

      for (j=0; j<ist.n_electrons; j++) {
     
        /*Determine If Orbital Already Stored*/
        store_orbital = trial_wf_electrons_temp[(i+1)*ist.n_electrons+j]-1; 
        check_orbital = 2; 
        for (k=0; k<n_max_trial_orbital; k++) {
 
          if ( store_orbital == list_occupied_orbitals[k] ) {
            check_orbital = 1; 
            k = n_max_trial_orbital; 
          }

        }
        if ( n_max_trial_orbital == 0 ) {
          check_orbital = 2 ; 
        }
  
        /*Store Orbital Number If Not Already Stored*/   
        if ( check_orbital == 2 ) { 
          list_occupied_orbitals[n_max_trial_orbital] = store_orbital; 
          n_max_trial_orbital++; 
        }

        if ( j < ist.n_up ) {
           l=i*ist.n_spatial_orbitals_n_up+(trial_wf_electrons_temp[(i+1)*ist.n_electrons+j]-1)*ist.n_up+j; 
           trial_wf_up_energy[l].real = 1.0; 
        }
        else {    
           l=i*ist.n_spatial_orbitals_n_down+(trial_wf_electrons_temp[(i+1)*ist.n_electrons+j]-1)*ist.n_down+(j-ist.n_up); 
           trial_wf_down_energy[l].real = 1.0; 
        } 
    
      }
     } /*If Determinants Bigger Than Phaseless*/
     else {
       trial_determinant_coefficients_energy[i] = trial_determinant_coefficients_phaseless[i] = trial_determinant_coefficients_energy_temp[i+1]; 
       fprintf(output, "%f\t %f\n", trial_determinant_coefficients_energy[i].real, trial_determinant_coefficients_energy[i].imag);      
 
       for (j=0; j<ist.n_electrons; j++) {

        /*Determine If Orbital Already Stored*/      
        store_orbital = trial_wf_electrons_temp[(i+1)*ist.n_electrons+j]-1;
        check_orbital = 2;
        for (k=0; k<n_max_trial_orbital; k++) {
          if ( store_orbital == list_occupied_orbitals[k] ) {
            check_orbital = 1;
            k = n_max_trial_orbital;
          }
        }
        if ( n_max_trial_orbital == 0 ) {
           check_orbital = 2; 
        }

        /*Store Orbital Number If Not Already Stored*/  
        if ( check_orbital == 2 ) {
          list_occupied_orbitals[n_max_trial_orbital] = store_orbital;
          n_max_trial_orbital++;
        }

          if ( j < ist.n_up ) {
             l= i*ist.n_spatial_orbitals_n_up+(trial_wf_electrons_temp[(i+1)*ist.n_electrons+j]-1)*ist.n_up+j;
             trial_wf_up_phaseless[l].real = trial_wf_up_energy[l].real = 1.0; 
          }
          else {   
             l=i*ist.n_spatial_orbitals_n_down+(trial_wf_electrons_temp[(i+1)*ist.n_electrons+j]-1)*ist.n_down+(j-ist.n_up);
             trial_wf_down_phaseless[l].real = trial_wf_down_energy[l].real = 1.0;    
         } 
    
       } /*j*/     
     } /*else*/
    }

   /*Now Print Out Determinants*/ 
   for (i=0; i<ist.n_determinants_trial_energy; i++) {
     for (j=0; j<ist.n_spatial_orbitals; j++) {
      for (k=0; k<ist.n_up; k++) {
        fprintf(output, "%f\t", trial_wf_up_energy[i*ist.n_spatial_orbitals_n_up+j*ist.n_up+k].real); 
      }
      fprintf(output, "\n"); 
     }
     fprintf(output, "\n\n"); fflush(output); 

     for (j=0; j<ist.n_spatial_orbitals; j++) {
      for (k=0; k<ist.n_down; k++) {
        fprintf(output, "%f\t", trial_wf_down_energy[i*ist.n_spatial_orbitals_n_down+j*ist.n_down+k].real);
      }
      fprintf(output, "\n");
     }
     fprintf(output, "\n\n"); fflush(output);  
   
   }

free(trial_wf_electrons_temp); 
free(trial_determinant_coefficients_energy_temp); 
free(temp_electron_config_storage); 
return(n_max_trial_orbital); 
}

/**************************************************************************************/

int trial_chemistry_restricted_readin(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,int *list_occupied_orbitals, int_st ist) {

   int j, k, orbital, m;
   int check_orbital, n_max_orbitals=0;  
   double coefficient; 
   FILE *input = fopen("trialwavefunction_readin.par", "r");
   FILE *output = fopen("trialwavefunction.dat", "a+");

   czero_vec(trial_wf_up,ist.n_spatial_orbitals_n_up);
   czero_vec(trial_wf_down,ist.n_spatial_orbitals_n_down);

   /*Read In Inputs*/
   fscanf(input, "%lf", &coefficient); 

   /*Get Orbitals From Which to Make WF*/
   for (j=0; j<ist.n_electrons; j++) {
      fscanf(input, "%d", &orbital);
      orbital -= 1;

      /*Store New Orbital Numbers*/
      check_orbital = 0; 
      for (m=0; m<n_max_orbitals; m++) {
        if ( orbital == list_occupied_orbitals[m] ) {
          check_orbital = 1; 
        }
      }
      if ( n_max_orbitals == 0 ) {
          check_orbital = 1; 
      }

      /*Store New Orbitals*/
      if ( check_orbital == 1 ) {
        list_occupied_orbitals[n_max_orbitals] = orbital; 
        n_max_orbitals++; 
      }  

      /*If an Up Orbital, Scan In and Make Up*/
      if ( j < ist.n_up ) {
        trial_wf_up[orbital*ist.n_up+j].real = 1.0;
      }
      else {
        trial_wf_down[orbital*ist.n_down+(j-ist.n_up)].real = 1.0;
      }

   } /*Run Through Electrons*/
   fprintf(output, "\n\n"); fflush(output);

   /*Now Print Out Determinants*/
   for (j=0; j<ist.n_spatial_orbitals; j++) {
      for (k=0; k<ist.n_up; k++) {
        fprintf(output, "%f\t", trial_wf_up[j*ist.n_up+k].real);
      }
      fprintf(output, "\n");
     }
     fprintf(output, "\n\n"); fflush(output);

     for (j=0; j<ist.n_spatial_orbitals; j++) {
      for (k=0; k<ist.n_down; k++) {
        fprintf(output, "%f\t", trial_wf_down[j*ist.n_down+k].real);
      }
      fprintf(output, "\n");
     }
    fflush(output); 

return(n_max_orbitals);
}

/**************************************************************************************/

void trial_fermions_readin(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,int_st ist) {

   int i, j, k;
   FILE *input = fopen("trialwavefunction_readin.par", "r");
   FILE *output = fopen("trialwavefunction.dat", "a+");

   /*Read In up*/
   for (i=0; i<ist.n_determinants_trial_energy; i++) {
    for (j=0; j<ist.n_sites; j++) {
     for (k=0; k<ist.n_up; k++) {
        fscanf(input, "%lf %lf", &trial_wf_up[i*ist.n_sites_n_up+j*ist.n_up+k].real, &trial_wf_up[i*ist.n_sites_n_up+j*ist.n_up+k].imag);
        fprintf(output, "%f+%fi\t", trial_wf_up[i*ist.n_sites_n_up+j*ist.n_up+k].real, trial_wf_up[i*ist.n_sites_n_up+j*ist.n_up+k].imag);
     }
     fprintf(output, "\n");
    }
    fprintf(output, "\n");
   }
   fprintf(output, "\n\n"); fflush(output);

  /*Read in Down*/
  for (i=0; i<ist.n_determinants_trial_energy; i++) {
    for (j=0; j<ist.n_sites; j++) {
      for (k=0; k<ist.n_down; k++) {
         fscanf(input, "%lf %lf", &trial_wf_down[i*ist.n_sites_n_down+j*ist.n_down+k].real, &trial_wf_down[i*ist.n_sites_n_down+j*ist.n_down+k].imag);
         fprintf(output, "%f+%fi\t", trial_wf_down[i*ist.n_sites_n_down+j*ist.n_down+k].real, trial_wf_down[i*ist.n_sites_n_down+j*ist.n_down+k].imag);
      }
      fprintf(output, "\n");
     }
     fprintf(output, "\n");
   }
   fprintf(output, "\n\n"); fflush(output);

return; 
}

/**************************************************************************************/

void trial_fermions_readin_up(MKL_Complex16 *trial_wf_up,int_st ist) {

   int i, j, k;
   FILE *input = fopen("trialwavefunction_readin.par", "r");
   FILE *output = fopen("trialwavefunction.dat", "a+");

   /*Read In up*/
   for (i=0; i<ist.n_determinants_trial_energy; i++) {
    for (j=0; j<ist.n_sites; j++) {
     for (k=0; k<ist.n_up; k++) {
        fscanf(input, "%lf %lf", &trial_wf_up[i*ist.n_sites_n_up+j*ist.n_up+k].real, &trial_wf_up[i*ist.n_sites_n_up+j*ist.n_up+k].imag);
        fprintf(output, "%f+%fi\t", trial_wf_up[i*ist.n_sites_n_up+j*ist.n_up+k].real, trial_wf_up[i*ist.n_sites_n_up+j*ist.n_up+k].imag);
     }
     fprintf(output, "\n");
    }
    fprintf(output, "\n");
   }
   fprintf(output, "\n\n"); fflush(output);

return; 
}

/**************************************************************************************/

void trial_fermions_readin_down(MKL_Complex16 *trial_wf_down,int_st ist) {

   int i, j, k;
   FILE *input = fopen("trialwavefunction_readin.par", "r");
   FILE *output = fopen("trialwavefunction.dat", "a+");

  /*Read in Down*/
  for (i=0; i<ist.n_determinants_trial_energy; i++) {
    for (j=0; j<ist.n_sites; j++) {
      for (k=0; k<ist.n_down; k++) {
         fscanf(input, "%lf %lf", &trial_wf_down[i*ist.n_sites_n_down+j*ist.n_down+k].real, &trial_wf_down[i*ist.n_sites_n_down+j*ist.n_down+k].imag);
         fprintf(output, "%f+%fi\t", trial_wf_down[i*ist.n_sites_n_down+j*ist.n_down+k].real, trial_wf_down[i*ist.n_sites_n_down+j*ist.n_down+k].imag);
      }
      fprintf(output, "\n");
     }
     fprintf(output, "\n");
   }
   fprintf(output, "\n\n"); fflush(output);

return; 
}

/*******************************************************************************************/

int trial_chemistry_restricted_noreadin(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *kinetic_matrix_original_fermions_up,int *list_occupied_orbitals,int_st ist) {

   /*Sets the Trial Wavefunction to the Restricted Hartree-Fock Wavefunction Using Coefficients from MolPro*/
   int count;
   int n_max_orbitals;  
   int flag_taken, hold_element; 
   int i, j, k; 
   int *sorted_up; 
   FILE *pf = fopen("trialwavefunction.dat", "a+"); 
   double current_element; 

   sorted_up = (int*)calloc(ist.n_up, sizeof(int)); 

   /*Set Sorted Up and Down to Random Numbers */
   for (i=0; i<ist.n_up; i++) {
     sorted_up[i] = 1000; 
    }

   /*Find Lowest Energy Orbitals*/
   count = 0; 
   while ( count < ist.n_up ) {

    current_element = 100000; 
    for (i=0; i<ist.n_spatial_orbitals; i++) {

     if ( kinetic_matrix_original_fermions_up[i*ist.n_spatial_orbitals+i].real < current_element ) {
        /*Ensure Element Not Taken*/
        flag_taken = 0; 
        for (j=0; j< ist.n_up; j++) {
         if ( sorted_up[j] == i ) {
           flag_taken = 1; 
         }
        }
        if ( flag_taken == 0 ) {
          hold_element = i; 
          current_element = kinetic_matrix_original_fermions_up[i*ist.n_spatial_orbitals+i].real; 
        }
     } 
     
   }
   sorted_up[count] = hold_element;
   count++;

  }

  /*Now Zero Trial WFs */
  czero_vec(trial_wf_up,ist.n_spatial_orbitals_n_up);
  czero_vec(trial_wf_down,ist.n_spatial_orbitals_n_down);

  /*Use Ordered Orbital Energy List to Initialize Trial WFs*/
  /*I Assume That All Trial Determinants Are the Same for Now - FIXXX!!!******************************************************************************************************/
  for (i=0; i<ist.n_up; i++) {
      trial_wf_up[sorted_up[i]*ist.n_up+i].real = 1.0;
  } 
  
  for (i=0; i<ist.n_down; i++) {
    trial_wf_down[sorted_up[i]*ist.n_down+i].real = 1.0; 
  }

  /*Store Max Number of Orbitals*/
  if ( ist.n_up > ist.n_down ) {
    n_max_orbitals = ist.n_up; 
    for (i=0; i<n_max_orbitals; i++) {
       list_occupied_orbitals[i] = i; 
    }   
  }
  else {
    n_max_orbitals = ist.n_down; 
    for (i=0; i<n_max_orbitals; i++) {
       list_occupied_orbitals[i] = i; 
    }
  }  

  /*Printing*/
   for (j=0; j<ist.n_spatial_orbitals; j++) {
     for (k=0; k<ist.n_up; k++) {
       fprintf(pf, "%f\t", trial_wf_up[j*ist.n_up+k].real); 
     }
     fprintf(pf, "\n"); 
    }
   fprintf(pf, "\n\n"); fflush(pf);  

  for (j = 0; j<ist.n_spatial_orbitals; j++) {
    for (k=0; k<ist.n_down; k++) {
      fprintf(pf, "%f\t", trial_wf_down[j*ist.n_down+k].real);
    }
    fprintf(pf, "\n");
  }
 fprintf(pf, "\n\n"); fflush(pf);
 

free(sorted_up); 
fclose(pf); 
return(n_max_orbitals);
}

/*******************************************************************************************************/

int trial_chemistry_unrestricted_noreadin(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *kinetic_matrix_original_fermions_up,MKL_Complex16 *kinetic_matrix_original_fermions_down,int *list_occupied_orbitals,int_st ist) {

   /*Sets the Trial Wavefunction to the Restricted Hartree-Fock Wavefunction Using Coefficients from MolPro*/
   int count;
   int flag_taken, hold_element;
   int n_max_orbitals; 
   int i, j, k;
   int *sorted_up, *sorted_down;
   FILE *pf = fopen("trialwavefunction.dat", "a+");
   double current_element;

   sorted_up = (int*)calloc(ist.n_up, sizeof(int));
   sorted_down = (int*)calloc(ist.n_down, sizeof(int)); 

   /*Set Sorted Up and Down to Random Numbers */
   for (i=0; i<ist.n_up; i++) {
     sorted_up[i] = 1000;
    }
   for (i=0; i<ist.n_down; i++) {
     sorted_down[i] = 1000; 
   } 

   /*Find Lowest Energy Orbitals*/
   count = 0;
   while ( count < ist.n_up ) {

    current_element = 100000;
    for (i=0; i<ist.n_spatial_orbitals; i++) {

     if ( kinetic_matrix_original_fermions_up[i*ist.n_spatial_orbitals+i].real < current_element ) {
        /*Ensure Element Not Taken*/
        flag_taken = 0;
        for (j=0; j< ist.n_up; j++) {
         if ( sorted_up[j] == i ) {
           flag_taken = 1;
         }
        }
        if ( flag_taken == 0 ) {
          hold_element = i;
          current_element = kinetic_matrix_original_fermions_up[i*ist.n_spatial_orbitals+i].real;
        }
     }

   }
   sorted_up[count] = hold_element;
   count++;
  }

   /*Find Lowest Energy Orbitals*/
   count = 0;
   while ( count < ist.n_down ) {

    current_element = 100000;
    for (i=0; i<ist.n_spatial_orbitals; i++) {

     if ( kinetic_matrix_original_fermions_down[i*ist.n_spatial_orbitals+i].real < current_element ) {
        /*Ensure Element Not Taken*/
        flag_taken = 0;
        for (j=0; j< ist.n_down; j++) {
         if ( sorted_down[j] == i ) {
           flag_taken = 1;
         }
        }
        if ( flag_taken == 0 ) {
          hold_element = i;
          current_element = kinetic_matrix_original_fermions_down[i*ist.n_spatial_orbitals+i].real;
        }
     }

   }
   sorted_down[count] = hold_element;
   count++;
  }


  /*Now Zero Trial WFs */
  czero_vec(trial_wf_up,ist.n_spatial_orbitals_n_up);
  czero_vec(trial_wf_down,ist.n_spatial_orbitals_n_down);

  /*Use Ordered Orbital Energy List to Initialize Trial WFs*/
  for (i=0; i<ist.n_up; i++) {
    trial_wf_up[sorted_up[i]*ist.n_up+i].real = 1.0;
  }

  for (i=0; i<ist.n_down; i++) {
    trial_wf_down[sorted_down[i]*ist.n_down+i].real = 1.0;
  }

  /*Make LIst of All Occupied Orbitals*/ 
  if ( ist.n_up > ist.n_down ) {
    n_max_orbitals = ist.n_up; 
    for (i=0; i<ist.n_up; i++) {
      list_occupied_orbitals[i] = i; 
     }
  }
  else {
     n_max_orbitals = ist.n_down; 
     for (i=0; i<ist.n_down; i++) {
       list_occupied_orbitals[i] = i; 
     }
  } 


  for (j=0; j<ist.n_spatial_orbitals; j++) {
   for (k=0; k<ist.n_up; k++) {
     fprintf(pf, "%f\t", trial_wf_up[j*ist.n_up+k].real);
   }
   fprintf(pf, "\n");
  }
  fprintf(pf, "\n\n"); fflush(pf);

  for (j = 0; j<ist.n_spatial_orbitals; j++) {
   for (k=0; k<ist.n_down; k++) {
     fprintf(pf, "%f\t", trial_wf_down[j*ist.n_down+k].real);
   }
   fprintf(pf, "\n");
  }
  fprintf(pf, "\n\n"); fflush(pf);


free(sorted_up);
free(sorted_down); 
fclose(pf);
return(n_max_orbitals); 
}

/***********************************************************************************/

void trial_free_fermions(MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,double *kinetic_eigs,double *kinetic_eigvecs,int_st ist,cns_st cns) {

   /*Sets the Trial Wavefunction to the Trial Kinetic WF*/        
   int i, j, k;
   double p;  
   FILE *pf2 = fopen("trial_wf.dat", "a+");

   /*Order Eigenvalues and Eigenvectors in Descending Order First*/
   for (i=1;i<ist.n_sites;i++) {
     k=i;
     p=kinetic_eigs[k-1];
     for (j=i+1;j<=ist.n_sites;j++) {
       if (kinetic_eigs[j-1] >= p) {k=j; p=kinetic_eigs[k-1]; };
     }
     if (k != i) {
        kinetic_eigs[k-1]=kinetic_eigs[i-1];
        kinetic_eigs[i-1]=p;
        for (j=1;j<=ist.n_sites;j++) {
           p=kinetic_eigvecs[(j-1)*ist.n_sites+i-1];
           kinetic_eigvecs[(j-1)*ist.n_sites+i-1]=kinetic_eigvecs[(j-1)*ist.n_sites+k-1];
           kinetic_eigvecs[(j-1)*ist.n_sites+k-1]=p;
         }
      }
    }

   /*Copy Lowest Up Down Eigenvectors Into Matrices*/
   for (i=0; i<ist.n_up; i++) {
      for (j=0; j<ist.n_sites; j++) {
        trial_wf_up[j*ist.n_up+i].real = kinetic_eigvecs[j*ist.n_sites+(ist.n_sites-i-1)]; 
        trial_wf_up[j*ist.n_up+i].imag = 0.0; 
      }
   }

   for (i=0; i<ist.n_down; i++) {
      for (j=0; j<ist.n_sites; j++) {
        trial_wf_down[j*ist.n_down+i].real = kinetic_eigvecs[j*ist.n_sites+(ist.n_sites-i-1)];
        trial_wf_up[j*ist.n_down+i].imag = 0.0;  
      }
   }  

   /*Save Trial Energy*/
   cns.trial_energy = kinetic_eigs[ist.n_sites-1]; 

   /*Print Trial Matrices*/
   print_cmat(trial_wf_up,ist.n_sites,ist.n_up,"trial_wf.dat");
   print_cmat(trial_wf_down,ist.n_sites,ist.n_down,"trial_wf.dat");
   fprintf(pf2, "Trial Energy: %g\n", cns.trial_energy); fflush(pf2);   

fclose(pf2); 
return; 
}  

/************************************************************************************/

void trial_free_fermions_up(MKL_Complex16 *trial_wf_up,double *kinetic_eigs,double *kinetic_eigvecs,int_st ist,cns_st cns) {

   /*Sets the Trial Wavefunction to the Trial Kinetic WF*/
   int i, j, k;
   double p;
   FILE *pf2 = fopen("trial_wf.dat", "a+");

   /*Order Eigenvalues and Eigenvectors in Descending Order First*/
   for (i=1;i<ist.n_sites;i++) {
     k=i;
     p=kinetic_eigs[k-1];
     for (j=i+1;j<=ist.n_sites;j++) {
       if (kinetic_eigs[j-1] >= p) {k=j; p=kinetic_eigs[k-1]; };
     }
     if (k != i) {
        kinetic_eigs[k-1]=kinetic_eigs[i-1];
        kinetic_eigs[i-1]=p;
        for (j=1;j<=ist.n_sites;j++) {
           p=kinetic_eigvecs[(j-1)*ist.n_sites+i-1];
           kinetic_eigvecs[(j-1)*ist.n_sites+i-1]=kinetic_eigvecs[(j-1)*ist.n_sites+k-1];
           kinetic_eigvecs[(j-1)*ist.n_sites+k-1]=p;
         }
      }
    }

   /*Copy Lowest Up Down Eigenvectors Into Matrices*/
   for (i=0; i<ist.n_up; i++) {
      for (j=0; j<ist.n_sites; j++) {
        trial_wf_up[j*ist.n_up+i].real = kinetic_eigvecs[j*ist.n_sites+(ist.n_sites-i-1)];
        trial_wf_up[j*ist.n_up+i].imag = 0.0;
      }
   }

   /*Save Trial Energy*/
   cns.trial_energy = kinetic_eigs[ist.n_sites-1];

   /*Print Trial Matrices*/
   print_cmat(trial_wf_up,ist.n_sites,ist.n_up,"trial_wf.dat");
   fprintf(pf2, "Trial Energy: %g\n", cns.trial_energy); fflush(pf2);

fclose(pf2);
return;
}

/************************************************************************************/

void trial_free_fermions_down(MKL_Complex16 *trial_wf_down,double *kinetic_eigs,double *kinetic_eigvecs,int_st ist,cns_st cns) {

   /*Sets the Trial Wavefunction to the Trial Kinetic WF*/
   int i, j, k;
   double p;
   FILE *pf2 = fopen("trial_wf.dat", "a+");

   /*Order Eigenvalues and Eigenvectors in Descending Order First*/
   for (i=1;i<ist.n_sites;i++) {
     k=i;
     p=kinetic_eigs[k-1];
     for (j=i+1;j<=ist.n_sites;j++) {
       if (kinetic_eigs[j-1] >= p) {k=j; p=kinetic_eigs[k-1]; };
     }
     if (k != i) {
        kinetic_eigs[k-1]=kinetic_eigs[i-1];
        kinetic_eigs[i-1]=p;
        for (j=1;j<=ist.n_sites;j++) {
           p=kinetic_eigvecs[(j-1)*ist.n_sites+i-1];
           kinetic_eigvecs[(j-1)*ist.n_sites+i-1]=kinetic_eigvecs[(j-1)*ist.n_sites+k-1];
           kinetic_eigvecs[(j-1)*ist.n_sites+k-1]=p;
         }
      }
    }

   /*Copy Lowest Up Down Eigenvectors Into Matrices*/
   for (i=0; i<ist.n_down; i++) {
      for (j=0; j<ist.n_sites; j++) {
        trial_wf_down[j*ist.n_down+i].real = kinetic_eigvecs[j*ist.n_sites+(ist.n_sites-i-1)];
        trial_wf_down[j*ist.n_down+i].imag = 0.0;
      }
   }

   /*Save Trial Energy*/
   cns.trial_energy = kinetic_eigs[ist.n_sites-1];

   /*Print Trial Matrices*/
   print_cmat(trial_wf_down,ist.n_sites,ist.n_down,"trial_wf.dat");
   fprintf(pf2, "Trial Energy: %g\n", cns.trial_energy); fflush(pf2);

fclose(pf2);
return;
}

/************************************************************************************/

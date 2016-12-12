#include "afqmc.h" 

/*******************************************************************************************/

int init_states(int *state_matrix,int_st ist) {

    /*Map The Different Possible Determinant States to Their Corresponding Configurations*/
    /*Right Now, Can Only Do a Max of 7 Electrons...Have to Think of How To Do Better*/
    int i=0, j=0, k=0, l=0, m=0, n=0, o=0;
    int count = 0;
    int state_number_1;
    int *config_vector_1;

    config_vector_1 = (int *)calloc(ist.n_electrons, sizeof(int));

    /*Map Numbers According to the Different Electron Numbers*/
    if ( ist.n_electrons == 1 ) { /*************************************************************/

       for (o=0; o<ist.n_spin_orbitals; o++) {
         state_number_1 = state_number(o,n,m,l,k,j,i,ist);

         /*Determine If State Is a Real State*/
         get_configuration(state_number_1,config_vector_1,ist);
         if ( check_electrons(config_vector_1,ist) ) {

             state_matrix[count] = state_number_1;
             count++;

            /*Exit If Count Greater Than Total States*/
            if ( count == ist.n_states ) {
              return(count);
            }
          }

        }
     }
     else if ( ist.n_electrons == 2 ) { /*******************************************************/

       for (o=0; o<ist.n_spin_orbitals; o++) {
        for (n=o+1; n<ist.n_spin_orbitals; n++) {
          state_number_1 = state_number(o,n,m,l,k,j,i,ist);

          /*Determine If State Is a Real State*/
          get_configuration(state_number_1,config_vector_1,ist);
          if ( check_electrons(config_vector_1,ist) ) {

             state_matrix[count] = state_number_1;
             count++;

            /*Exit If Count Greater Than Total States*/
            if ( count == ist.n_states ) {
              return(count);
            }
           }

         }
        }
     }
     else if ( ist.n_electrons == 3 ) { /******************************************************/

       for (o=0; o<ist.n_spin_orbitals; o++) {
        for (n=o+1; n<ist.n_spin_orbitals; n++) {
          for (m=n+1; m<ist.n_spin_orbitals; m++) {

            state_number_1 = state_number(o,n,m,l,k,j,i,ist);

            /*Determine If State Is a Real State*/
            get_configuration(state_number_1,config_vector_1,ist);
            if ( check_electrons(config_vector_1,ist) ) {
              state_matrix[count] = state_number_1;
count++;

              /*Exit If Count Greater Than Total States*/
              if ( count == ist.n_states ) {
                return(count);
              }
             }

          }
        }
       }

    }
    else if ( ist.n_electrons == 4 ) { /****************************************************/

      for (o=0; o<ist.n_spin_orbitals; o++) {
       for (n=o+1; n<ist.n_spin_orbitals; n++) {
        for (m=n+1; m<ist.n_spin_orbitals; m++) {
         for (l=m+1; l<ist.n_spin_orbitals; l++) {

            state_number_1 = state_number(o,n,m,l,k,j,i,ist);

            /*Determine If State Is a Real State*/
            get_configuration(state_number_1,config_vector_1,ist);
            if ( check_electrons(config_vector_1,ist) ) {

               state_matrix[count] = state_number_1;
               count++;

              /*Exit If Count Greater Than Total States*/
              if ( count == ist.n_states ) {
                 return(count);
              }
            }

         }
        }
       }
      }

    }
    else if ( ist.n_electrons == 5 ) { /****************************************************/

      for (o=0; o<ist.n_spin_orbitals; o++) {
       for (n=o+1; n<ist.n_spin_orbitals; n++) {
        for (m=n+1; m<ist.n_spin_orbitals; m++) {
         for (l=m+1; l<ist.n_spin_orbitals; l++) {
          for (k=l+1; k<ist.n_spin_orbitals; k++) {

            state_number_1 = state_number(o,n,m,l,k,j,i,ist);

            /*Determine If State Is a Real State*/
            get_configuration(state_number_1,config_vector_1,ist);
            if ( check_electrons(config_vector_1,ist) ) {

              state_matrix[count] = state_number_1;
              count++;

              /*Exit If Count Greater Than Total States*/
              if ( count == ist.n_states ) {
                return(count);
              }
            }

          }
         }
}
       }
      }

    }
    else if ( ist.n_electrons == 6 ) { /****************************************************/

      for (o=0; o<ist.n_spin_orbitals; o++) {
       for (n=o+1; n<ist.n_spin_orbitals; n++) {
        for (m=n+1; m<ist.n_spin_orbitals; m++) {
         for (l=m+1; l<ist.n_spin_orbitals; l++) {
          for (k=l+1; k<ist.n_spin_orbitals; k++) {
           for (j=k+1; j<ist.n_spin_orbitals; j++) {

             state_number_1 = state_number(o,n,m,l,k,j,i,ist);

             /*Determine If State Is a Real State*/
             get_configuration(state_number_1,config_vector_1,ist);
             if ( check_electrons(config_vector_1,ist) ) {

                state_matrix[count] = state_number_1;
                count++;

              /*Exit If Count Greater Than Total States*/
              if ( count == ist.n_states ) {
                return(count);
              }
            }


           }
          }
         }
        }
       }
      }

    }
    else if ( ist.n_electrons == 7 ){ /****************************************************/

      for (o=0; o<ist.n_spin_orbitals; o++) {
       for (n=o+1; n<ist.n_spin_orbitals; n++) {
        for (m=n+1; m<ist.n_spin_orbitals; m++) {
         for (l=m+1; l<ist.n_spin_orbitals; l++) {
          for (k=l+1; k<ist.n_spin_orbitals; k++) {
           for (j=k+1; j<ist.n_spin_orbitals; j++) {
            for (i=j+1; i<ist.n_spin_orbitals; i++) {

               state_number_1 = state_number(o,n,m,l,k,j,i,ist);

               /*Determine If State Is a Real State*/
               get_configuration(state_number_1,config_vector_1,ist);
               if ( check_electrons(config_vector_1,ist) ) {

                  state_matrix[count] = state_number_1;
                  count++;

                /*Exit If Count Greater Than Total States*/
                if ( count == ist.n_states ) {
                  return(count);
                }
              }

            }
           }
          }
 }
        }
       }
      }

    } /*Last Else If*/

free(config_vector_1);
return(count);
}

/***********************************************************************************************/

int state_number(int number_1,int number_2,int number_3,int number_4,int number_5,int number_6,int number_7,int_st ist) {

     /*Determines the Number of the State Given Its Orbitals*/
     int state_number;

     state_number = number_1 + number_2 * ist.n_spin_orbitals + number_3 * pow(ist.n_spin_orbitals, 2) + number_4 * pow(ist.n_spin_orbitals, 3)
                  + number_5 * pow(ist.n_spin_orbitals,4) + number_6 * pow(ist.n_spin_orbitals,5) + number_7 * pow(ist.n_spin_orbitals,6);

return(state_number);
}

/*******************************************************************************************************/

void get_configuration(int state_number,int *config_vector,int_st ist){

      /*Finds the Excited State Configuration Given the State Number*/
      int remainder, divided, i;

      /*Convert to Configuration*/
      for (i=ist.n_electrons-1; i>-1; i--) {
       divided=(int)(state_number/((double) pow(ist.n_spin_orbitals, i)));
       remainder=(int)(state_number-(int) divided*pow(ist.n_spin_orbitals, i));
       state_number=remainder;
       config_vector[i]=divided;
     }

return;
}

/*******************************************************************************************************/

int check_electrons(int *config_vector_1,int_st ist) {

    /*Determine If Correct Number Up or Down*/
    int i;
    int number_up = 0, number_down = 0;

    for (i=0; i<ist.n_electrons; i++) {
      if ( config_vector_1[i]%2 == 0 ) {
         number_up++;
      }
      else {
         number_down++;
      }
    }

   if ( number_up == ist.n_up && number_down == ist.n_down ) {
     return 1;
   }
   else {
     return 0;
   }

}

/*****************************************************************************************************/

MKL_Complex16 determine_overlap(MKL_Complex16 *hf_determinant,MKL_Complex16 *average_wave_function,int size1,int size2) {

   /*Determines The Overlap Between Two Determinants*/
   int j;
   int *indx;
   MKL_Complex16 *stored_product;
   MKL_Complex16 determinant;

   stored_product = (MKL_Complex16 *)calloc(size2*size2,sizeof(MKL_Complex16));
   indx = (int*)calloc(size2,sizeof(int));

   transpose_cmat_cmat(hf_determinant,average_wave_function,stored_product,size1,size2,size2);
   complex_ludmp(stored_product,size2,indx,&determinant);

   for (j=0; j<size2; j++) {
     determinant = Cmul(determinant, stored_product[j*size2+j]);
   }

free(indx);
free(stored_product);
return(determinant);
}

/*****************************************************************************************************/

void get_fci_expansion(MKL_Complex16 *average_wave_function_up,MKL_Complex16 *average_wave_function_down,int *state_matrix,int_st ist) {

   /*Determines and Prints the Expansion Coefficients that Come from the Overlap of the Average Wave Functions and Various HF+Excited State Determinants*/
   /*Can Be Used to Check Against FCI Code*/

   int i, j; 
   int j_orbital, j_spin, state_number;
   int number_up, number_down;
   MKL_Complex16 *hf_determinant_up, *hf_determinant_down;
   MKL_Complex16 *final_determinants;
   double normalization;
   int *config_vector_1;
   MKL_Complex16 determinant_up, determinant_down;
   FILE *pf = fopen("fci_coefficients.dat", "a+");

   hf_determinant_up = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_up,sizeof(MKL_Complex16));
   hf_determinant_down = (MKL_Complex16 *)calloc(ist.n_spatial_orbitals_n_down,sizeof(MKL_Complex16));
   final_determinants = (MKL_Complex16 *)calloc(ist.n_states,sizeof(MKL_Complex16));
   config_vector_1 = (int *)calloc(ist.n_electrons,sizeof(int));

   /*If Both Up and Down Electrons*/
   if ( ist.n_up > 0 && ist.n_down > 0 ) {

     for (i=0; i<ist.n_states; i++) {

       state_number = state_matrix[i];
       get_configuration(state_number,config_vector_1,ist);

       czero_vec(hf_determinant_up,ist.n_spatial_orbitals_n_up);
       czero_vec(hf_determinant_down,ist.n_spatial_orbitals_n_down);

       number_up = number_down = 0;
       for (j=0; j<ist.n_electrons; j++) {

         j_orbital = (int)(config_vector_1[j]/2.0);
         j_spin = config_vector_1[j]%2;
         if ( j_spin == 0 ) {
           hf_determinant_up[j_orbital*ist.n_up+number_up].real = 1.0 ;
           number_up++;
         }
         else {
           hf_determinant_down[j_orbital*ist.n_down+number_down].real = 1.0;
           number_down++;
          }
       }

       /*Now Find Determinants to Get Overlaps*/
       determinant_up = determine_overlap(hf_determinant_up,average_wave_function_up,ist.n_spatial_orbitals,ist.n_up);
       determinant_down = determine_overlap(hf_determinant_down,average_wave_function_down,ist.n_spatial_orbitals,ist.n_down);

      final_determinants[i] = Cmul(determinant_up, determinant_down);
     }

     /*Normalize the Final Components*/
    normalization = 0.0;
     for (i=0; i<ist.n_states; i++) {
       normalization += final_determinants[i].real * final_determinants[i].real; // - final_determinants[i].i * final_determinants[i].i;
     }
     normalization = 1.0/sqrt(normalization);

     for (i=0; i<ist.n_states; i++) {
       fprintf(pf, "%f+%fi\t", final_determinants[i].real*normalization, final_determinants[i].imag*normalization);
     }
     fprintf(pf, "\n"); fflush(pf);

   } /*If Both Occupied*/
   else if ( ist.n_up > 0 && ist.n_down == 0 ) {
     for (i=0; i<ist.n_states; i++) {

       state_number = state_matrix[i];
       get_configuration(state_number,config_vector_1,ist);
czero_vec(hf_determinant_up,ist.n_spatial_orbitals_n_up);

       number_up = 0;
       for (j=0; j<ist.n_electrons; j++) {

         j_orbital = (int)(config_vector_1[j]/2.0);
         hf_determinant_up[j_orbital*ist.n_up+number_up].real = 1.0 ;
         number_up++;
       }

       /*Now Find Determinants to Get Overlaps*/
       determinant_up = determine_overlap(hf_determinant_up,average_wave_function_up,ist.n_spatial_orbitals,ist.n_up);
       final_determinants[i] = determinant_up;

     }

     /*Normalize the Final Components*/
     normalization = 0.0;
     for (i=0; i<ist.n_states; i++) {
       normalization += final_determinants[i].real * final_determinants[i].real - final_determinants[i].imag * final_determinants[i].imag;
     }
     normalization = 1.0/sqrt(normalization);

     for (i=0; i<ist.n_states; i++) {
       fprintf(pf, "%f+%fi\t", final_determinants[i].real*normalization, final_determinants[i].imag*normalization);
     }
     fprintf(pf, "\n"); fflush(pf);

  }/*If Only Up*/
  else if ( ist.n_up == 0 && ist.n_down > 0 ) {

     for (i=0; i<ist.n_states; i++) {

       state_number = state_matrix[i];
       get_configuration(state_number,config_vector_1,ist);

       czero_vec(hf_determinant_down,ist.n_spatial_orbitals_n_down);

       number_down = 0;
       for (j=0; j<ist.n_electrons; j++) {
         j_orbital = (int)(config_vector_1[j]/2.0);
         hf_determinant_down[j_orbital*ist.n_down+number_down].real = 1.0;
         number_down++;
       }


       /*Now Find Determinants to Get Overlaps*/
       determinant_down = determine_overlap(hf_determinant_down,average_wave_function_down,ist.n_spatial_orbitals,ist.n_down);
       final_determinants[i] = determinant_down;
     }

     /*Normalize the Final Components*/
     normalization = 0.0;
     for (i=0; i<ist.n_states; i++) {
       normalization += final_determinants[i].real * final_determinants[i].real - final_determinants[i].imag * final_determinants[i].imag;
     }
     normalization = 1.0/sqrt(normalization);

     for (i=0; i<ist.n_states; i++) {
       fprintf(pf, "%f+%fi\t", final_determinants[i].real*normalization, final_determinants[i].imag*normalization);
     }
     fprintf(pf, "\n"); fflush(pf);

  }

free(config_vector_1);
free(hf_determinant_up);
free(hf_determinant_down);
free(final_determinants);
fclose(pf);
return;
}



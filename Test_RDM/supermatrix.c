#include "afqmc.h"

/**********************************************************************************/

void get_potential_matrix_elements(MKL_Complex16 *kinetic_matrix_original_fermions,MKL_Complex16 *potential_matrix_original_fermions,int_st ist,char *str,char *str2) {

      /*Scans File Containing Vijkl Elements and Stores Them in Memory*/
      int i; 
      FILE *pf = fopen(str, "r");
      FILE *pf2 = fopen(str2, "r");  

      for (i=0; i<ist.n_sites_fourth; i++) {
         fscanf(pf, "%lf %lf", &potential_matrix_original_fermions[i].real, &potential_matrix_original_fermions[i].imag); 
         if (i < ist.n_sites_sq) {
           fscanf(pf2, "%lf %lf", &kinetic_matrix_original_fermions[i].real, &kinetic_matrix_original_fermions[i].imag); 
         }
      }

fclose(pf);
fclose(pf2);    
return; 
}

/**********************************************************************************/

void get_potential_matrix_elements_molpro_restricted(MKL_Complex16 *kinetic_matrix_original_fermions,MKL_Complex16 *potential_matrix_original_fermions,MKL_Complex16 *potential_matrix_original_fermions_spin,MKL_Complex16 *energy_shifts,int *list_occupied_orbitals,int_st ist,cns_st cns,char *str) {

    /*Uses MolPro File with Vijkl Entries To Build Single Particle Kinetic Matrix and Two-Particle Potential Matrix*/
    int i, j, k, l; 
    int i2, j2, k2, l2; 
    int i2_even, j2_even, k2_even, l2_even; 
    int i2_odd, j2_odd, k2_odd, l2_odd; 
    int ij, in, jn, kn; 
    double matrix_element; 
    FILE *pf = fopen(str, "r"); 
    FILE *pf2 = fopen("checkmatrixelements.dat", "w+");   
    FILE *pf3 = fopen("errors.dat", "a+"); 

    czero_vec(potential_matrix_original_fermions,ist.n_spatial_orbitals_fourth);
    czero_vec(potential_matrix_original_fermions_spin,ist.n_spin_orbitals_fourth); 
    czero_vec(kinetic_matrix_original_fermions,ist.n_spatial_orbitals_sq);   
    czero_vec(energy_shifts,4); 

    /*Chemists' Notation Here!*/ 
    while ( fscanf(pf, "%lf %d %d %d %d", &matrix_element, &i, &j, &k, &l) != EOF ){
       i2 = i-1; j2 = j-1; k2 = k-1; l2 = l-1;  
 
       if ( k!=0 && l!=0 ) { /*Get Two-Particle Elements*/ 

         i2_even = i2*2; j2_even = j2*2; k2_even = k2*2; l2_even = l2*2;
         i2_odd = i2_even+1; j2_odd = j2_even+1; k2_odd = k2_even+1; l2_odd = l2_even+1;

         /*ijkl*/
         potential_matrix_original_fermions[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2].real = matrix_element; 
 

         /*Larger Matrix*/
         potential_matrix_original_fermions_spin[i2_even*ist.n_spin_orbitals_third+k2_even*ist.n_spin_orbitals_sq+j2_even*ist.n_spin_orbitals+l2_even].real = matrix_element; 
         potential_matrix_original_fermions_spin[i2_odd*ist.n_spin_orbitals_third+k2_odd*ist.n_spin_orbitals_sq+j2_odd*ist.n_spin_orbitals+l2_odd].real = matrix_element;   
         potential_matrix_original_fermions_spin[i2_even*ist.n_spin_orbitals_third+k2_odd*ist.n_spin_orbitals_sq+j2_even*ist.n_spin_orbitals+l2_odd].real = matrix_element; 
         potential_matrix_original_fermions_spin[i2_odd*ist.n_spin_orbitals_third+k2_even*ist.n_spin_orbitals_sq+j2_odd*ist.n_spin_orbitals+l2_even].real = matrix_element; 

         /*Should Add Part That Makes Assignments Faster If All Indices The Same, Etc.*/ 
         /*Get All Other Possible Matrix Elements Assuming Real*/
         /*jilk*/
         potential_matrix_original_fermions[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2].real = matrix_element;   
         
         /*Larger Matrix Elements*/
         potential_matrix_original_fermions_spin[j2_even*ist.n_spin_orbitals_third+l2_even*ist.n_spin_orbitals_sq+i2_even*ist.n_spin_orbitals+k2_even].real = matrix_element; 
         potential_matrix_original_fermions_spin[j2_odd*ist.n_spin_orbitals_third+l2_odd*ist.n_spin_orbitals_sq+i2_odd*ist.n_spin_orbitals+k2_odd].real = matrix_element;
         potential_matrix_original_fermions_spin[j2_even*ist.n_spin_orbitals_third + l2_odd*ist.n_spin_orbitals_sq+i2_even*ist.n_spin_orbitals+k2_odd].real = matrix_element; 
         potential_matrix_original_fermions_spin[j2_odd*ist.n_spin_orbitals_third+l2_even*ist.n_spin_orbitals_sq+i2_odd*ist.n_spin_orbitals+k2_even].real = matrix_element; 


         /*klij*/
         potential_matrix_original_fermions[k2*ist.n_spatial_orbitals_third+l2*ist.n_spatial_orbitals_sq+i2*ist.n_spatial_orbitals+j2].real = matrix_element;

         /*Larger Spin Matrix*/
         potential_matrix_original_fermions_spin[k2_even*ist.n_spin_orbitals_third+i2_even*ist.n_spin_orbitals_sq+l2_even*ist.n_spin_orbitals+j2_even].real = matrix_element; 
         potential_matrix_original_fermions_spin[k2_odd*ist.n_spin_orbitals_third+i2_odd*ist.n_spin_orbitals_sq+l2_odd*ist.n_spin_orbitals+j2_odd].real = matrix_element;
         potential_matrix_original_fermions_spin[k2_even*ist.n_spin_orbitals_third+i2_odd*ist.n_spin_orbitals_sq+l2_even*ist.n_spin_orbitals+j2_odd].real = matrix_element; 
         potential_matrix_original_fermions_spin[k2_odd*ist.n_spin_orbitals_third+i2_even*ist.n_spin_orbitals_sq+l2_odd*ist.n_spin_orbitals+j2_even].real = matrix_element; 


         /*ijlk*/
         potential_matrix_original_fermions[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2].real = matrix_element;

         /*Larger*/
         potential_matrix_original_fermions_spin[i2_even*ist.n_spin_orbitals_third+l2_even*ist.n_spin_orbitals_sq+j2_even*ist.n_spin_orbitals+k2_even].real = matrix_element; 
         potential_matrix_original_fermions_spin[i2_odd*ist.n_spin_orbitals_third+l2_odd*ist.n_spin_orbitals_sq+j2_odd*ist.n_spin_orbitals+k2_odd].real = matrix_element;    
         potential_matrix_original_fermions_spin[i2_even*ist.n_spin_orbitals_third+l2_odd*ist.n_spin_orbitals_sq+j2_even*ist.n_spin_orbitals+k2_odd].real = matrix_element; 
         potential_matrix_original_fermions_spin[i2_odd*ist.n_spin_orbitals_third+l2_even*ist.n_spin_orbitals_sq+j2_odd*ist.n_spin_orbitals+k2_even].real = matrix_element; 


         /*jikl*/
         potential_matrix_original_fermions[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2].real = matrix_element; 

         /*Larger*/
         potential_matrix_original_fermions_spin[j2_even*ist.n_spin_orbitals_third+k2_even*ist.n_spin_orbitals_sq+i2_even*ist.n_spin_orbitals+l2_even].real = matrix_element; 
         potential_matrix_original_fermions_spin[j2_odd*ist.n_spin_orbitals_third+k2_odd*ist.n_spin_orbitals_sq+i2_odd*ist.n_spin_orbitals+l2_odd].real = matrix_element;
         potential_matrix_original_fermions_spin[j2_even*ist.n_spin_orbitals_third+k2_odd*ist.n_spin_orbitals_sq+i2_even*ist.n_spin_orbitals+l2_odd].real = matrix_element; 
         potential_matrix_original_fermions_spin[j2_odd*ist.n_spin_orbitals_third+k2_even*ist.n_spin_orbitals_sq+i2_odd*ist.n_spin_orbitals+l2_even].real = matrix_element;       

         /*lkij*/
         potential_matrix_original_fermions[l2*ist.n_spatial_orbitals_third+k2*ist.n_spatial_orbitals_sq+i2*ist.n_spatial_orbitals+j2].real = matrix_element;

         /*Larger*/
         potential_matrix_original_fermions_spin[l2_even*ist.n_spin_orbitals_third+i2_even*ist.n_spin_orbitals_sq+k2_even*ist.n_spin_orbitals+j2_even].real = matrix_element; 
         potential_matrix_original_fermions_spin[l2_odd*ist.n_spin_orbitals_third+i2_odd*ist.n_spin_orbitals_sq+k2_odd*ist.n_spin_orbitals+j2_odd].real = matrix_element; 
         potential_matrix_original_fermions_spin[l2_even*ist.n_spin_orbitals_third+i2_odd*ist.n_spin_orbitals_sq+k2_even*ist.n_spin_orbitals+j2_odd].real = matrix_element; 
         potential_matrix_original_fermions_spin[l2_odd*ist.n_spin_orbitals_third+i2_even*ist.n_spin_orbitals_sq+k2_odd*ist.n_spin_orbitals+j2_even].real = matrix_element; 

         /*klji*/
         potential_matrix_original_fermions[k2*ist.n_spatial_orbitals_third+l2*ist.n_spatial_orbitals_sq+j2*ist.n_spatial_orbitals+i2].real = matrix_element; 

         /*Larger*/
         potential_matrix_original_fermions_spin[k2_even*ist.n_spin_orbitals_third+j2_even*ist.n_spin_orbitals_sq+l2_even*ist.n_spin_orbitals+i2_even].real = matrix_element; 
         potential_matrix_original_fermions_spin[k2_odd*ist.n_spin_orbitals_third+j2_odd*ist.n_spin_orbitals_sq+l2_odd*ist.n_spin_orbitals+i2_odd].real = matrix_element; 
         potential_matrix_original_fermions_spin[k2_even*ist.n_spin_orbitals_third+j2_odd*ist.n_spin_orbitals_sq+l2_even*ist.n_spin_orbitals+i2_odd].real = matrix_element; 
         potential_matrix_original_fermions_spin[k2_odd*ist.n_spin_orbitals_third+j2_even*ist.n_spin_orbitals_sq+l2_odd*ist.n_spin_orbitals+i2_even].real = matrix_element; 


         /*lkji*/
         potential_matrix_original_fermions[l2*ist.n_spatial_orbitals_third+k2*ist.n_spatial_orbitals_sq+j2*ist.n_spatial_orbitals+i2].real = matrix_element;  

         /*Larger*/
         potential_matrix_original_fermions_spin[l2_even*ist.n_spin_orbitals_third+j2_even*ist.n_spin_orbitals_sq+k2_even*ist.n_spin_orbitals+i2_even].real = matrix_element; 
         potential_matrix_original_fermions_spin[l2_odd*ist.n_spin_orbitals_third+j2_odd*ist.n_spin_orbitals_sq+k2_odd*ist.n_spin_orbitals+i2_odd].real = matrix_element;
         potential_matrix_original_fermions_spin[l2_even*ist.n_spin_orbitals_third+j2_odd*ist.n_spin_orbitals_sq+k2_even*ist.n_spin_orbitals+i2_odd].real = matrix_element; 
         potential_matrix_original_fermions_spin[l2_odd*ist.n_spin_orbitals_third+j2_even*ist.n_spin_orbitals_sq+k2_odd*ist.n_spin_orbitals+i2_even].real = matrix_element;  

       }
       else {
         /*Get One-Particle Elements*/ 
         if ( i!=0 && j!=0 ) { /*Ensure That This Is Not Just a Constant Shift*/
            kinetic_matrix_original_fermions[i2*ist.n_spatial_orbitals+j2].real = matrix_element;
            kinetic_matrix_original_fermions[j2*ist.n_spatial_orbitals+i2].real = matrix_element;  
         }
         else { /*Get Energy Shift*/
            energy_shifts[0].real = matrix_element;  
         }
       }
    }

    /*Print Out Checks*/
    for (i=0; i<ist.n_spatial_orbitals; i++) {
     in = i * ist.n_spatial_orbitals; 

     for (j=0; j<ist.n_spatial_orbitals; j++) {
      ij = in + j;
   
      fprintf(pf2, "%f+%fi\t", kinetic_matrix_original_fermions[ij].real, kinetic_matrix_original_fermions[ij].imag); 
     }
     fprintf(pf2, "\n"); 
    }
    fprintf(pf2, "\n\n"); fflush(pf2); 

    /*Print in Physics' Notation*/
    for (i=0; i<ist.n_spatial_orbitals; i++) {
     in = i * ist.n_spatial_orbitals_third; 

     for (k=0; k<ist.n_spatial_orbitals; k++) {
      kn = k * ist.n_spatial_orbitals; 

      for (j=0; j<ist.n_spatial_orbitals; j++) {
       jn = j * ist.n_spatial_orbitals_sq;   

       for (l=0; l<ist.n_spatial_orbitals; l++) {
        fprintf(pf2, "%d\t %d\t %d\t %d\t %f+%fi\n", i, j, k, l, potential_matrix_original_fermions[in+kn+jn+l].real, potential_matrix_original_fermions[in+kn+jn+l].imag); 
       }
      }
    }
   }
  fprintf(pf2, "\n\n"); 
  fprintf(pf2, "energy shift %f\n\n", energy_shifts[0].real); fflush(pf2); 
  

fclose(pf); 
fclose(pf2); 
return; 
}

/**********************************************************************************************************************************/

void get_number_matrix_elements_molpro_unrestricted(int *number_kinetic_sparse_up, int *number_kinetic_sparse_down, int *number_potential_sparse_up, int *number_potential_sparse_down, int *number_potential_sparse_updown, char *str) {

    /*Determines Which Mol Pro Elements are Non-Zero for Calculations*/ 
    int count = 0;  
    int i, j, k, l; 
    int number_kinetic_sparse_up_2 = 0, number_kinetic_sparse_down_2 = 0; 
    int number_potential_sparse_up_2 = 0, number_potential_sparse_down_2 = 0, number_potential_sparse_updown_2 = 0; 
    double matrix_element; 
    FILE *pf = fopen(str, "r"); 

    /*Scan Through Matrix Elements*/
    while ( fscanf(pf, "%lf %d %d %d %d", &matrix_element, &i, &j, &k, &l) != EOF ){

       /*If Matrix Element Sizeable*/
       if ( fabs(matrix_element) > .0000000000001 ) { 

        if ( k!=0 && l!=0 ) { /*Get Two-Particle Elements*/

           /*Count 8 Times for Up and Down and 4 Times for Updown*/
           if ( count == 0 ) { 
            number_potential_sparse_up_2+=8; 
           }
           else if ( count == 1 ) {
            number_potential_sparse_down_2+=8; 
           }
           else {
            number_potential_sparse_updown_2+=8; 
           }

        } /*k and l nonzero*/
        else { 

          /*Get One-Particle Elements*/
          if ( i!=0 && j!=0 ) { /*Ensure That This Is Not Just a Constant Shift*/

            /*Count Twice Because Need IJ and JI*/
            if ( count == 3 ) {
              if ( i!= j ) {
               number_kinetic_sparse_up_2+=2; 
              }
              else {
               number_kinetic_sparse_up_2+=1; 
              } 
            }
            else if ( count == 4 ) {
              if ( i!= j ) { 
                number_kinetic_sparse_down_2+=2; 
              }
              else {
                number_kinetic_sparse_down_2+=1; 
              }
            }

         }

       }/*If k and l zero*/
     } /*If Matrix Element NonZero*/   
     if ( i == 0 && j == 0 && k == 0 && l == 0 ) {
           count++;
     }


   } /*While Loop*/

   /*Record All Values*/
   *number_kinetic_sparse_up = number_kinetic_sparse_up_2; 
   *number_kinetic_sparse_down = number_kinetic_sparse_down_2; 
   *number_potential_sparse_up = number_potential_sparse_up_2; 
   *number_potential_sparse_down = number_potential_sparse_down_2; 
   *number_potential_sparse_updown = number_potential_sparse_updown_2;   

fclose(pf); 
}

/**********************************************************************************************************************************/

void get_potential_matrix_elements_molpro_unrestricted(MKL_Complex16 *kinetic_matrix_original_fermions_up,MKL_Complex16 *kinetic_matrix_original_fermions_down,MKL_Complex16 *potential_matrix_original_fermions_up,MKL_Complex16 *potential_matrix_original_fermions_down, MKL_Complex16 *potential_matrix_original_fermions_updown,MKL_Complex16 *potential_matrix_original_fermions_spin,double *kinetic_matrix_sparse_up, double *kinetic_matrix_sparse_down,double *potential_matrix_sparse_up,double *potential_matrix_sparse_down,double *potential_matrix_sparse_updown,int *kinetic_ij_sparse_up, int *kinetic_ij_sparse_down, int *potential_ij_sparse_up, int *potential_ij_sparse_down, int *potential_ij_sparse_updown, int *list_orbitals_occupied, MKL_Complex16 *energy_shifts,int *n_kinetic_sparse_up, int *n_kinetic_sparse_down, int *n_potential_sparse_up,int *n_potential_sparse_down, int *n_potential_sparse_updown, int_st ist,cns_st cns,char *str) {

    /*Uses MolPro File with Vijkl Entries To Build Single Particle Kinetic Matrix and Two-Particle Potential Matrix*/
    /*This Is If the System Is Unrestricted and Can Use Different Orbitals for Each Spin*/
    int i, j, k, l, m;
    int check_j, check_l;  
    int count = 0;
    int index_kinetic_sparse_up = 0, index_kinetic_sparse_down = 0;   
    int index_potential_sparse_up = 0, index_potential_sparse_down = 0, index_potential_sparse_updown = 0; 
    int f_index_up, f_index_down, f_index_updown;  
    int i2, j2, k2, l2;
    int i2_even, j2_even, k2_even, l2_even; 
    int i2_odd, j2_odd, k2_odd, l2_odd; 
    int ij, in, jn, kn;
    double matrix_element;
    FILE *pf = fopen(str, "r");
    FILE *pf2 = fopen("checkmatrixelements.dat", "w+");
    double *potential_matrix_original_fermions_up_2, *potential_matrix_original_fermions_down_2, *potential_matrix_original_fermions_updown_2;
    double *kinetic_matrix_original_fermions_up_2, *kinetic_matrix_original_fermions_down_2;  
     
    potential_matrix_original_fermions_up_2 = (double *)calloc(ist.n_spatial_orbitals_fourth,sizeof(double)); 
    potential_matrix_original_fermions_down_2 = (double *)calloc(ist.n_spatial_orbitals_fourth,sizeof(double));
    potential_matrix_original_fermions_updown_2 = (double *)calloc(ist.n_spatial_orbitals_fourth,sizeof(double)); 

    kinetic_matrix_original_fermions_up_2 = (double *)calloc(ist.n_spatial_orbitals_sq,sizeof(double)); 
    kinetic_matrix_original_fermions_down_2 = (double *)calloc(ist.n_spatial_orbitals_sq,sizeof(double)); 

  /*HERE*/
    czero_vec(potential_matrix_original_fermions_up,ist.n_spatial_orbitals_fourth);
    czero_vec(potential_matrix_original_fermions_down,ist.n_spatial_orbitals_fourth); 
    czero_vec(potential_matrix_original_fermions_updown,ist.n_spatial_orbitals_fourth);
    czero_vec(potential_matrix_original_fermions_spin,ist.n_spin_orbitals_fourth);  
    czero_vec(kinetic_matrix_original_fermions_up,ist.n_spatial_orbitals_sq);
    czero_vec(kinetic_matrix_original_fermions_down,ist.n_spatial_orbitals_sq);
    czero_vec(energy_shifts,4); 
 
    while ( fscanf(pf, "%lf %d %d %d %d", &matrix_element, &i, &j, &k, &l) != EOF ){
       i2 = i-1; j2 = j-1; k2 = k-1; l2 = l-1;

       if ( k!=0 && l!=0 ) { /*Get Two-Particle Elements*/
         i2_even = i2*2; j2_even = j2*2; k2_even = k2*2; l2_even = l2*2;
         i2_odd = i2_even+1; j2_odd = j2_even+1; k2_odd = k2_even+1; l2_odd = l2_even+1;

 
         if ( count == 0 ) { /*If All Spin Up*/  

           /*For Now I Assume That All Matrix Elements Are REal*/
           /*ijkl*/
           potential_matrix_original_fermions_up[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2].real = matrix_element;
           potential_matrix_original_fermions_spin[i2_even*ist.n_spin_orbitals_third+k2_even*ist.n_spin_orbitals_sq+j2_even*ist.n_spin_orbitals+l2_even].real = matrix_element;
           potential_matrix_original_fermions_up_2[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2] = matrix_element; 

           /*jilk*/
           potential_matrix_original_fermions_up[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2].real = matrix_element;
           potential_matrix_original_fermions_spin[j2_even*ist.n_spin_orbitals_third+l2_even*ist.n_spin_orbitals_sq+i2_even*ist.n_spin_orbitals+k2_even].real = matrix_element; 
           potential_matrix_original_fermions_up_2[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2] = matrix_element;

           /*klij*/
           potential_matrix_original_fermions_up[k2*ist.n_spatial_orbitals_third+l2*ist.n_spatial_orbitals_sq+i2*ist.n_spatial_orbitals+j2].real = matrix_element;
           potential_matrix_original_fermions_spin[k2_even*ist.n_spin_orbitals_third+i2_even*ist.n_spin_orbitals_sq+l2_even*ist.n_spin_orbitals+j2_even].real = matrix_element;
           potential_matrix_original_fermions_up_2[k2*ist.n_spatial_orbitals_third+l2*ist.n_spatial_orbitals_sq+i2*ist.n_spatial_orbitals+j2] = matrix_element;

           /*ijlk*/
           potential_matrix_original_fermions_up[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2].real = matrix_element;
           potential_matrix_original_fermions_spin[i2_even*ist.n_spin_orbitals_third+l2_even*ist.n_spin_orbitals_sq+j2_even*ist.n_spin_orbitals+k2_even].real = matrix_element;
           potential_matrix_original_fermions_up_2[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2] = matrix_element;

           /*jikl*/
           potential_matrix_original_fermions_up[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2].real = matrix_element;
           potential_matrix_original_fermions_spin[j2_even*ist.n_spin_orbitals_third+k2_even*ist.n_spin_orbitals_sq+i2_even*ist.n_spin_orbitals+l2_even].real = matrix_element;
           potential_matrix_original_fermions_up_2[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2] = matrix_element;

           /*lkij*/
           potential_matrix_original_fermions_up[l2*ist.n_spatial_orbitals_third+k2*ist.n_spatial_orbitals_sq+i2*ist.n_spatial_orbitals+j2].real = matrix_element;
           potential_matrix_original_fermions_spin[l2_even*ist.n_spin_orbitals_third+i2_even*ist.n_spin_orbitals_sq+k2_even*ist.n_spin_orbitals+j2_even].real = matrix_element;
           potential_matrix_original_fermions_up_2[l2*ist.n_spatial_orbitals_third+k2*ist.n_spatial_orbitals_sq+i2*ist.n_spatial_orbitals+j2] = matrix_element; 

           /*klji*/
           potential_matrix_original_fermions_up[k2*ist.n_spatial_orbitals_third+l2*ist.n_spatial_orbitals_sq+j2*ist.n_spatial_orbitals+i2].real = matrix_element;
           potential_matrix_original_fermions_spin[k2_even*ist.n_spin_orbitals_third+j2_even*ist.n_spin_orbitals_sq+l2_even*ist.n_spin_orbitals+i2_even].real = matrix_element; 
           potential_matrix_original_fermions_up_2[k2*ist.n_spatial_orbitals_third+l2*ist.n_spatial_orbitals_sq+j2*ist.n_spatial_orbitals+i2] = matrix_element;

           /*lkji*/
           potential_matrix_original_fermions_up[l2*ist.n_spatial_orbitals_third+k2*ist.n_spatial_orbitals_sq+j2*ist.n_spatial_orbitals+i2].real = matrix_element;
           potential_matrix_original_fermions_spin[l2_even*ist.n_spin_orbitals_third+j2_even*ist.n_spin_orbitals_sq+k2_even*ist.n_spin_orbitals+i2_even].real = matrix_element;
           potential_matrix_original_fermions_up_2[l2*ist.n_spatial_orbitals_third+k2*ist.n_spatial_orbitals_sq+j2*ist.n_spatial_orbitals+i2] = matrix_element;

         }
         else if ( count == 1 ) { /*If All Spin Down*/

           /*ijkl*/ 
           potential_matrix_original_fermions_down[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2].real = matrix_element;
           potential_matrix_original_fermions_spin[i2_odd*ist.n_spin_orbitals_third+k2_odd*ist.n_spin_orbitals_sq+j2_odd*ist.n_spin_orbitals+l2_odd].real = matrix_element;
           potential_matrix_original_fermions_down_2[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2] = matrix_element;


           /*jilk*/
           potential_matrix_original_fermions_down[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2].real = matrix_element;
           potential_matrix_original_fermions_spin[j2_odd*ist.n_spin_orbitals_third+l2_odd*ist.n_spin_orbitals_sq+i2_odd*ist.n_spin_orbitals+k2_odd].real = matrix_element;  
           potential_matrix_original_fermions_down_2[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2] = matrix_element; 

           /*klij*/
           potential_matrix_original_fermions_down[k2*ist.n_spatial_orbitals_third+l2*ist.n_spatial_orbitals_sq+i2*ist.n_spatial_orbitals+j2].real = matrix_element;
           potential_matrix_original_fermions_spin[k2_odd*ist.n_spin_orbitals_third+i2_odd*ist.n_spin_orbitals_sq+l2_odd*ist.n_spin_orbitals+j2_odd].real = matrix_element; 
           potential_matrix_original_fermions_down_2[k2*ist.n_spatial_orbitals_third+l2*ist.n_spatial_orbitals_sq+i2*ist.n_spatial_orbitals+j2] = matrix_element;

           /*ijlk*/
           potential_matrix_original_fermions_down[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2].real = matrix_element;
           potential_matrix_original_fermions_spin[i2_odd*ist.n_spin_orbitals_third+l2_odd*ist.n_spin_orbitals_sq+j2_odd*ist.n_spin_orbitals+k2_odd].real = matrix_element; 
           potential_matrix_original_fermions_down_2[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2] = matrix_element;

           /*jikl*/ 
           potential_matrix_original_fermions_down[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2].real = matrix_element;
           potential_matrix_original_fermions_spin[j2_odd*ist.n_spin_orbitals_third+k2_odd*ist.n_spin_orbitals_sq+i2_odd*ist.n_spin_orbitals+l2_odd].real = matrix_element; 
           potential_matrix_original_fermions_down_2[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2] = matrix_element;

           /*lkij*/
           potential_matrix_original_fermions_down[l2*ist.n_spatial_orbitals_third+k2*ist.n_spatial_orbitals_sq+i2*ist.n_spatial_orbitals+j2].real = matrix_element;
           potential_matrix_original_fermions_spin[l2_odd*ist.n_spin_orbitals_third+i2_odd*ist.n_spin_orbitals_sq+k2_odd*ist.n_spin_orbitals+j2_odd].real = matrix_element;  
           potential_matrix_original_fermions_down_2[l2*ist.n_spatial_orbitals_third+k2*ist.n_spatial_orbitals_sq+i2*ist.n_spatial_orbitals+j2] = matrix_element;

           /*klji*/
           potential_matrix_original_fermions_down[k2*ist.n_spatial_orbitals_third+l2*ist.n_spatial_orbitals_sq+j2*ist.n_spatial_orbitals+i2].real = matrix_element;
           potential_matrix_original_fermions_spin[k2_odd*ist.n_spin_orbitals_third+j2_odd*ist.n_spin_orbitals_sq+l2_odd*ist.n_spin_orbitals+i2_odd].real = matrix_element; 
           potential_matrix_original_fermions_down_2[k2*ist.n_spatial_orbitals_third+l2*ist.n_spatial_orbitals_sq+j2*ist.n_spatial_orbitals+i2] = matrix_element; 

           /*lkji*/
           potential_matrix_original_fermions_down[l2*ist.n_spatial_orbitals_third+k2*ist.n_spatial_orbitals_sq+j2*ist.n_spatial_orbitals+i2].real = matrix_element;
           potential_matrix_original_fermions_spin[l2_odd*ist.n_spin_orbitals_third+j2_odd*ist.n_spin_orbitals_sq+k2_odd*ist.n_spin_orbitals+i2_odd].real = matrix_element; 
           potential_matrix_original_fermions_down_2[l2*ist.n_spatial_orbitals_third+k2*ist.n_spatial_orbitals_sq+j2*ist.n_spatial_orbitals+i2] = matrix_element;

         }
         else { /*If Spin Up and Down*/

           /*ijkl*/
           potential_matrix_original_fermions_updown[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2].real = matrix_element;
           potential_matrix_original_fermions_updown_2[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2] = matrix_element;

           /*Larger*/
           potential_matrix_original_fermions_spin[i2_odd*ist.n_spin_orbitals_third+k2_even*ist.n_spin_orbitals_sq+j2_odd*ist.n_spin_orbitals+l2_even].real = matrix_element;
           potential_matrix_original_fermions_spin[i2_even*ist.n_spin_orbitals_third+k2_odd*ist.n_spin_orbitals_sq+j2_even*ist.n_spin_orbitals+l2_odd].real = matrix_element;  


           /*jilk*/ 
           potential_matrix_original_fermions_updown[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2].real = matrix_element;
           potential_matrix_original_fermions_updown_2[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2] = matrix_element;

           potential_matrix_original_fermions_spin[j2_odd*ist.n_spin_orbitals_third+l2_even*ist.n_spin_orbitals_sq+i2_odd*ist.n_spin_orbitals+k2_even].real = matrix_element;  
           potential_matrix_original_fermions_spin[j2_even*ist.n_spin_orbitals_third+l2_odd*ist.n_spin_orbitals_sq+i2_even*ist.n_spin_orbitals+k2_odd].real = matrix_element;       

           /*ijlk*/
           potential_matrix_original_fermions_updown[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2].real = matrix_element;
           potential_matrix_original_fermions_updown_2[i2*ist.n_spatial_orbitals_third+j2*ist.n_spatial_orbitals_sq+l2*ist.n_spatial_orbitals+k2] = matrix_element;

           potential_matrix_original_fermions_spin[i2_odd*ist.n_spin_orbitals_third+l2_even*ist.n_spin_orbitals_sq+j2_odd*ist.n_spin_orbitals+k2_even].real = matrix_element; 
           potential_matrix_original_fermions_spin[i2_even*ist.n_spin_orbitals_third+l2_odd*ist.n_spin_orbitals_sq+j2_even*ist.n_spin_orbitals+k2_odd].real = matrix_element; 

           /*jikl*/
           potential_matrix_original_fermions_updown[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2].real = matrix_element;
           potential_matrix_original_fermions_updown_2[j2*ist.n_spatial_orbitals_third+i2*ist.n_spatial_orbitals_sq+k2*ist.n_spatial_orbitals+l2] = matrix_element;

           potential_matrix_original_fermions_spin[j2_odd*ist.n_spin_orbitals_third+k2_even*ist.n_spin_orbitals_sq+i2_odd*ist.n_spin_orbitals+l2_even].real = matrix_element; 
           potential_matrix_original_fermions_spin[j2_even*ist.n_spin_orbitals_third+k2_odd*ist.n_spin_orbitals_sq+i2_even*ist.n_spin_orbitals+l2_odd].real = matrix_element; 

         } /*Up Down Potential*/
       }
       else {
         /*Get One-Particle Elements*/
         if ( i!=0 && j!=0 ) { /*Ensure That This Is Not Just a Constant Shift*/

           if ( count == 3 ) {
            kinetic_matrix_original_fermions_up[i2*ist.n_spatial_orbitals+j2].real = matrix_element;
            kinetic_matrix_original_fermions_up[j2*ist.n_spatial_orbitals+i2].real = matrix_element;

            kinetic_matrix_original_fermions_up_2[i2*ist.n_spatial_orbitals+j2] = matrix_element;
            kinetic_matrix_original_fermions_up_2[j2*ist.n_spatial_orbitals+i2] = matrix_element;

           }
           else if ( count == 4 ) {
            kinetic_matrix_original_fermions_down[i2*ist.n_spatial_orbitals+j2].real = matrix_element;
            kinetic_matrix_original_fermions_down[j2*ist.n_spatial_orbitals+i2].real = matrix_element;

            kinetic_matrix_original_fermions_down_2[i2*ist.n_spatial_orbitals+j2] = matrix_element;
            kinetic_matrix_original_fermions_down_2[j2*ist.n_spatial_orbitals+i2] = matrix_element;

           } /*count 4*/            

         }
         else { /*Get Energy Shift*/
            if ( count == 0 ) {  
              energy_shifts[0].real = matrix_element;
              count++; 
            }
            else if ( count == 1 ) {
              energy_shifts[1].real = matrix_element; 
              count++; 
            }
            else if ( count == 2 ) {
              energy_shifts[2].real = matrix_element; 
              count++; 
            }
            else if ( count == 3 ) {
              count++; 
            }
            else if ( count == 4 ) {
              count++; 
            }
           else if (count == 5 ) {
              energy_shifts[3].real = matrix_element; 
              count++; 
           }   
         }
       }
    }

    /*Print Out Checks*********************************************************************************************************************/
   /*Print Up Kinetic*/
   for (i=0; i<ist.n_kinetic_sparse_up; i++) {
    fprintf(pf2, "%f %d\n", kinetic_matrix_sparse_up[i], kinetic_ij_sparse_up[i]); fflush(pf2); 
   }
   fprintf(pf2, "\n\n"); fflush(pf2); 

   /*Print Down Kinetic*/
   for (i=0; i<ist.n_kinetic_sparse_down; i++) {
    fprintf(pf2, "%f %d\n", kinetic_matrix_sparse_down[i], kinetic_ij_sparse_down[i]); fflush(pf2);
   }
   fprintf(pf2, "\n\n"); fflush(pf2);


   /*Run Throguh Stored Elements and Restore in Sparse Matrices*/
   index_potential_sparse_up = index_potential_sparse_down = index_potential_sparse_updown = 0; 
   index_kinetic_sparse_up = index_kinetic_sparse_down = 0; 
   for (i=0; i<ist.n_spatial_orbitals; i++) {
 
    for (j=0; j<ist.n_spatial_orbitals; j++) {

      /*Check If J Is In Set of Orbitals in Trial Wavefunction*/
      check_j = 0; 
      for (m=0; m<ist.n_max_orbital_trial_energy; m++) {
        if ( j == list_orbitals_occupied[m] ) {
          check_j = 1; 
        }
       }

      /*Only Take Orbitals in List*/ 
      if ( check_j == 1 ) { 

        /*Obtain Kinetic Terms First*/
        if ( fabs(kinetic_matrix_original_fermions_up_2[i*ist.n_spatial_orbitals+j]) > .0000001 ) {
          kinetic_matrix_sparse_up[index_kinetic_sparse_up] = kinetic_matrix_original_fermions_up_2[i*ist.n_spatial_orbitals+j];
          kinetic_ij_sparse_up[index_kinetic_sparse_up] = i*ist.n_spatial_orbitals+j;
          index_kinetic_sparse_up++;
        }
        if ( fabs(kinetic_matrix_original_fermions_down_2[i*ist.n_spatial_orbitals+j]) > .0000001 ) {
          kinetic_matrix_sparse_down[index_kinetic_sparse_down] = kinetic_matrix_original_fermions_down_2[i*ist.n_spatial_orbitals+j];
          kinetic_ij_sparse_down[index_kinetic_sparse_down] = i*ist.n_spatial_orbitals+j;
          index_kinetic_sparse_down++;
        } 
        
       for (k=0; k<ist.n_spatial_orbitals; k++) {

        for (l=0; l<ist.n_spatial_orbitals; l++) {

          /*Check If L Is In Set of Orbitals in Trial Wavefunction*/
          check_l = 0; 
          for (m=0; m<ist.n_max_orbital_trial_energy; m++) {
           if ( l == list_orbitals_occupied[m] ) {
            check_l = 1;
          }
         }

         /*Ensure That L Orbital Is In List*/
         if ( check_l == 1 ) { 

            /*Store in potential matrix sparse*/
            if ( fabs(potential_matrix_original_fermions_up_2[i*ist.n_spatial_orbitals_third+j*ist.n_spatial_orbitals_sq+k*ist.n_spatial_orbitals+l]) > .00000001 ) {

                potential_matrix_sparse_up[index_potential_sparse_up] = potential_matrix_original_fermions_up_2[i*ist.n_spatial_orbitals_third+j*ist.n_spatial_orbitals_sq+k*ist.n_spatial_orbitals+l]; 
                f_index_up = 4*index_potential_sparse_up;
                potential_ij_sparse_up[f_index_up] = i*ist.n_spatial_orbitals+j;
                potential_ij_sparse_up[f_index_up+1] = k*ist.n_spatial_orbitals+l;
                potential_ij_sparse_up[f_index_up+2] = i*ist.n_spatial_orbitals+l;
                potential_ij_sparse_up[f_index_up+3] = k*ist.n_spatial_orbitals+j;
                index_potential_sparse_up++;
            }
 
            if ( fabs(potential_matrix_original_fermions_down_2[i*ist.n_spatial_orbitals_third+j*ist.n_spatial_orbitals_sq+k*ist.n_spatial_orbitals+l]) > .00000001 ) {
                potential_matrix_sparse_down[index_potential_sparse_down] = potential_matrix_original_fermions_down_2[i*ist.n_spatial_orbitals_third+j*ist.n_spatial_orbitals_sq+k*ist.n_spatial_orbitals+l];
                f_index_down = 4*index_potential_sparse_down;
                potential_ij_sparse_down[f_index_down] = i*ist.n_spatial_orbitals+j;
                potential_ij_sparse_down[f_index_down+1] = k*ist.n_spatial_orbitals+l;
                potential_ij_sparse_down[f_index_down+2] = i*ist.n_spatial_orbitals+l;
                potential_ij_sparse_down[f_index_down+3] = k*ist.n_spatial_orbitals+j;
                index_potential_sparse_down++;
            }  
   
            if ( fabs(potential_matrix_original_fermions_updown_2[i*ist.n_spatial_orbitals_third+j*ist.n_spatial_orbitals_sq+k*ist.n_spatial_orbitals+l]) > .00000001 ) {
                f_index_updown = 2*index_potential_sparse_updown;

                potential_matrix_sparse_updown[f_index_updown] = potential_matrix_original_fermions_updown_2[i*ist.n_spatial_orbitals_third+j*ist.n_spatial_orbitals_sq+k*ist.n_spatial_orbitals+l];
                potential_matrix_sparse_updown[f_index_updown+1] = potential_matrix_original_fermions_updown_2[k*ist.n_spatial_orbitals_third+l*ist.n_spatial_orbitals_sq+i*ist.n_spatial_orbitals+j]; 
                potential_ij_sparse_updown[f_index_updown] = i*ist.n_spatial_orbitals+j;
                potential_ij_sparse_updown[f_index_updown+1] = k*ist.n_spatial_orbitals+l;
                index_potential_sparse_updown++;
            }
 
         } /*Check l Orbital*/

       } /*l*/
      } /*k*/

     } /*Check J Orbital*/

    } /*j*/
   } /*i*/

  /*Print Potential Matrix*/
/*  fprintf(pf2, "Sparse Up Potential\n"); fflush(pf2); 
  for (i=0; i<ist.n_potential_sparse_up; i++) {
    if ( fabs(potential_matrix_sparse_up[i]) > .000001 ) { 
    fprintf(pf2, "%f %d %d %d %d\n", potential_matrix_sparse_up[i], potential_ij_sparse_up[4*i], potential_ij_sparse_up[4*i+1], potential_ij_sparse_up[4*i+2], potential_ij_sparse_up[4*i+3]); fflush(pf2); 
     }
  }
  fprintf(pf2, "\n\n"); fflush(pf2); 

 fprintf(pf2, "Sparse Down Potential\n"); fflush(pf2);
  for (i=0; i<ist.n_potential_sparse_down; i++) {
    if ( fabs(potential_matrix_sparse_down[i]) > .000001 ) { 
      fprintf(pf2, "%f %d %d %d %d\n", potential_matrix_sparse_down[i], potential_ij_sparse_down[4*i], potential_ij_sparse_down[4*i+1], potential_ij_sparse_down[4*i+2], potential_ij_sparse_down[4*i+3]); fflush(pf2); 
    }
  }
  fprintf(pf2, "\n\n"); fflush(pf2); 

  fprintf(pf2, "Sparse Up Down Potential\n"); fflush(pf2); 
  for (i=0; i<ist.n_potential_sparse_updown; i++) {
    if ( fabs(potential_matrix_sparse_updown[i]) > .000001 ) {
       fprintf(pf2, "%f %f %d %d\n", potential_matrix_sparse_updown[2*i], potential_matrix_sparse_updown[2*i+1], potential_ij_sparse_updown[2*i], potential_ij_sparse_updown[2*i+1]); fflush(pf2); 
    }
  }
  fprintf(pf2, "\n"); fflush(pf2);   
  */


  /*Print Energy Shift*/
  fprintf(pf2, "energy shift up-up %f\n", energy_shifts[0].real); 
  fprintf(pf2, "energy shift down-down %f\n", energy_shifts[1].real); 
  fprintf(pf2, "energy shift up-down %f\n\n", energy_shifts[2].real);  
  fprintf(pf2, "energy shift pseudopotential %f\n\n", energy_shifts[3].real); fflush(pf2); 

  /*Return Indices*/
  *n_kinetic_sparse_up = index_kinetic_sparse_up; 
  *n_kinetic_sparse_down = index_kinetic_sparse_down; 
  *n_potential_sparse_up = index_potential_sparse_up; 
  *n_potential_sparse_down = index_potential_sparse_down; 
  *n_potential_sparse_updown = index_potential_sparse_updown; 

fclose(pf);
fclose(pf2);

free(potential_matrix_original_fermions_up_2); 
free(potential_matrix_original_fermions_down_2); 
free(potential_matrix_original_fermions_updown_2); 
free(kinetic_matrix_original_fermions_up_2); 
free(kinetic_matrix_original_fermions_down_2); 

return; 
}

/**********************************************************************************/

void init_potential_chemistry_restricted(MKL_Complex16 *potential_matrix_original_fermions,MKL_Complex16 *potential_matrix_original_fermions_spin,MKL_Complex16 *potential_supermatrix_fermions,MKL_Complex16 *potential_eigs_fermions,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,double *potential_one_body_matrix,int_st ist,cns_st cns) {

     /*Forms the Vijkl Supermatrix and Diagonalizes It*/  
     int i, j, k, l;
     int info;  
     int i_orbital, j_orbital, k_orbital, l_orbital; 
     int i_spin, j_spin, k_spin, l_spin; 
     int il_total, jk_total; 
     MKL_Complex16 sum; sum.real = sum.imag = 0.0;
     double *potential_matrix_3, *potential_matrix_2, *eigs_real; 
     FILE *pf2 = fopen("finalpotentialmatrix.dat", "a+"); 
     FILE *pf = fopen("errors.dat", "a+"); 

  /*HERE*/
     potential_matrix_3 = (double *)calloc(ist.n_spin_orbitals_fourth,sizeof(double)); 
     eigs_real = (double *)calloc(ist.n_spin_orbitals_sq,sizeof(double)); 

     /*First Restore Vijkl Matrix in Supermatrix Form*/
     czero_vec(potential_supermatrix_fermions,ist.n_spin_orbitals_fourth);
     dzero_vec(potential_matrix_3,ist.n_spin_orbitals_fourth);
     dzero_vec(potential_one_body_matrix,ist.n_spatial_orbitals_sq);

     for (i=0; i<ist.n_spin_orbitals; i++) {
       i_spin = i%2; 
       i_orbital = (int)(i/2.0); 

       for (l=0; l<ist.n_spin_orbitals; l++) {
          l_spin = l%2; 
          l_orbital = (int)(l/2.0); 

          il_total = i*ist.n_spin_orbitals+l;

              for (j=0; j<ist.n_spin_orbitals; j++) {
                 j_spin = j%2; 
                 j_orbital = (int)(j/2.0); 
                
                  for (k=0; k<ist.n_spin_orbitals; k++) {

                    k_spin = k%2; 
                    k_orbital = (int)(k/2.0);

                    //Ensure That Proper Vijkls Are Zero
                    if ( k_spin == j_spin ) {

                       if ( i_spin == l_spin ) {

                        jk_total = k*ist.n_spin_orbitals+j;

                        potential_supermatrix_fermions[il_total*ist.n_spin_orbitals_sq + jk_total] = potential_matrix_original_fermions[i_orbital * ist.n_spatial_orbitals_third + l_orbital * ist.n_spatial_orbitals_sq + j_orbital* ist.n_spatial_orbitals + k_orbital];
                        potential_matrix_3[il_total*ist.n_spin_orbitals_sq+jk_total] = potential_supermatrix_fermions[il_total*ist.n_spin_orbitals_sq+jk_total].real; 
                        sum.real += 1; 

                       if ( i_spin == k_spin ) {
                        if ( j == l && k_spin == 1 ) {
                          potential_one_body_matrix[i_orbital*ist.n_spatial_orbitals+k_orbital] -= potential_matrix_original_fermions[i_orbital * ist.n_spatial_orbitals_third + l_orbital * ist.n_spatial_orbitals_sq + j_orbital * ist.n_spatial_orbitals + k_orbital].real; fflush(pf);
                         }
                        }


                       } //k not l
                    } //i not j


             } //k
            } //j

        } //l
      } //i

     fprintf(pf2, "one body matrix\n"); 
     for (i=0; i<ist.n_spatial_orbitals; i++) {
      for (j=0; j<ist.n_spatial_orbitals; j++) {
       fprintf(pf2, "%f \t", potential_one_body_matrix[i*ist.n_spatial_orbitals+j]); 
      }
      fprintf(pf2, "\n"); 
     }
     fprintf(pf2, "\n\n"); fflush(pf2);

     fprintf(pf2, "non zero potential matrix terms\n"); 
     for (i=0; i<ist.n_spin_orbitals; i++) {
      for (j=0; j<ist.n_spin_orbitals; j++) {
       for (k=0; k<ist.n_spin_orbitals; k++) {
        for (l=0; l<ist.n_spin_orbitals; l++) {
           if ( potential_matrix_3[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l] != 0 ) {
             fprintf(pf2, "%d %d %d %d %f\n", i, j, k, l, potential_matrix_3[i*ist.n_spin_orbitals_third+j*ist.n_spin_orbitals_sq+k*ist.n_spin_orbitals+l]); 
           }
        }
       }
      }
     }
     fflush(pf2); 


      /*Find Eigenvalues and Vectors of Supermatrix*/  
      if ( sum.real != 0.0 ) { 
        info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',ist.n_spin_orbitals_sq,potential_matrix_3,ist.n_spin_orbitals_sq,eigs_real);

         /*Hack for Now - Should Reall Insert Complex Diagonalization Routine*/
         for (i=0; i<ist.n_spin_orbitals_sq; i++) {
          for (j=0; j<ist.n_spin_orbitals_sq; j++) {
           potential_eigvecs_fermions[i*ist.n_spin_orbitals_sq+j].real = potential_matrix_3[j*ist.n_spin_orbitals_sq+i]; 
          }
          potential_eigs_fermions[i].real = eigs_real[i]; 
         } 

      }
      else { 
        czero_vec(potential_eigs_fermions,ist.n_spin_orbitals_sq); 
        cunity_vec(potential_eigvecs_fermions,ist.n_spin_orbitals_sq); 
      }

      /*Now Form HS Transform Constants from Eigenvalues*/
      transform_eigenvalues(potential_eigs_fermions,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,ist,cns); 

free(potential_matrix_3);   
free(eigs_real); 
return; 
}

/*******************************************************************************************************************/

void init_potential_chemistry_unrestricted(MKL_Complex16 *potential_matrix_original_fermions_up,MKL_Complex16 *potential_matrix_original_fermions_down,MKL_Complex16 *potential_matrix_original_fermions_updown,MKL_Complex16 *potential_supermatrix_fermions,MKL_Complex16 *potential_eigs_fermions,MKL_Complex16 *potential_eigvecs_fermions,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,double *potential_one_body_matrix_up,double *potential_one_body_matrix_down,int_st ist,cns_st cns) {

     /*Forms the Vijkl Supermatrix and Diagonalizes It*/
     int i, j, k, l;
     int info; 
     int i_orbital, j_orbital, k_orbital, l_orbital;
     int i_spin, j_spin, k_spin, l_spin;
     int il_total, jk_total;
     MKL_Complex16 sum; sum.real = sum.imag = 0.0;
     double *potential_matrix_2, *eigs_2; 
     FILE *pf = fopen("errors.dat", "a+");  

     potential_matrix_2 = (double *)calloc(ist.n_spin_orbitals_fourth,sizeof(double));
     eigs_2 = (double *)calloc(ist.n_spin_orbitals_sq,sizeof(double));

     /*First Restore Vijkl Matrix in Supermatrix Form*/
     czero_vec(potential_supermatrix_fermions,ist.n_spin_orbitals_fourth);
     dzero_vec(potential_matrix_2,ist.n_spin_orbitals_fourth);
     for (i=0; i<ist.n_spin_orbitals; i++) {
       i_spin = i%2;
       i_orbital = (int)(i/2.0);

       for (l=0; l<ist.n_spin_orbitals; l++) {
          l_spin = l%2;
          l_orbital = (int)(l/2.0);

          il_total = i*ist.n_spin_orbitals+l;

          for (j=0; j<ist.n_spin_orbitals; j++) {
           j_spin = j%2;
           j_orbital = (int)(j/2.0);

           for (k=0; k<ist.n_spin_orbitals; k++) {
            k_spin = k%2;
            k_orbital = (int)(k/2.0);

            if ( k_spin == j_spin ) { /*If K and J Spin Equal*/
 
              if ( i_spin == l_spin ) {

               jk_total = k*ist.n_spin_orbitals+j;

               /*If Both Spin Up*/
               if ( i_spin == 0 ) {
                 if ( j_spin == 0 ) {
                    potential_supermatrix_fermions[il_total*ist.n_spin_orbitals_sq + jk_total] = potential_matrix_original_fermions_up[i_orbital * ist.n_spatial_orbitals_third + l_orbital * ist.n_spatial_orbitals_sq + j_orbital* ist.n_spatial_orbitals + k_orbital];
                 }
                 else {
                    potential_supermatrix_fermions[il_total*ist.n_spin_orbitals_sq + jk_total] = potential_matrix_original_fermions_updown[i_orbital * ist.n_spatial_orbitals_third + l_orbital * ist.n_spatial_orbitals_sq + j_orbital* ist.n_spatial_orbitals + k_orbital];
                 }
               }
               else {
                 if ( j_spin == 1 ) {
                    potential_supermatrix_fermions[il_total*ist.n_spin_orbitals_sq + jk_total] = potential_matrix_original_fermions_down[i_orbital * ist.n_spatial_orbitals_third + l_orbital * ist.n_spatial_orbitals_sq + j_orbital* ist.n_spatial_orbitals + k_orbital];
                 }
                 else {
                     potential_supermatrix_fermions[il_total*ist.n_spin_orbitals_sq + jk_total] = potential_matrix_original_fermions_updown[j_orbital * ist.n_spatial_orbitals_third + k_orbital * ist.n_spatial_orbitals_sq + i_orbital* ist.n_spatial_orbitals + l_orbital];
                  }
               }

               potential_matrix_2[il_total*ist.n_spin_orbitals_sq+jk_total]=potential_supermatrix_fermions[il_total*ist.n_spin_orbitals_sq+jk_total].real;
               sum.real += 1;

               /*Check Factor of 1/2 Or Not*/
               if ( i_spin == k_spin ) {
                if ( j == l ) {
                  if ( k_spin == 1 ) {
                    potential_one_body_matrix_down[i_orbital*ist.n_spatial_orbitals+k_orbital] -= .5 * potential_matrix_original_fermions_down[i_orbital * ist.n_spatial_orbitals_third + l_orbital * ist.n_spatial_orbitals_sq + j_orbital * ist.n_spatial_orbitals + k_orbital].real; 
                  }
                  else { 
                    potential_one_body_matrix_up[i_orbital*ist.n_spatial_orbitals+k_orbital] -= .5 * potential_matrix_original_fermions_up[i_orbital * ist.n_spatial_orbitals_third + l_orbital * ist.n_spatial_orbitals_sq + j_orbital * ist.n_spatial_orbitals + k_orbital].real; 
                  }
                }
               }

              } /*i equals l spin*/
            } /*j equals k spin*/

          } /*k*/
         } /*j*/
       } /*l*/
      } /*i*/

      /*Find Eigenvalues and Vectors of Supermatrix*/
      if ( sum.real != 0.0 ) {
        /*Assume All Real for Now*/
        info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',ist.n_spin_orbitals_sq,potential_matrix_2,ist.n_spin_orbitals_sq,eigs_2);

         for (i=0; i<ist.n_spin_orbitals_sq; i++) {
          for (j=0; j<ist.n_spin_orbitals_sq; j++) {
           potential_eigvecs_fermions[i*ist.n_spin_orbitals_sq+j].real = potential_matrix_2[j*ist.n_spin_orbitals_sq+i];
          }
          potential_eigs_fermions[i].real = eigs_2[i];
         }

      }
      else {
        czero_vec(potential_eigs_fermions,ist.n_spin_orbitals_sq);
        cunity_vec(potential_eigvecs_fermions,ist.n_spin_orbitals_sq);
      }

      /*Now Form HS Transform Constants from Eigenvalues*/
      transform_eigenvalues(potential_eigs_fermions,potential_eigs_fermions_first_transform,potential_eigs_fermions_second_transform,ist,cns);

free(potential_matrix_2);
free(eigs_2);
return;
}

/************************************************************************/

void transform_eigenvalues(MKL_Complex16 *potential_eigs_fermions,MKL_Complex16 *potential_eigs_fermions_first_transform,MKL_Complex16 *potential_eigs_fermions_second_transform,int_st ist,cns_st cns){

    /*Transforms the Eigenvalues Into The Constants in the Exponentials of HS Transformation*/

    int i;

    for (i=0; i<ist.n_spin_orbitals_sq; i++) {

      /*Check If Real or Complex - Assuming Eigs Are Real for Now, But Not Necessarily...Should Be Careful Here*/

      /*Now, Why This Ordering...of +, - in Transform...*/ 
      potential_eigs_fermions_first_transform[i] = Csqrt(RCmul(-0.25*cns.dtau, potential_eigs_fermions[i])); 
      potential_eigs_fermions_second_transform[i] = Csqrt(RCmul(0.25*cns.dtau, potential_eigs_fermions[i])); 

    }


return; 
}

/**************************************************************************************************/


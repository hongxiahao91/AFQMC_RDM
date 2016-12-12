#include "afqmc.h"

/******************************************************************/

int population_control_fermions(MKL_Complex16 *weights,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *new_walker_weights,MKL_Complex16 *new_wf_up,MKL_Complex16 *new_wf_down,MKL_Complex16 *new_overlap_inverse_up,MKL_Complex16 *new_overlap_inverse_down,MKL_Complex16 *new_overlap_up,MKL_Complex16 *new_overlap_down,int *walkers_big,int *walkers_small,int *walkers_copied,double *weights_rescaled,int_st ist,cns_st cns) {
    /*Population Control Scheme Based Upon Shiwei's Original Scheme*/
    /*Current Version May Be Sped Up By Copying New Walkers More Efficiently...Test Version Now*/ 

    int i, j, k;   
    int ip, np;
    int number_walkers_new;  
    int number_copies, number_bad_walkers, number_big_walkers, number_small_walkers; 
    double phase; 
    double population_control_factor;
    double abs_value;
    double adjusted_weight;  
    long tidum;

    Randomize(); tidum = -random();

    adjusted_weight = 1.0; 

    begin_population_control_fermions: {}; 

    number_walkers_new = ist.n_walkers; 
    number_big_walkers = 0; 
    number_small_walkers = 0; 
    number_bad_walkers = 0; 

    izero_vec(walkers_big,ist.n_walkers); 
    izero_vec(walkers_small,ist.n_walkers); 
    izero_vec(walkers_copied,ist.n_walkers); 

    /*Run Through Walkers Figuring Out if Too Big or Too Small*/
    for (i=0; i<ist.n_walkers; i++) {

        /*If Walker Weight Too Small, Then Mark To Get Rid Of*/ 
        if ( Cabs(weights[i]) < cns.min_walker_weight ) {

          /*Store Walker Info*/
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          number_bad_walkers++;  
          if ( number_copies == 0 ) {
            walkers_small[i] = 1;
            number_small_walkers++; 
            number_walkers_new -= 1;      
          }

        }  /*If Walkers Weights Too Larger, Then Split Into Multiple*/
        else if ( Cabs(weights[i]) > cns.max_walker_weight) {

          /*Store Walkers Info*/
          walkers_big[number_big_walkers] = i;
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          walkers_copied[number_big_walkers] = number_copies - 1;  
          number_big_walkers++; 
          number_bad_walkers++; 
          number_walkers_new += number_copies - 1; 

        } 
    }


    /*End Population Control If There Are No Bad Walkers*/
    if ( number_bad_walkers == 0) { 
      goto end_population_control_fermions; 
    }
    else {    /*If Bad Walkers*/

      /*If Between Counter Values*****************/
      if ( number_walkers_new >= cns.min_number_walkers && number_walkers_new <= cns.max_number_walkers ) {

         /*Reset the weights of All Walkers That Are Too Big or Small*/
         number_walkers_new = 0; 
         for (i=0; i<ist.n_walkers; i++) {
         
            abs_value = Cabs(weights[i]);  
            if ( ( abs_value < cns.min_walker_weight ) || ( abs_value > cns.max_walker_weight )  ) {
                phase = atan(weights[i].imag/weights[i].real);  
                weights[i].real = cos(phase); 
                weights[i].imag = sin(phase); 
            } 
  
            /*Copy Walkers and Matrices*/
            if ( walkers_small[i] != 1 ) {
             new_walker_weights[number_walkers_new] = weights[i]; 
             new_overlap_up[number_walkers_new] = overlap_up[i]; 
             new_overlap_down[number_walkers_new] = overlap_down[i];  

               for (j=0; j<ist.n_sites_n_up; j++) {
                 np = number_walkers_new * ist.n_sites_n_up;
                 ip = i * ist.n_sites_n_up;  
 
                 new_wf_up[np+j] = wf_up[ip+j]; 
               }

               for (j=0; j<ist.n_sites_n_down; j++) {
                 np = number_walkers_new * ist.n_sites_n_down; 
                 ip = i * ist.n_sites_n_down;      
 
                 new_wf_down[np+j] = wf_down[ip+j]; 
               }        
               for (j=0; j<ist.n_up_sq; j++) {
                 np = number_walkers_new * ist.n_up_sq; 
                 ip = i * ist.n_up_sq; 
  
                 new_overlap_inverse_up[np+j] = overlap_inverse_up[ip+j]; 
               }
               for (j=0; j<ist.n_down_sq; j++) {
                 np = number_walkers_new * ist.n_down_sq;   
                 ip = i * ist.n_down_sq;   

                 new_overlap_inverse_down[np+j] = overlap_inverse_down[ip+j];
               }
 
              number_walkers_new++;
            }
        } /*Run Through Walkers*/
 
          /*Now Add New Walkers At the End of the LIst*/
          for (i=0; i<number_big_walkers; i++) {
 
           for (j=0; j< walkers_copied[i]; j++) {
             new_walker_weights[number_walkers_new] = weights[walkers_big[i]];
             new_overlap_up[number_walkers_new] = overlap_up[walkers_big[i]]; 
             new_overlap_down[number_walkers_new] = overlap_down[walkers_big[i]];  

             for (k=0; k<ist.n_sites_n_up; k++) {
               np = number_walkers_new*ist.n_sites_n_up; 
               ip = walkers_big[i]*ist.n_sites_n_up;
 
               new_wf_up[np+k] = wf_up[ip+k]; 
             }
             for (k=0; k<ist.n_sites_n_down; k++) {
               np = number_walkers_new*ist.n_sites_n_down;
               ip = walkers_big[i]*ist.n_sites_n_down; 

               new_wf_down[np+k] = wf_down[ip+k]; 
             }
             for (k=0; k<ist.n_up_sq; k++) {     
               np = number_walkers_new*ist.n_up_sq;
               ip = walkers_big[i]*ist.n_up_sq; 

               new_overlap_inverse_up[np+k] = overlap_inverse_up[ip+k]; 
             }
             for (k=0; k<ist.n_down_sq; k++) {
               np = number_walkers_new*ist.n_down_sq;
               ip = walkers_big[i]*ist.n_down_sq;

               new_overlap_inverse_down[np+k] = overlap_inverse_down[ip+k];
             }

             number_walkers_new++; 
           }
         }
  

        /*Last Thing Is To Now Copy Over Walkers*/
        cblas_zcopy(number_walkers_new,new_walker_weights,1,weights,1); 
        cblas_zcopy(number_walkers_new,new_overlap_up,1,overlap_up,1);
        cblas_zcopy(number_walkers_new,new_overlap_down,1,overlap_down,1);   
        //copy_cmat(weights,new_walker_weights,number_walkers_new); 
        //copy_cmat(overlap_up,new_overlap_up,number_walkers_new); 
        //copy_cmat(overlap_down,new_overlap_down,number_walkers_new);  

        cblas_zcopy(number_walkers_new*ist.n_sites_n_up,new_wf_up,1,wf_up,1); 
        cblas_zcopy(number_walkers_new*ist.n_sites_n_down,new_wf_down,1,wf_down,1); 
        cblas_zcopy(number_walkers_new*ist.n_up_sq,new_overlap_inverse_up,1,overlap_inverse_up,1); 
        cblas_zcopy(number_walkers_new*ist.n_down_sq,new_overlap_inverse_down,1,overlap_inverse_down,1); 
        //copy_cmat(wf_up,new_wf_up,number_walkers_new*ist.n_sites_n_up);
        //copy_cmat(wf_down,new_wf_down,number_walkers_new*ist.n_sites_n_down); 
        //copy_cmat(overlap_inverse_up,new_overlap_inverse_up,number_walkers_new*ist.n_up_sq); 
        //copy_cmat(overlap_inverse_down,new_overlap_inverse_down,number_walkers_new*ist.n_down_sq); 

    } /*If Pop Changes Within Bounds*/
    else {   /*If Pop Change Is Outside Bounds*/ 

        /*Increase Or Decrease Number of Walkers Appropriately*/        
        if ( number_walkers_new < cns.min_number_walkers ) {
           population_control_factor=pow(cns.population_control_factor, -1);
        }
        else {
           population_control_factor=cns.population_control_factor;
        }

        for (i=0; i<ist.n_walkers; i++) {
           weights[i]=RCmul(population_control_factor, weights[i]);
         }

         adjusted_weight /= population_control_factor;

         goto begin_population_control_fermions;
    }

  } /*else Nothing to change*/

   end_population_control_fermions: {};  

  /*Storing Rescaled Weights*/
  for (i=8; i>-1; i--) {
   weights_rescaled[i+1] = weights_rescaled[i]; 
  }
  weights_rescaled[0] = adjusted_weight; 

return(number_walkers_new);  
}

/******************************************************************/

int population_control_chemistry(MKL_Complex16 *weights,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *new_walker_weights,MKL_Complex16 *new_wf_up,MKL_Complex16 *new_wf_down,MKL_Complex16 *new_overlap_inverse_up,MKL_Complex16 *new_overlap_inverse_down,MKL_Complex16 *new_overlap_up,MKL_Complex16 *new_overlap_down,int *walkers_big,int *walkers_small,int *walkers_copied,double *weights_rescaled,int_st ist,cns_st cns) {
    /*Population Control Scheme Based Upon Shiwei's Original Scheme*/
    /*Current Version May Be Sped Up By Copying New Walkers More Efficiently...Test Version Now*/

    int i, j, k;
    int ip, np;
    int number_walkers_new;
    int number_copies, number_bad_walkers, number_big_walkers, number_small_walkers;
    double phase;
    double population_control_factor;
    double abs_value;
    double adjusted_weight;
    MKL_Complex16 total_weight_before, total_weight_after; 
    long tidum;
    FILE *pf = fopen("errors.dat", "a+"); 

    total_weight_before.real = total_weight_before.imag = total_weight_after.real = total_weight_after.imag = 0.0; 

    Randomize(); tidum = -random();

    adjusted_weight = 1.0;

    begin_population_control_chemistry: {};

    number_walkers_new = ist.n_walkers;
    number_big_walkers = 0;
    number_small_walkers = 0;
    number_bad_walkers = 0;

    izero_vec(walkers_big,ist.n_walkers);
    izero_vec(walkers_small,ist.n_walkers);
    izero_vec(walkers_copied,ist.n_walkers);

    /*Run Through Walkers Figuring Out if Too Big or Too Small*/
    for (i=0; i<ist.n_walkers; i++) {

        /*If Walker Weight Too Small, Then Mark To Get Rid Of*/
        if ( Cabs(weights[i]) < cns.min_walker_weight ) {

          /*Store Walker Info*/
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          number_bad_walkers++;
          if ( number_copies == 0 ) {
            walkers_small[i] = 1;
            number_small_walkers++;
            number_walkers_new -= 1;
          }

        }  /*If Walkers Weights Too Larger, Then Split Into Multiple*/
        else if ( Cabs(weights[i]) > cns.max_walker_weight) {

          /*Store Walkers Info*/
          walkers_big[number_big_walkers] = i;
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          walkers_copied[number_big_walkers] = number_copies - 1;
          number_big_walkers++;
          number_bad_walkers++;
          number_walkers_new += number_copies - 1;

        }
    }


    /*End Population Control If There Are No Bad Walkers*/
    if ( number_bad_walkers == 0) {
      goto end_population_control_chemistry;
    }
    else {    /*If Bad Walkers*/

      /*If Between Counter Values*****************/
      if ( number_walkers_new >= cns.min_number_walkers && number_walkers_new <= cns.max_number_walkers ) {


         /*Reset the weights of All Walkers That Are Too Big or Small*/
         number_walkers_new = 0;
         for (i=0; i<ist.n_walkers; i++) {

            abs_value = Cabs(weights[i]);
            if ( ( abs_value < cns.min_walker_weight ) || ( abs_value > cns.max_walker_weight )  ) {
                phase = atan(weights[i].imag/weights[i].real);
                weights[i].real = cos(phase);
                weights[i].imag = sin(phase);
            }

            /*Copy Walkers and Matrices*/
            if ( walkers_small[i] != 1 ) {
             new_walker_weights[number_walkers_new] = weights[i];
             new_overlap_up[number_walkers_new] = overlap_up[i];
             new_overlap_down[number_walkers_new] = overlap_down[i];

               for (j=0; j<ist.n_spatial_orbitals_n_up; j++) {
                 np = number_walkers_new * ist.n_spatial_orbitals_n_up;
                 ip = i * ist.n_spatial_orbitals_n_up;

                 new_wf_up[np+j] = wf_up[ip+j];
               }

               for (j=0; j<ist.n_spatial_orbitals_n_down; j++) {
                 np = number_walkers_new * ist.n_spatial_orbitals_n_down;
                 ip = i * ist.n_spatial_orbitals_n_down;

                 new_wf_down[np+j] = wf_down[ip+j];
               }
               for (j=0; j<ist.n_up_sq; j++) {
                 np = number_walkers_new * ist.n_up_sq;
                 ip = i * ist.n_up_sq;

                 new_overlap_inverse_up[np+j] = overlap_inverse_up[ip+j];
               }
               for (j=0; j<ist.n_down_sq; j++) {
                 np = number_walkers_new * ist.n_down_sq;
                 ip = i * ist.n_down_sq;

                 new_overlap_inverse_down[np+j] = overlap_inverse_down[ip+j];
               }

              number_walkers_new++;
            }

        } /*Run Through Walkers*/


          /*Now Add New Walkers At the End of the LIst*/
          for (i=0; i<number_big_walkers; i++) {

           for (j=0; j< walkers_copied[i]; j++) {
             new_walker_weights[number_walkers_new] = weights[walkers_big[i]];
             new_overlap_up[number_walkers_new] = overlap_up[walkers_big[i]];
             new_overlap_down[number_walkers_new] = overlap_down[walkers_big[i]];

             for (k=0; k<ist.n_spatial_orbitals_n_up; k++) {
               np = number_walkers_new*ist.n_spatial_orbitals_n_up;
               ip = walkers_big[i]*ist.n_spatial_orbitals_n_up;

               new_wf_up[np+k] = wf_up[ip+k];
             }
             for (k=0; k<ist.n_spatial_orbitals_n_down; k++) {
               np = number_walkers_new*ist.n_spatial_orbitals_n_down;
               ip = walkers_big[i]*ist.n_spatial_orbitals_n_down;

               new_wf_down[np+k] = wf_down[ip+k];
             }
             for (k=0; k<ist.n_up_sq; k++) {
               np = number_walkers_new*ist.n_up_sq;
               ip = walkers_big[i]*ist.n_up_sq;

               new_overlap_inverse_up[np+k] = overlap_inverse_up[ip+k];
             }
             for (k=0; k<ist.n_down_sq; k++) {
               np = number_walkers_new*ist.n_down_sq;
               ip = walkers_big[i]*ist.n_down_sq;

               new_overlap_inverse_down[np+k] = overlap_inverse_down[ip+k];
             }

             number_walkers_new++;
           }
         }


        /*Last Thing Is To Now Copy Over Walkers*/
        cblas_zcopy(number_walkers_new,new_walker_weights,1,weights,1); 
        cblas_zcopy(number_walkers_new,new_overlap_up,1,overlap_up,1); 
        cblas_zcopy(number_walkers_new,new_overlap_down,1,overlap_down,1); 

        cblas_zcopy(number_walkers_new*ist.n_spatial_orbitals_n_up,new_wf_up,1,wf_up,1); 
        cblas_zcopy(number_walkers_new*ist.n_spatial_orbitals_n_down,new_wf_down,1,wf_down,1); 
        cblas_zcopy(number_walkers_new*ist.n_up_sq,new_overlap_inverse_up,1,overlap_inverse_up,1); 
        cblas_zcopy(number_walkers_new*ist.n_down_sq,new_overlap_inverse_down,1,overlap_inverse_down,1);  

    } /*If Pop Changes Within Bounds*/
    else {   /*If Pop Change Is Outside Bounds*/

        /*Increase Or Decrease Number of Walkers Appropriately*/
        if ( number_walkers_new < cns.min_number_walkers ) {
           population_control_factor=pow(cns.population_control_factor, -1);
        }
        else {
           population_control_factor=cns.population_control_factor;
        }

        for (i=0; i<ist.n_walkers; i++) {
           weights[i]=RCmul(population_control_factor, weights[i]);
         }

         adjusted_weight /= population_control_factor;

         goto begin_population_control_chemistry;
    }

  } /*else Nothing to change*/

   end_population_control_chemistry: {};

  /*Storing Rescaled Weights*/
  for (i=8; i>-1; i--) {
   weights_rescaled[i+1] = weights_rescaled[i];
  }
  weights_rescaled[0] = adjusted_weight;

return(number_walkers_new);
}

/*******************************************************************/

int population_control_chemistry_multi(MKL_Complex16 *weights,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_total,MKL_Complex16 *det_overlap_up,MKL_Complex16 *det_overlap_down,MKL_Complex16 *new_walker_weights,MKL_Complex16 *new_wf_up,MKL_Complex16 *new_wf_down,MKL_Complex16 *new_overlap_inverse_up,MKL_Complex16 *new_overlap_inverse_down,MKL_Complex16 *new_overlap_total,MKL_Complex16 *new_det_overlap_up,MKL_Complex16 *new_det_overlap_down,int *walkers_big,int *walkers_small,int *walkers_copied,double *weights_rescaled,int_st ist,cns_st cns) {
    /*Population Control Scheme Based Upon Shiwei's Original Scheme*/
    /*Current Version May Be Sped Up By Copying New Walkers More Efficiently...Test Version Now*/

    int i, j, k;
    int ip, np;
    int number_walkers_new;
    int number_copies, number_bad_walkers, number_big_walkers, number_small_walkers;
    double phase;
    double population_control_factor;
    double abs_value;
    double adjusted_weight;
    MKL_Complex16 total_weight_before, total_weight_after;
    long tidum;

    total_weight_before.real = total_weight_before.imag = total_weight_after.real = total_weight_after.imag = 0.0;

    Randomize(); tidum = -random();

    adjusted_weight = 1.0;

    begin_population_control_chemistry: {};

    number_walkers_new = ist.n_walkers;
    number_big_walkers = 0;
    number_small_walkers = 0;
    number_bad_walkers = 0;

    izero_vec(walkers_big,ist.n_walkers);
    izero_vec(walkers_small,ist.n_walkers);
    izero_vec(walkers_copied,ist.n_walkers);

    /*Run Through Walkers Figuring Out if Too Big or Too Small*/
    for (i=0; i<ist.n_walkers; i++) {

        /*If Walker Weight Too Small, Then Mark To Get Rid Of*/
        if ( Cabs(weights[i]) < cns.min_walker_weight ) {

          /*Store Walker Info*/
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          number_bad_walkers++;
          if ( number_copies == 0 ) {
            walkers_small[i] = 1;
            number_small_walkers++;
            number_walkers_new -= 1;
          }

        }  /*If Walkers Weights Too Larger, Then Split Into Multiple*/
        else if ( Cabs(weights[i]) > cns.max_walker_weight ) {

          /*Store Walkers Info*/
          walkers_big[number_big_walkers] = i;
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          walkers_copied[number_big_walkers] = number_copies - 1;
          number_big_walkers++;
          number_bad_walkers++;
          number_walkers_new += number_copies - 1;

        }
    }


    /*End Population Control If There Are No Bad Walkers*/
    if ( number_bad_walkers == 0) {
      goto end_population_control_chemistry;
    }
    else {    /*If Bad Walkers*/

      /*If Between Counter Values*****************/
      if ( number_walkers_new >= cns.min_number_walkers && number_walkers_new <= cns.max_number_walkers ) {

         /*Reset the weights of All Walkers That Are Too Big or Small*/
         number_walkers_new = 0;
         for (i=0; i<ist.n_walkers; i++) {

            abs_value = Cabs(weights[i]);
            if ( ( abs_value < cns.min_walker_weight ) || ( abs_value > cns.max_walker_weight )  ) {
                phase = atan(weights[i].imag/weights[i].real);
                weights[i].real = cos(phase);
                weights[i].imag = sin(phase);
            }

            /*Copy Walkers and Matrices*/
            if ( walkers_small[i] != 1 ) {
             new_walker_weights[number_walkers_new] = weights[i];
             new_overlap_total[number_walkers_new] = overlap_total[i];

             for (j=0; j<ist.n_determinants_trial_energy; j++) {
                np = number_walkers_new * ist.n_determinants_trial_energy; 
                ip = i * ist.n_determinants_trial_energy; 

                new_det_overlap_up[np+j] = det_overlap_up[ip + j]; 
                new_det_overlap_down[np+j] = det_overlap_down[ip+j]; 
             }               


              for (j=0; j<ist.n_spatial_orbitals_n_up; j++) {
                 np = number_walkers_new * ist.n_spatial_orbitals_n_up;
                 ip = i * ist.n_spatial_orbitals_n_up;

                 new_wf_up[np+j] = wf_up[ip+j];
               }

               for (j=0; j<ist.n_spatial_orbitals_n_down; j++) {
                 np = number_walkers_new * ist.n_spatial_orbitals_n_down;
                 ip = i * ist.n_spatial_orbitals_n_down;

                 new_wf_down[np+j] = wf_down[ip+j];
               }
               for (j=0; j<ist.n_up_sq; j++) {
                 np = number_walkers_new * ist.n_determinants_energy_n_up_sq;
                 ip = i * ist.n_determinants_energy_n_up_sq;

                 new_overlap_inverse_up[np+j] = overlap_inverse_up[ip+j];
               }
               for (j=0; j<ist.n_down_sq; j++) {
                 np = number_walkers_new * ist.n_determinants_energy_n_down_sq;
                 ip = i * ist.n_determinants_energy_n_down_sq;

                 new_overlap_inverse_down[np+j] = overlap_inverse_down[ip+j];
               }

              number_walkers_new++;
            }

        } /*Run Through Walkers*/


          /*Now Add New Walkers At the End of the LIst*/
          for (i=0; i<number_big_walkers; i++) {

           for (j=0; j< walkers_copied[i]; j++) {
             new_walker_weights[number_walkers_new] = weights[walkers_big[i]];
             new_overlap_total[number_walkers_new] = overlap_total[walkers_big[i]];

             for (k=0; k<ist.n_determinants_trial_energy; k++) {
                np = number_walkers_new * ist.n_determinants_trial_energy;
                ip = walkers_big[i] * ist.n_determinants_trial_energy;

                new_det_overlap_up[np+k] = det_overlap_up[ip + k];
                new_det_overlap_down[np+k] = det_overlap_down[ip+k];
             }

             for (k=0; k<ist.n_spatial_orbitals_n_up; k++) {
               np = number_walkers_new*ist.n_spatial_orbitals_n_up;
               ip = walkers_big[i]*ist.n_spatial_orbitals_n_up;

               new_wf_up[np+k] = wf_up[ip+k];
             }
             for (k=0; k<ist.n_spatial_orbitals_n_down; k++) {
               np = number_walkers_new*ist.n_spatial_orbitals_n_down;
               ip = walkers_big[i]*ist.n_spatial_orbitals_n_down;

               new_wf_down[np+k] = wf_down[ip+k];
             }
             for (k=0; k<ist.n_up_sq; k++) {
               np = number_walkers_new*ist.n_determinants_energy_n_up_sq;
               ip = walkers_big[i]*ist.n_determinants_energy_n_up_sq;

               new_overlap_inverse_up[np+k] = overlap_inverse_up[ip+k];
             }
             for (k=0; k<ist.n_down_sq; k++) {
               np = number_walkers_new*ist.n_determinants_energy_n_down_sq;
               ip = walkers_big[i]*ist.n_determinants_energy_n_down_sq;

               new_overlap_inverse_down[np+k] = overlap_inverse_down[ip+k];
             }

             number_walkers_new++;
           }
         }


        /*Last Thing Is To Now Copy Over Walkers*/
        cblas_zcopy(number_walkers_new,new_walker_weights,1,weights,1);
        cblas_zcopy(number_walkers_new,new_overlap_total,1,overlap_total,1);
       
        cblas_zcopy(number_walkers_new*ist.n_determinants_trial_energy,new_det_overlap_up,1,det_overlap_up,1); 
        cblas_zcopy(number_walkers_new*ist.n_determinants_trial_energy,new_det_overlap_down,1,det_overlap_down,1); 

        cblas_zcopy(number_walkers_new*ist.n_spatial_orbitals_n_up,new_wf_up,1,wf_up,1);
        cblas_zcopy(number_walkers_new*ist.n_spatial_orbitals_n_down,new_wf_down,1,wf_down,1);
        cblas_zcopy(number_walkers_new*ist.n_determinants_energy_n_up_sq,new_overlap_inverse_up,1,overlap_inverse_up,1);
        cblas_zcopy(number_walkers_new*ist.n_determinants_energy_n_down_sq,new_overlap_inverse_down,1,overlap_inverse_down,1);


    } /*If Pop Changes Within Bounds*/
    else {   /*If Pop Change Is Outside Bounds*/

        /*Increase Or Decrease Number of Walkers Appropriately*/
        if ( number_walkers_new < cns.min_number_walkers ) {
           population_control_factor=pow(cns.population_control_factor, -1);
        }
        else {
           population_control_factor=cns.population_control_factor;
        }

        for (i=0; i<ist.n_walkers; i++) {
           weights[i]=RCmul(population_control_factor, weights[i]);
         }


         adjusted_weight /= population_control_factor;

         goto begin_population_control_chemistry;
    }

  } /*else Nothing to change*/

   end_population_control_chemistry: {};

  /*Storing Rescaled Weights*/
  for (i=8; i>-1; i--) {
   weights_rescaled[i+1] = weights_rescaled[i];
  }
  weights_rescaled[0] = adjusted_weight;

return(number_walkers_new);
}

/******************************************************************/

int population_control_bosons(MKL_Complex16 *weights,MKL_Complex16 *wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *new_walker_weights,MKL_Complex16 *new_wf_bosons,MKL_Complex16 *new_overlap_bosons,int *walkers_big,int *walkers_small,int *walkers_copied,double *weights_rescaled,int_st ist,cns_st cns) {
    /*Population Control Scheme Based Upon Shiwei's Original Scheme*/
    /*Current Version May Be Sped Up By Copying New Walkers More Efficiently...Test Version Now*/

    int i, j, k;
    int ip, np;
    int number_walkers_new;
    int number_copies, number_bad_walkers, number_big_walkers, number_small_walkers;
    double phase;
    double population_control_factor;
    double abs_value;
    double adjusted_weight;
    long tidum;

    Randomize(); tidum = -random();

    adjusted_weight = 1.0;

    begin_population_control_bosons: {};

    number_walkers_new = ist.n_walkers;
    number_big_walkers = 0;
    number_small_walkers = 0;
    number_bad_walkers = 0;

    izero_vec(walkers_big,ist.n_walkers);
    izero_vec(walkers_small,ist.n_walkers);
    izero_vec(walkers_copied,ist.n_walkers);

    /*Run Through Walkers Figuring Out if Too Big or Too Small*/
    for (i=0; i<ist.n_walkers; i++) {

        /*If Walker Weight Too Small, Then Mark To Get Rid Of*/
        if ( Cabs(weights[i]) < cns.min_walker_weight ) {

          /*Store Walker Info*/
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          number_bad_walkers++;
          if ( number_copies == 0 ) {
            walkers_small[i] = 1;
            number_small_walkers++;
            number_walkers_new -= 1;
          }

        }  /*If Walkers Weights Too Larger, Then Split Into Multiple*/
        else if ( Cabs(weights[i]) > cns.max_walker_weight) {

          /*Store Walkers Info*/
          walkers_big[number_big_walkers] = i;
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          walkers_copied[number_big_walkers] = number_copies - 1;
          number_big_walkers++;
          number_bad_walkers++;
          number_walkers_new += number_copies - 1;

        }
    }


    /*End Population Control If There Are No Bad Walkers*/
    if ( number_bad_walkers == 0) {
      goto end_population_control_bosons;
    }
    else {    /*If Bad Walkers*/


      /*If Between Counter Values*****************/
      if ( number_walkers_new > cns.min_number_walkers && number_walkers_new < cns.max_number_walkers ) {

         /*Reset the weights of All Walkers That Are Too Big or Small*/
         number_walkers_new = 0;
         for (i=0; i<ist.n_walkers; i++) {

            abs_value = Cabs(weights[i]);
            if ( ( abs_value < cns.min_walker_weight ) || ( abs_value > cns.max_walker_weight )  ) {
                phase = atan(weights[i].imag/weights[i].real);
                weights[i].real = cos(phase);
                weights[i].imag = sin(phase);
            }

            /*Copy Walkers and Matrices*/
            if ( walkers_small[i] != 1 ) {
             new_walker_weights[number_walkers_new] = weights[i];
             new_overlap_bosons[number_walkers_new] = overlap_bosons[i];

               for (j=0; j<ist.n_sites; j++) {
                 np = number_walkers_new * ist.n_sites; 
                 ip = i * ist.n_sites;

                 new_wf_bosons[np+j] = wf_bosons[ip+j];
               }

              number_walkers_new++;
            }
        } /*Run Through Walkers*/

          /*Now Add New Walkers At the End of the LIst*/
          for (i=0; i<number_big_walkers; i++) {

           for (j=0; j< walkers_copied[i]; j++) {
             new_walker_weights[number_walkers_new] = weights[walkers_big[i]];
             new_overlap_bosons[number_walkers_new] = overlap_bosons[walkers_big[i]];

             for (k=0; k<ist.n_sites; k++) {
               np = number_walkers_new*ist.n_sites;
               ip = walkers_big[i]*ist.n_sites;

               new_wf_bosons[np+k] = wf_bosons[ip+k];
             }

             number_walkers_new++;
           }
         }


        /*Last Thing Is To Now Copy Over Walkers*/
        cblas_zcopy(number_walkers_new,new_walker_weights,1,weights,1); 
        cblas_zcopy(number_walkers_new,new_overlap_bosons,1,overlap_bosons,1); 
        cblas_zcopy(number_walkers_new*ist.n_sites,new_wf_bosons,1,wf_bosons,1); 
        //copy_cmat(weights,new_walker_weights,number_walkers_new);
        //copy_cmat(overlap_bosons,new_overlap_bosons,number_walkers_new);
        //copy_cmat(wf_bosons,new_wf_bosons,number_walkers_new*ist.n_sites);

    } /*If Pop Changes Within Bounds*/
    else {   /*If Pop Change Is Outside Bounds*/

        /*Increase Or Decrease Number of Walkers Appropriately*/
        if ( number_walkers_new < cns.min_number_walkers ) {
           population_control_factor=pow(cns.population_control_factor, -1);
        }
        else {
           population_control_factor=cns.population_control_factor;
        }

        for (i=0; i<ist.n_walkers; i++) {
           weights[i]=RCmul(population_control_factor, weights[i]);
         }

         adjusted_weight /= population_control_factor;
         goto begin_population_control_bosons;
    }

  } /*else Nothing to change*/

   end_population_control_bosons: {};

  /*Storing Rescaled Weights*/
  for (i=8; i>-1; i--) {
   weights_rescaled[i+1] = weights_rescaled[i];
  }
  weights_rescaled[0] = adjusted_weight;

return(number_walkers_new);
}

/*********************************************************************************************************************************/

int population_control_bose_fermi(MKL_Complex16 *weights,MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *wf_bosons,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *overlap_bosons,MKL_Complex16 *new_walker_weights,MKL_Complex16 *new_wf_up,MKL_Complex16 *new_wf_down,MKL_Complex16 *new_wf_bosons,MKL_Complex16 *new_overlap_inverse_up,MKL_Complex16 *new_overlap_inverse_down,MKL_Complex16 *new_overlap_up,MKL_Complex16 *new_overlap_down,MKL_Complex16 *new_overlap_bosons,int *walkers_big,int *walkers_small,int *walkers_copied,double *weights_rescaled,int_st ist,cns_st cns) {
    /*Population Control Scheme Based Upon Shiwei's Original Scheme*/
    /*Current Version May Be Sped Up By Copying New Walkers More Efficiently...Test Version Now*/

    int i, j, k;
    int ip, np;
    int number_walkers_new;
    int number_copies, number_bad_walkers, number_big_walkers, number_small_walkers;
    double phase;
    double population_control_factor;
    double abs_value;
    double adjusted_weight;
    long tidum;

    Randomize(); tidum = -random();

    adjusted_weight = 1.0;

    begin_population_control_bose_fermi: {};

    number_walkers_new = ist.n_walkers;
    number_big_walkers = 0;
    number_small_walkers = 0;
    number_bad_walkers = 0;

    izero_vec(walkers_big,ist.n_walkers);
    izero_vec(walkers_small,ist.n_walkers);
    izero_vec(walkers_copied,ist.n_walkers);

    /*Run Through Walkers Figuring Out if Too Big or Too Small*/
    for (i=0; i<ist.n_walkers; i++) {

        /*If Walker Weight Too Small, Then Mark To Get Rid Of*/
        if ( Cabs(weights[i]) < cns.min_walker_weight ) {

          /*Store Walker Info*/
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          number_bad_walkers++;
          if ( number_copies == 0 ) {
            walkers_small[i] = 1;
            number_small_walkers++;
            number_walkers_new -= 1;
          }

        }  /*If Walkers Weights Too Larger, Then Split Into Multiple*/
        else if ( Cabs(weights[i]) > cns.max_walker_weight) {

          /*Store Walkers Info*/
          walkers_big[number_big_walkers] = i;
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          walkers_copied[number_big_walkers] = number_copies - 1;
          number_big_walkers++;
          number_bad_walkers++;
          number_walkers_new += number_copies - 1;

        }
    }


    /*End Population Control If There Are No Bad Walkers*/
    if ( number_bad_walkers == 0) {
      goto end_population_control_bose_fermi;
    }
    else {    /*If Bad Walkers*/

      /*If Between Counter Values*****************/
      if ( number_walkers_new > cns.min_number_walkers && number_walkers_new < cns.max_number_walkers ) {

         /*Reset the weights of All Walkers That Are Too Big or Small*/
         number_walkers_new = 0;
         for (i=0; i<ist.n_walkers; i++) {

            abs_value = Cabs(weights[i]);
            if ( ( abs_value < cns.min_walker_weight ) || ( abs_value > cns.max_walker_weight )  ) {
                phase = atan(weights[i].imag/weights[i].real);
                weights[i].real = cos(phase);
                weights[i].imag = sin(phase);
            }

            /*Copy Walkers and Matrices*/
            if ( walkers_small[i] != 1 ) {
             new_walker_weights[number_walkers_new] = weights[i];
             new_overlap_up[number_walkers_new] = overlap_up[i];
             new_overlap_down[number_walkers_new] = overlap_down[i];
             new_overlap_bosons[number_walkers_new] = overlap_bosons[i]; 

               for (j=0; j<ist.n_sites_n_up; j++) {
                 np = number_walkers_new * ist.n_sites_n_up;
                 ip = i * ist.n_sites_n_up;

                 new_wf_up[np+j] = wf_up[ip+j];
               }

               for (j=0; j<ist.n_sites_n_down; j++) {
                 np = number_walkers_new * ist.n_sites_n_down;
                 ip = i * ist.n_sites_n_down;

                 new_wf_down[np+j] = wf_down[ip+j];
               }
               
               for (j=0; j<ist.n_sites; j++) {
                 np = number_walkers_new * ist.n_sites;
                 ip = i * ist.n_sites;

                 new_wf_bosons[np+j] = wf_bosons[ip+j];
               }

               for (j=0; j<ist.n_up_sq; j++) {
                 np = number_walkers_new * ist.n_up_sq;
                 ip = i * ist.n_up_sq;

                 new_overlap_inverse_up[np+j] = overlap_inverse_up[ip+j];
               }
               for (j=0; j<ist.n_down_sq; j++) {
                 np = number_walkers_new * ist.n_down_sq;
                 ip = i * ist.n_down_sq;

                 new_overlap_inverse_down[np+j] = overlap_inverse_down[ip+j];
               }

              number_walkers_new++;
            }
        } /*Run Through Walkers*/

          /*Now Add New Walkers At the End of the LIst*/
          for (i=0; i<number_big_walkers; i++) {

           for (j=0; j< walkers_copied[i]; j++) {
             new_walker_weights[number_walkers_new] = weights[walkers_big[i]];
             new_overlap_up[number_walkers_new] = overlap_up[walkers_big[i]];
             new_overlap_down[number_walkers_new] = overlap_down[walkers_big[i]];
             new_overlap_bosons[number_walkers_new] = overlap_bosons[walkers_big[i]]; 

             for (k=0; k<ist.n_sites_n_up; k++) {
               np = number_walkers_new*ist.n_sites_n_up;
               ip = walkers_big[i]*ist.n_sites_n_up;

               new_wf_up[np+k] = wf_up[ip+k];
             }
             for (k=0; k<ist.n_sites_n_down; k++) {
               np = number_walkers_new*ist.n_sites_n_down;
               ip = walkers_big[i]*ist.n_sites_n_down;

               new_wf_down[np+k] = wf_down[ip+k];
             }
             for (k=0; k<ist.n_sites; k++) {
               np = number_walkers_new*ist.n_sites;
               ip = walkers_big[i]*ist.n_sites;

               new_wf_bosons[np+k] = wf_bosons[ip+k];
             }
             for (k=0; k<ist.n_up_sq; k++) {
               np = number_walkers_new*ist.n_up_sq;
               ip = walkers_big[i]*ist.n_up_sq;

               new_overlap_inverse_up[np+k] = overlap_inverse_up[ip+k];
             }
             for (k=0; k<ist.n_down_sq; k++) {
               np = number_walkers_new*ist.n_down_sq;
               ip = walkers_big[i]*ist.n_down_sq;

               new_overlap_inverse_down[np+k] = overlap_inverse_down[ip+k];
             }


             number_walkers_new++;
           }
         }


        /*Last Thing Is To Now Copy Over Walkers*/
        cblas_zcopy(number_walkers_new,new_walker_weights,1,weights,1); 
        cblas_zcopy(number_walkers_new,new_overlap_up,1,overlap_up,1); 
        cblas_zcopy(number_walkers_new,new_overlap_down,1,overlap_down,1); 
        cblas_zcopy(number_walkers_new,new_overlap_bosons,1,overlap_bosons,1); 
        //copy_cmat(weights,new_walker_weights,number_walkers_new);
        //copy_cmat(overlap_up,new_overlap_up,number_walkers_new);
        //copy_cmat(overlap_down,new_overlap_down,number_walkers_new);
        //copy_cmat(overlap_bosons,new_overlap_bosons,number_walkers_new); 

        cblas_zcopy(number_walkers_new*ist.n_sites_n_up,new_wf_up,1,wf_up,1); 
        cblas_zcopy(number_walkers_new*ist.n_sites_n_down,new_wf_down,1,wf_down,1); 
        cblas_zcopy(number_walkers_new*ist.n_sites,new_wf_bosons,1,wf_bosons,1); 
        cblas_zcopy(number_walkers_new*ist.n_up_sq,new_overlap_inverse_up,1,overlap_inverse_up,1); 
        cblas_zcopy(number_walkers_new*ist.n_down_sq,new_overlap_inverse_down,1,overlap_inverse_down,1); 
        /*copy_cmat(wf_up,new_wf_up,number_walkers_new*ist.n_sites_n_up);
        copy_cmat(wf_down,new_wf_down,number_walkers_new*ist.n_sites_n_down);
        copy_cmat(wf_bosons,new_wf_bosons,number_walkers_new*ist.n_sites); 
        copy_cmat(overlap_inverse_up,new_overlap_inverse_up,number_walkers_new*ist.n_up_sq);
        copy_cmat(overlap_inverse_down,new_overlap_inverse_down,number_walkers_new*ist.n_down_sq);
        */

    } /*If Pop Changes Within Bounds*/
    else {   /*If Pop Change Is Outside Bounds*/

        /*Increase Or Decrease Number of Walkers Appropriately*/
        if ( number_walkers_new < cns.min_number_walkers ) {
           population_control_factor=pow(cns.population_control_factor, -1);
        }
        else {
           population_control_factor=cns.population_control_factor;
        }

        for (i=0; i<ist.n_walkers; i++) {
           weights[i]=RCmul(population_control_factor, weights[i]);
         }

         adjusted_weight /= population_control_factor;

         goto begin_population_control_bose_fermi;
    }

  } /*else Nothing to change*/

   end_population_control_bose_fermi: {};

  /*Storing Rescaled Weights*/
  for (i=8; i>-1; i--) {
   weights_rescaled[i+1] = weights_rescaled[i];
  }
  weights_rescaled[0] = adjusted_weight;

return(number_walkers_new);
}

/*************************************************************************************************************/

int population_control_bose_fermi_spin_polarized(MKL_Complex16 *weights,MKL_Complex16 *wf_up,MKL_Complex16 *wf_bosons,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_bosons,MKL_Complex16 *new_walker_weights,MKL_Complex16 *new_wf_up,MKL_Complex16 *new_wf_bosons,MKL_Complex16 *new_overlap_inverse_up,MKL_Complex16 *new_overlap_up,MKL_Complex16 *new_overlap_bosons,int *walkers_big,int *walkers_small,int *walkers_copied,double *weights_rescaled,int_st ist,cns_st cns,int step) {
    /*Population Control Scheme Based Upon Shiwei's Original Scheme*/
    /*Current Version May Be Sped Up By Copying New Walkers More Efficiently...Test Version Now*/

    int i, j, k;
    int ip, np;
    int number_walkers_new;
    int number_copies, number_bad_walkers, number_big_walkers, number_small_walkers;
    double phase;
    double population_control_factor;
    double abs_value;
    double adjusted_weight;
    long tidum;

    Randomize(); tidum = -random();

    adjusted_weight = 1.0;

    begin_population_control_bose_fermi_spin_polarized: {};

    number_walkers_new = ist.n_walkers;
    number_big_walkers = 0;
    number_small_walkers = 0;
    number_bad_walkers = 0;

    izero_vec(walkers_big,ist.n_walkers);
    izero_vec(walkers_small,ist.n_walkers);
    izero_vec(walkers_copied,ist.n_walkers);

    /*Run Through Walkers Figuring Out if Too Big or Too Small*/
    for (i=0; i<ist.n_walkers; i++) {

        /*If Walker Weight Too Small, Then Mark To Get Rid Of*/
        if ( Cabs(weights[i]) < cns.min_walker_weight ) {

          /*Store Walker Info*/
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          number_bad_walkers++;
          if ( number_copies == 0 ) {
            walkers_small[i] = 1;
            number_small_walkers++;
            number_walkers_new -= 1;
          }
        }  /*If Walkers Weights Too Larger, Then Split Into Multiple*/
        else if ( Cabs(weights[i]) > cns.max_walker_weight) {

          /*Store Walkers Info*/
          walkers_big[number_big_walkers] = i;
          number_copies = (int)(ran_nrc(&tidum)+Cabs(weights[i]));
          walkers_copied[number_big_walkers] = number_copies - 1;
          number_big_walkers++;
          number_bad_walkers++;
          number_walkers_new += number_copies - 1;

        }
    }

    /*End Population Control If There Are No Bad Walkers*/
    if ( number_bad_walkers == 0) {
      goto end_population_control_bose_fermi_spin_polarized;
    }
    else {    /*If Bad Walkers*/

      /*If Between Counter Values*****************/
      if ( number_walkers_new > cns.min_number_walkers && number_walkers_new < cns.max_number_walkers ) {

         /*Reset the weights of All Walkers That Are Too Big or Small*/
         number_walkers_new = 0;
         for (i=0; i<ist.n_walkers; i++) {

            abs_value = Cabs(weights[i]);
            if ( ( abs_value < cns.min_walker_weight ) || ( abs_value > cns.max_walker_weight )  ) {
                phase = atan(weights[i].imag/weights[i].real);
                weights[i].real = cos(phase);
                weights[i].imag = sin(phase);
            }

            /*Copy Walkers and Matrices*/
            if ( walkers_small[i] != 1 ) {
             new_walker_weights[number_walkers_new] = weights[i];
             new_overlap_bosons[number_walkers_new] = overlap_bosons[i];
             new_overlap_up[number_walkers_new] = overlap_up[i];

               for (j=0; j<ist.n_sites; j++) {
                 np = number_walkers_new * ist.n_sites;
                 ip = i * ist.n_sites;
                 new_wf_bosons[np+j] = wf_bosons[ip+j];
               }
               for (j=0; j<ist.n_sites_n_up; j++) {
                 np = number_walkers_new * ist.n_sites_n_up;
                 ip = i * ist.n_sites_n_up;

                 new_wf_up[np+j] = wf_up[ip+j];
               }
               for (j=0; j<ist.n_up_sq; j++) {
                 np = number_walkers_new * ist.n_up_sq;
                 ip = i * ist.n_up_sq;

                 new_overlap_inverse_up[np+j] = overlap_inverse_up[ip+j];
               } 

              number_walkers_new++;
            }
        } /*Run Through Walkers*/

          /*Now Add New Walkers At the End of the LIst*/
          for (i=0; i<number_big_walkers; i++) {

           for (j=0; j< walkers_copied[i]; j++) {
             new_walker_weights[number_walkers_new] = weights[walkers_big[i]];
             new_overlap_bosons[number_walkers_new] = overlap_bosons[walkers_big[i]];
             new_overlap_up[number_walkers_new] = overlap_up[walkers_big[i]];

             for (k=0; k<ist.n_sites; k++) {
               np = number_walkers_new*ist.n_sites;
               ip = walkers_big[i]*ist.n_sites;

               new_wf_bosons[np+k] = wf_bosons[ip+k];
             }
             for (k=0; k<ist.n_sites_n_up; k++) {
               np = number_walkers_new*ist.n_sites_n_up;
               ip = walkers_big[i]*ist.n_sites_n_up;

               new_wf_up[np+k] = wf_up[ip+k];
             }
             for (k=0; k<ist.n_up_sq; k++) {
               np = number_walkers_new*ist.n_up_sq;
               ip = walkers_big[i]*ist.n_up_sq;

               new_overlap_inverse_up[np+k] = overlap_inverse_up[ip+k];
             } 

             number_walkers_new++;
           }
         }
   
        /*Last Thing Is To Now Copy Over Walkers*/
        cblas_zcopy(number_walkers_new,new_walker_weights,1,weights,1); 
        cblas_zcopy(number_walkers_new,new_overlap_bosons,1,overlap_bosons,1); 
        cblas_zcopy(number_walkers_new,new_overlap_up,1,overlap_up,1);   
        //copy_cmat(weights,new_walker_weights,number_walkers_new);
        //copy_cmat(overlap_bosons,new_overlap_bosons,number_walkers_new);
        //copy_cmat(overlap_up,new_overlap_up,number_walkers_new);

        cblas_zcopy(number_walkers_new*ist.n_sites,new_wf_bosons,1,wf_bosons,1); 
        cblas_zcopy(number_walkers_new*ist.n_sites_n_up,new_wf_up,1,wf_up,1); 
        cblas_zcopy(number_walkers_new*ist.n_up_sq,new_overlap_inverse_up,1,overlap_inverse_up,1); 
        //copy_cmat(wf_bosons,new_wf_bosons,number_walkers_new*ist.n_sites);
        //copy_cmat(wf_up,new_wf_up,number_walkers_new*ist.n_sites_n_up); 
        //copy_cmat(overlap_inverse_up,new_overlap_inverse_up,number_walkers_new*ist.n_up_sq);

    } /*If Pop Changes Within Bounds*/
    else {   /*If Pop Change Is Outside Bounds*/

        /*Increase Or Decrease Number of Walkers Appropriately*/
        if ( number_walkers_new < cns.min_number_walkers ) {
           population_control_factor=pow(cns.population_control_factor, -1);
        }
        else {
           population_control_factor=cns.population_control_factor;
        }

        for (i=0; i<ist.n_walkers; i++) {
           weights[i]=RCmul(population_control_factor, weights[i]);
         }
         adjusted_weight /= population_control_factor;
         goto begin_population_control_bose_fermi_spin_polarized;
    }

  } /*else Nothing to change*/

   end_population_control_bose_fermi_spin_polarized: {};

  /*Storing Rescaled Weights*/
  for (i=8; i>-1; i--) {
   weights_rescaled[i+1] = weights_rescaled[i];
  }
  weights_rescaled[0] = adjusted_weight;

return(number_walkers_new);
}

/********************************************************************************************************/

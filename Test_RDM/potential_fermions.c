#include "afqmc.h"

/************************************************************************************/

double find_gamma(double U,double dtau) { 

   /*Finds Gamma Given U*********************************/
   int i;
   double gamma;  
   double gamma_interval = 0.00000001; 
   double cosh_error = .000001; 
   double exponential_value = 2.0 * exp(dtau * U * .5);

   for (i=0; i<100000000000; i++) {
      gamma = i * gamma_interval;
      if ( fabs(exp(gamma) + exp(-1.0 * gamma) - exponential_value) < cosh_error ) {
         return(gamma); 
      }
   }    

return(0); 
} 

/************************************************************************************/

void propagate_forwards_potential_discrete_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   int flag_neg_overlap_up, flag_neg_overlap_down; 
   MKL_Complex16 probability_field_up, probability_field_down; 
   double spin_up = 1.0, spin_down = -1.0; 
   double field_up = 1.0, field_down = -1.0; 
   double field_selected; 
   double random_number; 
   MKL_Complex16 total_probability; 
   long tidum; 

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {
    
        /*Run Through Sites*/
        for (sites=0; sites<ist.n_sites; sites++) { 
 
          /*Check WEights*/
          if ( Cabs(weights[walker]) != 0 ) { 
 
            /*Ensure That Overlap Flags Are Positive*/
            flag_neg_overlap_up = flag_neg_overlap_down = 0;  

            /*Determine Probabilities for Each Scenario*/
            probability_field_up = determine_overlap_ratio_discrete_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,spin_up,field_up,sites);                 
            probability_field_up = Cmul(probability_field_up, determine_overlap_ratio_discrete_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,spin_down,field_up,sites));
            probability_field_up = RCmul(.5, probability_field_up); 


            probability_field_down = determine_overlap_ratio_discrete_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,spin_up,field_down,sites);
            probability_field_down = Cmul(probability_field_down, determine_overlap_ratio_discrete_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,spin_down,field_down,sites));
            probability_field_down = RCmul(.5, probability_field_down);

            /*If Constrained Path Approximation Applied - With Mirror Correction To Check If Next Timeslice Will Cause Overlap and to Reduce Weights If Next Step Does*/
            /*if ( ist.flag_cp == 1 ) {
              if ( probability_field_up < 0 ) {
                 flag_neg_overlap_up = 1; 
                 probability_field_up = 0;
                 weights[walker] /= (1-probability_field_up*2.0);  
              }
              if ( probability_field_down < 0 ) {
                 flag_neg_overlap_down = 1;
                 probability_field_down = 0;  
                 weights[walker] /= (1-probability_field_down*2.0); 
              }
            }*/

            total_probability = Cadd(probability_field_down, probability_field_up); 
 
            /*Now Propagate If Weight Is Not Negative************************************************/
            if ( total_probability.real <= 0 && ist.flag_cp == 1 ) {
                 weights[walker].real=0;
                 weights[walker].imag=0;  
             }
             else {  
 
                /*Perform Monte Carlo Decision for Which Field to Pick*/
                random_number = ran_nrc(&tidum);  
                if ( random_number < Cabs(Cdiv(probability_field_up, total_probability)) ) { /*If Spin Up Selected*/
                  weights[walker] = Cmul(weights[walker], total_probability); 
  
                  propagate_wave_functions_discrete_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],ist,cns,field_up,sites); 
                  update_overlaps_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],&overlap_up[walker],&overlap_down[walker],ist); 

  
                  field_selected = field_up; 
                }
                else { /*If Spin Down Selected*/ 
                  weights[walker] = Cmul(weights[walker], total_probability);  
  
                  propagate_wave_functions_discrete_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],ist,cns,field_down,sites);
                  update_overlaps_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],&overlap_up[walker],&overlap_down[walker],ist);
 
                  field_selected = field_down; 

               } 
            } /*If Not Traversing Boundary*/
            /*Done Propagating********************************************************************/

            /*Apply Mirror Correction Again Afterwards to Reduce Weight If Next Time Slice Will Result in Crossing*/
            /*if ( ist.flag_cp == 1 ) {

               if ( flag_neg_overlap_up == 0 && flag_neg_overlap_down == 0 ) {
                  if ( field_selected == field_up ) {
                      probability_field_up = determine_overlap_ratio(&wf_up[walker*n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*n_up_sq],ist,cns,spin_up,field_up,sites);    
                      probability_field_up *= determine_overlap_ratio(&wf_down[walker*n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*n_down_sq],ist,cns,spin_down,field_up,sites);
                      if ( probability_field_up < 0 ) {
                         weights[walker] /= (1-probability_field_up);  
                      }
                  }
                  else {
                      probability_field_down = determine_overlap_ratio(&wf_up[walker*n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*n_up_sq],ist,cns,spin_up,field_down,sites);
                      probability_field_down *= determine_overlap_ratio(&wf_down[walker*n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*n_down_sq],ist,cns,spin_down,field_down,sites);
                      if ( probability_field_down < 0 ) {
                         weights[walker] /= (1-probability_field_down); 
                      } 
                  } 
               }
                 
            }*/ /*Check Mirror Correction********************/  


         } /*Check Walker Weights*/

       } /*Run Through Sites*/
   } /*Run Through Walkers*/

return; 
}

/***************************************************************************************/

void propagate_forwards_potential_continuous_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant; 
   MKL_Complex16 ratio_determinant; 
   MKL_Complex16 exponential_field_up, exponential_field_down; 
   long tidum;

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Run Through Sites*/
        for (sites=0; sites<ist.n_sites; sites++) {

          /*Check WEights*/
          if ( Cabs(weights[walker]) != 0 ) {

            /*Select Random Field*/
            field = gasdev(&tidum);
   
            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * field * cns.propagation_exponent_fermions.real) * cns.U_exponent_fermions; 
            complex_constant = -1.0 * field * cns.propagation_exponent_fermions.imag;    
            exponential_field_up.real = real_constant * cos(complex_constant); 
            exponential_field_up.imag = real_constant * sin(complex_constant); 

            real_constant = exp(field * cns.propagation_exponent_fermions.real) * cns.U_exponent_fermions; 
            complex_constant = field * cns.propagation_exponent_fermions.imag; 
            exponential_field_down.real = real_constant * cos(complex_constant); 
            exponential_field_down.imag = real_constant * sin(complex_constant); 
           
            /*Determine Ratio of Determinants*/
            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_up,ist.n_up,sites);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,exponential_field_down,ist.n_down,sites));

             /*Propagate the WF and Get New Overlaps*/
             propagate_wave_functions_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],ist,cns,exponential_field_up,exponential_field_down,sites);
             update_overlaps_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],&overlap_up[walker],&overlap_down[walker],ist);

            /*Multiply Into Weights*/
            weights[walker] = Cmul(weights[walker], ratio_determinant); 

         } /*Check Walker Weights*/

       } /*Run Through Sites*/
   } /*Run Through Walkers*/

return;
}

/***************************************************************************************/

void propagate_forwards_potential_continuous_meanfield_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant;
   MKL_Complex16 ratio_determinant;
   MKL_Complex16 exponential_field_up, exponential_field_down;
   MKL_Complex16 mean_field_exponential; 
   long tidum;

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Run Through Sites*/
        for (sites=0; sites<ist.n_sites; sites++) {

          /*Check WEights*/
          if ( Cabs(weights[walker]) != 0 ) {

            /*Select Random Field*/
            field = gasdev(&tidum);

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * field * cns.propagation_exponent_fermions.real) * cns.U_exponent_fermions;
            complex_constant = -1.0 * field * cns.propagation_exponent_fermions.imag;
            exponential_field_up.real = real_constant * cos(complex_constant);
            exponential_field_up.imag = real_constant * sin(complex_constant);

            real_constant = exp(field * cns.propagation_exponent_fermions.real) * cns.U_exponent_fermions;
            complex_constant = field * cns.propagation_exponent_fermions.imag;
            exponential_field_down.real = real_constant * cos(complex_constant);
            exponential_field_down.imag = real_constant * sin(complex_constant);

            /*Determine Mean Field Exponential*/
            real_constant = exp(field * cns.propagation_exponent_fermions.real * cns.total_trial_density); 
            complex_constant = field * cns.propagation_exponent_fermions.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant); 
            mean_field_exponential.imag = real_constant * sin(complex_constant);    

            /*Modify Exponential Constant By Mean Field Term*/
            exponential_field_up = RCmul(cns.mean_field_exponential_fermions, exponential_field_up);
            exponential_field_down = RCmul(cns.mean_field_exponential_fermions, exponential_field_down);

            /*Determine Ratio of Determinants*/
            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_up,ist.n_up,sites);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,exponential_field_down,ist.n_down,sites));

            /*Multiply Mean Field Exponential Into Determinant*/
            ratio_determinant = Cmul(ratio_determinant, mean_field_exponential);   
           
            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],ist,cns,exponential_field_up,exponential_field_down,sites);
            update_overlaps_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],&overlap_up[walker],&overlap_down[walker],ist);
            
            /*Multiply Into Weights*/
            weights[walker] = Cmul(weights[walker], ratio_determinant); 

          } /*Check Walker Weights*/
        }/*Run Through Sites*/
    }

return; 
}

/********************************************************************************************************************************/

void propagate_forwards_potential_continuous_phaseless_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant;
   double phase, abs_value, projected_ratio_determinant; 
   MKL_Complex16 ratio_determinant; 
   MKL_Complex16 exponential_field_up, exponential_field_down; 
   long tidum;

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Run Through Sites*/
        for (sites=0; sites<ist.n_sites; sites++) {

          /*Check WEights*/
          if ( Cabs(weights[walker]) != 0 ) {

            /*Select Random Field*/
            field = gasdev(&tidum);

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * field * cns.propagation_exponent_fermions.real) * cns.U_exponent_fermions;
            complex_constant = -1.0 * field * cns.propagation_exponent_fermions.imag;
            exponential_field_up.real = real_constant * cos(complex_constant);
            exponential_field_up.imag = real_constant * sin(complex_constant);

            real_constant = exp(field * cns.propagation_exponent_fermions.real) * cns.U_exponent_fermions;
            complex_constant = field * cns.propagation_exponent_fermions.imag;
            exponential_field_down.real = real_constant * cos(complex_constant);
            exponential_field_down.imag = real_constant * sin(complex_constant);

            /*Determine Ratio of Determinants*/
            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_up,ist.n_up,sites);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,exponential_field_down,ist.n_down,sites));

             /*Propagate the WF and Get New Overlaps*/
             propagate_wave_functions_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],ist,cns,exponential_field_up,exponential_field_down,sites);
             update_overlaps_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],&overlap_up[walker],&overlap_down[walker],ist);

            /*Perform the Phaseless*/
            phase = atan(ratio_determinant.imag/ratio_determinant.real);
            abs_value = Cabs(ratio_determinant);
            if ( fabs(phase) < .5 * 3.1415 && ratio_determinant.real > 0 ) {
              projected_ratio_determinant = abs_value * cos(phase);
            }
            else {
              projected_ratio_determinant = 0.0;
            }

            /*Multiply Into Weights*/
            weights[walker] = RCmul(projected_ratio_determinant, weights[walker]); 
     
         } /*Check Walker Weights*/

       } /*Run Through Sites*/
   } /*Run Through Walkers*/

return;
}

/******************************************************************************************************************/

void propagate_forwards_potential_continuous_allshifts_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant;
   MKL_Complex16 shifted_field, shift;
   MKL_Complex16 shifted_field_exponential, shifted_field_exponent; 
   MKL_Complex16 prod_prop_shift;  
   MKL_Complex16 ratio_determinant;
   MKL_Complex16 exponential_field_up, exponential_field_down;
   MKL_Complex16 mean_field_exponential;
   long tidum;

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Run Through Sites*/
        for (sites=0; sites<ist.n_sites; sites++) {

          /*Check WEights*/
          if ( Cabs(weights[walker]) != 0 ) {

            /*Select Random Field*/
            field = gasdev(&tidum);
         
            /*Get Field Shift*/
            shift = compute_shift_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],&trial_wf_up[0],&trial_wf_down[0],&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],ist,cns);   
            shift.real = shift.imag = 0.0; 
 
            shifted_field.real = field - shift.real; 
            shifted_field.imag = -shift.imag;  

            shifted_field_exponent = Csub(RCmul(field, shift), RCmul(.5, Cmul(shift, shift)));
            shifted_field_exponential.real = exp(shifted_field_exponent.real) * cos(shifted_field_exponent.imag); 
            shifted_field_exponential.imag = exp(shifted_field_exponent.real) * sin(shifted_field_exponent.imag);  

            prod_prop_shift = Cmul(shifted_field, cns.propagation_exponent_fermions); 

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * prod_prop_shift.real) * cns.U_exponent_fermions;
            complex_constant = -1.0 * prod_prop_shift.imag;
            exponential_field_up.real = real_constant * cos(complex_constant);
            exponential_field_up.imag = real_constant * sin(complex_constant);

            real_constant = exp(prod_prop_shift.real) * cns.U_exponent_fermions;
            complex_constant = prod_prop_shift.imag;
            exponential_field_down.real = real_constant * cos(complex_constant);
            exponential_field_down.imag = real_constant * sin(complex_constant);

            /*Modify Exponential Constant By Mean Field Term*/
            exponential_field_up = RCmul(cns.mean_field_exponential_fermions, exponential_field_up);
            exponential_field_down = RCmul(cns.mean_field_exponential_fermions, exponential_field_down);

            /*Determine Mean Field Exponential*/
            real_constant = exp(prod_prop_shift.real * cns.total_trial_density);
            complex_constant = prod_prop_shift.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant);
            mean_field_exponential.imag = real_constant * sin(complex_constant);
            
            /*Determine Ratio of Determinants*/
            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_up,ist.n_up,sites);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,exponential_field_down,ist.n_down,sites));

            /*Multiply Mean Field Exponential Into Determinant*/
            ratio_determinant = Cmul(ratio_determinant, Cmul(mean_field_exponential, shifted_field_exponential));

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],ist,cns,exponential_field_up,exponential_field_down,sites);
            update_overlaps_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],&overlap_up[walker],&overlap_down[walker],ist);

            /*Multiply Into Weights*/
            weights[walker] = Cmul(weights[walker], ratio_determinant);

          } /*Check Walker Weights*/
        }/*Run Through Sites*/
    }

return;
}

/*******************************************************************************************************************************/

void propagate_forwards_potential_continuous_phaseless_meanfield_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant;
   double phase, abs_value, projected_ratio_determinant;
   MKL_Complex16 ratio_determinant; 
   MKL_Complex16 exponential_field_up, exponential_field_down;
   MKL_Complex16 mean_field_exponential; 
   long tidum;

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Run Through Sites*/
        for (sites=0; sites<ist.n_sites; sites++) {

          /*Check WEights*/
          if ( Cabs(weights[walker]) != 0 ) {

            /*Select Random Field*/
            field = gasdev(&tidum);

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * field * cns.propagation_exponent_fermions.real) * cns.U_exponent_fermions;
            complex_constant = -1.0 * field * cns.propagation_exponent_fermions.imag;
            exponential_field_up.real = real_constant * cos(complex_constant);
            exponential_field_up.imag = real_constant * sin(complex_constant);

            real_constant = exp(field * cns.propagation_exponent_fermions.real) * cns.U_exponent_fermions;
            complex_constant = field * cns.propagation_exponent_fermions.imag;
            exponential_field_down.real = real_constant * cos(complex_constant);
            exponential_field_down.imag = real_constant * sin(complex_constant);

            /*Determine Mean Field Exponential*/
            real_constant = exp(field * cns.propagation_exponent_fermions.real * cns.total_trial_density);
            complex_constant = field * cns.propagation_exponent_fermions.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant);
            mean_field_exponential.imag = real_constant * sin(complex_constant);

            /*Modify Exponential Constant By Mean Field Term*/
            exponential_field_up = RCmul(cns.mean_field_exponential_fermions, exponential_field_up);
            exponential_field_down = RCmul(cns.mean_field_exponential_fermions, exponential_field_down);

            /*Determine Ratio of Determinants*/
            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_up,ist.n_up,sites);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,exponential_field_down,ist.n_down,sites));

            /*Multiply Mean Field Exponential Into Determinant*/
            ratio_determinant = Cmul(ratio_determinant, mean_field_exponential);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],ist,cns,exponential_field_up,exponential_field_down,sites);
            update_overlaps_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],&overlap_up[walker],&overlap_down[walker],ist);

           /*Perform the Phaseless*/
           phase = atan(ratio_determinant.imag/ratio_determinant.real);
           abs_value = Cabs(ratio_determinant);
           if ( fabs(phase) < .5 * 3.1415 && ratio_determinant.real > 0 ) {
             projected_ratio_determinant = abs_value * cos(phase);
           }
           else {
             projected_ratio_determinant = 0.0;
           }

           /*Multiply Into Weights*/
           weights[walker] = RCmul(projected_ratio_determinant, weights[walker]);

     } /*Check Walker Weights*/

    } /*Run Through Sites*/
   } /*Run Through Walkers*/

return;
}

/**********************************************************************************************************************************/

void propagate_forwards_potential_continuous_phaseless_allshifts_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant;
   double phase, abs_value, projected_ratio_determinant; 
   MKL_Complex16 shifted_field, shift;
   MKL_Complex16 shifted_field_exponential, shifted_field_exponent;
   MKL_Complex16 prod_prop_shift;
   MKL_Complex16 ratio_determinant;
   MKL_Complex16 exponential_field_up, exponential_field_down;
   MKL_Complex16 mean_field_exponential;
   long tidum;

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Run Through Sites*/
        for (sites=0; sites<ist.n_sites; sites++) {

          /*Check WEights*/
          if ( Cabs(weights[walker]) != 0 ) {

            /*Select Random Field*/
            field = gasdev(&tidum);

            /*Get Field Shift*/
            //shift = compute_shift_fermions(&wf_up[walker*n_sites_n_up],&wf_down[walker*n_sites_n_down],&trial_wf_up[0],&trial_wf_down[0],&overlap_inverse_up[walker*n_up_sq],&overlap_inverse_down[walker*n_down_sq],ist,cns);
            shift.real = shift.imag = 0.0;

            shifted_field.real = field - shift.real;
            shifted_field.imag = -shift.imag;

            shifted_field_exponent = Csub(RCmul(field, shift), RCmul(.5, Cmul(shift, shift)));
            shifted_field_exponential.real = exp(shifted_field_exponent.real) * cos(shifted_field_exponent.imag);
            shifted_field_exponential.imag = exp(shifted_field_exponent.real) * sin(shifted_field_exponent.imag);

            prod_prop_shift = Cmul(shifted_field, cns.propagation_exponent_fermions);

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * prod_prop_shift.real) * cns.U_exponent_fermions;
            complex_constant = -1.0 * prod_prop_shift.imag;
            exponential_field_up.real = real_constant * cos(complex_constant);
            exponential_field_up.imag = real_constant * sin(complex_constant);

            real_constant = exp(prod_prop_shift.real) * cns.U_exponent_fermions;
            complex_constant = prod_prop_shift.imag;
            exponential_field_down.real = real_constant * cos(complex_constant);
            exponential_field_down.imag = real_constant * sin(complex_constant);

            /*Modify Exponential Constant By Mean Field Term*/
            exponential_field_up = RCmul(cns.mean_field_exponential_fermions, exponential_field_up);
            exponential_field_down = RCmul(cns.mean_field_exponential_fermions, exponential_field_down);

            /*Determine Mean Field Exponential*/
            real_constant = exp(prod_prop_shift.real * cns.total_trial_density);
            complex_constant = prod_prop_shift.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant);
            mean_field_exponential.imag = real_constant * sin(complex_constant);

            /*Determine Ratio of Determinants*/
            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_up,ist.n_up,sites);
            ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,exponential_field_down,ist.n_down,sites));

            /*Multiply Mean Field Exponential Into Determinant*/
            ratio_determinant = Cmul(ratio_determinant, Cmul(mean_field_exponential, shifted_field_exponential));


     /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],ist,cns,exponential_field_up,exponential_field_down,sites);
 
            update_overlaps_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],&overlap_up[walker],&overlap_down[walker],ist);


            /*Perform the Phaseless*/
            phase = atan(ratio_determinant.imag/ratio_determinant.real);
            abs_value = Cabs(ratio_determinant);
            if ( fabs(phase) < .5 * 3.1415 && ratio_determinant.real > 0 ) {
              projected_ratio_determinant = abs_value * cos(phase);
            }
            else {
              projected_ratio_determinant = 0.0;
            }

            /*Multiply Into Weights*/
            weights[walker] = RCmul(projected_ratio_determinant, weights[walker]); 

          } /*Check Walker Weights*/
        }/*Run Through Sites*/
    }

return;
}

/********************************************************************************************************************************/
 
MKL_Complex16 determine_overlap_ratio_discrete_fermions(MKL_Complex16 *wf,MKL_Complex16 *trial_wf,MKL_Complex16 *overlap_inverse,int_st ist,cns_st cns,double spin_up,double field_up,int site) {

      /*Determines Overlap Ratio*/
      MKL_Complex16 overlap_ratio;  
      MKL_Complex16 *stored_product1, *stored_product2; 
      MKL_Complex16 *One, *Zero;

      One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
      Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
      One[0].real = 1.0; One[0].imag = 0.0;
      Zero[0].real = 0.0; Zero[0].imag = 0.0;

      /*If An Up-Spin*/
      if (spin_up == 1.0 ) {
        stored_product1=(MKL_Complex16 *)calloc(ist.n_sites*ist.n_up,sizeof(MKL_Complex16)); 
        stored_product2=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));  

        /*Get Matrix For Green's Function*/
        cmat_transpose_cmat(overlap_inverse,trial_wf,stored_product1,ist.n_up,ist.n_sites,ist.n_up);
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_sites,ist.n_up,ist.n_sites,One,wf,ist.n_sites,stored_product1,ist.n_up,Zero,stored_product2,ist.n_up);
        //cmat_cmat(wf,stored_product1,stored_product2,ist.n_sites,ist.n_up,ist.n_sites);

        if ( field_up == 1.0 ) {
          overlap_ratio = RCmul(cns.U_exponent_fermions * cns.factor_spin_up_field_up - 1.0,stored_product2[site * ist.n_sites + site]); 
          overlap_ratio.real += 1.0; 
        }
        else { 
          overlap_ratio = RCmul(cns.U_exponent_fermions * cns.factor_spin_up_field_down - 1.0, stored_product2[site * ist.n_sites + site]);
          overlap_ratio.real += 1.0; 
        }  

      }
      else { /*If Spin-Down*/
        stored_product1=(MKL_Complex16 *)calloc(ist.n_sites*ist.n_down,sizeof(MKL_Complex16));
        stored_product2=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16)); 

        /*Get Matrix For Green's Function*/
        cmat_transpose_cmat(overlap_inverse,trial_wf,stored_product1,ist.n_down,ist.n_sites,ist.n_down);
        cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_sites,ist.n_down,ist.n_sites,One,wf,ist.n_sites,stored_product1,ist.n_down,Zero,stored_product2,ist.n_down);
        //cmat_cmat(wf,stored_product1,stored_product2,ist.n_sites,ist.n_down,ist.n_sites);

        if ( field_up == 1.0 ) {
          overlap_ratio = RCmul(cns.U_exponent_fermions * cns.factor_spin_up_field_down - 1.0,stored_product2[site * ist.n_sites + site]);
          overlap_ratio.real += 1.0; 
        }
        else {
           overlap_ratio = RCmul(cns.U_exponent_fermions * cns.factor_spin_up_field_up - 1.0,stored_product2[site * ist.n_sites + site]);
           overlap_ratio.real += 1.0; 
        }  

     }

free(stored_product1); 
free(stored_product2); 
free(One); 
free(Zero); 
return(overlap_ratio);  
}

/*****************************************************************************************/

MKL_Complex16 determine_overlap_ratio_continuous_fermions(MKL_Complex16 *wf,MKL_Complex16 *trial_wf,MKL_Complex16 *overlap_inverse,int_st ist,cns_st cns,MKL_Complex16 exponential,int size,int site) {

      /*Determines Overlap Ratio*/
      MKL_Complex16 One; One.real = 1.0; One.imag = 0.0; 
      MKL_Complex16 overlap_ratio;
      MKL_Complex16 *stored_product1, *stored_product2;
      MKL_Complex16 *One_2, *Zero;

      One_2 = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
      Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
      One_2[0].real = 1.0; One_2[0].imag = 0.0;
      Zero[0].real = 0.0; Zero[0].imag = 0.0; 
 
      /*If An Up-Spin*/
      stored_product1=(MKL_Complex16 *)calloc(ist.n_sites*size,sizeof(MKL_Complex16));
      stored_product2=(MKL_Complex16 *)calloc(ist.n_sites_sq,sizeof(MKL_Complex16));

      /*Get Matrix For Green's Function*/
      cmat_transpose_cmat(overlap_inverse,trial_wf,stored_product1,size,ist.n_sites,size);
      cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,ist.n_sites,size,ist.n_sites,One_2,wf,ist.n_sites,stored_product1,size,Zero,stored_product2,size); 
      //cmat_cmat(wf,stored_product1,stored_product2,ist.n_sites,size,ist.n_sites);
 
      /*Multiply Factor In*/
      overlap_ratio = Cmul(Csub(exponential, One), stored_product2[site * ist.n_sites + site]); 
      overlap_ratio = Cadd(overlap_ratio, One);   

free(stored_product1);
free(stored_product2);
free(One_2); 
free(Zero); 
return(overlap_ratio);
}

/****************************************************************************************/  

void propagate_wave_functions_discrete_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,int_st ist,cns_st cns,double field_up,int sites) {

     /*Propagates the Wave Function*/
     int i; 

     /*For Up Wave Function*/ 
     for (i=0; i<ist.n_up; i++) {
        if ( field_up == 1.0 ) {
          wf_up[sites*ist.n_up+i] = RCmul(cns.U_exponent_fermions * cns.factor_spin_up_field_up, wf_up[sites*ist.n_up+i]); 
        }
        else {
          wf_up[sites*ist.n_up+i] = RCmul(cns.U_exponent_fermions * cns.factor_spin_up_field_down, wf_up[sites*ist.n_up+i]);   
        }
     }

     /*For Down Wave Function*/
     for (i=0; i<ist.n_down; i++) {
        if ( field_up == 1.0 ) {
          wf_down[sites*ist.n_down+i] = RCmul(cns.U_exponent_fermions * cns.factor_spin_up_field_down, wf_down[sites*ist.n_down+i]);  
        }
        else {
          wf_down[sites*ist.n_down+i] = RCmul(cns.U_exponent_fermions * cns.factor_spin_up_field_up, wf_down[sites*ist.n_down+i]);  
        }
     }

return; 
}

/****************************************************************************************/

void propagate_wave_functions_continuous_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,int_st ist,cns_st cns,MKL_Complex16 exponential_up,MKL_Complex16 exponential_down,int sites) {

     /*Propagates the Wave Function*/
     int i;

     /*Multiply Factor Into Wavefunctions*/
     for (i=0; i<ist.n_up; i++) {
          wf_up[sites*ist.n_up+i] = Cmul(exponential_up, wf_up[sites*ist.n_up+i]);
     }

     /*Down*/
     for (i=0; i<ist.n_down; i++) {
          wf_down[sites*ist.n_down+i] = Cmul(exponential_down, wf_down[sites*ist.n_down+i]); 
     }

return;
}

/****************************************************************************************/

void propagate_wave_functions_continuous_fermions_up(MKL_Complex16 *wf_up,int_st ist,cns_st cns,MKL_Complex16 exponential_up,int sites) {

     /*Propagates the Wave Function*/
     int i;

     /*Multiply Factor Into Wavefunctions*/
     for (i=0; i<ist.n_up; i++) {
          wf_up[sites*ist.n_up+i] = Cmul(exponential_up, wf_up[sites*ist.n_up+i]);
     }

return;
}

/****************************************************************************************/

void propagate_wave_functions_continuous_fermions_down(MKL_Complex16 *wf_down,int_st ist,cns_st cns,MKL_Complex16 exponential_down,int sites) {

     /*Propagates the Wave Function*/
     int i;

     /*Down*/
     for (i=0; i<ist.n_down; i++) {
          wf_down[sites*ist.n_down+i] = Cmul(exponential_down, wf_down[sites*ist.n_down+i]);
     }

return;
}

/****************************************************************************************/

void propagate_wave_functions_continuous_fermions_spin_polarized(MKL_Complex16 *wf_up,int_st ist,cns_st cns,MKL_Complex16 exponential_up,int sites) {

     /*Propagates the Wave Function*/
     int i;

     /*Multiply Factor Into Wavefunctions*/
     for (i=0; i<ist.n_up; i++) {
          wf_up[sites*ist.n_up+i] = Cmul(exponential_up, wf_up[sites*ist.n_up+i]);
     }

return;
}

/********************************************************************************************/

void update_overlaps_fermions(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,int_st ist) {

    /*Updates the Values of the Overlaps To Reflect the Propagation*/ 
    MKL_Complex16 *stored_product_up; 
    MKL_Complex16 *stored_product_down; 
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_up=(MKL_Complex16 *)calloc(ist.n_up*ist.n_up,sizeof(MKL_Complex16)); 
    stored_product_down=(MKL_Complex16 *)calloc(ist.n_down*ist.n_down,sizeof(MKL_Complex16)); 

    /*Get Product of Trial and Actual*/
    cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,wf_up,ist.n_up,Zero,stored_product_up,ist.n_up);
    //transpose_cmat_cmat(trial_wf_up,wf_up,stored_product_up,ist.n_sites,ist.n_up,ist.n_up);
    (*overlap_up)=complex_inverse_det_fermions(stored_product_up,overlap_inverse_up,ist.n_up);

    cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_sites,One,trial_wf_down,ist.n_down,wf_down,ist.n_down,Zero,stored_product_down,ist.n_down);
    //transpose_cmat_cmat(trial_wf_down,wf_down,stored_product_down,ist.n_sites,ist.n_down,ist.n_down);
    (*overlap_down)=complex_inverse_det_fermions(stored_product_down,overlap_inverse_down,ist.n_down);

free(stored_product_up); 
free(stored_product_down); 
free(Zero); 
free(One); 
return; 
}

/********************************************************************************************/

void update_overlaps_fermions_up(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,int_st ist) {

    /*Updates the Values of the Overlaps To Reflect the Propagation*/
    MKL_Complex16 *stored_product_up;
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_up=(MKL_Complex16 *)calloc(ist.n_up*ist.n_up,sizeof(MKL_Complex16));

    /*Get Product of Trial and Actual*/
    cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,wf_up,ist.n_up,Zero,stored_product_up,ist.n_up);
    transpose_cmat_cmat(trial_wf_up,wf_up,stored_product_up,ist.n_sites,ist.n_up,ist.n_up);
    (*overlap_up)=complex_inverse_det_fermions(stored_product_up,overlap_inverse_up,ist.n_up);

free(stored_product_up);
free(Zero);
free(One);
return;
}

/********************************************************************************************/

void update_overlaps_fermions_down(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_down,int_st ist) {

    /*Updates the Values of the Overlaps To Reflect the Propagation*/
    MKL_Complex16 *stored_product_down;
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_down=(MKL_Complex16 *)calloc(ist.n_down_sq,sizeof(MKL_Complex16));

    /*Get Product of Trial and Actual*/
    cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_down,ist.n_down,ist.n_sites,One,trial_wf_down,ist.n_down,wf_down,ist.n_down,Zero,stored_product_down,ist.n_down);
    transpose_cmat_cmat(trial_wf_down,wf_down,stored_product_down,ist.n_sites,ist.n_down,ist.n_down);
    (*overlap_down)=complex_inverse_det_fermions(stored_product_down,overlap_inverse_down,ist.n_down);

free(stored_product_down);
free(Zero);
free(One);
return;
} 


/*****************************************************************************************/

void update_overlaps_fermions_spin_polarized(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,int_st ist) {

    /*Updates the Values of the Overlaps To Reflect the Propagation*/
    MKL_Complex16 *stored_product_up;
    MKL_Complex16 *One, *Zero;

    One = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    Zero = (MKL_Complex16 *)malloc(1*sizeof(MKL_Complex16));
    One[0].real = 1.0; One[0].imag = 0.0;
    Zero[0].real = 0.0; Zero[0].imag = 0.0;

    stored_product_up=(MKL_Complex16 *)calloc(ist.n_up*ist.n_up,sizeof(MKL_Complex16));

    /*Get Product of Trial and Actual*/
    cblas_zgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ist.n_up,ist.n_up,ist.n_sites,One,trial_wf_up,ist.n_up,wf_up,ist.n_up,Zero,stored_product_up,ist.n_up);
    //transpose_cmat_cmat(trial_wf_up,wf_up,stored_product_up,ist.n_sites,ist.n_up,ist.n_up);
    (*overlap_up)=complex_inverse_det_fermions(stored_product_up,overlap_inverse_up,ist.n_up);

free(stored_product_up);
free(One); 
free(Zero); 
return;
}

/*****************************************************************************************/

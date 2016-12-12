#include "afqmc.h"

/************************************************************************************/

void propagate_forwards_potential_continuous_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant; 
   MKL_Complex16 ratio_permanent; 
   MKL_Complex16 exponential_field_bosons; 
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
            real_constant = exp(-1.0 * field * cns.propagation_exponent_bosons.real) * cns.U_exponent_bosons; 
            complex_constant = -1.0 * field * cns.propagation_exponent_bosons.imag;    
            exponential_field_bosons.real = real_constant * cos(complex_constant); 
            exponential_field_bosons.imag = real_constant * sin(complex_constant); 

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,exponential_field_bosons,sites);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],exponential_field_bosons,sites); 
            update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);

            /*Multiply Into Weights*/
            weights[walker] = Cmul(weights[walker], ratio_permanent); 

         } /*Check Walker Weights*/

       } /*Run Through Sites*/
   } /*Run Through Walkers*/

return;
}

/***************************************************************************************/

void propagate_forwards_potential_continuous_meanfield_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant;
   MKL_Complex16 ratio_permanent;
   MKL_Complex16 exponential_field_bosons; 
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
            real_constant = exp(-1.0 * field * cns.propagation_exponent_bosons.real) * cns.U_exponent_bosons;
            complex_constant = -1.0 * field * cns.propagation_exponent_bosons.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            /*Determine Mean Field Exponential*/
            real_constant = exp(field * cns.propagation_exponent_bosons.real * cns.total_trial_density); 
            complex_constant = field * cns.propagation_exponent_bosons.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant); 
            mean_field_exponential.imag = real_constant * sin(complex_constant);    

            /*Modify Exponential Constant By Mean Field Term*/
            exponential_field_bosons = RCmul(cns.mean_field_exponential_bosons, exponential_field_bosons);

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,exponential_field_bosons,sites);

            /*Multiply Mean Field Exponential Into Determinant*/
            ratio_permanent = Cmul(ratio_permanent, mean_field_exponential);   
           
            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],exponential_field_bosons,sites); 
            update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);
            
            /*Multiply Into Weights*/
            weights[walker] = Cmul(weights[walker], ratio_permanent); 

          } /*Check Walker Weights*/
        }/*Run Through Sites*/
    }

return; 
}

/********************************************************************************************************************************/

void propagate_forwards_potential_continuous_phaseless_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant;
   double phase, abs_value, projected_ratio_permanent; 
   MKL_Complex16 ratio_permanent; 
   MKL_Complex16 exponential_field_bosons; 
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
            real_constant = exp(-1.0 * field * cns.propagation_exponent_bosons.real) * cns.U_exponent_bosons;
            complex_constant = -1.0 * field * cns.propagation_exponent_bosons.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,exponential_field_bosons,sites);

             /*Propagate the WF and Get New Overlaps*/
             propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],exponential_field_bosons,sites); 
             update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);

            /*Perform the Phaseless*/
            phase = atan(ratio_permanent.imag/ratio_permanent.real);
            abs_value = Cabs(ratio_permanent);
            if ( fabs(phase) < .5 * 3.1415 && ratio_permanent.real > 0 ) {
              projected_ratio_permanent = abs_value * cos(phase);
            }
            else {
              projected_ratio_permanent = 0.0;
            }

            /*Multiply Into Weights*/
            weights[walker] = RCmul(projected_ratio_permanent, weights[walker]); 
     
         } /*Check Walker Weights*/

       } /*Run Through Sites*/
   } /*Run Through Walkers*/

return;
}

/******************************************************************************************************************/

void propagate_forwards_potential_continuous_allshifts_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *weights,double *potential_matrix_bosons,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant;
   MKL_Complex16 shifted_field, shift;
   MKL_Complex16 shifted_field_exponential, shifted_field_exponent; 
   MKL_Complex16 prod_prop_shift;  
   MKL_Complex16 ratio_permanent;
   MKL_Complex16 exponential_field_bosons; 
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
            shift = compute_shift_bosons(&wf_bosons[walker*ist.n_sites],&trial_wf_bosons[0],potential_matrix_bosons,ist,cns,sites);   
            shift.real = shift.imag = 0.0; 
 
            shifted_field.real = field - shift.real; 
            shifted_field.imag = -shift.imag;  

            shifted_field_exponent = Csub(RCmul(field, shift), RCmul(.5, Cmul(shift, shift)));
            shifted_field_exponential.real = exp(shifted_field_exponent.real) * cos(shifted_field_exponent.imag); 
            shifted_field_exponential.imag = exp(shifted_field_exponent.real) * sin(shifted_field_exponent.imag);  

            prod_prop_shift = Cmul(shifted_field, cns.propagation_exponent_bosons); 

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * prod_prop_shift.real) * cns.U_exponent_bosons;
            complex_constant = -1.0 * prod_prop_shift.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            /*Modify Exponential Constant By Mean Field Term*/
            exponential_field_bosons = RCmul(cns.mean_field_exponential_bosons, exponential_field_bosons);

            /*Determine Mean Field Exponential*/
            real_constant = exp(prod_prop_shift.real * cns.total_trial_density);
            complex_constant = prod_prop_shift.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant);
            mean_field_exponential.imag = real_constant * sin(complex_constant);
            
            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,exponential_field_bosons,sites);

            /*Multiply Mean Field Exponential Into Determinant*/
            ratio_permanent = Cmul(ratio_permanent, Cmul(mean_field_exponential, shifted_field_exponential));

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],exponential_field_bosons,sites); 
            update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);

            /*Multiply Into Weights*/
            weights[walker] = Cmul(weights[walker], ratio_permanent);

          } /*Check Walker Weights*/
        }/*Run Through Sites*/
    }

return;
}

/*******************************************************************************************************************************/

void propagate_forwards_potential_continuous_phaseless_meanfield_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant;
   double phase, abs_value, projected_ratio_permanent;
   MKL_Complex16 ratio_permanent; 
   MKL_Complex16 exponential_field_bosons; 
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
            real_constant = exp(-1.0 * field * cns.propagation_exponent_bosons.real) * cns.U_exponent_bosons;
            complex_constant = -1.0 * field * cns.propagation_exponent_bosons.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            /*Determine Mean Field Exponential*/
            real_constant = exp(field * cns.propagation_exponent_bosons.real * cns.total_trial_density);
            complex_constant = field * cns.propagation_exponent_bosons.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant);
            mean_field_exponential.imag = real_constant * sin(complex_constant);

            /*Modify Exponential Constant By Mean Field Term*/
            exponential_field_bosons = RCmul(cns.mean_field_exponential_bosons, exponential_field_bosons);

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,exponential_field_bosons,sites);

            /*Multiply Mean Field Exponential Into Determinant*/
            ratio_permanent = Cmul(ratio_permanent, mean_field_exponential);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],exponential_field_bosons,sites); 
            update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);

           /*Perform the Phaseless*/
           phase = atan(ratio_permanent.imag/ratio_permanent.real);
           abs_value = Cabs(ratio_permanent);
           if ( fabs(phase) < .5 * 3.1415 && ratio_permanent.real > 0 ) {
             projected_ratio_permanent = abs_value * cos(phase);
           }
           else {
             projected_ratio_permanent = 0.0;
           }

           /*Multiply Into Weights*/
           weights[walker] = RCmul(projected_ratio_permanent, weights[walker]);

     } /*Check Walker Weights*/

    } /*Run Through Sites*/
   } /*Run Through Walkers*/

return;
}

/**********************************************************************************************************************************/

void propagate_forwards_potential_continuous_phaseless_allshifts_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,MKL_Complex16 *weights,double *potential_matrix_bosons,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field;
   double real_constant, complex_constant;
   double phase, abs_value, projected_ratio_permanent; 
   MKL_Complex16 shifted_field, shift;
   MKL_Complex16 shifted_field_exponential, shifted_field_exponent;
   MKL_Complex16 prod_prop_shift;
   MKL_Complex16 ratio_permanent;
   MKL_Complex16 exponential_field_bosons; 
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
            shift = compute_shift_bosons(&wf_bosons[walker*ist.n_sites],&trial_wf_bosons[0],potential_matrix_bosons,ist,cns,sites);
            shift.real = shift.imag = 0.0;

            shifted_field.real = field - shift.real;
            shifted_field.imag = -shift.imag;

            shifted_field_exponent = Csub(RCmul(field, shift), RCmul(.5, Cmul(shift, shift)));
            shifted_field_exponential.real = exp(shifted_field_exponent.real) * cos(shifted_field_exponent.imag);
            shifted_field_exponential.imag = exp(shifted_field_exponent.real) * sin(shifted_field_exponent.imag);

            prod_prop_shift = Cmul(shifted_field, cns.propagation_exponent_bosons);

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * prod_prop_shift.real) * cns.U_exponent_bosons;
            complex_constant = -1.0 * prod_prop_shift.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            /*Modify Exponential Constant By Mean Field Term*/
            exponential_field_bosons= RCmul(cns.mean_field_exponential_bosons, exponential_field_bosons);

            /*Determine Mean Field Exponential*/
            real_constant = exp(prod_prop_shift.real * cns.total_trial_density);
            complex_constant = prod_prop_shift.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant);
            mean_field_exponential.imag = real_constant * sin(complex_constant);

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,exponential_field_bosons,sites);

            /*Multiply Mean Field Exponential Into Determinant*/
            ratio_permanent = Cmul(ratio_permanent, Cmul(mean_field_exponential, shifted_field_exponential));

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],exponential_field_bosons,sites); 
            update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);

            /*Perform the Phaseless*/
            phase = atan(ratio_permanent.imag/ratio_permanent.real);
            abs_value = Cabs(ratio_permanent);
            if ( fabs(phase) < .5 * 3.1415 && ratio_permanent.real > 0 ) {
              projected_ratio_permanent = abs_value * cos(phase);
            }
            else {
              projected_ratio_permanent = 0.0;
            }

            /*Multiply Into Weights*/
            weights[walker] = RCmul(projected_ratio_permanent, weights[walker]); 

          } /*Check Walker Weights*/
        }/*Run Through Sites*/
    }

return;
}

/********************************************************************************************************************************/
 
MKL_Complex16 determine_overlap_ratio_continuous_bosons(MKL_Complex16 *wf,MKL_Complex16 *trial_wf,int_st ist,MKL_Complex16 exponential,int site) {

      /*Determines Overlap Ratio*/
      int i;
      MKL_Complex16 numerator, denominator;
      MKL_Complex16 overlap_ratio;

      /*Determine Numerator of Overlap Ratio Using Identical Orbital Representation*/
      numerator.real = numerator.imag = 0.0;
      for (i=0; i<ist.n_sites; i++) {
        if ( i!=site ) {
          numerator = Cadd(numerator, Cmul( conjugate(trial_wf[i]), wf[i]));
        }
        else {
          numerator = Cadd(numerator, Cmul( conjugate(trial_wf[i]), Cmul(exponential, wf[i])));
        }
      }

      /*Denominator Is Just Overlap*/
      denominator = perm_bosons(trial_wf,wf,ist);

      /*Run Trhough All Bosons*/
      for (i=0; i<ist.n_bosons-1; i++) {
        numerator = Cmul(numerator, numerator); 
        denominator = Cmul(denominator, denominator); 
      }

      /*Obtain Ratio*/
      overlap_ratio = Cdiv(numerator, denominator);

return(overlap_ratio);
}

/****************************************************************************************/  

void propagate_wave_functions_continuous_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 exponential_bosons,int sites) {

     /*Multiply Factor Into Wavefunctions*/
     wf_bosons[sites] = Cmul(exponential_bosons, wf_bosons[sites]);

return;
}

/****************************************************************************************/

void update_overlaps_bosons(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,int_st ist) {

    /*Updates the Values of the Overlaps To Reflect the Propagation*/ 
    int i; 
    MKL_Complex16 overlap; 

    /*Get Product of Trial and Actual*/
    overlap=perm_bosons(trial_wf_bosons,wf_bosons,ist); 

    /*Take to Nth Power*/
    for (i=0; i<ist.n_bosons-1; i++) {
     overlap = Cmul(overlap, overlap); 
    }
    (*overlap_bosons) = RCmul(factorial(ist.n_bosons), overlap); 
 
return; 
}

/*****************************************************************************************/

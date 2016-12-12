#include "afqmc.h"

/************************************************************************************/

void propagate_forwards_potential_continuous_bose_fermi_spin_polarized(MKL_Complex16 *wf_up, MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_bosons,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field_one, field_two; 
   double real_constant, complex_constant; 
   MKL_Complex16 ratio_permanent, ratio_determinant;  
   MKL_Complex16 exponential_field_bosons, exponential_field_fermions; 
   MKL_Complex16 total_exponential_field_bosons; 
   long tidum;

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Run Through Sites*/
        for (sites=0; sites<ist.n_sites; sites++) {

          /*Check WEights*/
          if ( Cabs(weights[walker]) != 0 ) {

             /*Select Random Field*/
            field_one = gasdev(&tidum);

            /*First Do Bose-Fermi Uncoupling*/
            real_constant = exp(-1.0 * field_one * cns.propagation_exponent_bose_fermi.real);
            complex_constant = -1.0 * field_one * cns.propagation_exponent_bose_fermi.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            real_constant = exp(field_one * cns.propagation_exponent_bose_fermi.real) * cns.U_exponent_fermions;
            complex_constant = field_one * cns.propagation_exponent_bose_fermi.imag;
            exponential_field_fermions.real = real_constant * cos(complex_constant);
            exponential_field_fermions.imag = real_constant * sin(complex_constant);

            total_exponential_field_bosons = exponential_field_bosons; 

            /*Now Propagate Using Boson Term*/
            field_two = gasdev(&tidum);

            /*Do Boson-Boson Decoupling*/
            real_constant = exp(-1.0 * field_two * cns.propagation_exponent_bosons.real);
            complex_constant = -1.0 * field_two * cns.propagation_exponent_bosons.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            total_exponential_field_bosons = Cmul(total_exponential_field_bosons, exponential_field_bosons);  

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,total_exponential_field_bosons,sites);

            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_fermions,ist.n_up,sites);
            ratio_permanent = Cmul(ratio_permanent, ratio_determinant);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_fermions_spin_polarized(&wf_up[walker*ist.n_sites_n_up],ist,cns,exponential_field_fermions,sites);
            update_overlaps_fermions_spin_polarized(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_up[walker],ist);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],total_exponential_field_bosons,sites);
            update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);

            /*Multiply Into Weights*/
            weights[walker] = Cmul(ratio_permanent, weights[walker]);           

         } /*Check Walker Weights*/

       } /*Run Through Sites*/
   } /*Run Through Walkers*/

return;
}

/***************************************************************************************/

void propagate_forwards_potential_continuous_meanfield_bose_fermi_spin_polarized(MKL_Complex16 *wf_up,MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_bosons,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field_one, field_two;
   double real_constant, complex_constant;
   MKL_Complex16 ratio_determinant, ratio_permanent;
   MKL_Complex16 exponential_field_bosons, exponential_field_fermions; 
   MKL_Complex16 total_exponential_field_bosons; 
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
            field_one = gasdev(&tidum);

            /*First Do Bose-Fermi Uncoupling*/
            real_constant = exp(-1.0 * field_one * cns.propagation_exponent_bose_fermi.real);
            complex_constant = -1.0 * field_one * cns.propagation_exponent_bose_fermi.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            real_constant = exp(field_one * cns.propagation_exponent_bose_fermi.real) * cns.U_exponent_fermions;
            complex_constant = field_one * cns.propagation_exponent_bose_fermi.imag;
            exponential_field_fermions.real = real_constant * cos(complex_constant);
            exponential_field_fermions.imag = real_constant * sin(complex_constant);

            total_exponential_field_bosons = exponential_field_bosons; 

            /*Now Propagate Using Boson Term*/
            field_two = gasdev(&tidum);

            /*Do Boson-Boson Decoupling*/
            real_constant = exp(-1.0 * field_two * cns.propagation_exponent_bosons.real);
            complex_constant = -1.0 * field_two * cns.propagation_exponent_bosons.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            exponential_field_bosons = RCmul(cns.mean_field_exponential_bosons, exponential_field_bosons);
            total_exponential_field_bosons = Cmul(total_exponential_field_bosons, exponential_field_bosons); 

            /*Determine Mean Field Exponential*/
            real_constant = exp(field_two * cns.propagation_exponent_bosons.real * cns.total_trial_density);
            complex_constant = field_two * cns.propagation_exponent_bosons.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant);
            mean_field_exponential.imag = real_constant * sin(complex_constant);

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,total_exponential_field_bosons,sites);

            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_fermions,ist.n_up,sites);
            ratio_permanent = Cmul(ratio_permanent, ratio_determinant);

            /*Multiply Ratio Permanent By Mean Field Term*/
            ratio_permanent = Cmul(ratio_permanent, mean_field_exponential);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_fermions_spin_polarized(&wf_up[walker*ist.n_sites_n_up],ist,cns,exponential_field_fermions,sites);
            update_overlaps_fermions_spin_polarized(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_up[walker],ist);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],total_exponential_field_bosons,sites);
            update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);

            /*Multiply Into Weights*/
            weights[walker] = Cmul(ratio_permanent, weights[walker]);

          } /*Check Walker Weights*/
        }/*Run Through Sites*/
    }

return; 
}

/********************************************************************************************************************************/

void propagate_forwards_potential_continuous_phaseless_bose_fermi_spin_polarized(MKL_Complex16 *wf_up,MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_bosons,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   int n_sites_n_up = ist.n_sites * ist.n_up; 
   int n_up_sq = ist.n_up * ist.n_up; 
   double field_one, field_two;
   double real_constant, complex_constant;
   double phase, abs_value, projected_ratio_permanent; 
   MKL_Complex16 ratio_permanent, ratio_determinant; 
   MKL_Complex16 exponential_field_bosons_one, exponential_field_bosons_two, exponential_field_bosons_total; 
   MKL_Complex16 exponential_field_fermions; 
   long tidum;

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Run Through Sites*/
        for (sites=0; sites<ist.n_sites; sites++) {

          /*Check WEights*/
          if ( Cabs(weights[walker]) != 0 ) {

            /*Select Random Field*/
            field_one = gasdev(&tidum);
            field_two = gasdev(&tidum);

            /*First Do Bose-Fermi Uncoupling*/
            real_constant = exp(-1.0 * field_one * cns.propagation_exponent_bose_fermi.real);
            complex_constant = -1.0 * field_one * cns.propagation_exponent_bose_fermi.imag;
            exponential_field_bosons_one.real = real_constant * cos(complex_constant);
            exponential_field_bosons_one.imag = real_constant * sin(complex_constant);

            real_constant = exp(field_one * cns.propagation_exponent_bose_fermi.real) * cns.U_exponent_fermions; 
            complex_constant = field_one * cns.propagation_exponent_bose_fermi.imag; 
            exponential_field_fermions.real = real_constant * cos(complex_constant); 
            exponential_field_fermions.imag = real_constant * sin(complex_constant);  

            /*Second Field*/
            real_constant = exp(-1.0 * field_two * cns.propagation_exponent_bosons.real);
            complex_constant = -1.0 * field_two * cns.propagation_exponent_bosons.imag;
            exponential_field_bosons_two.real = real_constant * cos(complex_constant);
            exponential_field_bosons_two.imag = real_constant * sin(complex_constant);

            exponential_field_bosons_total = Cmul(exponential_field_bosons_one, exponential_field_bosons_two); 

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,exponential_field_bosons_total,sites);

            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*n_up_sq],ist,cns,exponential_field_fermions,ist.n_up,sites);
            ratio_permanent = Cmul(ratio_permanent, ratio_determinant); 

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_fermions_spin_polarized(&wf_up[walker*n_sites_n_up],ist,cns,exponential_field_fermions,sites);
            update_overlaps_fermions_spin_polarized(&wf_up[walker*n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*n_up_sq],&overlap_up[walker],ist); 

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],exponential_field_bosons_total,sites); 
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

void propagate_forwards_potential_continuous_allshifts_bose_fermi_spin_polarized(MKL_Complex16 *wf_up,MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_bosons,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,double *potential_matrix_bosons,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field_one, field_two;
   double real_constant, complex_constant;
   MKL_Complex16 shifted_field_bose_fermi, shift_bose_fermi;
   MKL_Complex16 shifted_field_exponential_bose_fermi, shifted_field_exponent_bose_fermi;
   MKL_Complex16 shifted_field_bosons, shift_bosons; 
   MKL_Complex16 shifted_field_exponential_bosons, shifted_field_exponent_bosons;  
   MKL_Complex16 prod_prop_shift;  
   MKL_Complex16 ratio_permanent, ratio_determinant;
   MKL_Complex16 exponential_field_bosons, exponential_field_fermions; 
   MKL_Complex16 total_exponential_field_bosons; 
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
            field_one = gasdev(&tidum);

            /*Get shifts*/ 
            //shift_bose_fermi = compute_shift_bose_fermi(&wf_bosons[walker*ist.n_sites],&trial_wf_bosons[walker*ist.n_sites],potential_matrix_bosons,ist,cns,sites);
            shift_bose_fermi.real = shift_bose_fermi.imag = 0.0;

            shifted_field_bose_fermi.real = field_one - shift_bose_fermi.real;
            shifted_field_bose_fermi.imag = -shift_bose_fermi.imag;

            shifted_field_exponent_bose_fermi = Csub(RCmul(field_one, shift_bose_fermi), RCmul(.5, Cmul(shift_bose_fermi, shift_bose_fermi)));
            shifted_field_exponential_bose_fermi.real = exp(shifted_field_exponent_bose_fermi.real) * cos(shifted_field_exponent_bose_fermi.imag);
            shifted_field_exponential_bose_fermi.imag = exp(shifted_field_exponent_bose_fermi.real) * sin(shifted_field_exponent_bose_fermi.imag);

            prod_prop_shift = Cmul(shifted_field_bose_fermi, cns.propagation_exponent_bose_fermi);

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * prod_prop_shift.real); 
            complex_constant = -1.0 * prod_prop_shift.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            total_exponential_field_bosons = exponential_field_bosons; 

            /*Get Fermion Shift*/
            real_constant = exp(prod_prop_shift.real) * cns.U_exponent_fermions;
            complex_constant = prod_prop_shift.imag;
            exponential_field_fermions.real = real_constant * cos(complex_constant);
            exponential_field_fermions.imag = real_constant * sin(complex_constant);           

            /*Now Propagate Using Boson Term*/
            field_two = gasdev(&tidum);

            /*Get Shifted Fields*/
            /*shift_bosons = compute_shift_bosons(&wf_bosons[walker*ist.n_sites],&trial_wf_bosons[walker*ist.n_sites],potential_matrix_bosons,ist,cns,sites);*/
            shift_bosons.real = shift_bosons.imag = 0.0;

            shifted_field_bosons.real = field_two - shift_bosons.real;
            shifted_field_bosons.imag = -shift_bosons.imag;

            shifted_field_exponent_bosons = Csub(RCmul(field_two, shift_bosons), RCmul(.5, Cmul(shift_bosons, shift_bosons)));
            shifted_field_exponential_bosons.real = exp(shifted_field_exponent_bosons.real) * cos(shifted_field_exponent_bosons.imag);
            shifted_field_exponential_bosons.imag = exp(shifted_field_exponent_bosons.real) * sin(shifted_field_exponent_bosons.imag);

            prod_prop_shift = Cmul(shifted_field_bosons, cns.propagation_exponent_bosons);

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * prod_prop_shift.real);
            complex_constant = -1.0 * prod_prop_shift.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            /*Modify Exponential Constant By Mean Field Term*/
            exponential_field_bosons = RCmul(cns.mean_field_exponential_bosons, exponential_field_bosons);

            total_exponential_field_bosons = Cmul(total_exponential_field_bosons, exponential_field_bosons); 

            /*Mean Field Expnent*/
            real_constant = exp(prod_prop_shift.real * cns.total_trial_density);
            complex_constant = prod_prop_shift.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant);
            mean_field_exponential.imag = real_constant * sin(complex_constant);

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,total_exponential_field_bosons,sites);

            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_fermions,ist.n_up,sites);
            ratio_permanent = Cmul(ratio_permanent, ratio_determinant);

            /*Multiply Ratio Permanent By Mean Field Term*/
            ratio_permanent = Cmul(ratio_permanent, mean_field_exponential);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_fermions_spin_polarized(&wf_up[walker*ist.n_sites_n_up],ist,cns,exponential_field_fermions,sites);
            update_overlaps_fermions_spin_polarized(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_up[walker],ist);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],total_exponential_field_bosons,sites);
            update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);

            /*Multiply Into Weights*/
            weights[walker] = Cmul(ratio_permanent, weights[walker]);   

          } /*Check Walker Weights*/
        }/*Run Through Sites*/
    }

return;
}

/*******************************************************************************************************************************/

void propagate_forwards_potential_continuous_phaseless_meanfield_bose_fermi_spin_polarized(MKL_Complex16 *wf_up,MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_bosons,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){ 

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field_one, field_two; 
   double real_constant, complex_constant;
   double phase, abs_value, projected_ratio_permanent;
   MKL_Complex16 ratio_permanent, ratio_determinant;  
   MKL_Complex16 exponential_field_bosons, exponential_field_fermions;
   MKL_Complex16 total_exponential_field_bosons;  
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
            field_one = gasdev(&tidum);

            /*First Do Bose-Fermi Uncoupling*/
            real_constant = exp(-1.0 * field_one * cns.propagation_exponent_bose_fermi.real);
            complex_constant = -1.0 * field_one * cns.propagation_exponent_bose_fermi.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            real_constant = exp(field_one * cns.propagation_exponent_bose_fermi.real) * cns.U_exponent_fermions;
            complex_constant = field_one * cns.propagation_exponent_bose_fermi.imag;
            exponential_field_fermions.real = real_constant * cos(complex_constant);
            exponential_field_fermions.imag = real_constant * sin(complex_constant);

            total_exponential_field_bosons = exponential_field_bosons; 

            /*Now Propagate Using Boson Term*/
            field_two = gasdev(&tidum);

            /*Do Boson-Boson Decoupling*/
            real_constant = exp(-1.0 * field_two * cns.propagation_exponent_bosons.real);
            complex_constant = -1.0 * field_two * cns.propagation_exponent_bosons.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            exponential_field_bosons = RCmul(cns.mean_field_exponential_bosons, exponential_field_bosons);
            total_exponential_field_bosons = Cmul(total_exponential_field_bosons, exponential_field_bosons); 

            /*Determine Mean Field Exponential*/
            real_constant = exp(field_two * cns.propagation_exponent_bosons.real * cns.total_trial_density);
            complex_constant = field_two * cns.propagation_exponent_bosons.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant);
            mean_field_exponential.imag = real_constant * sin(complex_constant);

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,total_exponential_field_bosons,sites);

            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_fermions,ist.n_up,sites);
            ratio_permanent = Cmul(ratio_permanent, ratio_determinant);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_fermions_spin_polarized(&wf_up[walker*ist.n_sites_n_up],ist,cns,exponential_field_fermions,sites);
            update_overlaps_fermions_spin_polarized(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_up[walker],ist);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],total_exponential_field_bosons,sites);
            update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);

            /*Multiply Ratio Permanent By Mean Field Term*/
            ratio_permanent = Cmul(ratio_permanent, mean_field_exponential);

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

void propagate_forwards_potential_continuous_phaseless_allshifts_bose_fermi_spin_polarized(MKL_Complex16 *wf_up,MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_bosons,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *weights,double *potential_matrix_bosons,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites;
   double field_one, field_two;
   double real_constant, complex_constant;
   double phase, abs_value, projected_ratio_permanent; 
   MKL_Complex16 shifted_field_bosons, shift_bosons;
   MKL_Complex16 shifted_field_exponential_bosons, shifted_field_exponent_bosons;
   MKL_Complex16 shifted_field_bose_fermi, shift_bose_fermi; 
   MKL_Complex16 shifted_field_exponential_bose_fermi, shifted_field_exponent_bose_fermi; 
   MKL_Complex16 prod_prop_shift;
   MKL_Complex16 ratio_permanent, ratio_determinant;
   MKL_Complex16 exponential_field_bosons, exponential_field_fermions; 
   MKL_Complex16 total_exponential_field_bosons; 
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
            field_one = gasdev(&tidum);

            //shift_bose_fermi = compute_shift_bose_fermi(&wf_bosons[walker*ist.n_sites],&trial_wf_bosons[0],potential_matrix_bosons,ist,cns,sites);
            shift_bose_fermi.real = shift_bose_fermi.imag = 0.0;

            shifted_field_bose_fermi.real = field_one - shift_bose_fermi.real;
            shifted_field_bose_fermi.imag = -shift_bose_fermi.imag;

            shifted_field_exponent_bose_fermi = Csub(RCmul(field_one, shift_bose_fermi), RCmul(.5, Cmul(shift_bose_fermi, shift_bose_fermi)));
            shifted_field_exponential_bose_fermi.real = exp(shifted_field_exponent_bose_fermi.real) * cos(shifted_field_exponent_bose_fermi.imag);
            shifted_field_exponential_bose_fermi.imag = exp(shifted_field_exponent_bose_fermi.real) * sin(shifted_field_exponent_bose_fermi.imag);

            shifted_field_exponential_bose_fermi.real = shifted_field_exponential_bose_fermi.imag = 0.0; 

            prod_prop_shift = Cmul(shifted_field_bose_fermi, cns.propagation_exponent_bose_fermi);

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * prod_prop_shift.real);
            complex_constant = -1.0 * prod_prop_shift.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            total_exponential_field_bosons = exponential_field_bosons; 
 
            /*Get Fermion Shift*/
            real_constant = exp(prod_prop_shift.real) * cns.U_exponent_fermions;
            complex_constant = prod_prop_shift.imag;
            exponential_field_fermions.real = real_constant * cos(complex_constant);
            exponential_field_fermions.imag = real_constant * sin(complex_constant);

            /*Now Propagate Using Boson Term*/
            field_two = gasdev(&tidum);

            /*Get Shifted Fields*/
            shift_bosons = compute_shift_bosons(&wf_bosons[walker*ist.n_sites],&trial_wf_bosons[0],potential_matrix_bosons,ist,cns,sites);
            //shift_bosons.r = shift_bosons.i = 0.0;

            shifted_field_bosons.real = field_two - shift_bosons.real;
            shifted_field_bosons.imag = -shift_bosons.imag;

            shifted_field_exponent_bosons = Csub(RCmul(field_two, shift_bosons), RCmul(.5, Cmul(shift_bosons, shift_bosons)));
            shifted_field_exponential_bosons.real = exp(shifted_field_exponent_bosons.real) * cos(shifted_field_exponent_bosons.imag);
            shifted_field_exponential_bosons.imag = exp(shifted_field_exponent_bosons.real) * sin(shifted_field_exponent_bosons.imag);

            //shifted_field_exponential_bosons.r = shifted_field_exponential_bosons.i = 0.0; 

            prod_prop_shift = Cmul(shifted_field_bosons, cns.propagation_exponent_bosons);

            /*Get Multiplicative Constant*/
            real_constant = exp(-1.0 * prod_prop_shift.real);
            complex_constant = -1.0 * prod_prop_shift.imag;
            exponential_field_bosons.real = real_constant * cos(complex_constant);
            exponential_field_bosons.imag = real_constant * sin(complex_constant);

            /*Modify Exponential Constant By Mean Field Term*/
            exponential_field_bosons = RCmul(cns.mean_field_exponential_bosons, exponential_field_bosons);

            /*Mean Field Expnent*/
            real_constant = exp(prod_prop_shift.real * cns.total_trial_density);
            complex_constant = prod_prop_shift.imag * cns.total_trial_density;
            mean_field_exponential.real = real_constant * cos(complex_constant);
            mean_field_exponential.imag = real_constant * sin(complex_constant);

            total_exponential_field_bosons = Cmul(total_exponential_field_bosons, Cmul(exponential_field_bosons, mean_field_exponential));  

            /*Determine Ratio of Determinants*/
            ratio_permanent = determine_overlap_ratio_continuous_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,ist,total_exponential_field_bosons,sites);

            ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_fermions,ist.n_up,sites);
            ratio_permanent = Cmul(ratio_permanent, ratio_determinant);

            /*Multiply Ratio Permanent By Mean Field Term*/
            ratio_permanent = Cmul(ratio_permanent, mean_field_exponential);

            /*Propagate the WF and Get New Overlaps*/
            propagate_wave_functions_continuous_fermions_spin_polarized(&wf_up[walker*ist.n_sites_n_up],ist,cns,exponential_field_fermions,sites);
            update_overlaps_fermions_spin_polarized(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_up[walker],ist);

            propagate_wave_functions_continuous_bosons(&wf_bosons[walker*ist.n_sites],total_exponential_field_bosons,sites);
            update_overlaps_bosons(&wf_bosons[walker*ist.n_sites],trial_wf_bosons,&overlap_bosons[walker],ist);

            /*Multiply Into Weights*/
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
 
MKL_Complex16 determine_overlap_ratio_continuous_bose_fermi_spin_polarized(MKL_Complex16 *wf,MKL_Complex16 *trial_wf,int_st ist,MKL_Complex16 exponential,int site) {

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

void propagate_wave_functions_continuous_bose_fermi_spin_polarized(MKL_Complex16 *wf_bosons,MKL_Complex16 exponential_bosons,int sites) {

     /*Multiply Factor Into Wavefunctions*/
     wf_bosons[sites] = Cmul(exponential_bosons, wf_bosons[sites]);

return;
}

/****************************************************************************************/

void update_overlaps_bose_fermi_spin_polarized(MKL_Complex16 *wf_bosons,MKL_Complex16 *trial_wf_bosons,MKL_Complex16 *overlap_bosons,int_st ist) {

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

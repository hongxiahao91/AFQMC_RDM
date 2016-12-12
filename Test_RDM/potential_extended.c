#include "afqmc.h"

/***************************************************************************************/

void propagate_forwards_potential_continuous_extended(MKL_Complex16 *wf_up,MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_up,MKL_Complex16 *trial_wf_down,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_up,MKL_Complex16 *overlap_down,int *neighbors,int *number_neighbors,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites, i;
   double field_one, field_two, field_three, field_four, field_five; 
   double real_constant, complex_constant; 
   MKL_Complex16 exponential_field_up_site, exponential_field_up_nextsite, exponential_field_down_site, exponential_field_down_nextsite; 
   MKL_Complex16 sum_exponent; 
   MKL_Complex16 ratio_determinant;
   long tidum;
   FILE *pf = fopen("errors.dat", "a+"); 

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Check Weights*/
        if ( Cabs(weights[walker]) != 0 ) { 

          /*Run Through Sites*/
          for (sites=0; sites<ist.n_sites; sites++) {

            /*Run Through Neighbors for Each Site*/
            for (i=0; i<(int)(number_neighbors[sites]/2.0); i++) {

              /*Select Random Field*/
              field_one = gasdev(&tidum);
              field_two = gasdev(&tidum); 
              field_three = gasdev(&tidum); 
              field_four = gasdev(&tidum);   
              field_five = gasdev(&tidum); 

              /*************************************************************************************/
 
              /*Add Terms for Site and Up Spin Different Sites*/
              if ( i==0 ) {
                sum_exponent = Cadd(RCmul(field_one+field_two, cns.propagation_exponent_extended), RCmul(field_five, cns.propagation_exponent_fermions)); 
              }
              else {
                sum_exponent = RCmul(field_one+field_two, cns.propagation_exponent_extended); 
              }
              real_constant = exp(sum_exponent.real) * cns.U_exponent_fermions; 
              complex_constant = sum_exponent.imag; 
              exponential_field_up_site.real = real_constant * cos(complex_constant); 
              exponential_field_up_site.imag = real_constant * sin(complex_constant);
 
              /*Add Terms for Site and Down Spin Different Sites*/
              if ( i==0 ) {
                sum_exponent = Csub(RCmul(field_three+field_four, cns.propagation_exponent_extended), RCmul(field_five, cns.propagation_exponent_fermions)); 
              }
              else {
                sum_exponent = RCmul(field_three+field_four, cns.propagation_exponent_extended); 
              } 
              real_constant = exp(sum_exponent.real) * cns.U_exponent_fermions; 
              complex_constant = sum_exponent.imag; 
              exponential_field_down_site.real = real_constant * cos(complex_constant); 
              exponential_field_down_site.imag = real_constant * sin(complex_constant); 

              /*Now Update Based Upon the Particular Site Updates*/
              ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_up_site,ist.n_up,sites);
              ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,exponential_field_down_site,ist.n_down,sites));

              /*Propagate the WF and Get New Overlaps*/
              propagate_wave_functions_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],ist,cns,exponential_field_up_site,exponential_field_down_site,sites);
              update_overlaps_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],&overlap_up[walker],&overlap_down[walker],ist); 

              /*Multiply Into Weights*/
              weights[walker] = Cmul(weights[walker], ratio_determinant);

              /**************************************************************************************/

              /*Add Terms for the Neighboring Site and Spin Up*/
              sum_exponent = RCmul(-1*(field_one+field_three), cns.propagation_exponent_extended); 
              real_constant = exp(sum_exponent.real); 
              complex_constant = sum_exponent.imag; 
              exponential_field_up_nextsite.real = real_constant * cos(complex_constant); 
              exponential_field_up_nextsite.imag = real_constant * sin(complex_constant); 

              /*Add Terms for Neighboring Site and Spin Down*/
              sum_exponent = RCmul(-1*(field_two+field_four), cns.propagation_exponent_extended); 
              real_constant = exp(sum_exponent.real);
              complex_constant = sum_exponent.imag;
              exponential_field_down_nextsite.real = real_constant * cos(complex_constant);
              exponential_field_down_nextsite.imag = real_constant * sin(complex_constant);       

              /*Determine Ratio of Determinants*/
              ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_up_nextsite,ist.n_up,neighbors[4*sites+i]);
              ratio_determinant = Cmul(ratio_determinant, determine_overlap_ratio_continuous_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,exponential_field_down_nextsite,ist.n_down,neighbors[4*sites+i]));

               /*Propagate the WF and Get New Overlaps*/
               propagate_wave_functions_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],ist,cns,exponential_field_up_nextsite,exponential_field_down_nextsite,neighbors[4*sites+i]);
               update_overlaps_fermions(&wf_up[walker*ist.n_sites_n_up],&wf_down[walker*ist.n_sites_n_down],trial_wf_up,trial_wf_down,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_inverse_down[walker*ist.n_down_sq],&overlap_up[walker],&overlap_down[walker],ist);

              /*Multiply Into Weights*/
              weights[walker] = Cmul(weights[walker], ratio_determinant); 

            } /*Run Through Neighboring Sites*/

          } /*Run Through Sites*/

       } /*Check Walker Weights*/

   } /*Run Through Walkers*/


return;
}

/***************************************************************************************/

void propagate_forwards_potential_continuous_extended_up(MKL_Complex16 *wf_up,MKL_Complex16 *trial_wf_up,MKL_Complex16 *overlap_inverse_up,MKL_Complex16 *overlap_up,int *neighbors,int *number_neighbors,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites, i;
   double field_one; 
   double real_constant, complex_constant;
   MKL_Complex16 exponential_field_up_site, exponential_field_up_nextsite; 
   MKL_Complex16 sum_exponent;
   MKL_Complex16 ratio_determinant;
   long tidum;

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Check Weights*/
        if ( Cabs(weights[walker]) != 0 ) {

          /*Run Through Sites*/
          for (sites=0; sites<ist.n_sites; sites++) {

            /*Run Through Neighbors for Each Site*/
            for (i=0; i<(int)(number_neighbors[sites]/2.0); i++) {

              /*Select Random Field*/
              field_one = gasdev(&tidum);

              /*************************************************************************************/

              /*Add Terms for Site and Up Spin Different Sites*/
              sum_exponent = RCmul(field_one, cns.propagation_exponent_extended);
              real_constant = exp(sum_exponent.real) * cns.U_exponent_fermions;
              complex_constant = sum_exponent.imag;
              exponential_field_up_site.real = real_constant * cos(complex_constant);
              exponential_field_up_site.imag = real_constant * sin(complex_constant);

              /*Now Update Based Upon the Particular Site Updates*/
              ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_up_site,ist.n_up,sites);

              /*Propagate the WF and Get New Overlaps*/
              propagate_wave_functions_continuous_fermions_up(&wf_up[walker*ist.n_sites_n_up],ist,cns,exponential_field_up_site,sites);
              update_overlaps_fermions_up(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_up[walker],ist);

              /*Multiply Into Weights*/
              weights[walker] = Cmul(weights[walker], ratio_determinant);

              /**************************************************************************************/

              /*Add Terms for the Neighboring Site and Spin Up*/
              sum_exponent = RCmul(-1*field_one, cns.propagation_exponent_extended);  
              real_constant = exp(sum_exponent.real); 
              complex_constant = sum_exponent.imag; 
              exponential_field_up_nextsite.real = real_constant * cos(complex_constant); 
              exponential_field_up_nextsite.imag = real_constant * sin(complex_constant); 

              /*Determine Ratio of Determinants*/
              ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],ist,cns,exponential_field_up_nextsite,ist.n_up,neighbors[4*sites+i]);

              /*Propagate the WF and Get New Overlaps*/
              propagate_wave_functions_continuous_fermions_up(&wf_up[walker*ist.n_sites_n_up],ist,cns,exponential_field_up_nextsite,neighbors[4*sites+i]);
              update_overlaps_fermions_up(&wf_up[walker*ist.n_sites_n_up],trial_wf_up,&overlap_inverse_up[walker*ist.n_up_sq],&overlap_up[walker],ist);

              /*Multiply Into Weights*/
              weights[walker] = Cmul(weights[walker], ratio_determinant);

            } /*Run Through Neighboring Sites*/

          } /*Run Through Sites*/

       } /*Check Walker Weights*/

   } /*Run Through Walkers*/


return;
}

/***************************************************************************************/

void propagate_forwards_potential_continuous_extended_down(MKL_Complex16 *wf_down,MKL_Complex16 *trial_wf_down, MKL_Complex16 *overlap_inverse_down,MKL_Complex16 *overlap_down,int *neighbors,int *number_neighbors,MKL_Complex16 *weights,int_st ist,cns_st cns,long *idum){

   /*Propagates the Potential Forward for Each Walker*/
   /*Currently Updating Matrices the Stupid Way As a Check on the Formalism--Must Be Replaced By the Smart Way!!!*/
   int walker;
   int sites, i;
   double field_one; 
   double real_constant, complex_constant;
   MKL_Complex16 exponential_field_down_site, exponential_field_down_nextsite;
   MKL_Complex16 sum_exponent;
   MKL_Complex16 ratio_determinant;
   long tidum;

   Randomize(); tidum = -random();

   /*Run Through Walkers*/
   for (walker = 0; walker < ist.n_walkers; walker++) {

        /*Check Weights*/
        if ( Cabs(weights[walker]) != 0 ) {

          /*Run Through Sites*/
          for (sites=0; sites<ist.n_sites; sites++) {

            /*Run Through Neighbors for Each Site*/
            for (i=0; i<(int)(number_neighbors[sites]/2.0); i++) {

              /*Select Random Field*/
              field_one = gasdev(&tidum);

              /*************************************************************************************/

              /*Add Terms for Site and Up Spin Different Sites*/
              sum_exponent = RCmul(field_one, cns.propagation_exponent_extended);
              real_constant = exp(sum_exponent.real) * cns.U_exponent_fermions;
              complex_constant = sum_exponent.imag;
              exponential_field_down_site.real = real_constant * cos(complex_constant);
              exponential_field_down_site.imag = real_constant * sin(complex_constant);

              /*Now Update Based Upon the Particular Site Updates*/
              ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,exponential_field_down_site,ist.n_down,sites);

              /*Propagate the WF and Get New Overlaps*/
              propagate_wave_functions_continuous_fermions_down(&wf_down[walker*ist.n_sites_n_down],ist,cns,exponential_field_down_site,sites);
              update_overlaps_fermions_down(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],&overlap_down[walker],ist);

              /*Multiply Into Weights*/
              weights[walker] = Cmul(weights[walker], ratio_determinant);

              /**************************************************************************************/

              /*Add Terms for the Neighboring Site and Spin Up*/
              sum_exponent = RCmul(-1*field_one, cns.propagation_exponent_extended);
              real_constant = exp(sum_exponent.real);
              complex_constant = sum_exponent.imag;
              exponential_field_down_nextsite.real = real_constant * cos(complex_constant);
              exponential_field_down_nextsite.imag = real_constant * sin(complex_constant);

              /*Determine Ratio of Determinants*/
              ratio_determinant = determine_overlap_ratio_continuous_fermions(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],ist,cns,exponential_field_down_nextsite,ist.n_down,neighbors[4*sites+i]);

              /*Propagate the WF and Get New Overlaps*/
              propagate_wave_functions_continuous_fermions_down(&wf_down[walker*ist.n_sites_n_down],ist,cns,exponential_field_down_nextsite,neighbors[4*sites+i]);
              update_overlaps_fermions_down(&wf_down[walker*ist.n_sites_n_down],trial_wf_down,&overlap_inverse_down[walker*ist.n_down_sq],&overlap_down[walker],ist);
              
              /*Multiply Into Weights*/
              weights[walker] = Cmul(weights[walker], ratio_determinant);

            } /*Run Through Neighboring Sites*/

          } /*Run Through Sites*/

       } /*Check Walker Weights*/

   } /*Run Through Walkers*/


return;
}

/***************************************************************************************/

#include "afqmc.h"

/**************************************************************************************/

double error_real(double average, double stdev, double Number) {
  /*Finds the Error Based on the Average Over the Number of Runs and the Standard Deviation*/

  double average_runs, stdev_runs;
  double error_runs=0;
  double average_squared;
  double Holder;

  average_runs=1.0/Number*average;
  stdev_runs=1.0/Number*stdev;
  average_squared=average_runs*average_runs;

  Holder=stdev_runs-average_squared;

  if ( Holder > 0 ) {
    error_runs=sqrt(1.0/((double) Number-1)*Holder);
  }
  else {
   error_runs = 0.0; 
  }

return(error_runs);
}

/************************************************************************************************/

void accumulate_energy_blocks(MKL_Complex16 potential_energy, MKL_Complex16 kinetic_energy, MKL_Complex16 total_energy, MKL_Complex16 current_weight, MKL_Complex16 *total_potential_energy, MKL_Complex16 *total_kinetic_energy, MKL_Complex16 *total_energy_2, MKL_Complex16 *total_weight) {

  /*Accumulate the Average and Squared Average Energies To Calculate Statistics Later On*/
  /*Note That I Do Not Take Complex Parts Because They Confound Statistics*/
  total_potential_energy[0] = Cadd(total_potential_energy[0], potential_energy); 
  potential_energy = Cdiv(potential_energy, current_weight); 
  potential_energy.imag = 0; 
  total_potential_energy[1] = Cadd(total_potential_energy[1], potential_energy); 
  total_potential_energy[2] = Cadd(total_potential_energy[2], Cmul(potential_energy, potential_energy)); 

  total_kinetic_energy[0] = Cadd(total_kinetic_energy[0], kinetic_energy);
  kinetic_energy = Cdiv(kinetic_energy, current_weight);
  kinetic_energy.imag = 0; 
  total_kinetic_energy[1] = Cadd(total_kinetic_energy[1], kinetic_energy);
  total_kinetic_energy[2] = Cadd(total_kinetic_energy[2], Cmul(kinetic_energy, kinetic_energy));

  total_energy_2[0] = Cadd(total_energy_2[0], total_energy);
  total_energy = Cdiv(total_energy, current_weight);
  total_energy.imag = 0.0; 
  total_energy_2[1] = Cadd(total_energy_2[1], total_energy);
  total_energy_2[2] = Cadd(total_energy_2[2], Cmul(total_energy, total_energy));

  total_weight[0] = Cadd(total_weight[0], current_weight); 

return; 
}

/************************************************************************************/

int determine_production_runs(double *current_energy,double *previous_energies,int_st ist,cns_st cns) {

    /*Determine Whether The Energies Should Be Measured As In Production Runs*/
    int i; 
    int production_or_not = 0;
    int change_negative = 0, change_positive = 0;   
    double percent_change; 
    FILE *pf = fopen("checkdetermineproduction.dat", "a+"); 

    /*Move All Energies Stored Back By One*/
    for (i=8; i>-1; i--) {
       previous_energies[i+1]=previous_energies[i]; 
    }
    previous_energies[0] = *current_energy; 

    /*Calculate Percent Change in Energy*/
    for (i=0; i<9; i++) {
     if ( previous_energies[i+1] != 0 ) {
       percent_change = (previous_energies[i]-previous_energies[i+1])/previous_energies[i+1]; 
      fprintf(pf, "percent change %f i %d\n", percent_change, i); fflush(pf);  
     }
     else {
       percent_change = 0; 
       fprintf(pf, "percent change %f i %d\n", percent_change, i); fflush(pf); 
      }    
 
       if ( percent_change > 0 ) {
         change_positive++; 
       }
       else { 
         change_negative++; 
       }

    }

   /*Assume Equilibrated If More Than Three Changes are Positive*/
   if ( change_positive >= 3 ) {
     production_or_not = 1; 
   }

return(production_or_not); 
}

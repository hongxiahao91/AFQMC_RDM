/*This is a Code Designed to Determine the Best Block Length for Reblocking - 11/26/14*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int N_data_points = 2900; 
int N_start = 1200;  
int N_total_data_points = N_data_points-N_start; 

FILE *input = fopen("energy_only.dat", "r"); 
FILE *output = fopen("efficiency.dat", "a+"); 

int main() {

   int i, j, k;
   int max_power, max_number_blocks, number_blocks;  
   int block_size; 
   double time, current_energy; 
   double *energies, *energies_blocks;     
   double mean_run, mean_block; 
   double stdev_run, stdev_block; 

   max_number_blocks = (int)(N_total_data_points/2.0);

   energies = (double *)calloc(N_total_data_points,sizeof(double));  
   energies_blocks = (double *)calloc(max_number_blocks,sizeof(double)); 

   /*First Determine the Mean and Standard Deviation of the Entire Data Set*/   
   mean_run = 0; 
   for (i = 0; i<N_data_points; i++) {
     fscanf(input, "%lf %lf", &time, &current_energy); 
     if ( i >= N_start ) {
      energies[i-N_start] = current_energy; 
      mean_run += energies[i-N_start]; 
     }
   }
   mean_run /= ((double) N_total_data_points); 

   /*Now Determine the Standard Deviation*/
   stdev_run = 0; 
   for (i=0; i<N_total_data_points; i++) {
      stdev_run += (energies[i] - mean_run) * (energies[i] - mean_run); 
   }
   stdev_run /= ((double) N_total_data_points); 
   stdev_run = sqrt(stdev_run);

   /*Now Take STDEV of Mean*/
   stdev_run /= sqrt((double) N_total_data_points);  

   /*Find Max power*/
   max_power = (int) (log(N_total_data_points) / log(2));  
 
   /*Now Run Through Each Possible Power and Calculate St Dev*/  
   for (i=1; i<max_power; i++) {
     block_size = pow(2, i); 
     number_blocks = (int) (N_total_data_points/((double) block_size)); 

     /*Determine Block Averages*/
     for (j=0; j<number_blocks; j++) {
  
       mean_block = 0; 
       for (k=0; k<block_size; k++) {
         mean_block += energies[j*block_size+k]; 
       }
       mean_block /= ((double) block_size);          
       energies_blocks[j] = mean_block; 

     }

     /*Determine Block Standard Deviations*/
     stdev_block = 0.0; 
     for (j=0; j<number_blocks; j++) {
        stdev_block += (energies_blocks[j]-mean_run)*(energies_blocks[j]-mean_run);
     }
     stdev_block /= ((double) number_blocks); 

     fprintf(output, "%d\t %f\t %f\t %f\n", block_size, stdev_block, stdev_run, block_size*stdev_block/stdev_run); fflush(output); 

   }  

   fprintf(output, "%f\t %f\n", mean_run, stdev_run); fflush(output); 

return(0); 
}

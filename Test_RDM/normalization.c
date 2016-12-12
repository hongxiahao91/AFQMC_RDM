#include "afqmc.h"

/*********************************************************/

void normalize(MKL_Complex16 *weights,int_st ist,cns_st cns) {
 
   /*Normalize All of the Walker Weights by the Trial Energy*/
   int i; 

   for (i=0; i<ist.n_walkers; i++) {
    weights[i] = RCmul(cns.exp_trial_energy, weights[i]); 
   }

return; 
}

/***********************************************************/

#include "afqmc.h"

/*Determines the Neighbors of a Given Lattice Site*/

/***************************************************************************************/

void neighborsopenboundary(int site, int *neighbors, int *number_neighbors,int_st ist) {
  /*Obtains the 2*dimensions nearest-neighbors directly surrounding a site*/
  /*Should be Done More Generally, But This Works for the Moment!!!*/

  int nsites_squared_one_two=ist.n_sites_one*ist.n_sites_two; 
  int nsites_squared_one_one=ist.n_sites_one*ist.n_sites_one;  

   /*Ordering Such That Double Counting Can Be Avoided In Certain Cases*/
   if (ist.n_sites_two == 1) {
     if (site!=ist.n_sites_one-1 && site!=0) {
       neighbors[site*4]=site-1;
       neighbors[site*4+1]=site+1;
       number_neighbors[site]=2;
     }
     else if (site==ist.n_sites_one-1){
       neighbors[site*4]=site-1;
       number_neighbors[site]=1;
     }
     else if (site==0) {
       neighbors[site*4]=site+1;
       number_neighbors[site]=1;
     }
   }
   else { /*If Two Dimensions*/
     if (ist.n_sites_one==ist.n_sites_two) {
       if (site<ist.n_sites_one) {
        if (site!=0 && site!=ist.n_sites_one-1){
          neighbors[4*site]=site-1;
          neighbors[4*site+2]=site+1;
          neighbors[4*site+1]=site+ist.n_sites_one; 
          number_neighbors[site]=3;
        }
        else if (site==0) {
           neighbors[4*site]=site+1;
           neighbors[4*site+1]=site+ist.n_sites_one;
           number_neighbors[site]=2;
        }
        else if (site==ist.n_sites_one-1) {
           neighbors[4*site]=site-1;
           neighbors[4*site+1]=site+ist.n_sites_one;
           number_neighbors[site]=2;
        }
       }
       else if (site<nsites_squared_one_one && site>=nsites_squared_one_one-ist.n_sites_one) {
         if (site!=nsites_squared_one_one-ist.n_sites_one && site!=nsites_squared_one_one-1) {
          neighbors[4*site]=site-1;
          neighbors[4*site+2]=site+1;
          neighbors[4*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=3;
         }
         else if (site==nsites_squared_one_one-ist.n_sites_one) {
          neighbors[4*site]=site+1;
          neighbors[4*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=2;
        }
        else if (site==nsites_squared_one_one-1) {
          neighbors[4*site]=site-1;
          neighbors[4*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=2;
        }
      }
      else if (site%ist.n_sites_one==0) {
         if (site!=0 && site!=nsites_squared_one_one-ist.n_sites_one) {
           neighbors[4*site+2]=site+1;
           neighbors[4*site+1]=site+ist.n_sites_one;
           neighbors[4*site]=site-ist.n_sites_one;
           number_neighbors[site]=3;
         }
      }
      else if (site%ist.n_sites_one==ist.n_sites_one-1) {
         if (site!=ist.n_sites_one-1 && site!=nsites_squared_one_one-1) {
            neighbors[4*site]=site-1;
            neighbors[4*site+1]=site+ist.n_sites_one;
            neighbors[4*site+2]=site-ist.n_sites_one;
            number_neighbors[site]=3;
       }
      }
      else {
           neighbors[4*site]=site-1;
           neighbors[4*site+2]=site+1;
           neighbors[4*site+1]=site+ist.n_sites_one;
           neighbors[4*site+3]=site-ist.n_sites_one;
           number_neighbors[site]=4;
      }
     } /*If LengthSideOne==LengthSideTwo*/
     else {  /*I Assume Length Side One Is Horizontal and LengthSide Two Is Vertical*/
       if (site<ist.n_sites_one) {
        if (site!=0 && site!=ist.n_sites_one-1){
          neighbors[4*site]=site-1;
          neighbors[4*site+2]=site+1;
          neighbors[4*site+1]=site+ist.n_sites_one;
          number_neighbors[site]=3;
        }
        else if (site==0) {
           neighbors[4*site]=site+1;
           neighbors[4*site+1]=site+ist.n_sites_one;
           number_neighbors[site]=2;
        }
        else if (site==ist.n_sites_one-1) {
           neighbors[site*ist.n_sites_one]=site-1;
           neighbors[site*ist.n_sites_one+1]=site+ist.n_sites_one;
           number_neighbors[site]=2;
        }
       }
       else if (site<nsites_squared_one_two && site>=nsites_squared_one_two-ist.n_sites_one) {
        if (site!=nsites_squared_one_two-ist.n_sites_one && site!=nsites_squared_one_two-1) {
          neighbors[4*site]=site-1;
          neighbors[4*site+2]=site+1;
          neighbors[4*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=3;
        }
        else if (site==nsites_squared_one_two-ist.n_sites_one) {
          neighbors[4*site]=site+1;
          neighbors[4*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=2;
        }
        else if (site==nsites_squared_one_two-1) {
          neighbors[4*site]=site-1;
          neighbors[4*site+1]=site-ist.n_sites_one;
          number_neighbors[site]=2;
        }
       }
       else if (site%ist.n_sites_one==0) {
         if (site!=0 && site!=nsites_squared_one_two-ist.n_sites_one) {
           neighbors[4*site+2]=site+1;
           neighbors[4*site+1]=site+ist.n_sites_one;
           neighbors[4*site]=site-ist.n_sites_one;
           number_neighbors[site]=3;
         }
       }
       else if (site%ist.n_sites_one==ist.n_sites_one-1) {
         if (site!=ist.n_sites_one-1 && site!=nsites_squared_one_two-1) {
            neighbors[4*site]=site-1;
            neighbors[4*site+1]=site+ist.n_sites_one;
            neighbors[4*site+2]=site-ist.n_sites_one;
            number_neighbors[site]=3;
        }
       }
       else {
           neighbors[4*site]=site-1;
           neighbors[4*site+2]=site+1;
           neighbors[4*site+1]=site+ist.n_sites_one;
           neighbors[4*site+3]=site-ist.n_sites_one;
           number_neighbors[site]=4;
       }
    } /*If Sides Are Unequal LengthSideOne!=LengthSideTwo*/

  }/*else*/

return; 
}

/*********************************************************************************************************/

void neighborsperiodicboundary(int site, int *neighbors, int *number_neighbors,int_st ist) {
  /*Obtains the 2*dimensions nearest-neighbors directly surrounding a site*/
  /*Should be Done More Generally, But This Works for the Moment!!!*/

  int nsites_squared_one_one = ist.n_sites_one * ist.n_sites_one; 
  int nsites_squared_one_two = ist.n_sites_one * ist.n_sites_two; 

   /*Ordering Such That Double Counting Can Be Avoided In Certain Cases*/
  if ( ist.n_sites_two == 1 ) {
    if (site!=ist.n_sites_one-1 && site!=0) {
      neighbors[site*4]=site-1;
      neighbors[site*4+1]=site+1;
      number_neighbors[site]=2;
    }
    else if (site==ist.n_sites_one-1){
      neighbors[site*4]=site-1;
      neighbors[site*4+1]=0;
      number_neighbors[site]=2;
    }
    else if (site==0) {
      neighbors[site*4]=ist.n_sites_one-1;
      neighbors[site*4+1]=site+1;
      number_neighbors[site]=2;
    }
  }
  else { /*If 2 Dimensions*/ 

    if (ist.n_sites_one==ist.n_sites_two) {
      if (site<ist.n_sites_one) {
       if (site!=0 && site!=ist.n_sites_one-1){
         neighbors[4*site]=site-1;
         neighbors[4*site+2]=site+1;
         neighbors[4*site+1]=site+ist.n_sites_one;
         neighbors[4*site+3]=site+ist.n_sites_one*(ist.n_sites_one-1);
         number_neighbors[site]=4;
       }
       else if (site==0) {
          neighbors[4*site]=ist.n_sites_one-1;
          neighbors[4*site+2]=site+1;
          neighbors[4*site+1]=site+ist.n_sites_one;
          neighbors[4*site+3]=site+ist.n_sites_one*(ist.n_sites_one-1);
          number_neighbors[site]=4;
       }
       else if (site==ist.n_sites_one-1) {
          neighbors[4*site]=site-1;
          neighbors[4*site+2]=0;
          neighbors[4*site+1]=site+ist.n_sites_one;
          neighbors[4*site+3]=site+ist.n_sites_one*(ist.n_sites_one-1);
          number_neighbors[site]=4;
       }
      }
      else if (site<nsites_squared_one_one && site>=nsites_squared_one_one-ist.n_sites_one) {
        if (site!=nsites_squared_one_one-ist.n_sites_one && site!=nsites_squared_one_one-1) {
         neighbors[4*site]=site-1;
         neighbors[4*site+2]=site+1;
         neighbors[4*site+1]=site-ist.n_sites_one*(ist.n_sites_one-1);
         neighbors[4*site+3]=site-ist.n_sites_one;
         number_neighbors[site]=4;
        }
        else if (site==nsites_squared_one_one-ist.n_sites_one) {
         neighbors[4*site]=nsites_squared_one_one-1;
         neighbors[4*site+2]=site+1;
         neighbors[4*site+1]=site-ist.n_sites_one*(ist.n_sites_one-1);
         neighbors[4*site+3]=site-ist.n_sites_one;
         number_neighbors[site]=4;
       }
       else if (site==nsites_squared_one_one-1) {
         neighbors[4*site]=site-1;
         neighbors[4*site+2]=nsites_squared_one_one-ist.n_sites_one;
         neighbors[4*site+1]=site-ist.n_sites_one*(ist.n_sites_one-1);
         neighbors[4*site+3]=site-ist.n_sites_one;
         number_neighbors[site]=4;
       }
     }
     else if (site%ist.n_sites_one==0) {
        if (site!=0 && site!=nsites_squared_one_one-ist.n_sites_one) {
          neighbors[4*site]=site+(ist.n_sites_one-1);
          neighbors[4*site+2]=site+1;
          neighbors[4*site+1]=site+ist.n_sites_one;
          neighbors[4*site+3]=site-ist.n_sites_one;
          number_neighbors[site]=4;
        }
     }
     else if (site%ist.n_sites_one==ist.n_sites_one-1) {
        if (site!=ist.n_sites_one-1 && site!=nsites_squared_one_one-1) {
           neighbors[4*site]=site-1;
           neighbors[4*site+2]=site-(ist.n_sites_one-1);
           neighbors[4*site+1]=site+ist.n_sites_one;
           neighbors[4*site+3]=site-ist.n_sites_one;
           number_neighbors[site]=4;
      }
     }
     else {
          neighbors[4*site]=site-1;
          neighbors[4*site+2]=site+1;
          neighbors[4*site+1]=site+ist.n_sites_one;
          neighbors[4*site+3]=site-ist.n_sites_one;
          number_neighbors[site]=4;
     }
    } /*If LengthSideOne==LengthSideTwo*/
    else {  /*I Assume Length Side One Is Horizontal and LengthSide Two Is Vertical*/
       if (site<ist.n_sites_one) {
        if (site!=0 && site!=ist.n_sites_one-1){
         neighbors[4*site]=site-1;
         neighbors[4*site+2]=site+1;
         neighbors[4*site+1]=site+ist.n_sites_one;
         neighbors[4*site+3]=site+ist.n_sites_one*(ist.n_sites_two-1);
         number_neighbors[site]=4;
       }
       else if (site==0) {
          neighbors[4*site]=ist.n_sites_one-1;
          neighbors[4*site+2]=site+1;
          neighbors[4*site+1]=site+ist.n_sites_one;
          neighbors[4*site+3]=site+ist.n_sites_one*(ist.n_sites_two-1);
          number_neighbors[site]=4;
       }
       else if (site==ist.n_sites_one-1) {
          neighbors[4*site]=site-1;
          neighbors[4*site+2]=0;
          neighbors[4*site+1]=site+ist.n_sites_one;
          neighbors[4*site+3]=site+ist.n_sites_one*(ist.n_sites_two-1);
          number_neighbors[site]=4;
       }
      }
      else if (site<nsites_squared_one_two && site>=nsites_squared_one_two-ist.n_sites_one) {
       if (site!=nsites_squared_one_two-ist.n_sites_one && site!=nsites_squared_one_two-1) {
         neighbors[4*site]=site-1;
         neighbors[4*site+2]=site+1;
         neighbors[4*site+1]=site-ist.n_sites_one*(ist.n_sites_two-1);
         neighbors[4*site+3]=site-ist.n_sites_one;
         number_neighbors[site]=4;
       }
       else if (site==nsites_squared_one_two-ist.n_sites_one) {
         neighbors[4*site]=nsites_squared_one_two-1;
         neighbors[4*site+2]=site+1;
         neighbors[4*site+1]=site-ist.n_sites_one*(ist.n_sites_two-1);
         neighbors[4*site+3]=site-ist.n_sites_one; 
         number_neighbors[site]=4;
       }
       else if (site==nsites_squared_one_two-1) {
         neighbors[4*site]=site-1;
         neighbors[4*site+2]=nsites_squared_one_two-ist.n_sites_one;
         neighbors[4*site+1]=site-ist.n_sites_one*(ist.n_sites_two-1);
         neighbors[4*site+3]=site-ist.n_sites_one;
         number_neighbors[site]=4;
       }
      }
      else if (site%ist.n_sites_one==0) {
        if (site!=0 && site!=nsites_squared_one_two-ist.n_sites_one) {
          neighbors[4*site]=site+(ist.n_sites_one-1);
          neighbors[4*site+2]=site+1;
          neighbors[4*site+1]=site+ist.n_sites_one;
          neighbors[4*site+3]=site-ist.n_sites_one;
          number_neighbors[site]=4;
        }
      }
      else if (site%ist.n_sites_one==ist.n_sites_one-1) {
        if (site!=ist.n_sites_one-1 && site!=nsites_squared_one_two-1) {
           neighbors[4*site]=site-1;
           neighbors[4*site+2]=site-(ist.n_sites_one-1);
           neighbors[4*site+1]=site+ist.n_sites_one;
           neighbors[4*site+3]=site-ist.n_sites_one;
           number_neighbors[site]=4;
       }
      }
      else {
          neighbors[4*site]=site-1;
          neighbors[4*site+2]=site+1;
          neighbors[4*site+1]=site+ist.n_sites_one;
          neighbors[4*site+3]=site-ist.n_sites_one;
          number_neighbors[site]=4;
      }
    } /*If Sides Are Unequal LengthSideOne!=LengthSideTwo*/
 
  } /*If 2 Dimensions*/ 

return; 
}

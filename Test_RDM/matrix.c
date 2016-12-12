#include "afqmc.h"

/********************************************************************************************************/

void jacobireal(double *matrix, double *eigenvalues, double *eigenvectors, int size, int_st ist) {
     /*Eigenvalues Matrix*/

     int j,iq,ip,i;
     double tresh,theta,tau,t,sm,s,h,g,c; 
     double *b,*z;
 
     b=(double *)calloc(size,sizeof(double)); 
     z=(double *)calloc(size,sizeof(double)); 

     dzero_vec(eigenvalues,size); 
     dzero_vec(eigenvectors,size*size); 

     /*Numerical Recipes Jacobi Routine*/
     for (ip=1;ip<=size;ip++) {
                for (iq=1;iq<=size;iq++) eigenvectors[(ip-1)*size+(iq-1)]=0.0;
                eigenvectors[(ip-1)*size+(ip-1)]=1.0;
        }
        for (ip=1;ip<=size;ip++) {
                b[ip-1]=eigenvalues[ip-1]=matrix[(ip-1)*size+(ip-1)];
                z[ip-1]=0.0;
        }
        for (i=1;i<200;i++) {
                sm=0.0;
                for (ip=1;ip<=size-1;ip++) {
                        for (iq=ip+1;iq<=size;iq++)
                                sm += fabs(matrix[(ip-1)*size+(iq-1)]);
                }
                if (sm == 0.0) {
                        free(b); 
                        free(z);
                        return;
                }
                if (i < 4)
                        tresh=0.2*sm/((double) size*size);
                else
                        tresh=0.0;
                for (ip=1;ip<=size-1;ip++) {
                        for (iq=ip+1;iq<=size;iq++) {
                                g=100.0*fabs(matrix[(ip-1)*size+iq-1]);
                                if (i > 4 && (double)(fabs(eigenvalues[ip-1])+g) == (double)fabs(eigenvalues[ip-1])
                                        && (double)(fabs(eigenvalues[iq-1])+g) == (double)fabs(eigenvalues[iq-1]))
                                        matrix[(ip-1)*size+iq-1]=0.0;
                                else if (fabs(matrix[(ip-1)*size+iq-1]) > tresh) {
                                        h=eigenvalues[iq-1]-eigenvalues[ip-1];
                                        if ((double)(fabs(h)+g) == (double)fabs(h))
                                                t=(matrix[(ip-1)*size+(iq-1)])/h;
                                        else {
                                                theta=0.5*h/(matrix[(ip-1)*size+iq-1]);
                                                t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                                                if (theta < 0.0) t = -t;
                                        }
                                        c=1.0/sqrt(1+t*t);
                                        s=t*c;
                                        tau=s/(1.0+c);
                                        h=t*matrix[(ip-1)*size+iq-1];
                                        z[ip-1] -= h;
                                        z[iq-1] += h;
                                        eigenvalues[ip-1] -= h;
                                        eigenvalues[iq-1] += h;
                                        matrix[(ip-1)*size+iq-1]=0.0;
                                        for (j=1;j<=ip-1;j++) {
                                           g=matrix[(j-1)*size+ip-1];
                                           h=matrix[(j-1)*size+iq-1];
                                           matrix[(j-1)*size+ip-1]=g-s*(h+g*tau);
                                           matrix[(j-1)*size+iq-1]=h+s*(g-h*tau);
                                        }
                                        for (j=ip+1;j<=iq-1;j++) {
                                                g=matrix[(ip-1)*size+j-1];
                                                h=matrix[(j-1)*size+iq-1];
                                                matrix[(ip-1)*size+j-1]=g-s*(h+g*tau);
                                                matrix[(j-1)*size+iq-1]=h+s*(g-h*tau);
                                        }
                                        for (j=iq+1;j<=size;j++) {
                                                g=matrix[(ip-1)*size+j-1];
                                                h=matrix[(iq-1)*size+j-1];
                                                matrix[(ip-1)*size+j-1]=g-s*(h+g*tau);
                                                matrix[(iq-1)*size+j-1]=h+s*(g-h*tau);
                                        }
                                        for (j=1;j<=size;j++) {
                                                g=eigenvectors[(j-1)*size+ip-1];
                                                h=eigenvectors[(j-1)*size+iq-1];
                                                eigenvectors[(j-1)*size+ip-1]=g-s*(h+g*tau);
                                                eigenvectors[(j-1)*size+iq-1]=h+s*(g-h*tau);
                                        }
                                }
                        }
                }
                for (ip=1;ip<=size;ip++) {
                        b[ip-1] += z[ip-1];
                        eigenvalues[ip-1]=b[ip-1];
                        z[ip-1]=0.0;
                }
        }


free(b); 
free(z); 
return; 
}

/***************************************************************************************************************/

void complex_eigenvalues(MKL_Complex16 *matrix,MKL_Complex16 *eigenvalues,MKL_Complex16 *eigenvectors,int size,int_st ist) {
  /*Routine to Calculate Complex Eigenvalues of Complex Matrix*/
  /*From Wilkinson and Reinsch Book*/

  int *integer;

  integer=(int *)calloc((1+size),sizeof(int ));

  /*Comhes*/
  comhes(matrix,integer,size,ist);

  /*Then Run Through Comlr2*/
  comlr2(matrix,eigenvalues,eigenvectors,integer,size,ist);

free(integer);
return;
}

/**************************************************************************************************************/

void comlr2(MKL_Complex16 *hess, MKL_Complex16 *eigs, MKL_Complex16 *eigvecs, int *integer, int size, int_st ist) {
/*Determines Eigs and Eigenvectors from Hessenberg Mat*/

 int i, j, k, m, its, en;
 double sr, si, tr, ti, xr, xi, yr, yi, zr, zi, norm;
 MKL_Complex16 value1, value2;
 MKL_Complex16 x;
 MKL_Complex16 hold;

  tr=ti=0;

  cunity_vec(eigvecs,size);

  for (i=size-1; i>=2; i--) {
    j=integer[i-1];
    for (k=i+1; k<=size; k++) {
      eigvecs[(k-1)*size+i-1].real = hess[(k-1)*size+i-2].real;
      eigvecs[(k-1)*size+i-1].imag = hess[(k-1)*size+i-2].imag;
    }
    if (i!=j) {
      for (k=i; k<=size; k++) {
        eigvecs[(i-1)*size+k-1].real = eigvecs[(j-1)*size+k-1].real;
        eigvecs[(i-1)*size+k-1].imag = eigvecs[(j-1)*size+k-1].imag;
 
        eigvecs[(j-1)*size+k-1].real = eigvecs[(j-1)*size+k-1].imag = 0;
      }
      eigvecs[(j-1)*size+i-1].real = 1.0;
    }
  }

 en=size;
nextw:
 if (en<1) {
      goto fin;
   }
   its=0;


 nextit:
   for (k=en; k>=2; k--) {
     if ( fabs(hess[(k-1)*size+k-2].real) + fabs(hess[(k-1)*size+k-2].imag) <= EPS*( fabs(hess[(k-2)*size+k-2].real) + fabs(hess[(k-2)*size+k-2].imag) + fabs(hess[(k-1)*size+k-1].real) + fabs(hess[(k-1)*size+k-1].imag) )) {
         goto cont1;
     }
    }
    k=1;

 cont1:
    if (k==en) {
      goto root;
    }
    if (its==30) {
    }

    if (its==10 || its==20) {
      sr=fabs(hess[(en-1)*size+en-2].real)+fabs(hess[(en-2)*size+en-3].real);
      si=fabs(hess[(en-1)*size+en-2].imag)+fabs(hess[(en-2)*size+en-3].imag);
    }
    else {
      sr=hess[(en-1)*size+en-1].real;
      si=hess[(en-1)*size+en-1].imag;
      xr=hess[(en-2)*size+en-1].real*hess[(en-1)*size+en-2].real-hess[(en-2)*size+en-1].imag*hess[(en-1)*size+en-2].imag;
      xi=hess[(en-2)*size+en-1].real*hess[(en-1)*size+en-2].imag+hess[(en-2)*size+en-1].imag*hess[(en-1)*size+en-2].real;
      if ( xr!=0 || xi!=0) {
         yr= (hess[(en-2)*size+en-2].real-sr)/2.0;
         yi=(hess[(en-2)*size+en-2].imag-si)/2.0;
         value1.real=yr*yr-yi*yi+xr;
         value1.imag=2*yr*yi+xi;
         value2=Csqrt(value1);
         zr=value2.real;
         zi=value2.imag;
         if( yr*zr+yi*zi<0) {
            zr*=-1;
            zi*=-1;
         }
         value1.real=yr+zr;
         value1.imag=yi+zi;
         x.real=xr;
         x.imag=xi;
         value2=Cdiv(x, value1);
         xr=value2.real;
         xi=value2.imag;
         sr=sr-xr;
         si=si-xi;
      }
    }
    for(i=1; i<=en; i++) {
      hess[(i-1)*size+i-1].real=hess[(i-1)*size+i-1].real-sr;
      hess[(i-1)*size+i-1].imag=hess[(i-1)*size+i-1].imag-si;
    }

    tr=tr+sr;
    ti=ti+si;
    its++;
    j=k+1;

    xr=fabs(hess[(en-2)*size+en-2].real)+fabs(hess[(en-2)*size+en-2].imag);
    yr=fabs(hess[(en-1)*size+en-2].real)+fabs(hess[(en-1)*size+en-2].imag);
    zr=fabs(hess[(en-1)*size+en-1].real)+fabs(hess[(en-1)*size+en-1].imag);

for (m=en-1; m>=j; m--) {
       yi=yr;
       yr=fabs(hess[(m-1)*size+m-2].real)+fabs(hess[(m-1)*size+m-2].imag);
       xi=zr;
       zr=xr;
       xr=fabs(hess[(m-2)*size+m-2].real)+fabs(hess[(m-2)*size+m-2].imag);
       if ( yr<=EPS*zr/yi*(zr+xr+xi)) {
         goto cont2;
       }
    }
    m=k;

    cont2:
    for (i=m+1; i<=en; i++) {
          xr=hess[(i-2)*size+i-2].real;
          xi=hess[(i-2)*size+i-2].imag;
          yr=hess[(i-1)*size+i-2].real;
          yi=hess[(i-1)*size+i-2].imag;
          if ( (fabs(xr)+fabs(xi) ) < ( fabs(yr) + fabs(yi) ) ) {
             for (j=i-1; j<=size; j++) {
              zr=hess[(i-2)*size+j-1].real;
              hess[(i-2)*size+j-1].real=hess[(i-1)*size+j-1].real;
              hess[(i-1)*size+j-1].real=zr;
              zi=hess[(i-2)*size+j-1].imag;
              hess[(i-2)*size+j-1].imag=hess[(i-1)*size+j-1].imag;
              hess[(i-1)*size+j-1].imag=zi;
             }
           value1.real=xr;
           value1.imag=xi;
           value2.real=yr;
           value2.imag=yi;
           x=Cdiv(value1, value2);
           zr=x.real;
           zi=x.imag;

           eigs[i-1].real=1;
         }
         else {
           value1.real=yr;
           value1.imag=yi;
           value2.real=xr;
           value2.imag=xi;
           x=Cdiv(value1, value2);
           zr=x.real;
           zi=x.imag;
           eigs[i-1].real=-1;
         }
         hess[(i-1)*size+i-2].real=zr;
         hess[(i-1)*size+i-2].imag=zi;

         for (j=i; j<=size; j++) {
            hess[(i-1)*size+j-1].real=hess[(i-1)*size+j-1].real-zr*hess[(i-2)*size+j-1].real+zi*hess[(i-2)*size+j-1].imag;
            hess[(i-1)*size+j-1].imag=hess[(i-1)*size+j-1].imag-zr*hess[(i-2)*size+j-1].imag-zi*hess[(i-2)*size+j-1].real;
         }
       } /*i*/


       for (j=m+1; j<=en; j++) {
          xr=hess[(j-1)*size+j-2].real;
          xi=hess[(j-1)*size+j-2].imag;
          hess[(j-1)*size+j-2].real=hess[(j-1)*size+j-2].imag=0;

          if (eigs[j-1].real>0) {
             for (i=1; i<=j; i++) {
                zr=hess[(i-1)*size+j-2].real;
                hess[(i-1)*size+j-2].real=hess[(i-1)*size+j-1].real;
                hess[(i-1)*size+j-1].real=zr;

                zi=hess[(i-1)*size+j-2].imag;

                hess[(i-1)*size+j-2].imag=hess[(i-1)*size+j-1].imag;
                hess[(i-1)*size+j-1].imag=zi;
             }

             for (i=1; i<=size; i++) {
                zr=eigvecs[(i-1)*size+j-2].real;
                eigvecs[(i-1)*size+j-2].real=eigvecs[(i-1)*size+j-1].real;
                eigvecs[(i-1)*size+j-1].real=zr;

                zi=eigvecs[(i-1)*size+j-2].imag;
                eigvecs[(i-1)*size+j-2].imag=eigvecs[(i-1)*size+j-1].imag;
                eigvecs[(i-1)*size+j-1].imag=zi;
             }
          }

         for (i=1; i<=j; i++) {
           hess[(i-1)*size+j-2].real=hess[(i-1)*size+j-2].real+xr*hess[(i-1)*size+j-1].real-xi*hess[(i-1)*size+j-1].imag;
           hess[(i-1)*size+j-2].imag=hess[(i-1)*size+j-2].imag+xr*hess[(i-1)*size+j-1].imag+xi*hess[(i-1)*size+j-1].real;
         }

         for (i=1; i<=size; i++) {
           eigvecs[(i-1)*size+j-2].real=eigvecs[(i-1)*size+j-2].real+xr*eigvecs[(i-1)*size+j-1].real-xi*eigvecs[(i-1)*size+j-1].imag;
           eigvecs[(i-1)*size+j-2].imag=eigvecs[(i-1)*size+j-2].imag+xr*eigvecs[(i-1)*size+j-1].imag+xi*eigvecs[(i-1)*size+j-1].real;
         }

       } /*j*/

       goto nextit;

root:

      eigs[en-1].real=hess[(en-1)*size+en-1].real+tr;
      eigs[en-1].imag=hess[(en-1)*size+en-1].imag+ti;
      en=en-1;
      goto nextw;
fin:

     norm=0.0;
     for (i=1; i<=size; i++) {
        norm=norm+fabs(eigs[i-1].real)+fabs(eigs[i-1].imag);
        for (j=i+1; j<=size; j++) {
          norm=norm+fabs(hess[(i-1)*size+j-1].real)+fabs(hess[(i-1)*size+j-1].imag);
         }
      }

      for (en=size; en>=2; en--) {

         xr=eigs[en-1].real;
         xi=eigs[en-1].imag;
         for (i=en-1; i>=1; i--) {
           zr=hess[(i-1)*size+en-1].real;
           zi=hess[(i-1)*size+en-1].imag;
           for(j=i+1; j<=en-1; j++) {
            zr=zr+hess[(i-1)*size+j-1].real*hess[(j-1)*size+en-1].real-hess[(i-1)*size+j-1].imag*hess[(j-1)*size+en-1].imag;
            zi=zi+hess[(i-1)*size+j-1].real*hess[(j-1)*size+en-1].imag+hess[(i-1)*size+j-1].imag*hess[(j-1)*size+en-1].real;
           }
           yr=xr-eigs[i-1].real;
           yi=xi-eigs[i-1].imag;
           if ( yr==0.0 && yi==0.0 ) {
              yr=EPS*norm;
           }
           value1.real=zr;
           value1.imag=zi;
           value2.real=yr;
           value2.imag=yi;
           hold = Cdiv(value1,value2);
           hess[(i-1)*size+en-1].real = hold.real;
           hess[(i-1)*size+en-1].imag = hold.imag;
         }
       }
       for (j=size; j>=1; j--) {
        for (i=1; i<=size; i++) {

          zr=eigvecs[(i-1)*size+j-1].real;
          zi=eigvecs[(i-1)*size+j-1].imag;
          if (size<j) {
            m=size;
          }
          else {
            m=j-1;
          }

          for (k=1; k<=m; k++) {
              zr=zr+eigvecs[(i-1)*size+k-1].real*hess[(k-1)*size+j-1].real-eigvecs[(i-1)*size+k-1].imag*hess[(k-1)*size+j-1].imag;
              zi=zi+eigvecs[(i-1)*size+k-1].real*hess[(k-1)*size+j-1].imag+eigvecs[(i-1)*size+k-1].imag*hess[(k-1)*size+j-1].real;
           }
           eigvecs[(i-1)*size+j-1].real=zr;
           eigvecs[(i-1)*size+j-1].imag=zi;
        }
 }

return;
}

/***************************************************************************************************************/

void comhes(MKL_Complex16 *matrix,int *integer,int size,int_st ist) {
 /*Reduces the Matrix to Upper Hessenberg Form*/

 int i, j, la, m;
 double xr, xi, yr, yi;
 MKL_Complex16 x, y;

 la=size-1;
 for (m=2; m<=la; m++) {
    i=m;
    xr=xi=0;
    for (j=m; j<=size; j++) {
      if ( (fabs(matrix[(j-1)*size+m-2].real) + fabs(matrix[(j-1)*size+m-2].imag)) > (fabs(xr) + fabs(xi)) ) {
         xr=matrix[(j-1)*size+m-2].real;
         xi=matrix[(j-1)*size+m-2].imag;
         i=j;
      }
     }
     integer[m-1]=i;

      if (i!=m) {
         for (j=m-1; j<=size; j++) {
             yr=matrix[(i-1)*size+j-1].real;
             matrix[(i-1)*size+j-1].real=matrix[(m-1)*size+j-1].real;
             matrix[(m-1)*size+j-1].real=yr;

             yi=matrix[(i-1)*size+j-1].imag;
             matrix[(i-1)*size+j-1].imag=matrix[(m-1)*size+j-1].imag;
             matrix[(m-1)*size+j-1].imag=yi;
         }
         for(j=1; j<=size; j++) {
            yr=matrix[(j-1)*size+i-1].real;
            matrix[(j-1)*size+i-1].real=matrix[(j-1)*size+m-1].real;
            matrix[(j-1)*size+m-1].real=yr;
            yi=matrix[(j-1)*size+i-1].imag;
            matrix[(j-1)*size+i-1].imag=matrix[(j-1)*size+m-1].imag;
            matrix[(j-1)*size+m-1].imag=yi;
         }
       }

       if ( xr!=0 || xi!=0) {
          for (i=m+1; i<=size; i++) {
             yr=matrix[(i-1)*size+m-2].real;
             yi=matrix[(i-1)*size+m-2].imag;

             y.real=yr;
             y.imag=yi;
             x.real=xr;
             x.imag=xi;
             y=Cdiv(y, x);
             yr=y.real;
             yi=y.imag;
             matrix[(i-1)*size+m-2].real=yr;
             matrix[(i-1)*size+m-2].imag=yi;
             for (j=m; j<=size; j++) {
                  matrix[(i-1)*size+j-1].real=matrix[(i-1)*size+j-1].real-yr*matrix[(m-1)*size+j-1].real+yi*matrix[(m-1)*size+j-1].imag;
                  matrix[(i-1)*size+j-1].imag=matrix[(i-1)*size+j-1].imag-yr*matrix[(m-1)*size+j-1].imag-yi*matrix[(m-1)*size+j-1].real;
                }
                for (j=1; j<=size; j++) {
                   matrix[(j-1)*size+m-1].real+=yr*matrix[(j-1)*size+i-1].real-yi*matrix[(j-1)*size+i-1].imag;
                   matrix[(j-1)*size+m-1].imag+=yr*matrix[(j-1)*size+i-1].imag+yi*matrix[(j-1)*size+i-1].real;
                }
             }

        }
   }

return;
}

/***************************************************************************************************************/

void dmat_dmat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3) {

    /*Multiplies Two Matrices Together, Where One is of Size size_1xsize_2 and the other Is of Size size_2xsize_3*/

    int i, j, k; 
    dzero_vec(product,size_1*size_3);   

    for (i=0; i<size_1; i++) {
      for (j=0; j<size_3; j++) {
        product[i*size_3+j]=0.0; 
        
        for (k=0; k<size_2; k++) {
          product[i*size_3+j] += mat_1[i*size_2+k] * mat_2[k*size_3+j]; 
        }

      }
    }

return; 
}

/***************************************************************************************************************/

void cmat_cmat(MKL_Complex16 *mat_1, MKL_Complex16 *mat_2, MKL_Complex16 *product, int size_1, int size_2, int size_3) {

    /*Multiplies Two Matrices Together, Where One is of Size size_1xsize_2 and the other Is of Size size_2xsize_3*/

    int i, j, k;
    MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0; 
    czero_vec(product,size_1*size_3);

    for (i=0; i<size_1; i++) {
      for (j=0; j<size_3; j++) {
        product[i*size_3+j]=Zero;

        for (k=0; k<size_2; k++) {
          product[i*size_3+j] = Cadd(product[i*size_3+j], Cmul(mat_1[i*size_2+k], mat_2[k*size_3+j])); 
        }

      }
    }

return;
}

/***************************************************************************************************************/

void dmat_cmat(double *mat_1, MKL_Complex16 *mat_2, MKL_Complex16 *product, int size_1, int size_2, int size_3) {

    /*Multiplies Two Matrices Together, Where One is of Size size_1xsize_2 and the other Is of Size size_2xsize_3*/

    int i, j, k;
    MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0; 
    czero_vec(product,size_1*size_3);

    for (i=0; i<size_1; i++) {
      for (j=0; j<size_3; j++) {
        product[i*size_3+j]=Zero;

        for (k=0; k<size_2; k++) {
          product[i*size_3+j] = Cadd(product[i*size_3+j], RCmul(mat_1[i*size_2+k], mat_2[k*size_3+j]));   
        }

      }
    }

return;
}

/***************************************************************************************************************/

void mat_transpose_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3) {

    /*Multiplies Two Matrices Together, Where One is of Size size_1xsize_2 and the other Is of Size size_2xsize_3*/
    /*Note That I am Not Taking the Size of the Second Matrix As the Second Size and Third Size*/ 

    int i, j, k;
    dzero_vec(product,size_1*size_3);

    for (i=0; i<size_1; i++) {
      for (j=0; j<size_3; j++) {
        product[i*size_3+j]=0.0;

        for (k=0; k<size_2; k++) {
          product[i*size_3+j] += mat_1[i*size_1+k] * mat_2[j*size_2+k];
        }

      }
    }

return;
}

/***************************************************************************************************************/

void cmat_transpose_cmat(MKL_Complex16 *mat_1, MKL_Complex16 *mat_2, MKL_Complex16 *product, int size_1, int size_2, int size_3) {

    /*Multiplies Two Matrices Together, Where One is of Size size_1xsize_2 and the other Is of Size size_2xsize_3*/
    /*Note That I am Not Taking the Size of the Second Matrix As the Second Size and Third Size*/

    int i, j, k;
    MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;
    czero_vec(product,size_1*size_3);

    for (i=0; i<size_1; i++) {
      for (j=0; j<size_3; j++) {
        product[i*size_3+j]=Zero;

        for (k=0; k<size_2; k++) {
          product[i*size_3+j] = Cadd(product[i*size_3+j], Cmul(mat_1[i*size_2+k], mat_2[j*size_2+k]));
        }

      }
    }

return;
}


/***************************************************************************************************************/

void transpose_dmat_dmat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3) {

     /*Multiplies the Transpose of the First Matrix with the Second, As Needed to Compute Overlaps*/

     int i, j, k; 

     for (i=0; i<size_2; i++) {
      for (j=0; j<size_3; j++) {
        product[i*size_3+j]=0.0;

        for (k=0; k<size_1; k++) {
          product[i*size_3+j] += mat_1[k*size_2+i] * mat_2[k*size_3+j];
        }
      }
    }

return; 
}

/***************************************************************************************************************/

void transpose_cmat_cmat(MKL_Complex16 *mat_1, MKL_Complex16 *mat_2, MKL_Complex16 *product, int size_1, int size_2, int size_3) {

     /*Multiplies the Transpose of the First Matrix with the Second, As Needed to Compute Overlaps*/

     int i, j, k;
     MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0;

     for (i=0; i<size_2; i++) {
      for (j=0; j<size_3; j++) {
        product[i*size_3+j]=Zero;

        for (k=0; k<size_1; k++) {
          product[i*size_3+j] = Cadd(product[i*size_3+j], Cmul(mat_1[k*size_2+i], mat_2[k*size_3+j]));
        }
      }
    }

return;
}

/***************************************************************************************************************/

void copy_mat(double *mat_1, double *mat_2, int size) {

   int i; 

   for (i=0; i<size; i++) {
      mat_1[i] = mat_2[i]; 
   }

return; 
}

/***************************************************************************************************************/

void copy_cmat(MKL_Complex16 *mat_1, MKL_Complex16 *mat_2, int size) {

   int i;

   for (i=0; i<size; i++) {
      mat_1[i] = mat_2[i];
   }

return;
}

/***************************************************************************************************************/

void det(double *matrix, int size, double *determinant) {

    int i; 
    int *indx;   
   
    indx=(int *)calloc(size,sizeof(int)); 

    /*First Perform LU Decomposition*/
    ludmp(matrix,size,indx,determinant); 
 
    /*Now Multiply By Diagaonl*/
    for (i=0; i<size; i++) {
      (*determinant) *= matrix[i*size+i]; 
    }

free(indx); 
return; 
}

/***************************************************************************************************************/

void complex_det(MKL_Complex16 *matrix, int size, MKL_Complex16 *determinant) {

    int i;
    int *indx;

    indx=(int *)calloc(size,sizeof(int));

    /*First Perform LU Decomposition*/
    complex_ludmp(matrix,size,indx,determinant);
 

    /*Now Multiply By Diagaonl*/
    for (i=0; i<size; i++) {
      (*determinant) = Cmul( (*determinant), matrix[i*size+i]);
    }

free(indx);
return;
}

/***************************************************************************************************************/

void inverse(double *matrix, double *inverse, int size) {

     int i,j; 
     int *indx; 
     double *col; 
     double d;    

     indx=(int *)calloc(size,sizeof(int));
     col=(double *)calloc(size,sizeof(double));

     ludmp(matrix,size,indx,&d); 

     for (j=0; j<size; j++) {
        for (i=0; i<size; i++) {
           col[i]=0.0; 
        }
        col[j]=1.0; 

        lubksb(matrix,size,indx,col); 
        for(i=0; i<size; i++) {
          inverse[i*size+j]=col[i]; 
        }
    } 

free(col); 
free(indx); 
return; 
}

/***************************************************************************************************************/

double inverse_det(double *matrix, double *inverse, int size) {

  int i,j;
  int *indx;
  double *col;
  double determinant;

  /*Gets Inverse and Determinant*/ 

  indx=(int *)calloc(size,sizeof(int));
  col=(double *)calloc(size,sizeof(double));

  ludmp(matrix,size,indx,&determinant);
  for (j=0; j<size; j++) {
    determinant *= matrix[j*size+j];  

    for (i=0; i<size; i++) {
        col[i]=0.0;
     }
     col[j]=1.0;
     lubksb(matrix,size,indx,col);

     for(i=0; i<size; i++) {
       inverse[i*size+j]=col[i];
      }
    }

free(col); 
free(indx); 
return(determinant);  
}

/***************************************************************************************************************/

MKL_Complex16 complex_inverse_det_fermions(MKL_Complex16 *matrix, MKL_Complex16 *inverse, int size) {
   /*Should Totally be Able To Combine These Two Operations!!!*/

  int j;
  int *indx;
  MKL_Complex16 determinant;

  /*Gets Inverse and Determinant*/
  indx=(int *)calloc(size,sizeof(int));

  /*Get Inverse First*/
  complex_matrix_inverse(matrix,inverse,size);

  /*Then Get Det*/ 
  complex_ludmp(matrix,size,indx,&determinant);

  for (j=0; j<size; j++) {
    determinant = Cmul(determinant, matrix[j*size+j]);
  }
 
free(indx);
return(determinant);
}

/***************************************************************************************************************/

MKL_Complex16 perm_bosons(MKL_Complex16 *vector_1, MKL_Complex16 *vector_2, int_st ist) {
  /*Compute the Permanent for Two Row Vectors*/

  int j;
  MKL_Complex16 permanent; 

  permanent.real = permanent.imag = 0;  
  for (j=0; j<ist.n_sites; j++) {
    permanent = Cadd(permanent, Cmul(conjugate(vector_1[j]), vector_2[j])); 
  }

return(permanent);
}

/***************************************************************************************************************/

void ludmp(double *matrix, int n, int *indx, double *d) { 

   /*Performs an LU Decomposition Based on Numerical Recipes*/

   int i,imax=1,j,k; 
   double big,dum,sum,temp; 
   double *vv; 
 
   vv=(double*)calloc(n,sizeof(double)); 

   *d=1.0; 
   for (i=1; i<=n; i++){
     big=0.0;
     for (j=1; j<=n; j++) {
        if ((temp=fabs(matrix[(i-1)*n+j-1])) > big) big=temp; 
     }
     vv[i-1]=1.0/big; 
    }

    for (j=1; j<=n; j++) {

      for (i=1; i<j; i++) {
         sum=matrix[(i-1)*n+j-1]; 
         for (k=1; k<i; k++) { 
           sum -= matrix[(i-1)*n+k-1] * matrix[(k-1)*n+j-1]; 
         }
         matrix[(i-1)*n+j-1]=sum; 
       }

    big=0.0; 
    for (i=j; i<=n; i++) {
       sum=matrix[(i-1)*n+j-1]; 
       for (k=1; k<j; k++) {
         sum -= matrix[(i-1)*n+k-1]*matrix[(k-1)*n+j-1]; 
       }
       matrix[(i-1)*n+j-1]=sum; 
      
       if ( (dum=vv[i-1]*fabs(sum)) >= big) {
          big = dum; 
          imax = i; 
       }
     }

     if ( j!= imax ) {
        for (k=1; k<=n; k++) {
           dum = matrix[(imax-1)*n+k-1]; 
           matrix[(imax-1)*n+k-1] = matrix[(j-1)*n+k-1]; 
           matrix[(j-1)*n+k-1]=dum; 
         }
         (*d) = -(*d); 
         vv[imax-1]=vv[j-1]; 
     }
     indx[j-1]=imax; 

     if (matrix[(j-1)*n+j-1] == 0.0 ) { 
         matrix[(j-1)*n+j-1] = TINY; 
     }

     if (j != n ) {
        dum=1.0/matrix[(j-1)*n+j-1]; 
        for (i=j+1; i<=n; i++) {
          matrix[(i-1)*n+j-1] *= dum; 
        }
      }

   } 

free(vv); 
return; 
}

/********************************************************************************/

void complex_ludmp(MKL_Complex16 *matrix, int n, int *indx, MKL_Complex16 *d) {

   /*Performs an LU Decomposition Based on Numerical Recipes*/

   int i,imax,j,k;
   double big,dum,temp;
   double *vv;
   MKL_Complex16 dum2;
   MKL_Complex16 sum;
   MKL_Complex16 One, Negative_One;

   One.real = 1.0;
   One.imag = 0.0;

   Negative_One.real = -1.0;
   Negative_One.imag = 0.0;

   vv=(double*)calloc(n,sizeof(double));

   *d=One;
   for (i=1; i<=n; i++){
     big=0.0;
     for (j=1; j<=n; j++) {
        if ((temp=Cabs(matrix[(i-1)*n+j-1])) > big) big=temp;
     }
     vv[i-1]=1.0/big;
    }

    for (j=1; j<=n; j++) {
      for (i=1; i<j; i++) {
         sum=matrix[(i-1)*n+j-1];
         for (k=1; k<i; k++) {
           sum = Csub(sum, Cmul(matrix[(i-1)*n+k-1], matrix[(k-1)*n+j-1]));
         }
         matrix[(i-1)*n+j-1]=sum;
       }


    big=0.0;
    for (i=j; i<=n; i++) {
       sum=matrix[(i-1)*n+j-1];
       for (k=1; k<j; k++) {
         sum = Csub(sum,Cmul(matrix[(i-1)*n+k-1],matrix[(k-1)*n+j-1]));
       }
       matrix[(i-1)*n+j-1]=sum;
  if ( (dum=vv[i-1]*Cabs(sum)) >= big) {
          big = dum;
          imax = i;
       }
     }

     if ( j!= imax ) {
        for (k=1; k<=n; k++) {
           dum2 = matrix[(imax-1)*n+k-1];
           matrix[(imax-1)*n+k-1] = matrix[(j-1)*n+k-1];
           matrix[(j-1)*n+k-1]=dum2;
         }
         (*d) = Cmul(Negative_One, (*d));
         vv[imax-1]=vv[j-1];
     }
     indx[j-1]=imax;

     if (Cabs(matrix[(j-1)*n+j-1]) == 0 ){
         matrix[(j-1)*n+j-1].real = TINY;
         matrix[(j-1)*n+j-1].imag = 0;
     }

     if (j != n ) {
        dum2=Cdiv(One,matrix[(j-1)*n+j-1]);
        for (i=j+1; i<=n; i++) {
          matrix[(i-1)*n+j-1] =Cmul(matrix[(i-1)*n+j-1],dum2);
      }
    }
   }

free(vv);
return;
}

/***************************************************************************************************************/

void lubksb(double *matrix, int n, int *indx, double *b) {

    /*Peforms LU Forward and Backward Substitution From NRC*/
    
    int i,ii=0,ip,j; 
    double sum; 

    for (i=1; i<=n; i++) {
       ip=indx[i-1]; 
       sum=b[ip-1]; 
       b[ip-1]=b[i-1]; 
      
        if (ii) 
              for (j=ii; j<=i-1; j++) sum -= matrix[(i-1)*n+j-1] * b[j-1]; 
        else if (sum) ii=i; 
        b[i-1]=sum; 
     }

     for (i=n; i>=1; i--) {
        sum=b[i-1]; 
        for (j=i+1; j<=n; j++) sum -= matrix[(i-1)*n+j-1]*b[j-1]; 
        b[i-1]=sum/matrix[(i-1)*n+i-1]; 
     }

return; 
}

/***************************************************************************************************************/

void modified_gram_schmidt(double *q, double *r, int size1, int size2) {

   /*The Modified Gram-Schmidt Routine Used to Stabilize Matrix Products*/
   int i, j, k;   
   double temporary, anorm; 
   double *d; 

   d=(double *)calloc(size2,sizeof(double)); 
   dzero_vec(r,size2*size2); 

   for (i=1; i<=size2; i++) {

     temporary = 0.0;
     for (j=1; j<=size1; j++) {
       temporary += q[(j-1)*size2+i-1]*q[(j-1)*size2+i-1]; 
     }
     d[i-1]=sqrt(temporary); 
     anorm = 1.0/d[i-1];   

     for (j=1; j<=size1; j++) {
       q[(j-1)*size2+i-1] *= anorm; 
     }

     for ( j=i+1; j<=size2; j++) {
       temporary = 0.0; 
       for (k=1; k<=size1; k++) {
         temporary += q[(k-1)*size2+i-1] * q[(k-1)*size2+j-1]; 
       }
 
       for (k=1; k<=size1; k++) {
         q[(k-1)*size2+j-1]-= temporary * q[(k-1)*size2+i-1]; 
       }
       r[(i-1)*size2+j-1]=temporary*anorm;  
      }

   } /*i loop*/ 

   /*Now Make V Into R Matrix*/
   for (i=0; i<size2; i++) {
    r[i*size2+i]=d[i]; 
   } 
   
   free(d); 

return; 
}

/*********************************************************************************************************************/

void complex_modified_gram_schmidt(MKL_Complex16 *q, MKL_Complex16 *r, int size1, int size2) {

   /*The Modified Gram-Schmidt Routine Used to Stabilize Matrix Products*/
   int i, j, k;
   MKL_Complex16 Zero; Zero.real = Zero.imag = 0.0; 
   MKL_Complex16 One; One.real = 1.0; One.imag = 0.0; 
   MKL_Complex16 temporary, anorm;
   MKL_Complex16 *d;

   d=(MKL_Complex16 *)calloc(size2,sizeof(MKL_Complex16));
   czero_vec(r,size2*size2);
 
   for (i=1; i<=size2; i++) {

     temporary = Zero;
     for (j=1; j<=size1; j++) {
       temporary = Cadd(temporary, Cmul(q[(j-1)*size2+i-1], q[(j-1)*size2+i-1])); 
     }
     d[i-1]=Csqrt(temporary);
     anorm = Cdiv(One, d[i-1]);

     for (j=1; j<=size1; j++) {
       q[(j-1)*size2+i-1] = Cmul(q[(j-1)*size2+i-1], anorm);
     }

     for ( j=i+1; j<=size2; j++) {
       temporary = Zero;
       for (k=1; k<=size1; k++) {
         temporary = Cadd(temporary, Cmul(q[(k-1)*size2+i-1],q[(k-1)*size2+j-1]));
       }

       for (k=1; k<=size1; k++) {
         q[(k-1)*size2+j-1] = Csub(q[(k-1)*size2+j-1], Cmul(temporary, q[(k-1)*size2+i-1])); 
       }
       r[(i-1)*size2+j-1]=Cmul(temporary, anorm); 
      }

   } /*i loop*/

   /*Now Make V Into R Matrix*/
   for (i=0; i<size2; i++) {
    r[i*size2+i]=d[i];
   }

free(d);
return;
}

/***********************************************************************************************************/

void complex_matrix_inverse(MKL_Complex16 *input_matrix, MKL_Complex16 *output_matrix, int size) {

   /*Inverse Complex Matrix Using Only Real Inversion*/
   int size_squared = size * size;   
   double *real_matrix, *real_matrix_2, *complex_matrix;
   double *real_matrix_inverse;
   double *r0, *y11, *y11_inverse, *y10_positive, *y10_negative, *cr0_prod;

   real_matrix = (double*)calloc(size_squared,sizeof(double));
   real_matrix_2 = (double*)calloc(size_squared,sizeof(double));
   complex_matrix = (double*)calloc(size_squared,sizeof(double));
   real_matrix_inverse = (double*)calloc(size_squared,sizeof(double));
   r0 = (double*)calloc(size_squared,sizeof(double));
   y11 = (double*)calloc(size_squared,sizeof(double));
   y11_inverse  = (double*)calloc(size_squared,sizeof(double));
   y10_positive = (double*)calloc(size_squared,sizeof(double));
   y10_negative = (double*)calloc(size_squared,sizeof(double));
   cr0_prod = (double*)calloc(size_squared,sizeof(double));

   /*Store Real and Complex Parts of Original Matrix*/
   store_matrix_parts(input_matrix, real_matrix, complex_matrix,size,size); 
   double_copy_mat(real_matrix_2,real_matrix,size_squared);

   /*Get Inverse of Real Matrix*/
   inverse(real_matrix,real_matrix_inverse,size);

   /*Get r0 By Doing Mat Mat*/
   d_mat_mat(real_matrix_inverse,complex_matrix,r0,size,size,size);

   /*Get c r0 Prod*/
   d_mat_mat(complex_matrix,r0,cr0_prod,size,size,size);

   /*Add Two Matrices*/
   d_mat_add(cr0_prod,real_matrix_2,y11_inverse,size,size);

   /*Now Take Inverse to Get y11*/
   inverse(y11_inverse,y11,size);

   /*Now Take Product To Get y10*/
   d_mat_mat(r0,y11,y10_positive,size,size,size);

   /*Now Get Negative and Restore*/
   d_take_negative(y10_positive,y10_negative,size,size);

   /*Restore Matrix*/
   restore_matrix_parts(output_matrix,y11,y10_negative,size,size);

free(real_matrix);
free(real_matrix_2);
free(complex_matrix);
free(real_matrix_inverse);
free(r0);
free(y11);
free(y11_inverse);
free(y10_positive);
free(y10_negative);
free(cr0_prod);
return;
}

/****************************************************************************************************************/

void store_matrix_parts(MKL_Complex16 *input_matrix, double *real_matrix,double *complex_matrix,int size1,int size2){

    /*Store a Complex Matrix's Parts Into Real and Imaginary Matrices*/
    int i, j;
    int ip;

    for (i=0; i<size1; i++) {
      ip = i * size2;
      for (j=0; j<size2; j++) {
        real_matrix[ip+j] = input_matrix[ip+j].real;
        complex_matrix[ip+j] = input_matrix[ip+j].imag;
      }
    }

return;
}

/*************************************************************************************************************/

void double_copy_mat(double *mat_1, double *mat_2, int size) {

   int i;

   for (i=0; i<size; i++) {
      mat_1[i] = mat_2[i];
   }

return;
}

/*************************************************************************************************************/

void d_mat_mat(double *mat_1, double *mat_2, double *product, int size_1, int size_2, int size_3) {

    /*Multiplies Two Matrices Together, Where One is of Size size_1xsize_2 and the other Is of Size size_2xsize_3*/

    int i, j, k;
    dzero_vec(product,size_1*size_3);

    for (i=0; i<size_1; i++) {
      for (j=0; j<size_3; j++) {
        product[i*size_3+j]=0.0;

        for (k=0; k<size_2; k++) {
          product[i*size_3+j] += mat_1[i*size_2+k] * mat_2[k*size_3+j];
        }

      }
    }

return;
}

/************************************************************************************************************/ 

void d_mat_add(double *matrix_1, double *matrix_2, double *sum_matrix,int size1,int size2) {

     /*Take the Negative*/
    int i, j;
    int ip;

    /*Zero Sum Matrix*/
    dzero_vec(sum_matrix,size1*size2);

    for (i=0; i<size1; i++) {
      ip = i * size2;
      for (j=0; j<size2; j++) {
         sum_matrix[ip+j] = matrix_1[ip+j] + matrix_2[ip+j];
      }
    }

return;
}

/***********************************************************************************************************/

void d_take_negative(double *positive, double *negative,int size1,int size2) {

    /*Take the Negative*/
    int i, j;
    int ip;

    for (i=0; i<size1; i++) {
      ip = i * size2;
      for (j=0; j<size2; j++) {
        negative[ip+j] = -1 * positive[ip+j];
      }
    }

return;
}

/************************************************************************************************************/

void restore_matrix_parts(MKL_Complex16 *output_matrix, double *real_matrix, double *complex_matrix,int size1,int size2) {

    /*Stores Real and Complex Matrices Into A Full Complex Matrix*/
    int i, j;
    int ip;

    for (i=0; i<size1; i++) {
      ip = i * size2;
      for (j=0; j<size2; j++) {
        output_matrix[ip+j].real = real_matrix[ip+j];
        output_matrix[ip+j].imag = complex_matrix[ip+j];
      }
    }

return;
}

/***************************************************************************************************************/

MKL_Complex16 conjugate(MKL_Complex16 complex_number) {

    MKL_Complex16 conjugate;

    conjugate.real = complex_number.real;
    conjugate.imag = -1 * complex_number.imag;


return(conjugate);
}

/***************************************************************************************************************/

int factorial(int number) {

   int factorial = 1; 

   if ( number > 0 ) {
     while (number > 0 ) {
      factorial *= number; 
      number -= 1; 
     }
   }
   else if ( factorial == 0 ) { 
   }

return(factorial); 
}

/******************************************************************************************************************/

void complex_eigvecs_2(MKL_Complex16 *matrix, MKL_Complex16 *eigs, MKL_Complex16 *eigvecs, int size) {

   /*A Complex Eig Vec Routine Taken from Online*/

   int p, q, j; 
   int sweep; 
   double red, off, thresh; 
   MKL_Complex16 One; One.real=1.0; One.imag = 0.0; 
   MKL_Complex16 delta, t, s, invc, sx, sy, tx, ty; 
   MKL_Complex16 x, y; 
   MKL_Complex16 **ev; 

   ev = (MKL_Complex16 **)calloc(2, sizeof(MKL_Complex16*)); 
   for (j=0; j<2; j++) {
    ev[j]=(MKL_Complex16 *)calloc(size, sizeof(MKL_Complex16)); 
   } 

   for (p=0; p< size; p++) {
    ev[0][p].real = ev[0][p].imag = 0.0; 
    ev[1][p] = matrix[p*size+p]; 
    eigs[p] = ev[1][p]; 
   }

   for (p=0; p<size; p++) {
    for (q=0; q<size; q++) {
      eigvecs[p*size+q].real=eigvecs[p*size+q].imag=0; 
    }
    eigvecs[p*size+p].real = 1.0; 
   }  

   red = 0.01/((double)pow(size,4)); 

   for (sweep =0; sweep < 200; sweep++ ) {
 
      off = 0.0; 
      for (q=1; q<size; q++) {
        for (p=0; p<q; p++) {
          off += matrix[p*size+q].real * matrix[p*size+q].real + matrix[p*size+q].imag * matrix[p*size+q].imag + matrix[q*size+p].real * matrix[q*size+p].real + matrix[q*size+p].imag * matrix[q*size+p].imag;  
        }
      }
 
    if ( off < EPS ) {
      return; 
    } 

    thresh = 0.0; 
    if ( sweep < 4 ) { thresh = off*red; }

    for (q=1; q<size; q++) {
     for (p=0; p<q; p++) {
      off = matrix[p*size+q].real * matrix[p*size+q].real + matrix[p*size+q].imag * matrix[p*size+q].imag + matrix[q*size+p].real * matrix[q*size+p].real + matrix[q*size+p].imag * matrix[q*size+p].imag; 
      if ( sweep > 4 && off < EPS*(ev[1][p].real * ev[1][p].real + ev[1][p].imag * ev[1][p].imag + ev[1][q].real * ev[1][q].real + ev[1][q].imag * ev[1][q].imag )) {
          matrix[p*size+q].real = matrix[p*size+q].imag = 0.0; 
          matrix[q*size+p].real = matrix[q*size+p].imag = 0.0; 
      }
     else if ( off > thresh ) {

        delta = Cmul(matrix[p*size+q],matrix[q*size+p]);

        x = RCmul(.5, Csub(ev[1][p],ev[1][q])); 
        y= Csqrt(Cadd(Cmul(x,x),delta)); 
        t=Csub(x,y); 
        s=Cadd(x,y);
 
        if ( (t.real*t.real+t.imag*t.imag) < (s.real*s.real+s.imag*s.imag) ) { t=s; }
        t=Cdiv(One, t); 
        delta = Cmul(delta, t); 
        ev[0][p]= Cadd(ev[0][p], delta); 
        ev[1][p] = Cadd(eigs[p], ev[0][p]); 
        ev[0][q] = Csub(ev[0][q], delta); 
        ev[1][q] = Cadd(eigs[q], ev[0][q]);

        invc = Csqrt(Cadd(Cmul(delta, t), One));
        s=Cdiv(t,invc);
        t=Cdiv(t, Cadd(invc, One)); 
        sx=Cmul(s, matrix[p*size+q]); 
        ty = Cmul(t, matrix[p*size+q]); 
        sy = Cmul(s, matrix[q*size+p]); 
        tx = Cmul(t, matrix[q*size+p]); 

        for (j=0; j<size; j++) {
          x=matrix[j*size+p]; 
          y=matrix[j*size+q]; 
          matrix[j*size+p] = Cadd(x, Cmul(sy, Csub(y, Cmul(ty,x)))); 
          matrix[j*size+q] = Csub(y, Cmul(sx, Cadd(x, Cmul(tx,y)))); 
          x=matrix[p*size+j]; 
          y=matrix[q*size+j]; 
          matrix[p*size+j] = Cadd(x, Cmul(sx,Csub(y, Cmul(tx,x)))); 
          matrix[q*size+j] = Csub(y, Cmul(sy, Cadd(x, Cmul(ty,y)))); 
        }
   
        matrix[p*size+q].real = matrix[p*size+q].imag = 0.0; 
        matrix[q*size+p].real = matrix[q*size+p].imag = 0.0; 

        for (j=0; j<size; j++) {
          x=eigvecs[p*size+j]; 
          y=eigvecs[q*size+j]; 
          eigvecs[p*size+j] = Cadd(x, Cmul(sx, Csub(y, Cmul(tx,x)))); 
          eigvecs[q*size+j] = Csub(y, Cmul(sy, Cadd(x, Cmul(ty,y)))); 
        }
     } /*elseif*/
 
   }
  }

  for (p=0; p<size; p++) {
     ev[0][p].real = ev[0][p].imag = 0.0; 
     eigs[p] = ev[1][p]; 
  }      

 } /*sweeps*/ 
         
free(ev); 
return; 
}

/******************************************************************************************************************/

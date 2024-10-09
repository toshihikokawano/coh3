/***********************************************************/
/*      F.D.Becchetti, G.W.Greenlees                       */
/*      Phys. Rev. 182, 1190 (1969)                        */
/*      Ann. Rep. J.H.Williams Lab. p.116 (1969)           */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPBecchettiGreenlees(int at, int zt, int ai, int zi, Optical *omp)
{
   double eps;

   eps = 1.0-2.0*(double)zt/(double)at;

   if(ai == 1){ /* neutron or proton */
     if(zi==0){
       omp->v1  = 56.3-24.0*eps;
       omp->wv1 = -1.56        ;
       omp->ws1 = 13.0-12.0*eps;
       omp->rs  =  1.26 ;
       omp->rv  =  1.26 ;
       omp->rc  =  0.0  ;
       omp->as  =  0.58 ;
       omp->av  =  0.58 ;
     }
     else{
       omp->v1  = 54.0+24.0*eps +0.4*(double)zt/pow((double)at,1.0/3.0);
       omp->wv1 = -2.7         ;
       omp->ws1 = 11.8+12.0*eps;
       omp->rs  =  1.32 ;
       omp->rv  =  1.32 ;
       omp->rc  =  1.25 ;
       omp->as  =  0.51+0.7*eps;
       omp->av  =  0.51+0.7*eps;
     }
     omp->v2  = -0.32 ;
     omp->wv2 =  0.22 ;
     omp->ws2 = -0.25 ;
     omp->vso1=  6.2  ;
     omp->r0  =  1.17 ;
     omp->rvso=  1.01 ;
     omp->a0  =  0.75 ;
     omp->avso=  0.75 ;

     return(0x0013);
   }
   else if(ai == 3 && zi == 1){ /* triton*/
     omp->v1  = 151.9+50.0*eps;
     omp->v2  = -0.17 ;
     omp->wv1 =  41.7+44.0*eps;
     omp->wv2 = -0.33 ;
     omp->vso1=  2.5  ;

     omp->r0  =  1.20 ;
     omp->rv  =  1.40 ;
     omp->rvso=  1.20 ;
     omp->rc  =  1.3  ;

     omp->a0  =  0.72 ;
     omp->av  =  0.88 ;
     omp->avso=  0.72 ;

     return(0x0011);
   }
   else if(ai == 3 && zi == 2){ /* helion*/
     omp->v1  = 165.9-  6.4*eps;
     omp->v2  = -0.17 ;
     omp->wv1 =  46.0-110.0*eps;
     omp->wv2 = -0.33 ;
     omp->vso1=  2.5  ;

     omp->r0  =  1.20 ;
     omp->rv  =  1.40 ;
     omp->rvso=  1.20 ;
     omp->rc  =  1.3  ;

     omp->a0  =  0.72 ;
     omp->av  =  0.84 ;
     omp->avso=  0.72 ;

     return(0x0011);
   }
   return(0);
}



/***********************************************************/
/*      F.G.Perey                                          */
/*      Phys. Rev. 131, 745(1963)                          */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPPerey(int at, int zt, int ai, int zi, Optical *omp)
{
   double eps,a13;

   if(ai != 1 || zi != 1) return(0);

   eps = 1.0-2.0*(double)zt/(double)at;
   a13 = pow((double)at,1.0/3.0);

   omp->v1  = 53.3+27.0*eps + 0.4*(double)zt/a13;
   omp->v2  = -0.55 ;
   omp->v3  =  0.0  ;

   omp->wv1 =  0.0  ;
   omp->wv2 =  0.0  ;
   omp->wv3 =  0.0  ;

/* omp->ws1 =  3.0*a13 ; */
   omp->ws1 = 13.5  ;
   omp->ws2 =  0.0  ;
   omp->ws3 =  0.0  ;

   omp->vso1=  7.5  ;
   omp->vso2=  0.0  ;
   omp->vso3=  0.0  ;

   omp->r0  =  1.25 ;
   omp->rs  =  1.25 ;
   omp->rv  =  0.0  ;
   omp->rvso=  1.25 ;
   omp->rc  =  1.25 ;

   omp->a0  =  0.65 ;
   omp->as  =  0.47 ;
   omp->av  =  0.0  ;
   omp->avso=  0.47 ;

   return(0x0012);
}



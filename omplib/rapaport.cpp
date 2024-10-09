/***********************************************************/
/*      J.Rapaport, V.Kulkarni, R.W.Finlay                 */
/*      Nucl. Phys A330, 15(1979)                          */
/***********************************************************/

#include "omplib.h"

unsigned int OMPRapaport(double e,int at, int zt, int ai, int zi, Optical *omp)
{
   double eps;

   if(ai != 1 || zi != 0) return(0);

   eps = 1.0-2.0*(double)zt/(double)at;

   omp->v1  = 54.19-22.7*eps;
   omp->v2  = -0.33+0.19*eps;
   omp->v3  =  0.0;
   if(e <= 15.0){
       omp->wv1 =  0.0;           omp->wv2 =  0.0 ;
       omp->ws1 =  4.28-12.8*eps; omp->ws2 =  0.4 ;
   }
   else{
       omp->wv1 = -4.3;           omp->wv2 =  0.38;
       omp->ws1 = 14.0 -10.4*eps; omp->ws2 = -0.39;
   }

   omp->vso1=  6.2  ;

   omp->r0  =  1.198;
   omp->rs  =  1.295;
   omp->rv  =  1.295;
   omp->rvso=  1.01 ;

   omp->a0  =  0.663;
   omp->as  =  0.59 ;
   omp->av  =  0.59 ;
   omp->avso=  0.75 ;

   return(0x0013);
}



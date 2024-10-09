/***********************************************************/
/*      R.L.Varner, et al.                                 */
/*      Phys. Rep. 201, 57(1991)                           */
/*      CH89: Chapel Hill 1989                             */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPChapelHill89(double e,int at, int zt, int ai, int zi, Optical *omp)
{
   double eps,ec,a13;

   if(ai != 1) return(0);

   eps = (1.0-2.0*(double)zt/(double)at) * ((zi==0) ? -1.0 : 1.0);
   a13 = pow((double)at,1.0/3.0);

   omp->v1  = 52.9  +13.1 *eps;
   omp->v2  = -0.299;
   omp->v3  =  0.0  ;
   ec       =  0.0  ;
   if(zi == 1){
     omp->rc = 1.238+0.116/a13;
     ec=1.73*zt/(omp->rc*a13);
     omp->v1+= omp->v2*(-ec);
   }

   omp->wv1 =  7.8          /(1.0+exp((35.0-e+ec)/16.0));
   omp->wv2 =  0.0  ;
   omp->wv3 =  0.0  ;
   omp->ws1 =(10.0+18.0*eps)/(1.0+exp((e-ec-36.0)/37.0));
   omp->ws2 =  0.0  ;
   omp->ws3 =  0.0  ;

   omp->vso1=  5.9  ;

   omp->r0  =  1.250-0.225/a13;
   omp->rs  =  1.33 -0.42 /a13;
   omp->rv  =  1.33 -0.42 /a13;
   omp->rvso=  1.34 -1.2  /a13;

   omp->a0  =  0.690;
   omp->as  =  0.69 ;
   omp->av  =  0.69 ;
   omp->avso=  0.63 ;

   return(0x0013);
}



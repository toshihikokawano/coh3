/***********************************************************/
/*      J. Bojowald et al. for deuteron                    */
/*      Phys. Rev. C 38, 1153 (1988), Eq.(6)               */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPBojowald(double e, int at, int zt, int ai, int zi, Optical *omp)
{
   double a13;

   if(ai != 2 || zi != 1) return(0);

   a13 = pow((double)at,1.0/3.0);

   omp->v1  = 81.32 + 1.43 * zt/a13;  omp->v2  = -0.240;  omp->v3  =  0.0;
   omp->wv1 = 0.132*(e-45);
   if(omp->wv1<0.0) omp->wv1 = 0.0;

   omp->ws1 =  7.80 + 1.04 * a13 - 0.712*omp->wv1;
   if(omp->ws1 < 0.0) omp->ws1 = 0.0;

   omp->wv2 =  0.0;    omp->wv3 =  0.0;
   omp->ws2 =  0.0;    omp->ws3 =  0.0;
   omp->vso1=  6.0;    omp->vso2=  0.0;    omp->vso3=  0.0;

   omp->r0  =  1.180;  omp->a0  =  0.636 + 0.035 *a13;
   omp->rv  =  1.270;  omp->av  =  0.768 + 0.021 *a13;

   omp->rs  =  omp->rv;
   omp->as  =  omp->av;

   omp->rvso = omp->avso = 0.78 + 0.038*a13;

   omp->rc  =  1.3;

   return(0x0013);
}



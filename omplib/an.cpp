/***********************************************************/
/*      Haixia An and Chonghai Cai for deuteron            */
/*      Phys. Rev. C73 054605  (2006)                      */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPAnHaixia(int at, int zt, int ai, int zi, Optical *omp)
{
   double a13;

   if(ai != 2 || zi != 1) return(0);

   a13 = pow((double)at,1.0/3.0);


   omp->v1 = 91.85 + 0.642 * zt/a13; omp->v2 = -0.249; omp->v3 = 0.000116;

   omp->wv1 = 1.104; omp->wv2 = 0.0622; omp->wv3 = 0.0;

   omp->ws1 = 10.83; omp->ws2 = -0.0306;

   omp->vso1 = 7.114; omp->vso2 = 0.0; omp->vso3 = 0.0;

   omp->r0  =  1.152 - 0.00776/a13;
   omp->a0  =  0.719 + 0.0126 *a13;

   omp->rv  =  1.305 + 0.0997/a13;
   omp->av  =  0.855 - 0.1*a13;

   omp->rs  =  1.334 + 0.152/a13;
   omp->as  =  0.531 + 0.062*a13;

   omp->rvso = 0.972;
   omp->avso = 1.011;

   omp->rc  =  1.303;

   return(0x0013);
}


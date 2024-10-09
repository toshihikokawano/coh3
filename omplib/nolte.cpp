/***********************************************************/
/*      M.Nolte, H.Machner, and J.Bojowald                 */
/*      Phys. Rev. C36, 1312 (1987)                        */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPNolte(int at, int zt, int ai, int zi, Optical *omp)
{
   double a13;

   if(ai != 4 || zi != 2) return(0);

   a13 = pow((double)at,1.0/3.0);

   omp->v1  =101.1 + 6.051*zt/a13;    omp->v2  = -0.248;  omp->v3  =  0.0;
   omp->wv1 = 26.82- 1.706*a13   ;    omp->wv2 =  0.006;  omp->wv3 =  0.0;

   omp->ws1 =  0.0;    omp->ws2 =  0.0;    omp->ws3 =  0.0;
   omp->vso1=  0.0;    omp->vso2=  0.0;    omp->vso3=  0.0;

   omp->r0  =  1.245;  omp->a0  =  0.817 - 0.0085*a13;
   omp->rv  =  1.570;  omp->av  =  0.692 - 0.020 *a13;
   omp->rs  =  0.0;    omp->as  =  0.0;

   omp->rc  =  1.3;
/*
    This Coulomb radius was taken from L.W.Put, Nucl.Phys. A291, 93 (1977)
*/
   return(0x0011);
}



/***********************************************************/
/*  Yinlu Han, Yuyang Shi and Qingbiao Shen for deuteron   */
/*      Phys. Rev. C74 044615  (2006)                      */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPHanYinlu(int at, int zt, int ai, int zi, Optical *omp)
{
   double a13;
   double sym;

   if(ai != 2 || zi != 1) return(0);

   a13 = pow((double)at,1.0/3.0);
   sym = (double)(at-2*zt)/(double)at;


   omp->v1 = 82.18 -34.811*sym + 1.058*zt/a13; omp->v2 = -0.148; omp->v3 = -0.000886;

   omp->wv1 = -4.916 + 35.0*sym; omp->wv2 = 0.0555; omp->wv3 = 0.0000442;

   omp->ws1 = 20.968 - 43.398*sym; omp->ws2 = -0.0794;

   omp->vso1 = 7.406; omp->vso2 = 0.0; omp->vso3 = 0.0;

   omp->wso1 = -0.206; omp->wso2 = 0.0; omp->wso3 = 0.0;

   omp->r0  =  1.174;
   omp->a0  =  0.809;

   omp->rv  =  1.563;
   omp->av  =  0.700 + 0.045*a13;

   omp->rs  =  1.328;
   omp->as  =  0.465 + 0.045*a13;

   omp->rvso = 1.234;
   omp->avso = 0.813;

   omp->rc  =  1.698;

   return(0x0013);
}


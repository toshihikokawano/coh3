/***********************************************************/
/*     D.G. Madland and Schwandt                           */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPMadlandSchwandt(double e, int at, int zt, int ai, int zi, Optical *omp)
{
   double eps,a13,e1,e2;

   if(ai != 1 || zi != 1) return(0);

   eps = 1.0-2.0*(double)zt/(double)at;
   a13 = pow((double)at,1.0/3.0);
   e1 = e-80;
   e2 = log(e);

   omp->v1  =105.5*(1-0.1625*e2) + 16.5*eps;
   omp->v2  =  0.0  ;
   omp->v3  =  0.0  ;

   omp->wv1 =  6.6 +0.0273*e1 + 3.87E-06 * e1*e1*e1 ;
   omp->wv2 =  0.0  ;
   omp->wv3 =  0.0  ;

   omp->ws1 =  0.0  ;
   omp->ws2 =  0.0  ;
   omp->ws3 =  0.0  ;

   omp->vso1= 19.0*(1-0.166 *e2) - 3.75*eps;
   omp->vso2=  0.0  ;
   omp->vso3=  0.0  ;

   omp->wso1=  7.5*(1-0.248 *e2);
   omp->wso2=  0.0  ;
   omp->wso3=  0.0  ;

   omp->r0  =  1.125 + 0.001 *e;
   omp->rs  =  0.0  ;
   omp->rv  =  1.65  - 0.0024*e;
   omp->rvso=  0.92  + 0.0305*a13;
   omp->rwso=  0.877 + 0.0360*a13;
   omp->rc  =  1.25 ; /* ? */

   omp->a0  =  0.675 + 0.00031*e;
   omp->as  =  0.0  ;
   omp->av  =  0.328 + 0.00244*e;
   omp->avso=  0.768 - 0.0012 *e;
   omp->awso=  0.620;

   return(0x0011);
}


/*************************************************/
/*    Madland - Young                            */
/*        for Actinides                          */
/*************************************************/
unsigned int OMPMadlandYoung(double e, int at, int zt, Optical *omp)
{
  double eps = 1.0-2.0*(double)zt/(double)at;

  omp->v1  = 50.378 - 27.073 *eps;
  omp->v2  = -0.354;
  omp->ws1 =  9.265 - 12.666 *eps;
  omp->ws2 = -0.232;
  omp->ws3 =  0.03318;
  omp->vso1=  6.2;

  omp->r0  =  1.264;
  omp->rs  =  1.256;
  omp->rvso=  1.100;

  omp->a0  =  0.612;
  omp->as  =  0.553 + 0.0144*e;
  omp->avso=  0.750;

  return(0x0012);
}



/***********************************************************/
/*      D.Wilmore, P.E.Hodgson                             */
/*      Nucl. Phys. 55, 673 (1964)                         */
/***********************************************************/

#include "omplib.h"

unsigned int OMPWilmoreHodgson(int at, int ai, int zi, Optical *omp)
{
   double a;

   if(ai != 1 || zi != 0) return(0);

   a=(double)at;
   omp->v1  = 47.01 ;  omp->v2  = -0.267;  omp->v3  = -0.0018 ;
   omp->ws1 =  9.52 ;  omp->ws2 = -0.053;  omp->ws3 =  0.0    ;
   omp->r0  =  1.322-7.6e-04*a+4.0e-06*a*a-8.0e-09*a*a*a;
   omp->rs  =  1.266-3.7e-04*a+2.0e-06*a*a-4.0e-09*a*a*a;
   omp->a0  =  0.66 ;
   omp->as  =  0.48 ;
/****************  NOTE  ******************/
/* Original             Perey*            */
/* omp->v3  = -0.00118 omp->v3  = -0.0018 */
/* omp->ws2 = -0.53    omp->ws2 = -0.053  */
/*                                        */
/* *)Atomic Data and Nuclear Data Tables  */
/*   13, 293(1974)                        */
/******************************************/

   return(0x0012);
}



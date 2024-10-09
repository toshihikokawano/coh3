/***********************************************************/
/*      O.F.Lemos                                          */
/*      Orsay Report, Series A, No.136 (1972)              */
/***********************************************************/

#include "omplib.h"

unsigned int OMPLemos(int ai, int zi, Optical *omp)
{
   if(ai != 4 || zi != 2) return(0);

   omp->v1  =193.0;    omp->v2  = -0.15;   omp->v3  =  0.0;
   omp->wv1 = 21.0;    omp->wv2 =  0.25;   omp->wv3 =  0.0;
   omp->ws1 =  0.0;    omp->ws2 =  0.0;    omp->ws3 =  0.0;
   omp->vso1=  0.0;    omp->vso2=  0.0;    omp->vso3=  0.0;
   omp->r0  =  1.37;   omp->rv  =  1.37;   omp->rs  =  0.0;
   omp->rc  =  1.4;
   omp->a0  =  0.56;   omp->av  =  0.56;   omp->as  =  0.0;

   return(0x0011);
}



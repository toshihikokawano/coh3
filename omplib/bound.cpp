/***********************************************************/
/*      Shell Model Potential for Bound Nucleon            */
/*      Ross-Mark-Lawson Potential                         */
/***********************************************************/

#include "omplib.h"

unsigned int OMPBoundState(int zi, Optical *omp)
{
   omp->v1  = 42.8;
   omp->vso1= 39.5/4.0;

   omp->r0  =  1.300;
   omp->rvso=  1.300;

   omp->a0  =  0.690;
   omp->avso=  0.690;
   omp->rc  = (zi == 1) ? 1.25: 0.0;

   return(0x0010);
}



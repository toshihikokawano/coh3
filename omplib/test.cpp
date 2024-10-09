/***********************************************************/
/*    Test Parameters                                      */
/***********************************************************/

#include "omplib.h"

unsigned int OMPtest(int zi, Optical *omp)
{
  omp->v1  = 50.0;
  omp->wv1 =  0.0;
  omp->ws1 =  1.0;

  omp->vso1=  5.0;


  omp->r0  =  1.2;
  omp->rv  =  1.2;
  omp->rs  =  1.2;
  omp->rvso=  1.2;

  omp->a0  =  0.6;
  omp->av  =  0.6;
  omp->as  =  0.6;
  omp->avso=  0.4;

  omp->rc  = (zi == 1) ? 1.25: 0.0;

  return(0x0013);
}

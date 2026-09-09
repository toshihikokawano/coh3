/*****************************************************/
/*         User defined OMPs                         */
/*****************************************************/

#include "omplib.h"

unsigned int OMPscratch(double e,int at,int zt,int ai,int zi,Optical *omp)
{
  unsigned int pf = OMPKoningDelaroche(e,at,zt,ai,zi,omp);

  omp->wv1 = 0.0;

  return(pf);
}

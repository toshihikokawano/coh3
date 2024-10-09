/*****************************************************/
/*         User defined OMPs for unknown particles   */
/*****************************************************/

#include "omplib.h"

unsigned int OMPuserdef(double e,int at,int zt,int ai,int zi,Optical *omp)
{
  unsigned int potform = 0;

  // tetraneutron Z = 0, A = 4
  if(zi == 0 && ai == 4){
    Optical ompn;
    OMPKoningDelaroche(e/4.0,at,zt,1,0,&ompn);

    *omp = ompn;

    omp->v1   = 4 * ompn.v1;
    omp->ws1  = 4 * ompn.ws1;
    omp->wv1  = 4 * ompn.wv1;

    omp->vso1 = 0.0;
    omp->wso1 = 0.0;
    potform = 0x0013;
  }
  // hexaneutron Z = 0, A = 6
  else if(zi == 0 && ai == 6){
    Optical ompn;
    OMPKoningDelaroche(e/6.0,at,zt,1,0,&ompn);

    *omp = ompn;

    omp->v1   = 6 * ompn.v1;
    omp->ws1  = 6 * ompn.ws1;
    omp->wv1  = 6 * ompn.wv1;

    omp->vso1 = 0.0;
    omp->wso1 = 0.0;
    potform = 0x0013;
  }

  return(potform);
}

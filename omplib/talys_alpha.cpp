/***********************************************************/
/*      A. J. Koning                                       */
/*      TALYS folding method                               */
/***********************************************************/

#include "omplib.h"

unsigned int OMPTALYS_alpha(double e, int at, int zt, int ai, int zi, Optical *omp)
{
   Optical nopt,popt;

   if(ai != 4 || zi != 2) return(0);

   OMPKoningDelaroche(0.25*e, at,  zt, 1, 0, &nopt);
   OMPKoningDelaroche(0.25*e, at,  zt, 1, 1, &popt);

   omp->v1   = 2*nopt.v1  + 2*popt.v1;    omp->v2  =  0.0  ;  omp->v3  =  0.0;
   omp->wv1  = 2*nopt.wv1 + 2*popt.wv1;   omp->wv2 =  0.0  ;  omp->wv3 =  0.0;
   omp->ws1  = 2*nopt.ws1 + 2*popt.ws1;   omp->ws2 =  0.0  ;  omp->ws3 =  0.0;

   omp->wso1=  0.0;    omp->wso2=  0.0;    omp->wso3=  0.0;
   omp->vso1=  0.0;    omp->vso2=  0.0;    omp->vso3=  0.0;


   omp->r0  =  (nopt.r0 + popt.r0)/2.0;
   omp->a0  =  (nopt.a0 + popt.a0)/2.0;

   omp->rv  = omp->rs  =  omp->r0;
   omp->av  = omp->as  =  omp->a0;

   omp->rc  =  popt.rc;

   return(0x0013);
}


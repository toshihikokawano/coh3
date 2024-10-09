/***********************************************************/
/*      A.Koning, and J.-P.Delaroche                       */
/*      Nucl.Phys. A (2002)                                */
/*      KD2001                                             */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPKoningDelaroche(double e,int at, int zt, int ai, int zi, Optical *omp)
{
   double eps,ef,a13,ex,v1,v2,v3,d1,d2,w1,w2,vc;

   if(ai != 1) return(0);

   a13 = pow((double)at,1.0/3.0);
   eps = 1.0-2.0*(double)zt/(double)at;

/* common geometry part */
   omp->r0   = 1.3039 - 0.4054    /a13;
   omp->rs   = 1.3424 - 0.01585   *a13;
   omp->rv   = omp->r0;
   omp->rvso = 1.1854 - 0.647     /a13;
   omp->rwso = omp->rvso;
   omp->rc   = (zi==0) ? 0.0:
                         1.198 
                      +  0.697*pow((double)at, -2.0/3.0) 
                      + 12.994*pow((double)at, -5.0/3.0);

   omp->a0   = 0.6778 - 1.487e-04 *at;
   omp->av   = omp->a0;
   omp->avso = 0.59;
   omp->awso = omp->avso;

 /* for neutron */
   if(zi == 0){
     ef = -11.2814             + 0.02646 *at;
     v1 = 59.30     - 21.0*eps - 0.024   *at;
     v2 = 7.228e-03            - 1.48e-06*at;
     v3 = 1.994e-05            - 2.00e-08*at;
     w1 = 12.195               + 0.0167  *at;
     d1 = 16.0      - 16.0*eps;
     vc = 0.0;

     omp->as = 0.5446 - 1.656e-04 *at;
   }
 /* for proton */
   else{
     ef = -8.4075              + 0.01378 *at;
     v1 = 59.30     + 21.0*eps - 0.024   *at;
     v2 = 7.067e-03            + 4.36e-06*at;
     v3 = 1.747e-05            + 1.50e-08*at;
     w1 = 14.336               + 0.0189  *at;
     d1 = 14.3      + 20.0*eps;
     vc = 0.42*(double)zt/a13;

     omp->as = 0.5413 + 3.963e-04 *at;
   }

   w2 = 73.55   + 0.0795*at;
   d2 =  0.0180 + 3.802e-03/(1.0+exp((at-156.0)/8.0));

   ex = e-ef;

   omp->v1  =  v1*(1.0 - ex*(v2 - ex*(v3 - ex*7.0e-09))) + vc;
   omp->v2  =  0.0;
   omp->v3  =  0.0;

   omp->ws1 =  d1*exp(-d2*ex)*ex*ex/(ex*ex+132.25);
   omp->ws2 =  0.0;
   omp->ws3 =  0.0;

   omp->wv1 =  w1            *ex*ex/(ex*ex+w2*w2);
   omp->wv2 =  0.0;
   omp->wv3 =  0.0;

   omp->vso1= (5.922 + 0.0030*at)*exp(-0.0040*ex);
   omp->vso2=  0.0;

   omp->wso1= -3.1           *ex*ex/(ex*ex+2.56e+04);
   omp->wso2=  0.0;

   return(0x0013);
}



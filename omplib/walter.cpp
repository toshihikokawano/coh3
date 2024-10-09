/***********************************************************/
/*      R.L.Walter, P.P.Guss                               */
/*      Proc. Int. Conf. Santa Fe 1079 (1985)              */
/*      Rad. Effect, 95, 73 (1986)                         */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPWalterGuss(double e,int at, int zt, int ai, int zi, Optical *omp)
{
   double eps,a13;

   if(ai != 1) return(0);

   eps = 1.0-2.0*(double)zt/(double)at;
   a13 = pow((double)at,1.0/3.0);

   if(e <= 40.0){
     if(zi==0){  omp->v1  = 52.56 -16.5  *eps;
                 omp->v2  = -0.310+ 0.081*eps;
     }
     else{       omp->v1  = 52.56 +16.5  *eps +0.4*(double)zt/a13;
                 omp->v2  = -0.310- 0.081*eps;
     }
   }
   else{
     if(zi == 0){omp->v1  = 52.56 - 0.310*40.0*(1.0+log(e/40.0))
                          -(16.50 - 0.081*40.0*(1.0+log(e/40.0)))*eps;
     }
     else{       omp->v1  = 52.56 - 0.310*40.0*(1.0+log(e/40.0))
                          +(16.50 - 0.081*40.0*(1.0+log(e/40.0)))*eps
                          +  0.4*(double)zt/a13 * 678.0/(e*e);
     }
     omp->v2  =  0.0;
   }

   if(e <= 39.4){omp->wv1 = -0.963;
                 omp->wv2 =  0.153;
   }
   else{         omp->wv1 = -0.963+ 0.153*e*(1.0-0.33*log(e/39.4));
                 omp->wv2 =  0.0;
   }

   if(zi == 0){  omp->ws1 = 10.85 -14.94 *eps;
   }
   else{         omp->ws1 = 10.85 +14.94 *eps-1.30;
   }
   omp->ws2 = -0.157;

   omp->vso1=  5.767+ 2.0*eps;
   omp->vso2= -0.015;
   omp->wso1=  0.791;
   omp->wso2= -0.018;

   omp->r0  =  1.219;
   omp->rs  =  1.282; omp->rv  = 1.38 +3.76/(double)at;
   omp->rvso=  1.103; omp->rwso= 1.364;
   omp->rc  = (zi == 0) ? 0.0 : 1.219;

   omp->a0  =  0.688;
   omp->as  =  0.512; omp->av  = 0.557-0.462/sqrt((double)at);
   omp->avso=  0.560; omp->awso= 0.632;
/****************  NOTE  ******************/
/* In original paper                      */
/* omp->av  = 0.557-0.462*sqrt((double)at)*/
/******************************************/

   return(0x0013);
}


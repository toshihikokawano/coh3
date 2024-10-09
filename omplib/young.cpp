/***********************************************************/
/*     P.G. Young Potentials for Actinides                 */
/***********************************************************/

#include <cmath>

#include "omplib.h"

/*************************************************/
/* Am OMP                                        */
/*************************************************/
unsigned int OMPYoung_Am(double e, int at, int zt, Optical *omp)
{
   double eps;

   eps = 1.0-2.0*(double)zt/(double)at;

/*    P.G.Young for Am241 */

   omp->v1  = 52.102 - 27.75*eps;
   omp->v2  = -0.30;
   omp->ws1 = (e<= 8.0) ? 5.3243-9.4995*eps :  8.9243-9.4995*eps + 0.046*8.0;
   omp->ws2 = (e<= 8.0) ? 0.45            : -0.046;
   omp->wv1 = (e<= 8.0) ? 0.0             : -1.6;
   omp->wv2 = (e<= 8.0) ? 0.0             :  0.2;

   omp->vso1=  6.2;

   omp->r0  =  1.25;
   omp->rs  =  1.24;
   omp->rv  =  1.24;
   omp->rvso=  1.01;

   omp->a0  =  0.60;
   omp->as  =  0.55;
   omp->av  =  0.55;
   omp->avso=  0.75;

/* adjusted to S0 */

   omp->ws1 = 2.3;
   omp->ws2 = 0.577;

   return(0x0013);
}


/*************************************************/
/* Pu239 OMP                                     */
/*************************************************/
unsigned int OMPYoung_Pu(double e, int at, int zt, Optical *omp)
{

  if(at != 239 || zt != 94) return 0;

/*    P.G.Young for Pu239 */

   omp->v1  = 46.2;
   omp->v2  = -0.30;
   omp->ws1 = (e<= 8.0) ? 3.3   :  7.268;
   omp->ws2 = (e<= 8.0) ? 0.45  : -0.046;
   omp->wv1 = (e<= 7.0) ? 0.0   : -0.7;
   omp->wv2 = (e<= 7.0) ? 0.0   :  0.1;

   omp->vso1=  6.2;

   omp->r0  =  1.26;
   omp->rs  =  1.24;
   omp->rv  =  1.26;
   omp->rvso=  1.12;

   omp->a0  =  0.63;
   omp->as  =  0.50;
   omp->av  =  0.63;
   omp->avso=  0.47;

/*    Lagrange for Pu240 */

/*
   omp->v1  = 47.5;
   omp->v2  = -0.30;
   omp->ws1 = (e<= 10.0) ? 2.7   :  6.700;
   omp->ws2 = (e<= 10.0) ? 0.40  :  0.000;

   omp->vso1=  7.5;

   omp->r0  =  1.24;
   omp->rs  =  1.26;
   omp->rvso=  1.24;

   omp->a0  =  0.62;
   omp->as  =  0.58;
   omp->avso=  0.62;

*/
   return(0x0013);
}


/*************************************************/
/* Re OMP                                        */
/*************************************************/
unsigned int OMPYoung_Re(double e, int at, int zt, Optical *omp)
{
  if(zt != 75) return(0);

  double eps = 1.0-2*75.0/at;


  omp->v1  = 49.8 - 16.0*eps;
  omp->v2  = -0.30;

  omp->ws1 = ((e<= 9.0) ? 4.020 : 11.220) - 8.0*eps;
  omp->ws2 =  (e<= 9.0) ? 0.75  : -0.05;

  omp->wv1 = -1.8;
  omp->wv2 =  0.2;

  omp->vso1=  7.5;

  omp->r0  =  1.26;
  omp->rs  =  1.26;
  omp->rv  =  1.26;
  omp->rvso=  1.26;

  omp->a0  =  0.61;
  omp->as  =  0.47;
  omp->av  =  0.61;
  omp->avso=  0.61;

  return(0x0013);
}

/*************************************************/
/*    Flap 2.2 by F.Dietrich(LLNL) for Actinides */
/*    deformation and coupling scheme are not    */
/*    provided but use of some standard values   */
/*    are recommended                            */
/*************************************************/

#include "omplib.h"

static inline double flapint(double v1, double v2, double e1, double e2, double x){
  return( (v2-v1)/(e2-e1)*(x-e1) + v1 );
}

unsigned int OMPFlap22(double e, int at, int zt, Optical *omp)
{
  double eps;

  eps = 1.0-2.0*(double)zt/(double)at;

  if(e <= 1.0){
    omp->v1  = flapint(52.000,52.000, 0, 1,e) - flapint(26.000,26.000, 0, 1,e)*eps;
    omp->wv1 = flapint( 0.000, 0.000, 0, 1,e) - flapint( 0.000, 0.000, 0, 1,e)*eps;
    omp->ws1 = flapint( 3.080, 3.480, 0, 1,e) - flapint( 1.540, 1.740, 0, 1,e)*eps;
    omp->r0  = flapint( 1.250, 1.249, 0, 1,e);
  }
  else if(e <= 5.0){
    omp->v1  = flapint(52.000,51.661, 1, 5,e) - flapint(26.000,25.830, 1, 5,e)*eps;
    omp->wv1 = flapint( 0.000, 0.000, 1, 5,e) - flapint( 0.000, 0.000, 1, 5,e)*eps;
    omp->ws1 = flapint( 3.480, 4.737, 1, 5,e) - flapint( 1.740, 2.369, 1, 5,e)*eps;
    omp->r0  = flapint( 1.249, 1.245, 1, 5,e);
  }
  else if(e <= 10.0){
    omp->v1  = flapint(51.661,49.856, 5,10,e) - flapint(25.830,24.928, 5,10,e)*eps;
    omp->wv1 = flapint( 0.000, 0.338, 5,10,e) - flapint( 0.000, 0.169, 5,10,e)*eps;
    omp->ws1 = flapint( 4.737, 6.768, 5,10,e) - flapint( 2.369, 3.384, 5,10,e)*eps;
    omp->r0  = flapint( 1.245, 1.240, 5,10,e);
  }
  else if(e <= 20.0){
    omp->v1  = flapint(49.856,46.810,10,20,e) - flapint(24.928,23.405,10,20,e)*eps;
    omp->wv1 = flapint( 0.338, 2.143,10,20,e) - flapint( 0.169, 1.072,10,20,e)*eps;
    omp->ws1 = flapint( 6.768, 6.768,10,20,e) - flapint( 3.384, 3.384,10,20,e)*eps;
    omp->r0  = flapint( 1.240, 1.230,10,20,e);
  }
  else if(e <= 50.0){
    omp->v1  = flapint(46.810,38.351,20,50,e) - flapint(23.405,19.175,20,50,e)*eps;
    omp->wv1 = flapint( 2.143, 7.557,20,50,e) - flapint( 1.072, 3.779,20,50,e)*eps;
    omp->ws1 = flapint( 6.768, 1.354,20,50,e) - flapint( 3.384, 0.677,20,50,e)*eps;
    omp->r0  = flapint( 1.230, 1.210,20,50,e);
  }

  omp->vso1=  6.2;

  omp->rv  =  1.27;
  omp->rs  =  1.27;
  omp->rvso=  1.15;

  omp->a0  =  0.63;
  omp->av  =  0.62;
  omp->as  =  0.62;
  omp->avso=  0.75;

  return(0x0013);
}



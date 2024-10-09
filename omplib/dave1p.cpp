/***********************************************************/
/*    J.H. Dave and C.R. Gould OMP for 1p shell            */
/*    Phys. Rev. C 28, 2212 (1983)                         */
/***********************************************************/

#include "omplib.h"

unsigned int OMPDave1p(int at, int zt, Optical *omp)
{
  /*** Li6 */
  if( (zt == 3) && (at == 6) ){
    omp->v1  = 39.620;    omp->v2  = -0.027;
    omp->r0  =  1.507;
    omp->a0  =  0.663;

    omp->ws1 =  9.298;    omp->ws2 =  0.646;
    omp->rs  =  1.616;
    omp->as  =  0.196;
  }

  /*** Li7 */
  else if( (zt == 3) && (at == 7) ){
    omp->v1  = 37.740;    omp->v2  = -0.001;
    omp->r0  =  1.504;
    omp->a0  =  0.565;

    omp->ws1 =  0.165;    omp->ws2 =  1.113;
    omp->rs  =  1.512;
    omp->as  =  0.185;
  }

  /*** Be9 */
  else if( (zt == 4) && (at == 9) ){
    omp->v1  = 38.500;    omp->v2  = -0.145;
    omp->r0  =  1.447;
    omp->a0  =  0.387;

    omp->ws1 =  1.666;    omp->ws2 =  0.365;
    omp->rs  =  1.368;
    omp->as  =  0.424;
  }

  /*** B10 */
  else if( (zt == 5) && (at == 10) ){
    omp->v1  = 47.910;    omp->v2  = -0.346;
    omp->r0  =  1.387;
    omp->a0  =  0.464;

    omp->ws1 =  0.675;    omp->ws2 =  0.810;
    omp->rs  =  1.336;
    omp->as  =  0.278;
  }

  /*** B11 */
  else if( (zt == 5) && (at == 11) ){
    omp->v1  = 45.660;    omp->v2  = -0.004;
    omp->r0  =  1.337;
    omp->a0  =  0.548;

    omp->ws1 =  0.003;    omp->ws2 =  1.219;
    omp->rs  =  1.438;
    omp->as  =  0.182;
  }

  /*** C12 */
  else if( (zt == 6) && (at == 12) ){
    omp->v1  = 58.860;    omp->v2  = -0.663;
    omp->r0  =  1.211;
    omp->a0  =  0.434;

    omp->ws1 = 12.650;    omp->ws2 =  0.045;
    omp->rs  =  1.387;
    omp->as  =  0.163;
  }

  /*** C13 */
  else if( (zt == 6) && (at == 13) ){
    omp->v1  = 61.230;    omp->v2  = -0.765;
    omp->r0  =  1.131;
    omp->a0  =  0.561;

    omp->ws1 = 16.400;    omp->ws2 =  0.136;
    omp->rs  =  1.412;
    omp->as  =  0.112;
  }

  /*** N14 */
  else if( (zt == 7) && (at == 14) ){
    omp->v1  = 50.540;    omp->v2  = -0.006;
    omp->r0  =  1.209;
    omp->a0  =  0.573;

    omp->ws1 = 12.070;    omp->ws2 =  0.703;
    omp->rs  =  1.415;
    omp->as  =  0.105;
  }

  /*** O16 */
  else if( (zt == 8) && (at == 16) ){
    omp->v1  = 48.250;    omp->v2  = -0.053;
    omp->r0  =  1.255;
    omp->a0  =  0.536;

    omp->ws1 =  4.418;    omp->ws2 =  0.556;
    omp->rs  =  1.352;
    omp->as  =  0.205;
  }

  /*** global fit */
  else{
    double eps = (double)(at - 2*zt)/(double)at;
    omp->v1  = 45.140 - 23.48*eps;  omp->v2 = -0.020;
    omp->r0  =  1.508 - 0.0133*at;
    omp->a0  =  0.500;

    omp->ws1 = 11.320 - 16.08*eps;  omp->ws2 =  0.237;
    omp->rs  =  1.353;
    omp->as  =  0.200;
  }

  omp->vso1=  5.500;
  omp->rvso=  1.150;
  omp->avso=  0.500;

  return(0x0013);
}

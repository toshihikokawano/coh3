/***********************************************************/
/*      A.B. Smith, ANL/NDM and Nucl Phys. A               */
/***********************************************************/

#include <cmath>

#include "omplib.h"

static double jvol_to_dep (double, double, double, double);
static double jsrf_to_dep (double, double, double, double);

/*************************************************/
/*    A.B.Smith ANL/NDM-120                      */
/*        for Ni58 (and modified for Ni60,Ni61   */
/*************************************************/
unsigned int OMPSmithA120(double e, int at, int zt, Optical *omp)
{
  if(zt != 28) return(0);

  double eps   = 1.0-2.0*(double)zt/(double)at;
  double eps58 = 1.0-2.0*28/58.0;

  double jv = 512.0 - 6.02*e;
  double jw = (e < 10.0) ? 133.0 - 3.9 *e : 94.0;

  /*** Asymmetry correction, Rapaport, Phys. Rep. 87, 25 (1982) */
  jv = jv - 310.0*(eps - eps58);
  jw = jw - 210.0*(eps - eps58);

  omp->vso1=  5.5;

  omp->r0  =  1.305 - 0.006*e;
  omp->a0  =  0.6461;
  omp->rs  =  (e<4.7) ? 1.50  - 0.07 *e : 1.160 + 0.002*e;
  omp->as  =  0.260 + 0.02*e;
  omp->rvso=  1.000;
  omp->avso=  0.650;

  omp->v1  = jvol_to_dep((double)at,jv,omp->r0,omp->a0);
  omp->ws1 = jsrf_to_dep((double)at,jw,omp->rs,omp->as);

  return(0x0012);
}


/*************************************************/
/*    A.B.Smith ANL/NDM-136                      */
/*        for Fe56                               */
/*************************************************/
unsigned int OMPSmithA136(double e, int at, int zt, Optical *omp)
{
  if(zt != 26 || at != 56) return(0);

  double jv  = 506.0 - 6.6837*e + 0.059882*e*e;

  double jws = (e<  7.5) ? 150.35  - 7.3301  *e :  95.975 - 0.080151*e;
  double jwv = (e< 12.0) ?   0.0                : -16.15  + 1.35    *e;

  double ef  = (e+9.225)*(e+9.225);

  omp->vso1=  5.9099;
  omp->vso2=  -0.015;


  omp->r0  =  1.300 - 0.0038*e;
  omp->a0  =  0.6338;
  omp->rv  =  omp->rs  =  1.260 - 0.0024 *e;
  omp->av  =  omp->as  =  0.70*ef/(ef + 12.65*12.65);

  omp->rvso=  1.103;
  omp->avso=  0.560;

  omp->v1  = jvol_to_dep((double)at,jv ,omp->r0,omp->a0);
  omp->ws1 = jsrf_to_dep((double)at,jws,omp->rs,omp->as);
  omp->wv1 = jvol_to_dep((double)at,jwv,omp->rv,omp->av);

  return(0x0013);
}


/*************************************************/
/*    A.B.Smith Nucl. Phys.A415                  */
/*        regional parmaeter for Z=39 to 51      */
/*************************************************/
unsigned int OMPSmithA415(int at, int zt, Optical *omp)
{
  double eps = 1.0-2.0*(double)zt/(double)at;

  omp->v1  = 52.58 - 30.0 *eps;
  omp->v2  = -0.30;
  omp->ws1 = 11.7 - 25.0 *eps - 1.8 * cos(2*3.141592*(at-90)/29.0);
  omp->vso1= 6.0;

  omp->r0  = 1.131+0.00107*at;
  omp->rs  = 2.028-0.00683*at;
  omp->rvso=  omp->r0;

  omp->a0  =  1.203-0.00511*at;
  omp->as  = -0.1061+0.005551*at;
  omp->avso=  omp->a0;

  return(0x0013);
}


/*************************************************/
/*    Convert from Volume Integral J to depth    */
/*************************************************/
double jvol_to_dep(double anum, double jv, double r, double a)
{
  if(jv == 0.0) return(0.0);

  double rx = r*pow(anum,1.0/3.0);
  double x1 = 4.0*PI/3.0*r*r*r;
  double x2 = PI*a/rx;
  double x3 = x1*(x2*x2+1.0);

  return(jv/x3);
}


double jsrf_to_dep(double anum, double jw, double r, double a)
{
  if(jw == 0.0) return(0.0);

  double rx = r*pow(anum,1.0/3.0);
  double x1 = 16.0*PI/anum*rx*rx*a;
  double x2 = PI*a/rx;
  double x3 = x1*(x2*x2/3.0+1.0);

  return(jw/x3);
}


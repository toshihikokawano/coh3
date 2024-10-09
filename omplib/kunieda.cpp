/*************************************************/
/*     S. Kunieda, S. Chiba, K. Shibata,         */
/*     A. Ichihara, E. Soukhovitskii             */
/*     J. Nucl. Sci. Technol, 44, 838 (2007)     */
/*                                               */
/*     The original FORTRAN code written by      */
/*     S. Kunieda, translated into C++           */
/*************************************************/

#include <cmath>

#include "omplib.h"

static inline double KOM_vr0   (const int, const int);
static inline double KOM_lmdvr (const int, double);
static inline double KOM_wdbw  (const int);
static inline double KOM_wdwid (const int, const int);
static inline double KOM_wcwid (const int, const int);
static inline double KOM_rr    (const int);
static inline double KOM_aro   (const int, const int);

static double at2 = 0.0, at3 = 0.0;

unsigned int OMPKunieda(double e,int at, int zt, int ai, int zi, Optical *omp)
{
  const double pdis = 2.0;
  bool proton = (zi == 1) ? true : false;

  at2 = at * at;
  at3 = at2 * at;
  double sym = -(1.0 - 2.0*zt/(double)at) * (1.0 - 2.0*zi/(double)ai);

  double ef = OMPFermiEnergy(at,zt,zi);
  double ex = e - ef;
  double ep = pow(ex,pdis);

  double ex2 = ex*ex;
  double ex3 = ex2*ex;

  /*** relativistic correction factor */
  double xm = (proton) ? 938.790002 : 939.572354;
  double fact = 2.0 * (e + xm) / (e + 2.0 * xm);

  /*** real potential depth */
  double vr0   = KOM_vr0(at,zt);
  double vr1   = 0.027;
  double vr2   = 0.000120;
  double vr3   = 0.00000035;
  double vrla  = 94.88;
  double lmdvr = KOM_lmdvr(zt,vr0);
  double ciso  = 24.3;
  double ccoul = 0.9;

  double v = vr0 + vr1 * ex + vr2 * ex2 + vr3 * ex3 + vrla * exp(-lmdvr * ex);

  v  *= 1.0 + (ciso*sym) / (vr0 + vrla);

  if(proton){
    double p = (lmdvr * vrla * exp(-lmdvr * ex) - vr1 - 2.0*vr2*ex - 3.0*vr3*ex2)
              *(1.0 + (ciso*sym) / (vr0 + vrla));
    v += ccoul * zt * pow(at,-1.0/3.0) * p;
  }
  v  /= fact;
  omp->v1 = v;

  /*** imaginary surface potential */
  double wdbw  = KOM_wdbw(zt);
  double wdwid = KOM_wdwid(at,zt);
  double lmdwd = 0.0140;
  double wciso = 18.0;

  double wd = ( wdbw + wciso * sym ) * exp(-lmdwd * ex);
  wd *= ep /(ep + pow(wdwid,pdis)) / fact;
  omp->ws1 = wd;

  /*** imaginary volume potential */
  double wv = 17.0;
  double wcwid = KOM_wcwid(at,zt);

  wv *= ep /( ep + pow(wcwid,pdis) )/ fact;
  omp->wv1 = wv;

  /*** real spin-orbit potential */
  double vs = 5.922 + 0.003 * at;
  double lmdso = 0.005;

  omp->vso1 = vs * exp(-lmdso * ex);

  /*** imaginary spin-orbint potential */
  double wsbw  = -3.1;
  double wswid = 160.0;

  omp->wso1 = wsbw * ( ep /( ep + pow(wswid,pdis) ) );

  /*** geometry parameters */
  omp->r0 = omp->rs = omp->rv = KOM_rr(zt);
  omp->a0 = omp->as = omp->av = KOM_aro(at,zt);

  omp->rvso = omp->rwso = 1.18 - 0.65 * pow((double)at,-1.0/3.0);;
  omp->avso = omp->awso = 0.59;

  omp->rc = (proton) ? 1.264 : 0.0;

  return(0x0013);
}


double KOM_vr0(const int at, const int zt)
{
  double vr0 = 0.0;

  switch(zt){
  case 12:
  case 13:
  case 14: vr0 = -39.0; break;
  case 20: vr0 = -37.5; break;
  case 26:
  case 28:
  case 29: vr0 = -37.0; break;
  case 39:
  case 40: vr0 = -37.5; break;
  case 41:
  case 42:
  case 50: vr0 = -38.0; break;
  case 58:
  case 73:
  case 74:
  case 79:
  case 82:
  case 83:
  case 90: vr0 = -38.5; break;
  case 92: vr0 = -39.0; break;
  default: vr0 = -34.40 - 6.18e-02 * at + 3.37e-04 * at2 - 6.52e-07 * at3; break;
  }
  return(vr0);
}


double KOM_lmdvr(const int zt, double vr0)
{
  double lmdvr = 0.0;

  switch(zt){
  case 13: lmdvr = 0.00430; break;
  case 20:
  case 26:
  case 28:
  case 29: lmdvr = 0.00450; break;
  case 73:
  case 74:
  case 79: lmdvr = 0.00423; break;
  case 82:
  case 83:
  case 90: lmdvr = 0.00425; break;
  case 92: lmdvr = 0.00417; break;
  default: lmdvr = 1.05e-02 + 1.63e-04 * vr0; break;
  }
  return(lmdvr);
}


double KOM_wdbw(const int zt)
{
  double wdbw = 0.0;

  switch(zt){
  case 12: wdbw = 11.5; break;
  case 13:
  case 14:
  case 20:
  case 26:
  case 28: wdbw = 11.0; break;
  case 29:
  case 39:
  case 40: wdbw = 11.5; break;
  case 41:
  case 42: wdbw = 12.0; break;
  case 50: wdbw = 12.8; break;
  case 58: wdbw = 12.5; break;
  case 73:
  case 74: wdbw = 13.0; break;
  case 79: wdbw = 13.5; break;
  case 82:
  case 83: wdbw = 14.0; break;
  case 90:
  case 92: wdbw = 15.0; break;
  default: wdbw = 11.08 + 6.83e-05 * at2; break;
  }
  return(wdbw);
}


double KOM_wdwid(const int at, const int zt)
{
  double wdwid = 0.0;

  switch(zt){
  case 12:
  case 13:
  case 14:
  case 20: wdwid  = 10.0; break;
  case 26:
  case 28:
  case 29:
  case 39:
  case 40:
  case 41:
  case 42:
  case 58:
  case 73:
  case 74:
  case 79: wdwid  = 12.0; break;
  case 82:
  case 83:
  case 90:
  case 92: wdwid  = 13.0; break;
  default: wdwid  = 12.72 - 1.54e-02 * at + 7.14e-05 * at2; break;
  }
  return(wdwid);
}


double KOM_wcwid(const int at, const int zt)
{
  double wcwid = 0.0;

  if(zt <= 14) wcwid = 90.0;
  else{
    switch(zt){
    case 20:
    case 26:
    case 28: wcwid =  90.0; break;
    case 29:
    case 39: wcwid =  95.0; break;
    case 40: wcwid =  95.0; break;
    case 41: wcwid =  98.0; break;
    case 42: wcwid =  95.0; break;
    case 50: wcwid = 103.0; break;
    case 58:
    case 73:
    case 74:
    case 79:
    case 82:
    case 83:
    case 90:
    case 92: wcwid = 105.0; break;
    default: wcwid = 105.05 - 14.13 / (1.0 + exp((at-100.21)/13.47)); break;
    }
  }
  return(wcwid);
}


double KOM_rr(const int zt)
{
  double rr = 1.210;
  if     (zt == 20) rr = 1.220;
  else if(zt == 73) rr = 1.197;
  return(rr);
}


double KOM_aro(const int at, const int zt)
{
  double aro = 0.0;

  switch(zt){
  case 12: aro = 0.560; break;
  case 13: aro = 0.520; break;
  case 14: aro = 0.540; break;
  case 20: aro = 0.610; break;
  case 26: aro = 0.610; break;
  case 28: aro = 0.620; break;
  case 29: aro = 0.620; break;
  case 39: aro = 0.650; break;
  case 40: aro = 0.635; break;
  case 41: aro = 0.650; break;
  case 42: aro = 0.650; break;
  case 50: aro = 0.655; break;
  case 58: aro = 0.650; break;
  case 73: aro = 0.660; break;
  case 74: aro = 0.660; break;
  case 79: aro = 0.660; break;
  case 82: aro = 0.671; break;
  case 83: aro = 0.670; break;
  case 90: aro = 0.685; break;
  case 92: aro = 0.685; break;
  default: aro = 0.490 + 3.22e-03*at - 2.07e-05*at2 + 4.47e-08*at3;
  }
  return(aro);
}


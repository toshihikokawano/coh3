/******************************************************************************/
/*  intwave.cpp                                                               */
/*        wave functions inside nucleus, Fox-Goodwin or Numerov method        */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "optical.h"
#include "etc.h"

static inline void            ccUcalc1 (const int, const int, const double, const double, std::complex<double> *, std::complex<double> *);
static inline std::complex<double> ccUcalc2 (const int, const int, std::complex<double> *, std::complex<double> *);


/**********************************************************/
/*      Solve Schroedinger Equation for Optical Potential */
/**********************************************************/
void omInternalFunction(
const int    integ,           // number of radial points
const double delta,           // radial integration bin width
const double wavesq,          // k^2
const double xl,              // L: orbital angular momentum
const double xs,              // S: projectile intrinsic spin
const double xj,              // J = L + S
Potential    *pot,            // optical potential
Wavefunc     *wfn)            // wave-functions (output)
{
  double deltasq = delta*delta/12.0;
  double xl_term = xl*(xl+1.0);
  double so_term = xj*(xj+1.0)-xl_term-xs*(xs+1.0);

/***************************************/
/*     Fox-Goodwin Method              */
/***************************************/
  std::complex<double> a1,a2,a3,q1,q2,q3,p1,p2,p3;

  wfn->internal[0] = p1 = std::complex<double>(0.0,0.0);
  double p2r = pow(delta,xl+1.0);
  double p2i = p2r * (-1.0e-10);
  wfn->internal[1] = p2 = std::complex<double>(p2r,p2i);

  q1 = std::complex<double>(0.0,0.0);
  q2 = xl_term*pot->r2inv[1] - wavesq - pot->mean_field[1] - so_term*pot->spin_orbit[1];
  
  for(int i=2 ; i<=integ ; i++){
    q3 = xl_term*pot->r2inv[i] - wavesq - pot->mean_field[i] - so_term*pot->spin_orbit[i];

    a1 = 1.0 - deltasq * q1;
    a2 = 2.0 + 10.0*deltasq * q2;
    p3 = a2 * p2 - a1 * p1;
    a3 = 1.0/(1.0 - deltasq * q3);

    wfn->internal[i] = a3*p3;

    q1 = q2; q2 = q3;
    p1 = p2; p2 = wfn->internal[i];
  }
}


/**********************************************************/
/*      Irregular Solution at Origin                      */
/**********************************************************/
void omIrregularFunction(const int integ, const double delta, const double wavesq, const double xl, const double xs, const double xj, Potential *pot, Wavefunc *wfn)
{
  double deltasq = delta*delta/12.0;
  double xl_term = xl*(xl+1.0);
  double so_term = xj*(xj+1.0)-xl_term-xs*(xs+1.0);

/***************************************/
/*     Fox-Goodwin Method              */
/***************************************/
  std::complex<double> a1,a2,a3,q1,q2,q3,p1,p2,p3;

  p1 = wfn->internal[integ-2];
  p2 = wfn->internal[integ-3];

  q1 = xl_term*pot->r2inv[integ-2]-wavesq - pot->mean_field[integ-2] - so_term*pot->spin_orbit[integ-2];
  q2 = xl_term*pot->r2inv[integ-3]-wavesq - pot->mean_field[integ-3] - so_term*pot->spin_orbit[integ-3];

  for(int i=integ-4 ; i>0 ; i--){
    q3 = xl_term*pot->r2inv[i] - wavesq - pot->mean_field[i] - so_term*pot->spin_orbit[i];

    a1 = 1.0 - deltasq * q1;
    a2 = 2.0 + 10.0*deltasq * q2;
    p3 = a2 * p2 - a1 * p1;
    a3 = 1.0/(1.0 - deltasq * q3);

    wfn->internal[i] = a3*p3;

    q1 = q2; q2 = q3;
    p1 = p2; p2 = wfn->internal[i];
  }
  wfn->internal[0] = std::complex<double>(0.0,0.0);
}


/**********************************************************/
/*      Normalization Factor of  Wave Function            */
/**********************************************************/
std::complex<double> omNormalizationFactor(const int l, const int i, const CoulombData cdata, std::complex<double> s, Wavefunc *wfn)
{
  /*** W = F + C(G+iF) */
  std::complex<double> c = std::complex<double>(s.imag()/2.0, (1.0 - s.real())/2.0);

  double wr = wfn->external[l].real()*c.real() + wfn->external[l].imag()*(1-c.imag());
  double wi = wfn->external[l].imag()*c.real() + wfn->external[l].real()*   c.imag() ;
  double d = norm(wfn->internal[i]);

  double xr = (wr*wfn->internal[i].real() + wi*wfn->internal[i].imag())/d;
  double xi = (wi*wfn->internal[i].real() - wr*wfn->internal[i].imag())/d;

  std::complex<double> x(xr,xi);

  /*** W = (G - iF) - S (G + iF) */
/*
  std::complex<double> w = (conj(wfn->external[l]) - s * wfn->external[l]) * std::complex<double>(0.0,1.0) * 0.5;
  std::complex<double> x = w / wfn->internal[i];
*/

  if(cdata.eta != 0.0) x *= omCoulombPhaseFactor(l,cdata);

  return(x);
}


/**********************************************************/
/*      Coulomb Phase Factor                              */
/**********************************************************/
std::complex<double> omCoulombPhaseFactor(const int l, CoulombData cdata)
{
  for(int l0=1 ; l0<=l ; l0++) cdata.sigma0 += atan(cdata.eta/l0);
  std::complex<double> coul(cos(cdata.sigma0), sin(cdata.sigma0));
  return(coul);
}


/**********************************************************/
/*      Normalize Wave Function (Distorted Wave)          */
/**********************************************************/
void omNormalization(const int integ, const double beta, std::complex<double> norm, Potential *pot, Wavefunc *wfn)
{
  double alpha = beta*beta/4.0;
  double c = 1.0;
  for(int i=0 ; i<=integ ; i++){
    if(beta > 0.0){
      c = 1.0/sqrt(1.0 + alpha*(pot->mean_field[i].real()+pot->coulomb[i]));
    }
    std::complex<double> w = c * wfn->internal[i] * norm;

    wfn->internal[i] = w;
  }
}


/**********************************************************/
/*      Schroedinger Equation for Deformed Potential      */
/**********************************************************/
void ccInternalFunction(const int m, const int integ, const double delta, CCdata *cdt, Wavefunc *wfn, std::complex<double> **v, std::complex<double> *p[])
{
  std::complex<double> u,w[MAX_CCMATRIX];

  double deltasq = delta*delta;
  double delta12 = deltasq/12.0;

/***************************************/
/*     Modified Numerov Method         */
/***************************************/

  for(int i=0 ; i<m ; i++){
    int id = i*(i+1)/2+i;
    u = delta * sqrt(v[0][id] + std::complex<double>(cdt[i].lev->wavesq,0.0));
    u = pow(u,cdt[i].chn.l+1.0);

    double f  = 1.0; for(int l=3 ; l<=(2*cdt[i].chn.l+1) ; l+=2){ f *= l; }
    u /= f;

    for(int j=0 ; j<m ; j++){
      p[0][i*m+j] = p[1][i*m+j] = std::complex<double>(0.0,0.0);
      if(i == j)  p[1][i*m+j] = u;
    }
  }

  int idx = 0, k0, k1, k2, kk;
  for(int k=1 ; k<=integ ; k++){

    k0 = (k-1)%3;
    k1 =  k   %3;
    k2 = (k+1)%3;
    kk = k - integ + 6;

    for(int i=0 ; i<m ; i++){
      ccUcalc1(m,i,deltasq,delta12,v[k],w);
      for(int j=0 ; j<m ; j++){
        u = ccUcalc2(m,j,w,p[k1]);
        idx = i*m + j;
        p[k2][idx] = 2.0*p[k1][idx] - p[k0][idx] + u;
        if(kk >= 0) wfn[idx].internal[kk] = p[k1][idx] + u/12.0;
      }
    }
  }
}


inline void ccUcalc1(const int m, const int i, const double d, const double d12, std::complex<double> *v, std::complex<double> *u)
{
  for(int l=0 ; l<m ; l++){
    std::complex<double> w(0.0,0.0);
    for(int k=0 ; k<m ; k++){
      w += v[(k <= i) ?  i*(i+1)/2+k : k*(k+1)/2+i]
         * v[(l <= k) ?  k*(k+1)/2+l : l*(l+1)/2+k];
    }
    u[l] = (v[(l <= i) ?  i*(i+1)/2+l : l*(l+1)/2+i] + d12*w)*d;
  }
}


inline std::complex<double> ccUcalc2(const int m, const int j, std::complex<double> *z, std::complex<double> *x)
{
  std::complex<double> u(0.0,0.0);
  for(int l=0 ; l<m ; l++) u += z[l] * x[l*m+j];
  return(u);
}

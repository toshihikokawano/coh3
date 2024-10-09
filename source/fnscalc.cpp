/******************************************************************************/
/*  fncalc.cpp                                                                */
/*        Madland-Nix model main part                                         */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "omcalc.h"
#include "fns.h"

#define GAUSS_LEGENDRE32
#define GAUSS_LAGUERRE32
#include "gauss.h"

/*** energy grid for inverse reaction cross sections */
static const int nesig = 80;
static double esiggrid[nesig] = {
 1.0e-06, 2.0e-06, 5.0e-06, 1.0e-05, 2.0e-05, 5.0e-05, 1.0e-04, 2.0e-04, 5.0e-04, 1.0e-03,
 2.0e-03, 3.0e-03, 4.0e-03, 6.0e-03, 8.0e-03, 1.0e-02, 2.0e-02, 3.0e-02, 4.0e-02, 6.0e-02,
 8.0e-02, 1.0e-01, 1.5e-01, 2.0e-01, 3.0e-01, 4.0e-01, 5.0e-01, 6.0e-01, 8.0e-01, 1.0e+00,
 1.1e+00, 1.2e+00, 1.3e+00, 1.4e+00, 1.5e+00, 1.6e+00, 1.7e+00, 1.8e+00, 1.9e+00, 2.0e+00,
 2.2e+00, 2.4e+00, 2.6e+00, 2.8e+00, 3.0e+00, 3.5e+00, 4.0e+00, 4.5e+00, 5.0e+00, 5.5e+00,
 6.0e+00, 7.0e+00, 8.0e+00, 9.0e+00, 1.0e+01, 1.1e+01, 1.2e+01, 1.3e+01, 1.4e+01, 1.5e+01,
 1.6e+01, 1.8e+01, 2.0e+01, 2.2e+01, 2.4e+01, 2.6e+01, 2.8e+01, 3.0e+01, 3.2e+01, 3.4e+01,
 3.6e+01, 3.8e+01, 4.0e+01, 4.5e+01, 5.0e+01, 6.0e+01, 7.0e+01, 8.0e+01, 9.0e+01, 1.0e+02},
  sigL[nesig],sigH[nesig];

static double fnsLABSpec (const double, const double, const double, double *, double *, double *, const double, const double);
static double fnsCMSSpec (const double, const double, double *, double *, const double);
static void   fnsCMSNormalization (const double, const double, double *, double *);
static double fnsInterpolateSigma (const double, double *);
static double fnsAverageCMSSpectrumEnergy (const double, double *, double *, double *, const double);


/**********************************************************/
/*      Madland-Nix Model Calculation                     */
/*      ---                                               */
/*      2-triangular temperature distributions,           */
/*      F.-J. Hambsch, et al.,  ANE 32, 1032 (2005)       */
/**********************************************************/
void fnsMadlandNixModel(const int mcn, const int ne, double *energy, double **spec, FChance *fc)
{
  double cl1[MAX_GAUSSLEG*2],ch1[MAX_GAUSSLEG*2];
  double cl2[MAX_GAUSSLEG*2],ch2[MAX_GAUSSLEG*2];

  double tmh  = fc->hf.tmax;
  double tml  = fc->lf.tmax;
  double efl  = fc->hf.za.getA() /((double)fc->lf.za.getA() * mcn) * fc->tke;
  double efh  = fc->lf.za.getA() /((double)fc->hf.za.getA() * mcn) * fc->tke;
  double s    = fc->tps;
  double b    = fc->anisotropy;

  /*** k(T) for 2-region calc. */
  double f1 = 2.0  /(s+1.0);
  double f2 = 2.0*s/(s+1.0);
  fnsCMSNormalization(0.0,f1*tml,sigL,cl1);
  fnsCMSNormalization(0.0,f1*tmh,sigH,ch1);
  if(s > 1.0){
    fnsCMSNormalization(f1*tml,f2*tml,sigL,cl2);
    fnsCMSNormalization(f1*tmh,f2*tmh,sigH,ch2);
  }

  /*** at each outgoing energy, average spetra from light and heavy fragments
       weighted by nu-bars */
  double norml = 1.0/(2.0*sqrt(efl)*tml*tml);
  double normh = 1.0/(2.0*sqrt(efh)*tmh*tmh);
  for(int i=0 ; i<ne ; i++){
    spec[0][i] = norml*fnsLABSpec(energy[i],efl,tml,sigL,cl1,cl2,s,b);
    spec[1][i] = normh*fnsLABSpec(energy[i],efh,tmh,sigH,ch1,ch2,s,b);
    spec[2][i] = (fc->nuratio*spec[0][i] + spec[1][i])/(fc->nuratio + 1.0);
  }

  /*** re-normalize spetrum, due to the limited energy range, and average energy  */
  for(int i=0 ; i<3 ; i++){
    fnsNormalizeSpectrum(ne,energy,spec[i]);
  }

  fc->ecms[0] = fnsAverageCMSSpectrumEnergy(tml,sigL,cl1,cl2,s);
  fc->ecms[1] = fnsAverageCMSSpectrumEnergy(tmh,sigH,ch1,ch2,s);
  fc->ecms[2] = (fc->nuratio*fc->ecms[0] + fc->ecms[1])/(fc->nuratio + 1.0);
}


/**********************************************************/
/*      Temperature Dependent Normalization, k(T)         */
/**********************************************************/
void fnsCMSNormalization(const double b1, const double b2, double *sigma, double *c)
{
  for(int i=0 ; i<MAX_GAUSSLEG ; i++){
    double t1 = 0.5*(1-gaussleg_x[i])*(b2-b1) + b1;
    double t2 = 0.5*(1+gaussleg_x[i])*(b2-b1) + b1;
    
    double z1 = 0.0, z2 = 0.0;
    for(int j=0 ; j<MAX_GAUSSLAG ; j++){
      double x1 = gausslag_x[j] * t1;
      double x2 = gausslag_x[j] * t2;
      
      z1 += x1*t1*fnsInterpolateSigma(x1,sigma)*gausslag_a[j];
      z2 += x2*t2*fnsInterpolateSigma(x2,sigma)*gausslag_a[j];
    }
    c[i]                  = (z1>0.0) ? 1.0/z1 : 0.0;
    c[2*MAX_GAUSSLEG-i-1] = (z2>0.0) ? 1.0/z2 : 0.0;
  }
}


/**********************************************************/
/*      Integrate Over Available CMS Energy               */
/**********************************************************/
double fnsLABSpec(const double e, const double ef, const double tm, double *sigma, double *c1, double *c2, const double s, const double b)
{
  /*** integrate over CMS outgoing energy, the range defined [elo,ehi] */
  double elo = sqrt(e) - sqrt(ef); elo = elo*elo;
  double ehi = sqrt(e) + sqrt(ef); ehi = ehi*ehi;

  double z = 0.0;
  double c = 3.0/(3.0+b);

  for(int i=0 ; i<MAX_GAUSSLEG*2 ; i++){
    double x = (i<MAX_GAUSSLEG) ? 0.5*(1-gaussleg_x[i])
                                : 0.5*(1+gaussleg_x[2*MAX_GAUSSLEG-1-i]);
    double a = (i<MAX_GAUSSLEG) ? gaussleg_a[i] : gaussleg_a[2*MAX_GAUSSLEG-1-i];
    double eps = (ehi - elo)*x + elo;

    /*** integ over temperature from zero to Tm */
    z += sqrt(eps) * fnsInterpolateSigma(eps,sigma) * fnsCMSSpec(eps,tm,c1,c2,s) * a
       * (1.0 + b*(e-eps-ef)*(e-eps-ef)/(4.0*eps*ef));
  }
  z *= (ehi-elo)*0.5*c;

  return(z);
}


/**********************************************************/
/*      Integrate CMSpec Over Temperature Distribution    */
/**********************************************************/
double fnsCMSSpec(const double eps, const double tm, double *c1, double *c2, const double s)
{
  double f1 = (s+1.0)/(2.0*s);
  double f2 = tm/(s+1.0);

  double z1 = 0.0,  z2 = 0.0;
  for(int j=0 ; j<MAX_GAUSSLEG*2 ; j++){
    double x = (j<MAX_GAUSSLEG) ? 0.5*(1-gaussleg_x[j])
                                : 0.5*(1+gaussleg_x[2*MAX_GAUSSLEG-1-j]);
    double a = (j<MAX_GAUSSLEG) ? gaussleg_a[j] : gaussleg_a[2*MAX_GAUSSLEG-1-j];
    x = 2.0*x-1.0;

    /*** first part of integral */
    double t1 = (x+1.0)*f2;
    z1 += t1 * exp(-eps/t1) * c1[j] * a; 

    /*** second part of integral */
    if(s > 1.0){
      double t2 = x*(s-1.0)*f2 + tm;
      z2 += (-1.0*t2*f1 + tm) * exp(-eps/t2) * c2[j] * a;
    }
  }  

  double z = (z1*f1 + z2) * tm * 0.5;

  return(z);
}


/**********************************************************/
/*     Prepare Inverse Reaction Cross Section             */
/**********************************************************/
void fnsInverseCrossSection(Pdata *pn, ZAnumber *zal, ZAnumber *zah)
{
  CrossSection cx;
  double *tc = new double [3*MAX_J];

  ZAnumber lf = *zal - pn->za; // for inv. reaction, target should be FF mass - 1
  ZAnumber hf = *zah - pn->za;

  double mul = lf.getA() * pn->mass / (lf.getA() + pn->mass);
  double muh = hf.getA() * pn->mass / (hf.getA() + pn->mass);

  for(int i=0 ; i<nesig ; i++){
    omCalc(esiggrid[i], pn, &lf, mul, tc, &cx);  sigL[i] = cx.reaction;
    omCalc(esiggrid[i], pn, &hf, muh, tc, &cx);  sigH[i] = cx.reaction;
  }

  delete [] tc;
}


/**********************************************************/
/*      Average Energy for CMS Spectrum                   */
/**********************************************************/
double fnsAverageCMSSpectrumEnergy(const double tm, double *sigma, double *c1, double *c2, const double s)
{
  /*** integrate CMS spectrum from zero to MAX_SPECTRUM_ENERGY using Simpson method
       energy region divided by two, 0-100 keV and 100 keV - 30 MeV */

  int    ndiv = 100;
  double em   = 1.0;  // upper boundary of first integration
  double de   = em/(double)ndiv;

  double s1 = 0.0, e1 = 0.0;
  for(int i=1 ; i<=ndiv ; i++){
    double eps = i*de;
    double a   = (i%2 == 0) ? 2 : 4;
    if(i == 100) a = 1.0;
    double phi = 2.0*eps/(tm*tm)*fnsInterpolateSigma(eps,sigma)*fnsCMSSpec(eps,tm,c1,c2,s);
    s1 += phi*a;
    e1 += phi*a*eps;
  }
  s1 *= de/3.0;
  e1 *= de/3.0;

  double s2 = 0.0, e2 = 0.0;
  de = (MAX_SPECTRUM_ENERGY - em)/ndiv;
  for(int i=0 ; i<=ndiv ; i++){
    double eps = em + i*de;
    double a   = (i%2 == 0) ? 2 : 4;
    if( (i == 0) || (i == ndiv) ) a = 1.0;
    double phi = 2.0*eps/(tm*tm)*fnsInterpolateSigma(eps,sigma)*fnsCMSSpec(eps,tm,c1,c2,s);
    s2 += phi*a;
    e2 += phi*a*eps;
  }
  s2 *= de/3.0;
  e2 *= de/3.0;

  return( (e1+e2)/(s1+s2) );
}


/**********************************************************/
/*      Look Up Inverse Reaction Cross Section            */
/**********************************************************/
double fnsInterpolateSigma(const double e, double *sigma)
{
  double s = 0.0;

  if(e <= 0.0){
    s = 0.0;
  }
  /*** if lower than the first energy point, use 1/v */
  else if(e < esiggrid[0]){
    s = sigma[0] * sqrt(esiggrid[0] / e);
  }
  /*** if larger than the defined region, use the last point */
  else if(e > esiggrid[nesig-1]){
    s = sigma[nesig-1];
  }
  /*** linear interpolation */
  else{
    for(int i=1 ; i<nesig ; i++){
      if(e < esiggrid[i]){
        s = (sigma[i]-sigma[i-1])/(esiggrid[i]-esiggrid[i-1])*(e-esiggrid[i-1]) + sigma[i-1];
        break;
      }
    }
  }

  return (s);
}


/**********************************************************/
/*     Renormalize Spectrum                               */
/**********************************************************/
void fnsNormalizeSpectrum(const int n, double *x, double *y)
{
  double s = 0.0;
  for(int i=0 ; i<n-1 ; i++){
    s += (y[i+1]+y[i])*(x[i+1]-x[i])*0.5;
  }
  if(s > 0.0){
    s = 1.0/s;
    for(int i=0 ; i<n ; i++) y[i] *= s;
  }
}


/**********************************************************/
/*      Average Secondary Neutron Energy                  */
/**********************************************************/
double fnsAverageSpectrumEnergy(const int n, double *x, double *y)
{
  /*** since Y(i) and Y(i+1) vary linearly,
       int_a^b xY dx = (b-a)/6 [Y(i+1)(b+2a) + Y(i)(a+2b)] */
  double s1 = 0.0, s2 = 0.0;
  for(int i=0 ; i<n-1 ; i++){
    double d0 = x[i+1] - x[i];
    double d1 = x[i+1] + 2*x[i  ];
    double d2 = x[i  ] + 2*x[i+1];
    s1 += (y[i+1]    + y[i]   ) * d0 * 0.5;
    s2 += (y[i+1]*d1 + y[i]*d2) * d0 / 6.0;
  }
  double e = (s1 != 0.0) ? s2/s1 : 0.0;

  return(e);
}


/**********************************************************/
/*     Renormalize Excitation Energy Distribution         */
/**********************************************************/
void fnsNormalizeExcitationDist(const int n, double *x, double *y)
{
  double s = 0.0;
  for(int i=0 ; i<n-1 ; i++){
    s += y[i] * fabs(x[i+1]-x[i]);
  }
  if(s > 0.0){
    s = 1.0/s;
    for(int i=0 ; i<n ; i++) y[i] *= s;
  }
}



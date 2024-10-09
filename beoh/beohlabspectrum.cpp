/******************************************************************************/
/*  beohlabspectrum.cpp                                                       */
/*        convert CMS neutron spectrum into LAB                               */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "beoh.h"
#include "global.h"

#define GAUSS_LEGENDRE20
#include "gauss.h"

static double beohCmsSpectrumInterpolate(const int, const double, const double, double *);

//static const double eps = 1.0e-99;
//static inline double lowfilter (double x){if(fabs(x)<eps) return 0.0; else return(x);}

/**********************************************************/
/*      Integrate Over TKE Distribution                   */
/**********************************************************/
void beohIntegCmsSpectrum(const int n, const int nc, const unsigned int a, BData *bdt, double de, double *spc)
{
  double *y, *z;
  y = new double [n];
  z = new double [n];
  for(int k=0; k<n ; k++) y[k] = z[k] = 0.0;

  int k0 = (bdt[0].exmean - EXRANGE_FACTOR * bdt[0].exwidth)/de;
  int k1 = (bdt[0].exmean + EXRANGE_FACTOR * bdt[0].exwidth)/de;
  if(k0 <= 0) k0 = 0;

  double q = bdt[0].ekmean + bdt[0].exmean;

  /*** Gaussian distribution */
  double s = 0.0;
  for(int i=0 ; i<nc ; i++){
    for(int k=k0 ; k<=k1 ; k++){
      double x = (k+0.5) * de;
      double t = q - x;
      if(t < 0.0){ k1 = k-1; break; }

      y[k] += exp(-(x - bdt[i].exmean)*(x - bdt[i].exmean)/(2.0*bdt[i].exwidth*bdt[i].exwidth));
      s += y[k] * bdt[i].fraction;
    }
  }

  /*** normalize distribution */
  if(s > 0.0){
    for(int k=k0 ; k<=k1 ; k++) s += y[k];
    for(int k=k0 ; k<=k1 ; k++) y[k] /= s;
  }

  double etop = beohZeroCut(spc) * de;

  for(int k=k0 ; k<=k1 ; k++){
    double x = (k+0.5) * de;
    double ef = (q - x) / (double)a;
    double efsq = sqrt(ef);
    double emax = sqrt(etop) + efsq; emax = emax * emax;

    for(int i=0 ; i<n ; i++){
      double el = (i+0.5)*de;
      if(el > emax + 1.0) break;
      double sl = beohLabConvert(n,sqrt(el),efsq,etop,de,spc);
      z[i] += y[k] * sl;
    }
  }

  for(int k=0 ; k<n ; k++) spc[k] = z[k];

  delete [] y;
  delete [] z;
}


/**********************************************************/
/*      Convert CMS Spectrum into LAB Spectrum            */
/**********************************************************/
void beohLabSpectrum(const int n, const double ef, const int k0, const double de, double *spc, double *spl)
{
  if(ef <= 0.0) return;

  double efsq = sqrt(ef);
  double etop = de * k0;

  for(int k=0 ; k<n ; k++) spl[k] = beohLabConvert(n,sqrt((k+0.5)*de),efsq,etop,de,spc);
}


double beohLabConvert(const int n, const double elsq, const double efsq, const double emax, const double de, double *spc)
{
  double spl = 0.0;

  double e1 = elsq - efsq;   e1 = e1*e1;  // {sqrt(E) - sqrt(Ef)}^2
  double e2 = elsq + efsq;   e2 = e2*e2;  // {sqrt(E) + sqrt(Ef)}^2
  double e3 = (e2+e1)*0.5;
  double e4 = (e2-e1)*0.5;

  if(e1 > emax) return spl;

  /*** Gauss-Legendre Integration in [e1,e2] */
  for(int ig=0 ; ig<MAX_GAUSSLEG ; ig++){
    double x1 = e3 + e4*gaussleg_x[ig];
    double x2 = e3 - e4*gaussleg_x[ig];
    double s1 = beohCmsSpectrumInterpolate(n,x1,de,spc);
    double s2 = beohCmsSpectrumInterpolate(n,x2,de,spc);
    spl += (s1/sqrt(x1) + s2/sqrt(x2)) * gaussleg_a[ig];
  }
  spl *= e4*0.25/efsq;  // int phi(e)/(4*sqrt(e Ef)) de
  
  return spl;
}


/***********************************************************/
/*      Interpolate CMS Neutron Spectrum                   */
/***********************************************************/
double beohCmsSpectrumInterpolate(const int n, const double x, const double de, double *s)
{
  double y = 0.0;

  if(x <= 0.0)  return(y);

  double e1min = 0.0;
  double e1max = de*0.5;
  double e1    = (e1min + e1max)*0.5;

  double e2min = e1max;
  double e2max = e2min + de;
  double e2    = (e2min + e2max)*0.5;


  /*** if lower than the first point, interpolate by sqrt(E) */ 
  if( x < e1 ){
    y = s[0] * sqrt(x/e1);
    return(y);
  }

  /*** linear interpolation */
  for(int k=1 ; k<n ; k++){
    e1    = (e1min + e1max)*0.5;
    e2    = (e2min + e2max)*0.5;

    if( (e1 <=x ) && (x < e2) ){
      y = (x-e1)/(e2-e1) * (s[k]-s[k-1]) + s[k-1];
      break;
    }

    e1min = e1max;
    e1max = e2max;
    e2min = e2max;
    e2max = e2min + de;
  }

  return(y);
}

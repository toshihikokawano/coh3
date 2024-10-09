/******************************************************************************/
/*  FRDMshell.cpp                                                             */
/*        Strutinsky Shell Correction Energy                                  */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "FRDM.h"

static double ShCLambda (const int, const int, const double, double *);
static double ShCParticleNumber (const int, const double, const double, double *);
static double ShCIntegEnergy (const int, const double, const double, double *);
static double ShCSmoothDensity (const int, const double, const double, double *);

static double ShCHermite (const int m, const double);
static double ShCHermiteCoeff (const int m);

extern double *fact;

static const int p_order = 8;
static double lambda = 0.0;

/**********************************************************/
/*      Shell Correction Energy                           */
/**********************************************************/
double FRDMShellEnergy(MFTSystem *sys, const int np, SPEnergy *spe, SCParameter *scp)
{
  const double gamma = gCurvature / sys->a13 * sys->bs;

  /*** determine the highest level number, and double it */
  int nmax = spe->ne * 2;

  /*** arrange s.p. levels by taking account of degeneracy of two */
  double *epsn;
  epsn = new double [nmax];
  for(int i=0 ; i<nmax ; i++) epsn[i] = spe->energy[i/2];

  /*** lambda-bar */
  lambda = ShCLambda(np,nmax,gamma,epsn);

  /*** shell correction energy */
  double esmac = ShCIntegEnergy(nmax,lambda,gamma,epsn);

  /*** save some quantities for later BCS calculation */
  scp->chemical_potential = lambda;
  scp->rho = 0.5*ShCSmoothDensity(nmax,lambda,gamma,epsn);

  /*
  double de = 0.2;
  for(int k=0 ; ; k++){
    double e = k*de -50;
    double x = ShCParticleNumber(nmax,e,gamma,epsn);
    std::cout << x << " "  << e << std::endl;
    if(e >= 10.0) break;
  }
  std::cout << std::endl;
  std::cout << std::endl;
  for(int n=0 ; n<300 ; n++){
    std::cout << n <<" " <<epsn[n] << std::endl;
  }
  double s = 0.0, z = 0.0;
  for(int n=0 ; n<np ; n++){
    double e1 = ShCLambda(n,nmax,gamma,epsn);
    s += e1 - epsn[n];
    if(n == 0) z += epsn[n]*0.5;
    else z += epsn[n];
    std::cout << n <<" " <<epsn[n] << " "<<  e1 << " "  << s<< std::endl;
  }
  std::cout << esmac <<" " << z << " " << z - esmac << std::endl;
*/


  delete [] epsn;
  return(esmac);
}


/**********************************************************/
/*      Sum of Single Particle Energies                   */
/**********************************************************/
double FRDMSumSingleParticleEnergy(const int np, SPEnergy *spe)
{
  /*** sum of single particle energies */
  double es = 0.0;
  for(int i=0 ; i<np ; i++) es += spe->energy[i/2];

  return(es);
}


/**********************************************************/
/*      Lambda-bar by Integrating n(e)                    */
/**********************************************************/
double ShCLambda(const int np, const int nmax, const double g, double *ep)
{
  const double eps = 1e-06;
  const int kmax = 200;

  double lambda = 0.0;
  double e1 = ep[0];
  double e2 = ep[nmax-1];

  for(int k=0 ; k<kmax ; k++){
    lambda = 0.5*(e1 + e2);
    if(fabs(e1 - e2) < eps) break;

    double x = ShCParticleNumber(nmax,lambda,g,ep);
    if(x < (double)np){ e1 = lambda; }
    else{               e2 = lambda; }
  }

  return(lambda);
}


/**********************************************************/
/*      Function n-bar(e)                                 */
/*      Bolsterli et al. PRC 5, 1050 (1972), Eq.(7)       */
/**********************************************************/
double ShCParticleNumber(const int nmax, const double e, const double g, double *ep)
{
  double s = 0.0;
  for(int i=0 ; i<nmax ; i++){
    double u = (e - ep[i])/g;
    double x1 = 0.5*(1.0 + erf(u));
    double x2 = 0.0;
    for(int m=2 ; m<=p_order ; m+=2){
      x2 += ShCHermiteCoeff(m) * ShCHermite(m-1,u);
    }
    s += x1 - exp(-u*u)/SQRTPI * x2;
  }

  return(s);
}


/**********************************************************/
/*      Function Int eps(n) dn                            */
/*      Bolsterli et al. PRC 5, 1050 (1972), Eq.(9)       */
/**********************************************************/
double ShCIntegEnergy(const int nmax, const double ef, const double g, double *ep)
{
  double s = 0.0;
  for(int i=0 ; i<nmax ; i++){
    double u = (ef - ep[i])/g;
    double x1 = 0.5*ep[i]*(1.0 + erf(u));
    double x2 = exp(-u*u)/SQRTPI;
    double x3 = 0.0;
    for(int m=2 ; m<=p_order ; m+=2){
      double cm = ShCHermiteCoeff(m);
      double hm = 0.5*g * ShCHermite(m  ,u)
                 + ep[i]* ShCHermite(m-1,u)
                 + m*g  * ShCHermite(m-2,u);
      x3 += cm*hm;
    }
    s += x1 - 0.5*g*x2 - x2*x3;
  }

  return(s);
}


/**********************************************************/
/*      Smooth Part of Function g(e) for rho_bar          */
/*      Bolsterli et al. PRC 5, 1050 (1972), Eq.(6)       */
/**********************************************************/
double ShCSmoothDensity(const int nmax, const double ef, const double g, double *ep)
{
  double s = 0.0;
  for(int i=0 ; i<nmax ; i++){
    double u = (ef - ep[i])/g;
    double x1 = exp(-u*u);
    double x2 = 0.0;
    for(int m=0 ; m<=p_order ; m+=2){
      x2 += ShCHermiteCoeff(m) * ShCHermite(m,u);
    }
    s += x1 * x2;
  }
  s = s/(SQRTPI*g);

  return(s);
}


/**********************************************************/
/*      Hermite Polynomial and Coefficient                */
/**********************************************************/
double ShCHermite(const int n, const double x)
{
  double fh = 0.0;

  if(n == 0) fh = 1.0;
  else{
    int p = n%2;

    if((x == 0.0) && (p == 0)){
      fh = ((n/2)%2 == 0 ? 1.0 : -1.0) * exp( fact[n] - fact[n/2] );
    }
    else{
      double s = 0.0;
      for(int m=0 ; m<=(n-p)/2 ; m++){
        s += (m%2 == 0 ? 1.0 : -1.0) * pow(2.0*x,(n-2.0*m)) / exp(fact[m] + fact[n-2*m]);
      }
      fh = s * exp(fact[n]);
    }
  }
  return(fh);
}

double ShCHermiteCoeff(const int m)
{
  if(m%2 != 0) return(0.0);

  int m2 = m/2;
  double x1 = pow(2.0,(double)m);
  double x2 = exp(fact[m2]);
  return( ((m2%2 == 0) ? 1.0 : -1.0) / (x1 * x2) );
}

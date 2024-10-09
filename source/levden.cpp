/******************************************************************************/
/*  levden.cpp                                                                */
/*        Gilbert-Cameron level density calculation                           */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "levden.h"
#include "etc.h"

/*** lower threshold of spin distribtion in the level density
     set zero for no threshold */
static const double SpinDistributionCut = 1e-10;


/**********************************************************/
/*      A.Gilbert and A.G.W.Cameron                       */
/*      Can. J. Phys., 43, 1446(1965)                     */
/*      modified for Ignatyuk Shell Correction            */
/**********************************************************/
double ldLevelDensity(const double excitation_energy, const double mass, LevelDensity *ldp)
{
  double dens = 0.0;

  /*** if Em = 0, use const.temp all */
  if(ldp->match_energy == 0.0){
    dens = ldConstantTemperature(excitation_energy,ldp->E0,ldp->temperature);
  }
  else{
    double u  = excitation_energy - ldp->pairing_energy;
    double a  = ldDensityParameter(excitation_energy, mass, ldp);
    ldp->spin_cutoff = ldSpinCutoff(u,a,mass,ldp->sigma0); 
    dens = (excitation_energy < ldp->match_energy) ?
           ldConstantTemperature(excitation_energy,ldp->E0,ldp->temperature) :
           ldFermiGas(u,a,ldp->spin_cutoff);
  }

  return(dens);
}


/**********************************************************/
/*      Energy Dependent a-Parameter                      */
/**********************************************************/
double ldDensityParameter(const double excitation_energy, const double mass, LevelDensity *ldp)
{
  double a = 0.0;

  if(ldp->match_energy == 0.0) return(ldp->a);
  else{
    double u  = excitation_energy - ldp->pairing_energy;
    double um = ldp->match_energy - ldp->pairing_energy;
    a = (excitation_energy < ldp->match_energy) ?
        ldShellCorrection(um,ldp->a,ldp->shell_correct,mass) :
        ldShellCorrection(u ,ldp->a,ldp->shell_correct,mass);
  }

  return(a);
}


/**********************************************************/
/*      Gaussian Spin Distribution                        */
/**********************************************************/
double ldSpinDistribution(const double j, const double sc, const double sc0, const double f)
{
  double s = cfmax(sc,sc0) * f;
  double r = (j+0.5)*exp( -(j+0.5)*(j+0.5)/(2.0*s*s) )/(s*s);
  if( (SpinDistributionCut > 0.0) && (r < SpinDistributionCut) ) r = 0.0;

  return(r);
}


/**********************************************************/
/*      Parity Distribution                               */
/**********************************************************/
double ldParityDistribution(const int pt)
{
  double pd = 0.5;
  if(pt > 0) return(pd);
  else       return(1.0-pd);
}


double ldParityDistributionEdepend(const int pt, const double e, const double f)
{
  /*** parity unevenness is linearly interporlated up to Ex */
  double ex = 10.0;
  double pd0 = 0.5, pd = 0.5;
  if(e <= ex){
    if(pt > 0) pd = pd0       * (( 1.0 - f)/ex*e + f);
    else       pd = (1.0-pd0) * ((-1.0 + f)/ex*e + 2.0 - f);
  }
  else{
    if(pt > 0) pd = pd0;
    else       pd = 1.0-pd;
  }
  return(pd);
}


/**********************************************************/
/*      Spin Cut-Off from Discrete Levels                 */
/**********************************************************/
double ldLevelSpinCutoff(const int n, Level *lev)
{
  double sig2 = 0.0;

  if(n <= 1) return(0.0);
  for(int k=0 ; k<n ; k++){
    sig2 += (lev[k].spin+1.0)*lev[k].spin;
  }

  return( sqrt(sig2/(2.0*n)) );
}


/**********************************************************/
/*      Smooth Connection of Rho(T) and Rho(G)            */
/**********************************************************/
void ldTGconnect(const double mass, const double lve, const double lvn, LevelDensity *ldp)
{
  if(lvn <= 2.0) return;

  double lvu  = lve - ldp->pairing_energy;
  double a    = ldShellCorrection(lvu,ldp->a,ldp->shell_correct,mass);
  double ex1  = cfmax(ldp->pairing_energy,9.0/(4.0*a));
  double ex2  = 100.0;

  double t, e01, e02, ex, ux, eps = 0.0;

  int i=0;
  do{
    ex  = (ex1+ex2)/2.0;
    ux  = ex - ldp->pairing_energy;

    /*** if Ux is lower than max discrete level energy, a-val is constant */
    a   = ldShellCorrection(cfmax(lvu,ux),ldp->a,ldp->shell_correct,mass);
    t   = 1.0/(sqrt(a/ux)-3.0/(2.0*ux)); 

    /*** cannot connect const. temp and Fermi gas */
    if( (t < 0.0) || ((i++) > 100) ){
      t   = 0.0;
      e01 = 0.0;
      ex  = 0.0;
      break;
    }
    e01 = t*log( (exp(lve/t)-1.0) / lvn);
    e02 = ex-t*log(t*ldFermiGas(ux,a,ldSpinCutoff(ux,a,mass,ldp->sigma0)));

    if( (eps = e01-e02) < 0.0 ) ex1 = ex;
    else                        ex2 = ex;

  }while(fabs(eps) > 1e-6);

  ldp->temperature   = t;
  ldp->E0            = e01;
  ldp->match_energy  = ex;
}


/**********************************************************/
/*      Estimate a-parameter from Rho(T)                  */
/**********************************************************/
void ldTextrapolate(const double mass, const double lve, LevelDensity *ldp)
{
  double ex1  = cfmax(ldp->pairing_energy,lve);
  double ex2  = 100.0;

  double a=0.0, r1, r2, ex, ux, eps = 0.0;

  int i=0;
  do{
    ex = (ex1+ex2)/2.0;
    ux = ex-ldp->pairing_energy;
    a  = 1.0/ldp->temperature+1.5/ux;
    a  = a*a*ux;

    r1 = ldConstantTemperature(ex,ldp->E0,ldp->temperature);
    r2 = ldFermiGas(ux,a,ldSpinCutoff(ux,a,mass,ldp->sigma0));

    /*** cannot connect const. temp and Fermi gas */
    if((i++) > 100){
      a = 0.0;
      ldp->match_energy = 0.0;
      break;
    }

    if( (eps = r1-r2) < 0.0 ) ex2=ex;
    else                      ex1=ex;

  }while(fabs(eps) > 1e-6);

  if(a > 0.0){
    ldp->a = a/ldShellCorrection(ux,1.0,ldp->shell_correct,mass);
    ldp->match_energy = ux + ldp->pairing_energy;
  }
}


/**********************************************************/
/*      Extrapolate Rho(G) Down to Zero                   */
/**********************************************************/
void ldGextrapolate(const double mass, LevelDensity *ldp)
{
  double ux = 0.5;
  double a  = ldShellCorrection(ux,ldp->a,ldp->shell_correct,mass);
  double r  = ldFermiGas(ux,a,ldSpinCutoff(ux,a,mass,ldp->sigma0));

  ldp->match_energy = ldp->pairing_energy + ux;
  ldp->E0 = ldp->match_energy - ldp->temperature * log(ldp->temperature * r);
}


/**********************************************************/
/*      Fermi Gas Model                                   */
/**********************************************************/
double ldFermiGas(const double u, const double a, const double s)
{
  if(u < 0.0) return(0.0);
  else{
    double c = 1.0/(12.0*sqrt(2.0));
    return( c*exp(2.0*sqrt(a*u))/(pow(a,0.25)*pow(u,1.25))/s );
  }
}


/**********************************************************/
/*      Constant Temperature Model                        */
/**********************************************************/
double ldConstantTemperature(const double e, const double e0, const double temp) 
{
  return(exp((e-e0)/temp)/temp);
}


/**********************************************************/
/*      Spin Cut-Off Parameter                            */
/**********************************************************/
double ldSpinCutoff(double u, const double a, const double mass, const double sig0)
{
  double sig2 = 0.0;
  if(u < 0.5) u = 0.5; /* In order to avoid sqrt(negative). The value 0.5 is empirical */
//  sig2 = 0.146*sqrt(a*u)*pow(mass,(2.0/3.0)); // original Gilbert - Cameron
//  sig2 = 0.0034725*sqrt(u/a)*pow(mass,(5.0/3.0));  // jNST 43
    sig2 = 0.006945 *sqrt(u/a)*pow(mass,(5.0/3.0));

  sig2 = sqrt(sig2);

  return( cfmax(sig2,sig0) );
}


/**********************************************************/
/*      Shell Energy Correction                           */
/**********************************************************/
double ldShellCorrection(const double u, const double astar, const double dw, const double a)
{
  double gamma = 0.31*pow(a,-1/3.0);
  double ax    = astar;

  if(u > 0.0){
    ax = (u == 0.0) ? astar * (1.0 + gamma*dw) 
                    : astar * (1.0 + (1.0-exp(-gamma*u))*dw/u);
  }
  return(ax);
}

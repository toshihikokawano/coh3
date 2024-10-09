/******************************************************************************/
/*  fisdens.cpp                                                               */
/*        Fission level density calculation                                   */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "levden.h"
#include "parameter.h"

#define  FISDENS_COH
#undef  FISDENS_GNASH
#undef  FISDENS_CCONE

#ifdef FISDENS_COH
static double  fldDensityEnhancement (const int, const double, ZAnumber *, LevelDensity *);
#endif
#ifdef FISDENS_GNASH
static double  fldDensityEnhancementGNASH (const double, const double, const double);
#endif
#ifdef FISDENS_CCONE
static double  fldFissionDensityCCONE (const int, const double, const double, ZAnumber *, LevelDensity *);
#endif

static double  fldHighestTransitionState       (Barrier *);
//static inline double fldVibrationalEnhancement (double, double);
static inline double fldRotationalEnhancement (const int, const int, const double, const double);
static inline double fldSigmaParalell (const double, const double, const double, const double);
static inline double fldSigmaPerpendicular (const double, const double);


/**********************************************************/
/*      Fission Level Density                             */
/**********************************************************/
void statSetupFissionLevelDensity(Nucleus *n, Fission *fb)
{
  for(int i=0 ; i<MAX_HUMP ; i++){

    if(fb->barrier[i].height == 0.0) continue;

    double fl = 1.0;
    switch(i){
    case 0: fl = parmGetFactor(parmFL1 ,n->za); break;
    case 1: fl = parmGetFactor(parmFL2 ,n->za); break;
    case 2: fl = parmGetFactor(parmFL3 ,n->za); break;
    default: break;
    }

    /*** determine the highest discrete transition state */
    if(fb->barrier[i].elmax == 0.0) fb->barrier[i].elmax = fldHighestTransitionState(&fb->barrier[i]);

    /*** calculate fission level densities in each bin */
    for(int k=0 ; k<n->ntotal ; k++){
      for(int j=0 ; j<MAX_J ; j++){

        /*** in case fission density is already stored */
        double d1 = fb->barrier[i].density[k][j].even * fl;
        double d2 = fb->barrier[i].density[k][j].odd  * fl;

        /*** normal calculation */
#ifdef FISDENS_COH
        double x = fldDensityEnhancement(i,n->excitation[k],&n->za,&n->ldp) * fl;
        d1 = n->density[k][j].even * x;
        d2 = n->density[k][j].odd  * x;
#endif

        /*** GNASH compatible enhancement, parameger given by fl */
#ifdef FISDENS_GNASH
        double x = fldDensityEnhancementGNASH(n->excitation[k],fl,n->ldp.pairing_energy);
        d1 = n->density[k][j].even * x;
        d2 = n->density[k][j].odd  * x;
#endif

        /*** fission density in CCONE emulation */
#ifdef FISDENS_CCONE
        double j0 = halfint(n->lev[0].spin);
        d1 = fldFissionDensityCCONE(i,j0+j,n->excitation[k],&n->za,&n->ldp) * fl;
        d2 = d1;
#endif

        /*** replace fission density by calculated ones */
        fb->barrier[i].density[k][j].even = d1;
        fb->barrier[i].density[k][j].odd  = d2;
      }
    }
  }
}


/**********************************************************/
/*      Collective Enhancement Factor for Fission         */
/**********************************************************/
double fldDensityEnhancement(const int i, const double ex, ZAnumber *za, LevelDensity *ldp)
{
  /*** normal level density */
  double mass = (double)(za->getA());
  double eshell = 2.6;
  if(i > 0) eshell = 0.6 + 0.1*(za->getZ()-97.0) + 0.04*(za->getN() - 143.0);

  double pfac = 1.5; // pairing energy scaling factor at the suddle point
  double um   = 2.0;
  double pair = ldp->pairing_energy*pfac;
  double em   = um + pair;
  double ux   = ex - pair;

  if(ex < em)  ux = um;

  /*** level density parameters */
  double a    = ldShellCorrection(ux,ldp->a,eshell,mass);
  double sig  = ldSpinCutoff(ux,a,mass,ldp->sigma0);
  double t    = sqrt(ux/a);


  /*** sigma2 for paralell and perpendicular at barriers,
       beta_2  inner and outer barriers assumed */
  double b2g = 0.2;
  double b2f = (i == 0) ? 0.4 : 0.8;

  /*** rotational enhancement, devided by g.s. enhancement
       sqrt(2) from the factor 0.006945, instead of 0.0138 */
  double ef0 = fldSigmaPerpendicular(b2g,sig*sqrt(2.0));
  double ef1 = fldRotationalEnhancement(i,za->getN(),
                                       fldSigmaPerpendicular(b2f,sig),
                                       fldSigmaParalell(b2f,a,t,mass));
  double ef  = ef1/ef0;

  if((ex < pair) && (pair > 0.0)) ef = (ef - 2.0)*ex/pair + 2.0;

  /*** empirical damping factor, Fermi functional form */
  double df = 1.0 / (1.0 + exp( (ux - 10.0)/2.0 ) );
  ef *= df;

  if(ef < 1.0) ef = 1.0;

  return(ef);
}


/**********************************************************/
/*      Vibrational Enhancement Factor                    */
/**********************************************************/
// double fldVibrationalEnhancement(double mass, double t)
// {
//   double kvib = 1.0;
//   kvib = exp(0.0555 * pow(mass,2.0/3.0) * pow(t,4.0/3.0));
//   return(kvib);
// }


/**********************************************************/
/*      Rotational Enhancement Factor                     */
/**********************************************************/
double fldRotationalEnhancement(const int i, const int n, const double sperp, const double spara)
{
  double krot = 1.0;

  /*** inner barrier with N>144, axially asymmetric (gamma-deformation) */
  if( (i == 0) && (n > 144) ){
    krot = sqrt(8.0*PI*spara) * sperp;
  }
  /*** inner barrier with N<144, outer barrier, axially symmetric */
  else{
    krot = sperp;
    if(i > 0) krot *= 2.0;
  }
  return(krot);
}


/**********************************************************/
/*      Sigma2 for Paralell to the Axis                   */
/**********************************************************/
double fldSigmaParalell(const double b2, const double a, const double t, const double mass)
{
//b2 *= 1.0 / (1.0+exp( (t-2.0) / 0.5 ));
  double sig2 = 6.0/(PI*PI) * 0.24*pow(mass,2.0/3.0) * (1.0-2.0*b2/3.0) *a*t;
  return(sig2);
}


/**********************************************************/
/*      Sigma2 for Perpendicular to the Axis              */
/**********************************************************/
double fldSigmaPerpendicular(const double b2, const double s)
{
//b2 *= 1.0 / (1.0+exp( (t-2.0) / 0.5 ));
  double sig2 = s*s*( 1.0 + b2 * sqrt(5.0/(16.0*PI)) );
  return(sig2);
}


/**********************************************************/
/*      Highest Discrete Transition State                 */
/**********************************************************/
double fldHighestTransitionState(Barrier *b)
{
  const int jcut = 8;  // maximum J-value in the discrete region

  double elmax = 0.0;
  for(int k=0 ; k<b->nband ; k++){
    for(int jj=0 ; jj<=2*jcut ; jj+=2){
      double x = specRotationalBand(jj,&b->kband[k],b->inertia);
      if(x < 0.0) continue;
      if(x > elmax) elmax = x;
    }
  }
  return(elmax);
}


/**********************************************************/
/*      Collective Enhancement Factor in GNASH            */
/**********************************************************/
#ifdef FISDENS_GNASH
double fldDensityEnhancementGNASH(const double ex, const double f, const double pair)
{
  double ux = ex - pair;
  double ef = 1.0;

  if(ux < 0.001) ux = 0.001;

  if(ux < 0.0){
    ef = 2.0;
  }else if(ux < 1.0){
    ef = (f-2.0)*ux + 2.0;
  }
  else{
    ef = f * pow(ux,0.25);
  }

  /*** empirical damping factor, Fermi functional form */
  double df = 1.0 / (1.0 + exp( (ux - 10.0)/2.0 ) );

  return(ef*df);
}
#endif


/**********************************************************/
/*      O. Iwamoto, J.Nucl. Sci. Tech. 44, 687 (2007)     */
/**********************************************************/
#ifdef FISDENS_CCONE
double fldFissionDensityCCONE(const int i, const double j, const double ex, ZAnumber *za, LevelDensity *ldp)
{
  /*** normal level density */
  double mass = (double)(za->getA());
  double eshell = 2.6;
  if(i>0) eshell = 0.6 + 0.1*(za->getZ()-97.0) + 0.04*(za->getN() - 143.0);

  double pfac = 14.0/12.0; // pairing energy scaling factor at the saddle point
  double um   = 2.0;
  double em   = um + ldp->pairing_energy*pfac;
  double u    = ex - ldp->pairing_energy*pfac;
  double astr = (0.04468*mass + 0.2014*pow(mass,2/3.0))*1.02;

  double t,a,sig,dens,sperp,spara,krot;
  double b2f = (i == 0) ? 0.4 : 0.8;

  double ux = u;
  if(ex < em)  ux = um;

  a     = ldShellCorrection(ux,astr,eshell,mass);
  t     = sqrt(ux/a);

  sig   = 0.01389 * pow(mass,(5.0/3.0)) * (1.0+ 0.2*sqrt(5.0/(16.0*PI)))*t;
  sperp = 0.01389 * pow(mass,(5.0/3.0)) * (1.0+ b2f*sqrt(5.0/(16.0*PI)))*t;
  spara = fldSigmaParalell(b2f,mass/8.0,t,mass);

  if(i==0) krot = sqrt(8.0*PI*spara) * sperp;
  else     krot = sperp * 2.0;

  dens  = ldFermiGas(ux,a,ldp->spin_cutoff)*krot;

  if(ex < em){
    double e0   = em - ldp->temperature * log(ldp->temperature * dens);
    dens = ldConstantTemperature(ex,e0,ldp->temperature);
  }

  double sd = (j+0.5)*exp( -(j+0.5)*(j+0.5)/(2.0*sig) )/sig;

  return(dens*sd*0.5);
}
#endif

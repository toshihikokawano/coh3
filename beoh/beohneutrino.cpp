/******************************************************************************/
/*  beohneutrino.cpp                                                          */
/*        calculate electron and neutrino energy spectra                      */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "beoh.h"
#include "etc.h"
#include "terminate.h"

static void   beohAddElectronSpectrum (const double, const double, double *, double *, double *);
static double beohEnergySpectrum (const double);
static double beohFermiFunction (const double);

static int z,a;
static double endp = 0.0; // end-point energy


/***********************************************************/
/*      Electron and Neutrino Spectra                      */
/***********************************************************/
void beohElectronSpectrum(const double q, double *rc, double *rd, Nucleus *n, double *spc1, double *spc2)
{
  double *t = nullptr;
  try{
    t = new double [MAX_ENERGY_BIN];
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("beohElectronSpectrum");
  }

  z = n->za.getZ();
  a = n->za.getA();

  for(int k=0 ; k<MAX_ENERGY_BIN  ; k++) spc1[k] = spc2[k] = 0.0;

  /*** final state in continuum */
  for(int k=0 ; k<n->ncont ; k++){
    if(rc[k] == 0.0) continue;
    endp = q - n->excitation[k];
    beohAddElectronSpectrum(rc[k],n->de,t,spc1,spc2);
  }

  /*** final state in discrete levels */
  for(int i=0 ; i<n->ndisc ; i++){
    if(rd[i] == 0.0) continue;
    endp = q - n->lev[i].energy;
    beohAddElectronSpectrum(rd[i],n->de,t,spc1,spc2);
  }

  /*** convert into per-MeV */
  double x = 1.0/n->de;
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    spc1[k] *= x;
    spc2[k] *= x;
  }

  delete [] t;
}


/***********************************************************/
/*      Partial Spectra for Given Endpoint                 */
/***********************************************************/
void beohAddElectronSpectrum(const double r, const double de, double *t, double *spc1, double *spc2)
{
  /*** clear work array */
  for(int k=0 ; k<MAX_ENERGY_BIN  ; k++) t[k] = 0.0;

  /*** normalization constant */
  double s  = 0.0;
  int    km = 0;
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    double e0 = (k>0)  ? ((double)k-0.5)*de : 0;
    double e1 = ((double)k+0.5)*de;
    double e  =  (e0+e1)/2.0;
    s += (t[k] = beohEnergySpectrum(e));
    if(e>endp){
      km = k;
      break;
    }
  }

  /*** weighted sum of electron and neutrino spectra */
  if(s>0.0) s = r/s;
  for(int k=0 ; k<=km ; k++){
    spc1[   k] += t[k]*s;
    spc2[km-k] += t[k]*s;
  }
}


/***********************************************************/
/*      Fermi Function                                     */
/*      Y. Isozumi, F. Ishizuka,                           */
/*      H. Miyatake, T. Kato, M. Tosaki,                   */
/*      Nihon Houshasen Anzen Kanri, 4, page 71 (2005)     */
/***********************************************************/
double beohEnergySpectrum(const double e)
{
  double s = 0.0;
  if(e > endp) return(s);

  /*** this is unnormalized spectrum */
  double f = beohFermiFunction(e);

  if(e == 0.0){
    s = endp * endp * f;
  }
  else{
    s = sqrt(e*(e+2*EELECTRON)) * (endp-e)*(endp-e) * (e+EELECTRON) * f;
  }
  return(s);
}


double beohFermiFunction(const double e)
{
  double f = 0.0;
  double g = sqrt(1.0 - ALPHA*ALPHA*z*z);

  if(e == 0.0){
    f = PI2 * pow(ALPHA*z,2*g-1.0) * pow(EELECTRON,2*g);
  }
  else{
    double radius = ALPHA/2.0 * pow((double)a,1.0/3.0);
    double w  = e + EELECTRON;
    double y  = ALPHA * z * w / sqrt( e*(e+2*EELECTRON) );
    double p  = sqrt( e*(e+2*EELECTRON) ) / VLIGHT / HBAR;
    double r  = y/(1.0+g);
    double r2 = r*r;

    double v  = 2*(1-g) * ( 1.0 - log(y) + r/(24*y*(1+r2)) + 0.5*log(r2/(1+r2)) + r*atan(r) );
    double x1 = PI_4 * y*y*y/(g*g+y*y) * (4.0+1/r2) * 2.0/(exp(PI2*r)-exp(-PI2*r)) * exp(v);
    double x2 = gam(2*g+1.0); x2= x2*x2;

    f = 2.0 * (1+g) * pow(2*p*radius,2*g-2.0) * exp(PI*y) * x1/x2;
  }

  return(f);
}


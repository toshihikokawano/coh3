/******************************************************************************/
/*  omsetform.cpp                                                             */
/*        calculate radial part of optical potentials                         */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "optical.h"

static inline double   omWoodsSaxon             (const double, const double, const double);
static inline double   omSurfaceWoodsSaxon      (const double, const double, const double);
static inline double   omGaussian               (const double, const double, const double);
static inline double   omThomas                 (const double, const double, const double);
static inline double   omFirstDerivSurfaceWS    (const double, const double, const double);
static inline double   omFirstDerivGaussian     (const double, const double, const double);
static inline double   omSecondDerivSurfaceWS   (const double, const double, const double);
//static inline double   omSecondDerivGaussian    (const double, const double, const double);

static inline double   omFirstDerivWS           (const double, const double, const double);
static inline double   omSecondDerivWS          (const double, const double, const double);
static inline double   omThirdDerivWS           (const double, const double, const double);


/**********************************************************/
/*     Spherical Optical Potential Radial Form            */
/**********************************************************/
int omPotentialForm(const int zzprod, Optical *omp, LevelData *ldt, Potential *pot)
{
  double r, x[6], eps = log(CRIT_MATCHING);
  int    n;
  unsigned int potform = omp->getPotForm(); 

  /*** set lowest energy = 1 keV to avoid long integration */
  double energy =  (ldt->energy < 0.001) ? 0.001 : ldt->energy;

  x[0] = x[1] = x[2] = x[3] = x[4] = x[5] = 0.0;
  if(potform & 0x0010) x[0]=omp->R0+omp->a0   *(log(    omp->volume.real()  /energy)-eps);
  if(potform & 0x0020) x[1]=omp->R0s+omp->a0s *(log(4.0*fabs(omp->surface.real())/energy)-eps);
  if(potform & 0x0001) x[2]=omp->Rv+omp->av   *(log(    omp->volume.imag()  /energy)-eps);
  if(potform & 0x0002) x[3]=omp->Rs+omp->as   *(log(4.0*omp->surface.imag() /energy)-eps);
  if(potform & 0x0004) x[3]=omp->Rs+omp->as   *sqrt(log(4.0*omp->surface.imag() /energy)-eps);

  if(omp->spin_orbit.real() != 0.0) x[4]=omp->Rvso+omp->avso*(log(     omp->spin_orbit.real() *CSPO
                                                        /(omp->avso*omp->Rvso*energy))-eps);
  if(omp->spin_orbit.imag() != 0.0) x[5]=omp->Rwso+omp->awso*(log(fabs(omp->spin_orbit.imag())*CSPO
                                                        /(omp->awso*omp->Rwso*energy))-eps);

  pot->rad_match = x[0];
  for(int i=1 ; i<6 ; i++) if(x[i] > pot->rad_match) pot->rad_match=x[i];

  n = (int)(pot->rad_match/pot->width)+1;
  r = pot->width*n; if(r < pot->rad_match) n++;
  n += 3;

  if(n >= MAX_POINTS){
//  stdd::cerr << "integration number too large " << n << "  Ecms: " << energy << std::endl;
    n = MAX_POINTS-1;
  }

  pot->setMatching(n);
  pot->setRho(fabs(ldt->wave_number));

  omPotentialFixedLength(zzprod,omp,ldt->reduced_mass,pot);

  return(pot->n_match);
}


/**********************************************************/
/*     Radial Form for Fixed Matching Radius              */
/**********************************************************/
void omPotentialFixedLength(const int zzprod, Optical *omp, const double mu, Potential *pot)
{
  double r;
  double c1 = 2.0*mu*AMUNIT/VLIGHTSQ/HBARSQ;
  double c2 = c1*CSPO;
  double c3 = c1*zzprod*PERMITTIV*COULOMBSQ;

  for(int i=0 ; i<=pot->n_match ; i++){
    r = (double)i*pot->width;
    pot->r2inv[i]      =  (i == 0) ? 0.0 : 1.0 / (r*r);
    pot->mean_field[i] = c1 * omPotentialRadialMean(r,omp);
    pot->spin_orbit[i] = c2 * omPotentialRadialSpo(r,omp);
    pot->coulomb[i]    = (zzprod != 0) ? c3 * omPotentialRadialCoulomb(r,omp) : 0.0;

    /*** add Coulomb potential to the optical potential */
    pot->mean_field[i] -= pot->coulomb[i];

    // std::cout << std::setw(4) << i;
    // std::cout << std::setw(12) << i * pot->width;
    // std::cout << std::setw(12) << pot->r2inv[i];
    // std::cout << std::setw(12) << pot->mean_field[i].real() / c1;
    // std::cout << std::setw(12) << pot->mean_field[i].imag() / c1;
    // std::cout << std::setw(12) << pot->spin_orbit[i].real() / c2;
    // std::cout << std::setw(12) << pot->spin_orbit[i].imag() / c2;
    // if(zzprod != 0) std::cout << std::setw(12) << pot->coulomb[i] / c3;
    // std::cout << std::endl;
  }
}


/**********************************************************/
/*     spherical Coulomb potential                        */
/**********************************************************/
double omPotentialRadialCoulomb(const double r, Optical *omp)
{
  return( (r <= omp->Rc) ? (3.0-r*r/(omp->Rc*omp->Rc))/(2.0*omp->Rc) : 1.0/r );
}


/**********************************************************/
/*     Store one-dimensional potential radial form        */
/**********************************************************/
std::complex<double> omPotentialRadialMean(const double r, Optical *omp)
{
  double pr = 0.0, pi = 0.0;

  unsigned int real_pot = omp->getRealForm();
  unsigned int imag_pot = omp->getImagForm();

  switch(real_pot){
  case 1 : pr = omp->volume.real() *        omWoodsSaxon(r,omp->R0 ,omp->a0 );  break;
  case 2 : pr = omp->surface.real()* omSurfaceWoodsSaxon(r,omp->R0s,omp->a0s);  break;
  case 3 : pr = omp->volume.real() *        omWoodsSaxon(r,omp->R0 ,omp->a0 )
               +omp->surface.real()* omSurfaceWoodsSaxon(r,omp->R0s,omp->a0s);  break;
  default: pr = 0.0;
  }

  switch(imag_pot){
  case 1 : pi = omp->volume.imag() *        omWoodsSaxon(r,omp->Rv,omp->av);    break;
  case 2 : pi = omp->surface.imag()* omSurfaceWoodsSaxon(r,omp->Rs,omp->as);    break;
  case 3 : pi = omp->volume.imag() *        omWoodsSaxon(r,omp->Rv,omp->av)
               +omp->surface.imag()* omSurfaceWoodsSaxon(r,omp->Rs,omp->as);    break;
  case 4 : pi = omp->surface.imag()*          omGaussian(r,omp->Rs,omp->as);    break;
  case 5 : pi = omp->volume.imag() *        omWoodsSaxon(r,omp->Rv,omp->av)
               +omp->surface.imag()*          omGaussian(r,omp->Rs,omp->as);    break;
  default: pi = 0.0;
  }

  return(std::complex<double>(pr,pi));
}


std::complex<double> omPotentialRadialSpo(const double r, Optical *omp)
{
  double pr = 0.0, pi = 0.0;
  if(omp->spin_orbit.real() != 0.0) pr = omp->spin_orbit.real()*omThomas(r,omp->Rvso,omp->avso);
  if(omp->spin_orbit.imag() != 0.0) pi = omp->spin_orbit.imag()*omThomas(r,omp->Rwso,omp->awso);

  return(std::complex<double>(pr,pi));
}


/**********************************************************/
/*     DWBA Collective form factor                        */
/**********************************************************/
int dwCollectiveForm(const double beta, Optical *omp, Potential *pot)
{
  for(int i=0 ; i<=pot->n_match ; i++){
    double r = (double)i*(pot->width);

    pot->mean_field[i] = beta * omPotentialRadialFirstDeriv(r,omp);

    if(fabs(pot->mean_field[i].real()) < 1.0e-06){
      if(i < pot->n_match) pot->n_match = i;
      break;
    }
  }
  return(pot->n_match);
}


/**********************************************************/
/*     Collective form factor for coupled-channels        */
/**********************************************************/
std::complex<double> omPotentialRadialFirstDeriv(const double r, Optical *omp)
{
  double pr = 0.0, pi = 0.0;

  switch(omp->getRealForm()){
  case 1 : pr = omp->volume.real()  * omp->R0
              * omFirstDerivWS(r,omp->R0,omp->a0);          break;
  case 2 : pr = omp->surface.real() * omp->R0s
              * omFirstDerivSurfaceWS(r,omp->R0s,omp->a0s); break;
  case 3 : pr = omp->volume.real()  * omp->R0
              * omFirstDerivWS(r,omp->R0,omp->a0)
              + omp->surface.real() * omp->R0s
              * omFirstDerivSurfaceWS(r,omp->R0s,omp->a0s); break;
  default: pr = 0.0; break;
  }

  switch(omp->getImagForm()){
  case 1 : pi = omp->volume.imag()  * omp->Rv
              * omFirstDerivWS(r,omp->Rv,omp->av);          break;
  case 2 : pi = omp->surface.imag() * omp->Rs
              * omFirstDerivSurfaceWS(r,omp->Rs,omp->as);   break;
  case 3 : pi = omp->volume.imag()  * omp->Rv
              * omFirstDerivWS(r,omp->Rv,omp->av)
              + omp->surface.imag() * omp->Rs
              * omFirstDerivSurfaceWS(r,omp->Rs,omp->as);   break;
  case 4 : pi = omp->surface.imag() * omp->Rs
              * omFirstDerivGaussian(r,omp->Rs,omp->as);    break;
  case 5 : pi = omp->volume.imag()  * omp->Rv
              * omFirstDerivWS(r,omp->Rv,omp->av)
              + omp->surface.imag() * omp->Rs
              * omFirstDerivGaussian(r,omp->Rs,omp->as);    break;
  default: pi = 0.0; break;
  }

  return( std::complex<double>(-pr,-pi) );
}


/**********************************************************/
/*     2nd-Order Collective form factor                   */
/**********************************************************/
std::complex<double> omPotentialRadialSecondDeriv(const double r, Optical *omp)
{
  double pr = 0.0, pi = 0.0;

  switch(omp->getRealForm()){
  case 1 : pr = omp->volume.real()
              * 0.5 * omp->R0*omp->R0 
              * omSecondDerivWS(r,omp->R0,omp->a0);
    break;
  case 2 : pr = omp->surface.real()
              * 0.5 * omp->R0s*omp->R0s
              * omSecondDerivSurfaceWS(r,omp->R0s,omp->a0s);
    break;
  case 3 : pr = omp->volume.real()
              * 0.5 * omp->R0*omp->R0
              * omSecondDerivWS(r,omp->R0,omp->a0)
              + omp->surface.real()
              * 0.5 * omp->R0s*omp->R0s
              * omSecondDerivSurfaceWS(r,omp->R0s,omp->a0s);
    break;
  default: pr = 0.0;
    break;
  }

  switch(omp->getImagForm()){
  case 1 : pi = omp->volume.imag()
              * 0.5 * omp->Rv*omp->Rv
              * omSecondDerivWS(r,omp->Rv,omp->av);
    break;
  case 2 : pi = omp->surface.imag()
              * 0.5 * omp->Rs*omp->Rs
              * omSecondDerivSurfaceWS(r,omp->Rs,omp->as);
    break;
  case 3 : pi = omp->volume.imag()
              * 0.5 * omp->Rv*omp->Rv
              * omSecondDerivWS(r,omp->Rv,omp->av)
              + omp->surface.imag()
              * 0.5 * omp->Rs*omp->Rs
              * omSecondDerivSurfaceWS(r,omp->Rs,omp->as);
    break;
  default: pi = 0.0;
  }

  return(std::complex<double>(pr,pi));
}


/**********************************************************/
/*     Various potential radial forms                     */
/**********************************************************/
/*** Woods - Saxon */
inline double omWoodsSaxon(const double x, const double r, const double a)
{ return(  1/( 1+exp((x-r)/a) )  ); }

/*** Surface Woods - Saxon, scaled for 1 at r=R */
inline double omSurfaceWoodsSaxon(const double x, const double r, const double a)
{ double y = exp((x-r)/a);
  double z = 1.0 + y;
  return( 4*y/(z*z) ); }

/*** Gaussian shape */
inline double omGaussian(const double x, const double r, const double a)
{ double y = (x-r)*(x-r)/(a*a);
  return(  exp(-y)  ); }

/*** Thomas shape for spin-orbit force */
inline double omThomas(const double x, const double r, const double a)
{ if(x == 0.0) return(0.0);
  double y = exp((x-r)/a);
  double z = 1.0 + y;
  return(  y/(x*a*z*z)  ); }

/*** derivative Surface Woods -Saxon = omSecondDerivWS x (-4a) */
inline double omFirstDerivSurfaceWS(const double x, const double r, const double a)
{ return( -4.0*a * omSecondDerivWS(x,r,a) ); }

/*** second derivative Surface Woods -Saxon = omThirdDerivWS x (-4a) */
inline double omSecondDerivSurfaceWS(const double x, const double r, const double a)
{ return( -4.0*a * omThirdDerivWS(x,r,a) ); }

/*** first-order derivative Gaussian */
inline double omFirstDerivGaussian(const double x, const double r, const double a)
{ double y = (x-r)*(x-r)/(a*a);
  return( -2*(x-r)/(a*a)*exp(-y)  ); }

/*** second-order derivative Gaussian (not used) */
//inline double omSecondDerivGaussian(double x, double r, double a)
//{ double y = (x-r)*(x-r)/(a*a);
//  return(  (4*(x-r)*(x-r)/(a*a) - 2)/(a*a) * exp(-y)  ); }

/*** first-order derivative Woods - Saxon */
inline double omFirstDerivWS(const double x, const double r, const double a)
{ double y = exp((x-r)/a);
  double z = 1.0 + y;
  return( - 1.0/a * y/(z*z) ); }

/*** second-order derivative Woods -Saxon */
inline double omSecondDerivWS(const double x, const double r, const double a)
{ double y = exp((x-r)/a);
  double z = 1.0 + y;
  return( - 1.0/(a*a) * y*(1.0 - y)/(z*z*z) ); }

/*** third-order derivative Woods -Saxon */
inline double omThirdDerivWS(const double x, const double r, const double a)
{  double y = exp((x-r)/a);
   double z = 1+y;
   return( - 1.0/(a*a*a) * y*(1.0 - 4*y + y*y)/(z*z*z*z) ); }


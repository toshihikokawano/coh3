/******************************************************************************/
/*  omsetparm.cpp                                                             */
/*        setting up parameters in optical model calculations                 */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "structur.h"
#include "optical.h"
#include "parameter.h"

static CoulombData omSetCoulomb (const int, const double, const double);
static void omAdjustOmp (const int, const int, Optical *);

/**********************************************************/
/*      Setup Energy Dependent Parameters for OM Calc     */
/*      cms energy, wave number, Coulomb parameter of     */
/*      ejectile, excitation of residual nucleus          */
/**********************************************************/
void omSetEnergy(const double e, const int zz, const double mu, LevelData *lev)
{
  lev->energy       = e;
  lev->reduced_mass = mu;
  lev->wavesq       = lev->energy * mu *2.0*AMUNIT/VLIGHTSQ/HBARSQ;
  lev->wave_number  = sqrt(lev->wavesq);
  
  if(zz == 0){
    lev[0].coulomb.eta = lev->coulomb.sigma0 = 0.0;
  }
  else{
    lev->coulomb = omSetCoulomb(zz,mu,lev->energy);
  }
}


/**********************************************************/
/*      Setup Energy Dependent Parameters for CC Calc     */
/**********************************************************/
int ccSetEnergy(const double e, const int zz, const double mu, Collective *col)
{
  int open_channel = 0;
  
  for(int k=0 ; k<col->nlevel ; k++){

    col->lev[k].reduced_mass = mu;
    col->lev[k].energy = e - col->lev[k].excitation 
                           + col->lev[col->target_index].excitation;
    col->lev[k].wavesq = col->lev[k].energy*mu*2.0*AMUNIT/VLIGHTSQ/HBARSQ;

    /* Open Channel */
    if(col->lev[k].energy > 0.0){
      col->lev[k].wave_number = sqrt(col->lev[k].wavesq);

      if(zz == 0){
        col->lev[k].coulomb.eta = col->lev[k].coulomb.sigma0 = 0.0;
      }else{
        col->lev[k].coulomb = omSetCoulomb(zz,mu,col->lev[k].energy);
      }
      open_channel++;
    }

    /* Closed Channel, the wave number is set to negative value */
    else{
      col->lev[k].wave_number = -sqrt(-col->lev[k].wavesq);

      if(zz == 0){
        col->lev[k].coulomb.eta = col->lev[k].coulomb.sigma0 = 0.0;
      }else{
        col->lev[k].coulomb = omSetCoulomb(zz,mu,-col->lev[k].energy);
        col->lev[k].coulomb.sigma0 = 0.0;
      }
    }
  }

  return(open_channel);
}


/**********************************************************/
/*     Coulomb Parameters                                 */
/**********************************************************/
CoulombData omSetCoulomb(const int zz, const double mu, const double e)
{
  double y1,y2,y3,y4,eta,sig0;

  eta = y1 = PERMITTIV*COULOMBSQ * zz * sqrt(AMUNIT*mu/(2.0*e))/VLIGHT/HBAR;
  y2 = y1*y1;
  y3 = 16.0 + y2;
  y4 = y3*y3;
  sig0 = -y1+y1*log(y3)/2.+3.5*atan(y1/4.)-(atan(y1)+atan(y1/2.)+atan(y1/3.))
         -y1*(1.+(y2-48.)/(30.*y4)+(y2*y2-160.*y2+1280.)/(105.*y4*y4))
         /(12.*y3);
/*
  double w = 0.0;
  if(l > 0){
    for(int n=1 ; n<=l ; n++) w += atan(eta / (double)n);
  }
  sig0 = w;
*/
  CoulombData c = {eta,sig0};
  return(c);
}


/**********************************************************/
/*     Setup Optical Potential Parameters                 */
/**********************************************************/
unsigned int omSetOmp(const int ompidx, const double e, const int z0, const int a0, const int z1, const int a1, Optical *omp)
{
  /***  Clear optical potential */
  omp->init();

  /***  Retrieve potential parameters from library */
  unsigned int potfm = omp_library(ompidx, z0, a0, z1, a1, e, omp);

  /***  Energy dependent depths */
  omp->setDepth(e);

  if(omp->volume.real() <= 0.0){
    potfm = (potfm & 0x009f);
    omp->volume.real(0.0);
  }
  if(omp->volume.imag() <= 0.0){
    potfm = (potfm & 0x00fe);
    omp->volume.imag(0.0);
  }
  if(omp->surface.imag() <= 0.0){
    potfm = (potfm & 0x00f9);
    omp->surface.imag(0.0);
  }
  if(omp->spin_orbit.real() <= 0.0){
    omp->spin_orbit.real(0.0);
  }

  /***  Parameter adjustment */
  omAdjustOmp(z1,a1,omp);

  /***  Reduced parameters */
  omp->setReducedRadius(pow((double)a0,1.0/3.0));

  /***  Store potential form and ID */
  potfm = ((unsigned int)ompidx << 8) | potfm;
  omp->setIndex(potfm);

  return( potfm );
}


/**********************************************************/
/*     Adjustment Optical Potential Parameters            */
/**********************************************************/
void omAdjustOmp(const int z, const int a, Optical *omp)
{
  Particle p = unknown;

  if(     z == 0 && a == 1) p = neutron;
  else if(z == 1 && a == 1) p = proton;
  else if(z == 1 && a == 2) p = deuteron;
  else if(z == 1 && a == 3) p = triton;
  else if(z == 2 && a == 3) p = helion;
  else if(z == 2 && a == 4) p = alpha;

  double vs = parmGetFactor(parmOV,p);
  double ws = parmGetFactor(parmOW,p);

  double rvs = parmGetFactor(parmORV,p);
  double rws = parmGetFactor(parmORW,p);

  double avs = parmGetFactor(parmOAV,p);
  double aws = parmGetFactor(parmOAW,p);

  double rc  = parmGetFactor(parmOC,p);

  /***  Scale optical potential depth */
  if(vs != 1.0) omp->volume.real( vs * omp->volume.real() );
  if(ws != 1.0) omp->surface.imag( ws * omp->surface.imag() );

  /***  Scale optical potential radius */
  if(rvs != 1.0) omp->r0 *= rvs;
  if(rws != 1.0) omp->rs *= rws;

  /***  Scale optical potential diffuseness */
  if(avs != 1.0) omp->a0 *= avs;
  if(aws != 1.0) omp->as *= aws;

  /***  Scale Coulomb radius */
  if(rc != 1.0)  omp->rc *= rc;
}

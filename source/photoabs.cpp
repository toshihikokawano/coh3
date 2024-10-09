/******************************************************************************/
/*  photoabs.cpp                                                              */
/*        photon absorption cross section                                     */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "global.h"
#include "terminate.h"

static double photoQD    (const double, const int, const int);

//  Fraction of QD, which goes into Pre-equilibrium calc
static double qdratio = 0.0;
double photoQDratio()
{ return qdratio; }


/**********************************************************/
/*     Gamma-ray induced entrance channel calculation     */
/*         GDR model and quasi deuteron break-Up          */
/**********************************************************/
double photoAbsorption(const int targlev, const int targid, const double eg, const double coef, Nucleus *n, Transmission *tin)
{
  if(eg <= 0.0){
    message << "incident gamma-ray energy " << eg << " negative or zero";
    cohTerminateCode("photoAbsorption");
  }

  for(int j=0 ; j<3*MAX_J ; j++) tin->tran[j] = 0.0;

  /*** calculate photon transmission coefficients */
  double tg[MAX_MULTIPOL];
  statStoreDiscreteGammaTransmission(eg,tg,n);

  /*** the order of coefficients is special for the photon case,
       (E0,M0,-), (E1,M1,-), (E2,M2,-) ... */
  tin->tran[ 3] = tg[E1];
  tin->tran[ 4] = tg[M1];
  tin->tran[ 6] = tg[E2];
  tin->tran[ 7] = tg[M2];
  tin->tran[ 9] = tg[E3];
  tin->tran[10] = 0.0;
  tin->lmax = 3;

  /*** GDR absorption, calculate SigR from transmissions */
  double absgdr = statStoreInitPopulationPhoton(targlev,targid,coef,tin);

  /*** quasi-deutron absorption */
  double absqd  = photoQD(eg,n->za.getZ(),n->za.getA());

  /*** add QD to GDR, and renormalize, is it OK ? */
  qdratio = (absgdr > 0.0) ? absqd/absgdr : 0.0;
  double f = 1.0 + qdratio;
  qdratio = absqd / (absqd + absgdr);  if(qdratio > 1.0) qdratio = 1.0;

  /*** set max J = I + Lmax */
  n->jmax = n->lev[targlev].spin + tin->lmax;

  return(absgdr*f);
}


/**********************************************************/
/*      Quasi Deuteron Absorption                         */
/*      Chadwick et al., PRC44, 814 (1991)                */
/**********************************************************/
double photoQD(const double eg, const int z, const int a)
{
  const double deuteron_bind = 2.22452;
  const double levinger_const = 6.5;
  double sigd = 0.0, pbf = 0.0;

  double e2 = eg * eg;
  double e3 = e2 * eg;
  double e4 = e3 * eg;

  if(eg > deuteron_bind) sigd = 61.2*( pow((eg-2.224),1.5) )/e3;

  if(eg < 20.0)       pbf = exp(-73.3   /eg);
  else if(eg > 140.0) pbf = exp(-24.2348/eg);
  else                pbf = 0.083714 - 0.0098343*eg + 4.1222e-04*e2 - 3.4762e-06*e3 + 9.3537e-09*e4;

  double x = (double)(a-z)*(double)z/(double)a;

  return(levinger_const * x * sigd * pbf);
}


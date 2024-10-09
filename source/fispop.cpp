/******************************************************************************/
/*  fispop.cpp                                                                */
/*        fission of compound nucleus                                         */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "levden.h"
#include "nucleus.h"
#include "parameter.h"
#include "terminate.h"

static double specFissionTransmission (const int, const int, const int, const double, double *, Barrier *, Nucleus *);
static double specFissionDoubleHumpModel (const int, const int, const int, const double, double *, Nucleus *);
static inline double specHillWheeler (const double, const double, const double);

static double tjsave = 0.0, dfsave = 0.0;


/**********************************************************/
/*      Fission Transmission Coefficient                  */
/**********************************************************/
double specTransitionFission
( const Statcalcmode  mode,        // calc. control for sum T/Moldauer/HF
  const int           k0,          // index of parent excited state
  const int           p0,          // parent state parity (-1/+1)
  const int           j0,          // parent spate spin (doubled)
  const double        q,           // pop(parent)/Sigma T, if mode=HF
  Nucleus             *n)          // nucleus information
{
  double  ch[MAX_HUMP];
  Fission *fb = n->fission;

  /*** excitation energy of compound nucleus */
  double ex = n->excitation[k0];

  /*** if well parameter is given, check the first barrier */
  bool potwell = (fb->potwell[0].height != 0.0 || fb->potwell[0].curvature != 0.0) ? true : false;

  double tj = 0.0;
  double df = 0.0;
  for(int i=0 ; i<MAX_HUMP ; i++) ch[i] = 0.0;

  /*** calculate Tj once, then save the value */
  if(mode == sumall){

    if(potwell) tj = specFissionPotentialModel(k0,p0,j0,ex,ch,n);
    else        tj = specFissionDoubleHumpModel(k0,p0,j0,ex,ch,n);

    /*** fission transmission adjustment */
    tj *= parmGetFactor(parmTJ,n->za,fission);

    if(fb->fisenhance.energy > 0.0){
      double w2 = fb->fisenhance.width * fb->fisenhance.width;
      double e2 = ex - fb->barrier[0].height - fb->fisenhance.energy; e2 = e2 * e2;
      tj *= 1.0 + fb->fisenhance.peak * w2 / (e2 + w2);
    }

    /*** average number of channels - degree-of-freedom */
    if(potwell) df = ch[0];
    else{
      for(int i=0 ; i<MAX_HUMP ; i++) df += ch[i];
      df /= fb->hump;
    }

    tjsave = tj;
    dfsave = df;
  }
  else{
    tj = tjsave;
    df = dfsave;
  }

  if(mode == wfactor) statWidthFluctuation(tj,1.0,df);
  else if(mode == hauser || mode == fluctuation){
    double wfc = (mode==fluctuation) ? statMoldauer(-1,0,0,0,tj,df) : 1.0;

    /*** next-channel index for gamma contains CN itself */
    int idx = n->cdt[gammaray].next;
    tj *= q*wfc;

    crx.prod[idx].fiss += tj;
  }

  return(tj);
}


/**********************************************************/
/*      Double Hump Fission Barrier Model                 */
/**********************************************************/
double  specFissionDoubleHumpModel(const int k0, const int p0, const int j0, const double ex, double *ch, Nucleus *n)
{
  double tc[2],tf[MAX_HUMP];
  Fission *fb = n->fission;

  /*** Tf for the all barriers */
  for(int i=0 ; i<MAX_HUMP ; i++){
    if(fb->barrier[i].height == 0.0) continue;

    /*** fission transmission and number of effective fission channels */
    ch[i] = specFissionTransmission(k0,p0,j0,ex,tc,&fb->barrier[i],n);

    /*** add continuum and discrete transitions */
    tf[i] = tc[0] + tc[1];
  }

  /*** effective transmission coefficient */
  double tj = tf[0];
  if(fb->hump >= 2){
    tj = ((tf[0]+tf[1]) > 0.0) ? tf[0]*tf[1] / (tf[0] + tf[1]) : 0.0;

    if(fb->hump == 3){
      double tj2 = tj;
      tj = ((tj2+tf[2]) > 0.0) ? tj2*tf[2] / (tj2 + tf[2]) : 0.0;
    }
  }

  return(tj);
}


/**********************************************************/
/*      Fission Transmission for Single Hump              */
/**********************************************************/
double  specFissionTransmission(const int k0, const int p0, const int j0, const double ex, double *tf, Barrier *b, Nucleus *n)
{
  /*** index for spin J */
  int idx = (j0-(int)(2.0*halfint(n->lev[0].spin)))/2;

  tf[0] = tf[1] = 0.0;

  double nch = 0.0; // effective number of fission channels

  /*** for all K-bands on top of this barrier
       determine critical energy Ec for discrete transition */
  for(int k=0 ; k<b->nband ; k++){
    if(p0 != b->kband[k].parity) continue;
    if(j0 <  b->kband[k].k2)     continue;
 
    /*** energy of state J0, measured from the top of barrier */
    double el = specRotationalBand(j0,&b->kband[k],b->inertia);
    if(el < 0.0) continue;
    if(el > b->elmax) continue;

    /*** sum transmission coefficients for all possible J-pi states */
    tf[0] += specHillWheeler(b->height,b->curvature,ex-el);
    nch += 1.0;
  }

  /*** continuum states above E(critical) */
  for(int k1=k0 ; k1<n->ntotal ; k1++){
    double el = n->excitation[k1];
    if(el < b->elmax) break;

    /*** number of levels in a bin */
    double rho = ((p0 > 0) ? b->density[k1][idx].even
                           : b->density[k1][idx].odd ) * n->de;

    tf[1] += specHillWheeler(b->height,b->curvature,ex-el) * rho;
    nch += rho;
  }

  return(nch);
}


/**********************************************************/
/*      Hill-Wheeler Double Parabolic Expression          */
/**********************************************************/
double  specHillWheeler(const double barrier_height, const double curvature, const double energy)
{
  double trans = 0.0;
  if(curvature > 0.0){
    trans = 1.0 /(1.0+exp(PI2*(barrier_height-energy)/curvature));
  }

  return(trans);
}


/**********************************************************/
/*      Build Rotational Band Spectrum of K-Band          */
/**********************************************************/
double  specRotationalBand(const int jj, KBand *kb, double x)
{
  double e = -1.0;
  double j = (double)jj/2.0;
  double k = (double)kb->k2/2.0;

  if(j < k) return(e);

  /*** K=0+ case, no J-odd state */
  if( (kb->k2 == 0) && (kb->parity > 0) && ((jj%4) != 0) ) return(e);

  /*** K=0- case, no J-even state */
  if( (kb->k2 == 0) && (kb->parity < 0) && ((jj%4) == 0) ) return(e);

  e = (j*(j+1.0)-k*(k+1.0))*x + kb->excitation;
  return (e);
}



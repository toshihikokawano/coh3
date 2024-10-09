/******************************************************************************/
/*  FRDMpairing.cpp                                                           */
/*        Pairing Energy by Lipkin-Nogami BCS Model                           */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "FRDM.h"


/**********************************************************/
/*      Pairing Correlation and QuasiParticle Energy, Epc */
/**********************************************************/
double FRDMPairingEnergy(LNParameter *bcs, SPEnergy *spe)
{
  const int np = bcs->np;

  /*** determine sharp distribution factor, nk */
  int *nk;
  nk = new int [bcs->n2 + 1];

  if((np%2) != 0){ // odd target
    int ns = 1;
    for(int i=0 ; i<=bcs->n2 ; i++){
      if(ns == np){     nk[i] = 1; }
      else if(ns < np){ nk[i] = 2; }
      else{             nk[i] = 0; }
      ns += 2;
    }
  }else{           // even target
    int ns = 0;
    for(int i=0 ; i<=bcs->n2 ; i++){
      nk[i] = (ns < np) ? 2 : 0;
      ns += 2;
    }
  }

  double eodd = 0.0;
  /*** odd nucleon particle energy */
  if((np%2) != 0){
    double e = spe->energy[bcs->nf] + (4*bcs->lambda2 - bcs->strength)*spe->v2[bcs->nf];
    double x = e - bcs->chemical_potential;
    eodd = sqrt(x * x + bcs->pairing_gap * bcs->pairing_gap) + bcs->lambda2;
  }

  double e1 = bcs->pairing_gap * bcs->pairing_gap / bcs->strength;
  double e0 = 0.0, e2 = 0.0, e3 = 0.0;
  for(int i=bcs->n1 ; i<=bcs->n2 ; i++){
    e0 += (2.0*spe->v2[i] - nk[i]) * spe->energy[i];
    e2 += 2.0*spe->v2[i]*spe->v2[i] - nk[i];
    e3 += spe->v2[i] * (1.0 - spe->v2[i]);
  }
  double ep = e0 - e1 - e2 * 0.5 * bcs->strength - e3 * 4.0 * bcs->lambda2 + eodd;

  delete [] nk;

  return(ep);
}


/**********************************************************/
/*      Average Pairing Correlation Energy, E-bar(p.c.)   */
/**********************************************************/
double FRDMAveragePairingEnergy(const double rhobar, LNParameter *bcs)
{
  double deltag = bcs->pairing_eff;

  double g,x1,x2,y1,y2,z1,z2,z3,z4,z5;

  g  = bcs->strength;

  double xf = (int)(bcs->np/2);
  if((bcs->np % 2) == 0) xf -= 0.5;

  y1 = (bcs->n1 - xf)/rhobar;
  y2 = (bcs->n2 - xf)/rhobar;

  x1 = y1*y1 + deltag*deltag;
  x2 = y2*y2 + deltag*deltag;

  z1 = 2.0/(g * rhobar);
  z2 = log( sqrt(x2) / sqrt(x1) );
  z3 = rhobar * deltag / 4.0;
  z4 = atan(y2/deltag) - atan(y1/deltag);
  z5 = deltag * (y2/x2 - y1/x1);

  double a  = z3*z3 * (z1*z1 - z2*z2);
  double b  = z3*z3 * z4*z4;
  double c  = z3/8.0 * (z5 + z4);

  double lambda2 = 0.25*g * (a-c)/(b-c);

  double ep = 0.5*rhobar * ( (y2-g)*(y2-sqrt(x2)) + (y1-g)*(y1+sqrt(x1)) )
            + (0.25*g - lambda2) * rhobar * deltag * z4
            + bcs->pairing_mac;

  return(ep);
}


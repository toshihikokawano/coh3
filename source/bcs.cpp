/******************************************************************************/
/*  bcs.cpp                                                                   */
/*        BCS calculation, pairing energy systematics, solve gap equation     */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "optical.h"
#include "dsd.h"
#include "terminate.h"

static const int MAX_BCSITER = 100;

/**********************************************************/
/*     BCS Gap Energy Systematics from KTUY04             */
/**********************************************************/
double bcsEnergyGap(const Particle p, const int z, const int a)
{
  double delta, deltan, deltaz;
  const double
         c1n = 0.45, c2n = 1.60, c1p = 0.40, c2p = 1.58, alpha = 0.25,
         r0 = 1.04, rad = 0.1268, ri = 3.0, rq = 1.4996, crq = 1.0, aad = 2.0,
         coo = 0.1070, coo1 = 20.14;

  /*** KTUY04 pairing energy */
  int n    = a-z;
  double eps  = (n-z)/a;
  double pisq =  PI*PI;

  double rn = r0*pow((a+aad),(1/3.0))*(1+crq*eps*eps) + rad + ri*eps + rq*eps*eps;
  double rp = r0*pow((a+aad),(1/3.0))*(1+crq*eps*eps) + rad - ri*eps + rq*eps*eps;

  double vn = 4.0*PI/3.0 * rn*rn*rn;  vn = pow(vn,2/3.0);
  double vp = 4.0*PI/3.0 * rp*rp*rp;  vp = pow(vp,2/3.0);

  double dn = HBARSQ*VLIGHTSQ*pisq / (2*MNEUTRON*AMUNIT*pow((3*pisq*n),(1/3.0))*vn);
  double dp = HBARSQ*VLIGHTSQ*pisq / (2*MPROTON *AMUNIT*pow((3*pisq*z),(1/3.0))*vp);

  double mn = c1n*dn + c2n* pow(dn,alpha);
  double mz = c1p*dp + c2p* pow(dp,alpha);

  double f = (n==z) ? 0.0 : coo*(1.0+coo1*(pow(a,(-1.0/3.0)) - pow(a,(-2.0/3.0))));
  double mo = f*mz*mn/(mz+mn);

  deltaz = ((int)z%2 == 0) ? 1 : 0;
  deltan = ((int)n%2 == 0) ? 1 : 0;
  delta  =  mz*deltaz + mn*deltan - mo*deltaz*deltan;

  if(p == neutron) delta = mn;
  else             delta = mz;

  return(delta);
}


/**********************************************************/
/*     BCS U/V Factor                                     */
/**********************************************************/
double bcsVfactor(const double ef, const double eb, const double delta)
{
  double de = eb - ef;
  double ej = sqrt( de*de + delta*delta );
  return( sqrt( (1.0 - de/ej)*0.5 ) );
}


/**********************************************************/
/*     Solve BCS Gap Equation for Fixed Delta             */
/**********************************************************/
double bcsGapEquation(const int n, const double ef, const double delta, Bstate *bst)
{
  double d2  = delta*delta;
  double ef1 = ef*0.1;
  double ef2 = ef*2.0;
  int    k   = 0;

  for(k=0 ; k<MAX_BCSITER ; k++){
    if( fabs(ef1-ef2) < 1e-06) break;

    double ef = (ef1+ef2)/2.0;
    double g  = 0.0;
    double c  = 0.0;

    for(int i=0 ; i<MAX_SPLEVEL ; i++){
      if(bst[i].bind > 0.0) break;
      double y  = bst[i].bind-ef;
      double z  = sqrt( y*y+d2 );
      g += (bst[i].j2+1.0)/z;
      c += (bst[i].j2+1.0)*(1.0 - y/z);
    }
    g = 4.0/g;
    c = c/2.0;

    if(c < (double)n){
      ef2 = ef;
    }else{
      ef1 = ef;
    }
  }
  if(k == MAX_BCSITER){
    message << "gap equation cannot be solved";
    cohTerminateCode("bcsGapEquation");
  }

  return(ef);
}


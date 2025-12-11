/******************************************************************************/
/*  ddxkalbach.cpp                                                            */
/*        calculate Kalbach's systematics                                     */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "ddxkalbach.h"
#include "polysq.h"


static const int NDDX =   61;   // calculate angular points

static double ddxKalbachAfactor (const int, const int, const double, const double, const double);

static double  ecms = 0.0;      // center-of-mass energy
static double  sa = 0.0, sb[7]; // separation energies
static int     incid = 0;       // incident particle ID

/**********************************************************/
/*      Separation Energies Defined in Kalbach Syst.      */
/**********************************************************/
void ddxKalbachSetParm(const double ein, ZAnumber comp, ZAnumber proj)
{
  double   ib = 0.0;
  ZAnumber ejec(0,0),resd(0,0);

  sb[0] = 0.0;
  ecms  = ein;

  for(int c=1 ; c<7 ; c++){

    switch(c){
    case  1: ejec.setZA(0,1); ib =  0.0  ;  break;
    case  2: ejec.setZA(1,1); ib =  0.0  ;  break;
    case  3: ejec.setZA(2,4); ib = 28.296;  break;
    case  4: ejec.setZA(1,2); ib =  2.225;  break;
    case  5: ejec.setZA(1,3); ib =  7.718;  break;
    case  6: ejec.setZA(2,3); ib =  8.482;  break;
    default: ejec.setZA(0,0); ib =  0.0  ;  break;
    }

    resd = comp - ejec;
    if( (proj.getZ() == ejec.getZ()) && (proj.getA() == ejec.getA()) ) incid = c;

    double xc = (double)(comp.getN()-comp.getZ())*(double)(comp.getN()-comp.getZ());
    double xb = (double)(resd.getN()-resd.getZ())*(double)(resd.getN()-resd.getZ());
    double yc = pow((double)comp.getA(),1.0/3.0);
    double yb = pow((double)resd.getA(),1.0/3.0);
    double zc = (double)comp.getZ()*(double)comp.getZ();
    double zb = (double)resd.getZ()*(double)resd.getZ();

    sb[c] = 15.68 * (comp.getA()    - resd.getA()   )
          - 28.07 * (xc/comp.getA() - xb/resd.getA())
          - 18.56 * (yc*yc - yb*yb)
          + 33.22 * (xc/pow(yc,4.0) - xb/pow(yb,4.0))
          - 0.717 * (zc/yc - zb/yb)
          + 1.211 * (zc/comp.getA() - zb/resd.getA())
          - ib;
  }

  sa = sb[incid];
}


/**********************************************************/
/*      Kalbach Systematics for DDX                       */
/**********************************************************/
void ddxKalbach(
const int    c,      // channel number
const double ep,     // secondary particle energy
const double f,      // MSD fraction
const int    nleg,   // order of Legender 
double *cl)          // Legendre coefficients (output)
{
  double *ddx = new double [NDDX]; // double differential cross section, particle energies
  double *ang = new double [NDDX];

  double dt = 180.0/(NDDX-1);

  /*** Kalbach systematics a-parameter */
  double a = ddxKalbachAfactor(incid,c,ep,sa,sb[c]);
  double x = a/(2.0*sinh(a));

  /*** calculate angular distribution */
  for(int i=0 ; i<NDDX ; i++){
    ang[i] = (double)i*dt;
    double q = ang[i]/180.0 * PI;
    ddx[i] = x * (cosh(a*cos(q)) + f*sinh(a*cos(q)));
  }

  /*** calculate Legendre expansion coefficients by LSQ method */
  LSQLegendre(false,NDDX,nleg,ang,ddx,cl);

  delete [] ddx;
  delete [] ang;
}


/**********************************************************/
/*      a(E) Factor in Kalbach Systematics                */
/**********************************************************/
double ddxKalbachAfactor(const int a, const int b, const double eout, const double sa, const double sb)
{
  const double et1 = 130.0;
  const double et3 =  41.0;
  const double c1 = 0.04, c2 = 1.8e-06, c3 = 6.7e-07;

  double ea = ecms + sa;
  double eb = eout + sb;

  double e1 = fmin(ea,et1);
  double e3 = fmin(ea,et3);

  double x1 = e1*eb/ea;
  double x3 = e3*eb/ea;

  double ma = (a == 3) ? 0.0 : 1.0;
  double mb = (b == 3) ? 2.0 : ((b == 1) ? 0.5 : 1.0);

  return( c1*x1 + c2*pow(x1,3.0) + c3*pow(x3,4.0)*ma*mb );
}

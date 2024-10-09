/******************************************************************************/
/*  FRDMzeroenergy.cpp                                                        */
/*        Zero-Point Energy in FRDM                                           */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "FRDM.h"
#include "polysq.h"
#include "global.h"


/**********************************************************/
/*   Zero-Point Energy                                    */
/**********************************************************/
double FRDMZeropointEnergy(PolarGLPoint *gp, MFTSystem *sys, Basis *basis, HFInterface *spspec, FRDMEnergy *energy)
{
  const double kparm = 0.33;
  /*** surpress output */

  GlobalPrintMFT psave;
  psave = pmf;
  pmf.off();
  ctl.expansion = false;

  double a[4],x[5],y[5],e1,e2;
  double eps0 = sys->getEps(2);

  /*** perturb eps2 */
  double deps = 0.025;
  if(sys->getA() < 40) deps = 0.05;

  /*** calculate energies at eps +/- 2d and eps +/- d */
  for(int i=0 ; i<5 ; i++){
    x[i] = eps0 + deps*(i-2);
    if(i != 2){
      sys->setEps(2,x[i]);
      e1 = FRDMMacroEnergy(gp,sys);
      e2 = FRDMMicroEnergy(sys,basis,spspec);
    }
    else{
      e1 = energy->macro;
      e2 = energy->micro;
    }
    y[i] =  e1 + e2;
  }

  /*** fit by the third-order poloynomial, dE^2/d eps^2 = 6a eps + 2b */
  LSQPolynomial(5,4,x,y,a);
  double ce = 6.0*a[3]*eps0 + 2.0*a[2];

  /*** in case of convex, return zero */
  double ezero = 0.0;
  if(ce > 0.0){
    /*** nucleus mass devided by (h-bar c)^2 = 1/MeV fm^2 */
    double m0 = (ENEUTRON * sys->getN() + EPROTON * sys->getZ() + y[2])
               / (VLIGHTSQ * HBARSQ);
    double mr = m0 * sys->R0*sys->R0; // in [1/ MeV]

    /*** incompressible irrotational inertia */
    double x1 = 1.0 + 2.0/9.0 * eps0*eps0;
    double x2 = 1.0 - 2.0/3.0 * eps0;
    double x3 = 1.0 - 1.0/3.0 * eps0*eps0 - 2.0/27.0 * eps0*eps0*eps0;
    double be = 2.0*x1/(15.0*x2*x2) * pow(x3,-4.0/3.0) * mr;

    /*** zero-point energy, omega x hbar in MeV  (not sure...) */
    ezero = 0.5 * kparm * sqrt(ce / be) / (VLIGHT*HBAR);
  }

  /*** restore print control */
  pmf = psave;
  sys->setEps(2,eps0);

  return(ezero);
}

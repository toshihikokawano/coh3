/******************************************************************************/
/*  HFmultipole.cpp                                                           */
/*        multipole moment calculation                                        */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "HF.h"

static double renorm[2];
static double MPMNeckPosition(Basis *, GridData *, LaguerrePolynomial *, HermitePolynomial *);

/**********************************************************/
/*      Multipole Moment                                  */
/*      ----------------                                  */
/*      Calculate Q10, Q20, Q30, Q40                      */
/**********************************************************/
void HFMultipoleMoment(Basis *basis, Moment *moment, GridData *rho, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  const double eps = 1.0e-12;
  const double aqn = 1.0;

  /*** multipole moments,  Q10 (zcm), Q20, Q30 and Q40 */

  renorm[0] = renorm[1] = 0.0;
  double zcm = 0.0, q20 = 0.0, q30 = 0.0, q40 = 0.0, qn = 0.0;

  double c1 = PI/(basis->bz * basis->bp2);

  for(int i=0 ; i<fl->nr ; i++){
    double r  = sqrt(fl->x[i])/basis->bp;
    double r2 = r*r;

    for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;
      double z  = fh->x[j]/basis->bz;
      double z2 = z*z;

      double w = c1 * fh->G[j] * fl->G[i];

      renorm[0] += w * rho[0].p[i][j];
      renorm[1] += w * rho[1].p[i][j];

      double c2 = w * (rho[0].p[i][j] + rho[1].p[i][j]);

      zcm += c2 * z;
      q20 += c2 * (2.0*z2 - r2);
      q30 += c2 * z * (z2 - 1.5*r2);
      q40 += c2 * (z2*z2 - 3.0 * z2 * r2 + 3.0/8.0 * r2*r2);
    }
  }
  zcm /= renorm[0] + renorm[1];
  q30 *= sqrt(7.0/PI4);
  q40 *= sqrt(9.0/PI4);

  /*** neck coordinate QN */
  double zmin = MPMNeckPosition(basis,rho,fl,fh);

  for(int i=0 ; i<fl->nr ; i++){
    for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;

      double w = fh->G[j] * fl->G[i];
      double z = fh->x[j]/basis->bz;
      double y = (z - zmin)/aqn;

      qn += c1 * w * (rho[0].p[i][j] + rho[1].p[i][j]) * exp(-y*y);
    }
  }

  moment->zcm  = (fabs(zcm) > eps) ? zcm : 0.0;
  moment->q20  = (fabs(q20) > eps) ? q20 : 0.0;
  moment->q30  = (fabs(q30) > eps) ? q30 : 0.0;
  moment->q40  = (fabs(q40) > eps) ? q40 : 0.0;
  moment->qn   = (fabs(qn ) > eps) ? qn  : 0.0;
  moment->zmin = zmin;
/*
  std::cout << zcm << std::endl;
  std::cout << q20 << std::endl;
  std::cout << q30 << std::endl;
  std::cout << q40 << std::endl;
  std::cout << qn  << std::endl;
  std::cout << zmin << std::endl;
*/
}


/**********************************************************/
/*      Find Minumim Z                                    */
/**********************************************************/
double MPMNeckPosition(Basis *basis, GridData *rho, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  double *zbuf, *zdens;
  zbuf = new double [fh->getNsize()];
  zdens = &zbuf[(fh->getNsize()-1)/2];

  int jmin = 0;

  /*** integrate density in the Z-direction */
  for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j == 0) continue;
    zdens[j] = 0.0;
    for(int i=0 ; i<fl->nr ; i++){
      zdens[j] += fh->G[j] * fl->G[i] * (rho[0].p[i][j] + rho[1].p[i][j]);
    }
  }

  /*** find highest density location on the left side */
  int    j1 = 0;
  double d1 = 0.0;
  for(int j=-fh->nz ; j<0 ; j++){
    if(zdens[j] > d1){ j1 = j; d1 = zdens[j]; }
  }

  /*** on the right side */
  int    j2 = 0;
  double d2 = 0.0;
  for(int j=1 ; j<=fh->nz ; j++){
    if(zdens[j] > d2){ j2 = j; d2 = zdens[j]; }
  }

  /*** find density takes minimum between two highest points */
  if( (j1 == -1) && (j2 == 1)){
    /*** if center */
    jmin = 0;
  }
  else{
    double zdensmin = 1.0e+6;
    jmin = j1;
    for(int j = j1 ; j<=j2 ; j++){
      if(j == 0) continue;
      if(zdens[j] < zdensmin){
        jmin = j;
        zdensmin = zdens[j];
      }
    }
  }

  delete [] zbuf;

  double zmin = fh->x[jmin]/basis->bz;
  return(zmin);
}


/**********************************************************/
/*      Renormalization of Local Density                  */
/**********************************************************/
void HFRenormalizeLocalDensity(const int n, const int z, GridData *rho, GridData *tau, GridData *d2rho, GridData *divJ, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  for(int p=0 ; p<2 ; p++){
    double c = (p == 0) ? n/renorm[p] : z/renorm[p];

    for(int i=0 ; i<fl->nr ; i++){
      for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;

        rho[p].p[i][j]   *= c;
        tau[p].p[i][j]   *= c;
        d2rho[p].p[i][j] *= c;
        divJ[p].p[i][j]  *= c;
      }
    }
  }
}


/******************************************************************************/
/*  HFcoulomb.cpp                                                             */
/*        Coulomb potential calculation                                       */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "HF.h"

/**********************************************************/
/*      Coulomb Potential                                 */
/**********************************************************/
void HFCoulomb(Basis *basis, GridData *d2rho, GridData *vcoul, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  const double e2 = COULOMBSQ * PERMITTIV;

  double c = e2 /(basis->bz * basis->bp2);

  for(int i1=0 ; i1<fl->nr ; i1++){
    double r1 = sqrt(fl->x[i1])/basis->bp;

    for(int j1=-fh->nz ; j1<=fh->nz ; j1++){ if(j1==0) continue;
      double z1 = fh->x[j1]/basis->bz;

      vcoul->p[i1][j1] = 0.0;

      for(int i2=0 ; i2<fl->nr ; i2++){
        double r2 = sqrt(fl->x[i2])/basis->bp;
        double dr = r1 + r2;

        for(int j2=-fh->nz ; j2<=fh->nz ; j2++){ if(j2==0) continue;
          double z2 = fh->x[j2]/basis->bz;
          double dz = z1 - z2;
          double d2 = dz*dz + dr*dr;
          double d  = sqrt(d2);

          double arg = 1.0 - 4.0*r1*r2/d2;
          double w = c * fl->G[i2] * fh->G[j2];
          double elint = 1.0 + arg*(0.4630151 + 0.1077812 * arg);
          if(arg > 1.0e-6) elint -= arg*(0.2452727 + 0.0412496*arg)*log(arg);
          vcoul->p[i1][j1] += w * d2rho->p[i2][j2] * elint * d;
        }
      }
    }
  }
/*
  for(int i=0 ; i<fl->nr ; i++){
    double r = sqrt(fl->x[i])/basis->bp;

    for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j == 0) continue;
      double z = fh->x[j]/basis->bz;

      std::cout << std::setw(12) << r << std::setw(12) << z;
      std::cout << std::setw(12) << vcoul->p[i][j] << std::endl;
    }
    std::cout << std::endl;
  }
*/
}

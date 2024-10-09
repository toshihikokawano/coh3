/******************************************************************************/
/*  HFdensity.cpp                                                             */
/*        Hartree-Fock local density calculation                              */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#undef DEBUG

#include "HF.h"


/**********************************************************/
/*      Density Matrix, T diag(2v) T'                     */
/*      Vautherin PRC 7, 296 (1973), Eq.(5.4)             */ 
/**********************************************************/
void HFDensityMatrix(Basis *b, double *v2, double *matV, Operator *H)
{
  /*** clear matrix */
  for(int i=0 ; i<b->nv ; i++){
    for(int j=0 ; j<=i ; j++) matV[i*(i+1)/2+j] = 0.0;
  }

  /*** store T diag(2*v) T' into a working symmetric matrix matV */
  for(int ib=0 ; ib<b->nb ; ib++){
    int n = b->ndata[ib];

    for(int i=0 ; i<n ; i++){
      int i0 = b->index[ib] + i;
      
      for(int j=0 ; j<=i ; j++){
        int j0 = b->index[ib] + j;

        int ij = i0*(i0+1)/2+j0;
        matV[ij] = 0.0;

        for(int k=0 ; k<n ; k++){
          int k0 = b->index[ib] + k;

          matV[ij] += 2.0*v2[k0] * H[ib].vector[k][i] * H[ib].vector[k][j];
        }
      }
    }
  }

#ifdef DEBUG
  for(int i=0 ; i<b->nv ; i++){
    for(int j=0 ; j<=i ; j++){
      int ij = i*(i+1)/2+j;
      std::cout << std::setw(12) << matV[ij];
    }
    std::cout << std::endl;
  }
#endif
}


/**********************************************************/
/*      rho, d2rho, divJ, and tau                         */
/*      Vautherin PRC 7, 296 (1973), Eq.(5.5)             */ 
/**********************************************************/
void HFLocalDensity(Basis *basis, double *matV, GridData *rho, GridData *d2rho, GridData *tau, GridData *divJ, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  const double rhocut = 1.0e-07;

  rho->clear();
  d2rho->clear();
  tau->clear();
  divJ->clear();

  for(int kb=0 ; kb<basis->nb ; kb++){
    for(int k1=basis->index[kb] ; k1<basis->index[kb] + basis->ndata[kb] ; k1++){

      int nr1 = basis->state[k1].getNr();
      int nz1 = basis->state[k1].getNz();
      int l1  = basis->state[k1].getL();
      int s1  = basis->state[k1].getS();

      for(int k2=basis->index[kb] ; k2<=k1 ; k2++){

        int nr2 = basis->state[k2].getNr();
        int nz2 = basis->state[k2].getNz();
        int l2  = basis->state[k2].getL();
        int s2  = basis->state[k2].getS();

        if( (2*l1 + s1) != (2*l2 + s2) ) continue;

        double rho12 = matV[k1*(k1+1)/2 + k2];
        if( fabs(rho12) < rhocut ) continue;

        /*** Sigma(k1) != Sigma(k2) case */
        if( s1 != s2 ){
          for(int i=0 ; i<fl->nr ; i++){
            double eta = fl->x[i];
            double L00 = fl->P0[nr1][l1][i] * fl->P0[nr2][l2][i];
            double L01 = fl->P0[nr1][l1][i] * fl->P1[nr2][l2][i];
            double L10 = fl->P1[nr1][l1][i] * fl->P0[nr2][l2][i];

            for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;
              double H01 = fh->P0[nz1][j] * fh->P1[nz2][j];
              double H10 = fh->P1[nz1][j] * fh->P0[nz2][j];

              divJ->p[i][j] += rho12 * pow(eta, (l1+l2-1)/2.0)
                *  ( (s1-s2)/2.0 * (H01*L10 - H10*L01) - L00*(l1*H01 + l2*H10) );
            }
          }
        }
        /*** Sigma(k1) == Sigma(k2)
             then Lambda(k1) == Lambda(k2), because omega is the same */
        else{
          double f = ((k1 == k2) ? 1.0 : 2.0) * rho12;
          double xj = 0.5 * f * basis->bp / basis->bz * s1 * l1;
          double x1 = nz1 + nz2 + 1.0;
          double x2 = 2.0*(nr1 + nr2 + l1 + 1.0);

          for(int i=0 ; i<fl->nr ; i++){
            double eta = fl->x[i];
            double L00 = fl->P0[nr1][l1][i] * fl->P0[nr2][l2][i];
            double L01 = fl->P0[nr1][l1][i] * fl->P1[nr2][l2][i];
            double L10 = fl->P1[nr1][l1][i] * fl->P0[nr2][l2][i];
            double L11 = fl->P1[nr1][l1][i] * fl->P1[nr2][l2][i];

            for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;
              double zeta = fh->x[j];
              double H00 = fh->P0[nz1][j] * fh->P0[nz2][j];
              double H11 = fh->P1[nz1][j] * fh->P1[nz2][j];

              rho->p[i][j] += f * pow(eta,(double)l1) * H00 * L00;
              tau->p[i][j] += f * pow(eta,(double)(l1-1)) *
                             ( basis->bz2 * eta * H11 * L00
                             + basis->bp2 * H00 * (L11 + l1 * l2 * L00) );
              d2rho->p[i][j] += f * pow(eta,(double)l1) * H00 * L00
                             * (basis->bz2*(zeta*zeta - x1) + basis->bp2*(eta - x2));
              divJ->p[i][j] += xj * pow(eta,(double)(l1-1)) * H00 * (L01 + L10);
            }
          }
        }
      }
    }
  }

  /*** multiply beta_z beta_perp^2/pi x exp(-eta - xi^2) */
  double c1 = basis->bz * basis->bp2 / PI;
  double c2 = basis->bz2 * basis->bp2*basis->bp / PI;
  for(int i=0 ; i<fl->nr ; i++){
    for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;

      double c3 = c1 * exp(-fh->x[j]*fh->x[j] - fl->x[i]);
      double c4 = c2 * exp(-fh->x[j]*fh->x[j] - fl->x[i]);

      rho->p[i][j]  *= c3;
      tau->p[i][j]  *= c3;
      d2rho->p[i][j] = 2.0 * (c3 * d2rho->p[i][j] + tau->p[i][j]);
      divJ->p[i][j] *= 2.0 * c4;
    }
  }

#ifdef DEBUG
  for(int i=0 ; i<fl->nr ; i++){
    for(int j=-fh->nz ; j<=fh->nz ; j++){
      if(j == 0) continue;
      std::cout << std::setw(4) << i << std::setw(4) << j;
      std::cout << std::setw(12) << rho->p[i][j];
      std::cout << std::setw(12) << d2rho->p[i][j];
      std::cout << std::setw(12) << divJ->p[i][j];
      std::cout << std::setw(12) << tau->p[i][j] << std::endl;
    }
  }
#endif
}


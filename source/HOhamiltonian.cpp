/******************************************************************************/
/*  HOhamiltonian.cpp                                                         */
/*        Harmonic Oscillator Hamiltonian                                     */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "HObasis.h"


static inline double HOweight(const double, const double, const double);

#undef DEBUG

/**********************************************************/
/*      Calculate HO Hamiltonian and Store Elements in H  */
/*      see Damgaard et al. NPA 135, 432 (1969)           */
/*          Vautherin PRC 7, 296 (1973)                   */ 
/**********************************************************/
void HOHamiltonian(const int k, Basis *basis, OneBodyPotential *pot, Operator *H, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  double *buf, *L00, *L10, *L01, *L11, *H00, *H01, *H10, *H11;

  buf = new double [4*fl->nr + 4*(2*fh->nz + 1)];

  int idx = 0;
  L00 = &buf[idx];          idx += fl->nr;
  L01 = &buf[idx];          idx += fl->nr;
  L10 = &buf[idx];          idx += fl->nr;
  L11 = &buf[idx];          idx += fl->nr;
  H00 = &buf[idx + fh->nz]; idx += 2*fh->nz + 1;
  H01 = &buf[idx + fh->nz]; idx += 2*fh->nz + 1;
  H10 = &buf[idx + fh->nz]; idx += 2*fh->nz + 1;
  H11 = &buf[idx + fh->nz];

  /*** clear H */
  H->clear();

  /*** calculate matrix elements for Skyrme interaction */
  for(int i1=0 ; i1<basis->ndata[k] ; i1++){
    int j1  = i1 + basis->index[k];
    int l1  = basis->state[j1].getL();
    int s1  = basis->state[j1].getS();
    int nz1 = basis->state[j1].getNz();
    int nr1 = basis->state[j1].getNr();

    for(int i2=0 ; i2<=i1 ; i2++){
      int j2  = i2 + basis->index[k];
      int l2  = basis->state[j2].getL();
      int s2  = basis->state[j2].getS();
      int nz2 = basis->state[j2].getNz();
      int nr2 = basis->state[j2].getNr();

      double hkin = 0.0, hpot = 0.0, hsp0 = 0.0, hsp1 = 0.0;

      for(int i=0 ; i<fl->nr ; i++){
        L00[i] = fl->P0[nr1][l1][i] * fl->P0[nr2][l2][i];
        L11[i] = fl->P1[nr1][l1][i] * fl->P1[nr2][l2][i];
        L01[i] = fl->P0[nr1][l1][i] * fl->P1[nr2][l2][i];
        L10[i] = fl->P1[nr1][l1][i] * fl->P0[nr2][l2][i];
      }

      for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j == 0) continue;
        H00[j] = fh->P0[nz1][j] * fh->P0[nz2][j];
        H11[j] = fh->P1[nz1][j] * fh->P1[nz2][j];
        H01[j] = fh->P0[nz1][j] * fh->P1[nz2][j];
        H10[j] = fh->P1[nz1][j] * fh->P0[nz2][j];
      }

      for(int i=0 ; i<fl->nr ; i++){
        double wk = HOweight(fl->W[i],fl->x[i],(double)l1-1.0);
        double wp = HOweight(fl->W[i],fl->x[i],(double)l1);
        double ws = HOweight(fl->W[i],fl->x[i], 0.5*(l1+l2-1.0)) * basis->bp * basis->bz;

        for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j == 0) continue;

          /*** Kinetic Operator, V:Eq.(4.21) */
          if( (s1 == s2) && (l1 == l2) ){
            double x00 =  fl->x[i] * L00[i];
            double x11 =  L11[i] + l1*l2 * L00[i];
            double vt  = (basis->bz2 * H11[j] * x00 + basis->bp2 * H00[j] * x11);
            hkin += fh->W[j] * wk * vt * pot->effmass[i][j];
          }

          /*** Central Potential, V:Eq.(4.22) or D:Eq.(42) */
          if( (s1 == s2) && (l1 == l2) ){
            hpot += fh->W[j] * wp * H00[j] * L00[i] * pot->central[i][j];
          }

          /*** Spin-Orbit Potential, V:Eq.(4.23) or D:Eq.(55), (56) */
          /*** diagonal matrix elements */
          if( (s1 == s2) && (l1 == l2) ){
            double xsp = -2.0*basis->bp2 * l1 * (s1*0.5) * (L01[i] + L10[i]);
            hsp0 += fh->W[j] * wk * xsp * H00[j] * pot->spinorb[i][j];
          }

          /*** off-diagonal elements, sigma1 != sigma2 */
          if(s1 != s2){
            double x01 = (l2-l1) * L01[i] + l2 * L00[i];
            double x10 = (l1-l2) * L10[i] + l1 * L00[i];
            hsp1 += fh->W[j] * ws * (H10[j] * x01 + H01[j] * x10) * pot->spinorb[i][j];
          }
        }
      }
      double sum = hkin + hpot + hsp0 + hsp1;
      H->setVal(i1,i2,sum);
    }
  }

#ifdef DEBUG
  for(int i=0 ; i<basis->ndata[k] ; i++){
    for(int j=0 ; j<basis->ndata[k] ; j++){
      std::cout << std::setw(12) << H->getVal(i,j);
    }
    std::cout << std::endl;
  }
#endif
  delete [] buf;
}


static inline double HOweight(const double w, const double x, const double p)
{
  double z = 0.0;
  /*** avoid 0^0 */
  if( (x == 0.0) && (p == 0.0) ){
    z = w;
  }
  else if(x >= 0){
    z = w * pow(x,p);
  }
  return(z);
}


/******************************************************************************/
/*  HFsolvebcs.cpp                                                            */
/*        solve gap equation                                                  */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "HF.h"


static void BCSSingleParticleDensity (const double, Basis *, SPEnergy *, Operator *, GridData *, LaguerrePolynomial *, HermitePolynomial *);
static void BCSMatrixElement (Basis *, double *, BCSParameter *, SPEnergy *, GridData *, GridData *, LaguerrePolynomial *, HermitePolynomial *);
static void BCSConstantMatrixElement (const int, Basis *, double *, BCSParameter *);
static double BCSGapEquation (const bool, const int, double, Basis *, double *, SPEnergy *, double *);
static inline double BCSCutoffFactor1 (const double, const double);
static inline double BCSCutoffFactor2 (const double, const double);
static inline double BCSDenominator (const double, const double, const double, const double);


/**********************************************************/
/*      BCS Equation                                      */
/*      ------------                                      */
/*      Density dependent pairing model                   */
/**********************************************************/
void HFSolveBCS(const int np, Basis *basis, BCSParameter *bcs, SPEnergy *spe, double *v2, double *matV, Operator *H, GridData *rho, GridData *spd, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  double *d0;
  d0 = new double [basis->nv];

  double lambda  = (spe->energy[np/2-1] + spe->energy[np/2])*0.5;


  /*** density dependent delta interaction */
  bool delta = false;
  if(bcs->pairing == "delta"){
    /*** single-particle density sum_s |phi_i(z,r,s)|^2 */
    BCSSingleParticleDensity(lambda, basis, spe, H, spd, fl, fh);
    BCSMatrixElement(basis, matV, bcs, spe, rho, spd, fl, fh);
    delta = true;
  }
  /*** constant seniority G */
  else if(bcs->pairing == "seniority"){
    BCSConstantMatrixElement(np, basis, matV, bcs);
    delta = false;
  }

  /*** solve BCS equation */
  lambda = BCSGapEquation(delta, np, lambda, basis, matV, spe, d0);

  /*** occupation number and pairing energy */
  double epair = 0.0;
  double gap   = 0.0;
  double denom = 0.0;

  for(int k=0 ; k<basis->nv ; k++){
    double f = (delta) ? BCSCutoffFactor1(spe->energy[k],lambda)
                       : BCSCutoffFactor2(spe->energy[k],lambda);
    double e = BCSDenominator(spe->energy[k],lambda,f,d0[k]);
    v2[spe->index[k]] = 0.5 * (1.0 - (spe->energy[k] - lambda)/e);

    epair += f*d0[k]*d0[k]/e;
    double x = sqrt( (1.0 - v2[spe->index[k]]) * v2[spe->index[k]] );
    gap   += x*d0[k];
    denom += x;
  }

  bcs->chemical_potential = lambda;
  bcs->pairing_energy     = -0.5 * epair;
  bcs->average_gap        = (denom == 0.0) ? 0.0 : gap/denom;

  /*** copy v2 into spe in the same order */
  for(int k=0 ; k<basis->nv ; k++) spe->v2[k] = v2[spe->index[k]];
/*
  cout << "Chemical Potential " << setw(12) << bcs->chemical_potential << endl;
  cout << "Pairing Energy     " << setw(12) << bcs->pairing_energy << endl;
  cout << "Average Gap        " << setw(12) << bcs->average_gap << endl;
*/
  delete [] d0;
}


/**********************************************************/
/*      Singple Particle Density  |phi_k(i,j)|^2          */
/**********************************************************/
void BCSSingleParticleDensity(const double lambda, Basis *basis, SPEnergy *spe, Operator *H, GridData *spd, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  const double eps = 1.0e-10;

  /*** truncate if occupation probability is too small */
  for(int k=0 ; k<basis->nv ; k++){
    double f = BCSCutoffFactor1(spe->energy[k],lambda);
    spe->truncate[k] = (f < eps) ? 1 : 0;
  }

  /*** sum of s.p. densities for spin up and down */
  double c1 = basis->bz * basis->bp2 / PI;

  for(int k=0 ; k<basis->nv ; k++){
    if(spe->truncate[k]){
      for(int i=0 ; i<fl->nr ; i++){
        for(int j=-fh->nz ; j<=fh->nz ; j++) spd[k].p[i][j] = 0.0;
      }
      continue;
    }

    /*** indices before sorting */
    int k1 = spe->subindex[k];  // index of k-th eigenvalue inside the block
    int kb = spe->block[k];     // block index

    for(int i=0 ; i<fl->nr ; i++){
      for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j == 0) continue;

        double c2 = c1 * exp( -(fl->x[i] + fh->x[j]*fh->x[j]) );

        double sum1 = 0.0; // spin up sum
        double sum2 = 0.0; // spin down sum

        for(int k2 = 0 ; k2<basis->ndata[kb] ; k2++){
          int k0 = basis->index[kb] + k2;

          if(fabs(H[kb].vector[k1][k2]) < eps) continue;

          int nz = basis->state[k0].getNz();
          int nr = basis->state[k0].getNr();
          int nl = basis->state[k0].getL();
          int ns = basis->state[k0].getS();

          double p = 0.0;
          if(nl == 0){ p = 1.0; }
          else if(fl->x[i] >= 0.0){ p = pow(fl->x[i],(double)nl*0.5); }
          else continue;

          p *= H[kb].vector[k1][k2] * fh->P0[nz][j] * fl->P0[nr][nl][i];

          if(ns == 1){ sum1 += p; }
          else{        sum2 += p; }
        }

        spd[k].p[i][j] = c2 * (sum1*sum1 + sum2*sum2);
      }
    }
  }
}


/**********************************************************/
/*      Pairing Matrix Elements                           */
/**********************************************************/
void BCSMatrixElement(Basis *basis, double *matV, BCSParameter *bcs, SPEnergy *spe, GridData *rho, GridData *spd, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  const double eta   = 1.0;
  const double beta  = 1.0;
  const double rho_c = 0.16;

  for(int k1=0 ; k1<basis->nv ; k1++){
    for(int k2=0 ; k2<=k1 ; k2++) matV[k1*(k1+1)/2 + k2] = 0.0;
  }

  double c = bcs->strength * PI / (basis->bz * basis->bp2);

  for(int k1=0 ; k1<basis->nv ; k1++){
    if(spe->truncate[k1]){
      for(int k2=0 ; k2<=k1 ; k2++) matV[k1*(k1+1)/2 + k2] = 0.0;
      continue;
    }

    for(int k2=0 ; k2<=k1 ; k2++){
      if(spe->truncate[k2]){
        matV[k1*(k1+1)/2 + k2] = 0.0;
        continue;
      }

      double v = 0.0;
      for(int i=0 ; i<fl->nr ; i++){
        for(int j=-fh->nz; j<=fh->nz ; j++){ if(j == 0) continue;
          double x = rho[0].p[i][j] + rho[1].p[i][j];
          v += fh->G[j] * fl->G[i] * spd[k1].p[i][j] * spd[k2].p[i][j]
             * (1.0 - eta * pow(x/rho_c, beta));
        }
      }
      matV[k1*(k1+1)/2 + k2] = c * v;
    }
  }
/*
  for(int k1=0 ; k1<basis->nv ; k1++){
    for(int k2=0 ; k2<=k1 ; k2++) cout << setw(12) << matV[k1*(k1+1)/2 + k2];;
    cout << endl;
  }
*/
}


/**********************************************************/
/*      Constant Pairing Matrix Elements                  */
/**********************************************************/
void BCSConstantMatrixElement(const int np, Basis *basis, double *matV, BCSParameter *bcs)
{
  double v = - fabs(bcs->strength)/(11.0 + (double)np);
  for(int k1=0 ; k1<basis->nv ; k1++){
    for(int k2=0 ; k2<=k1 ; k2++) matV[k1*(k1+1)/2 + k2] = v;
  }
/*
  for(int k1=0 ; k1<basis->nv ; k1++){
    for(int k2=0 ; k2<=k1 ; k2++) cout << setw(12) << matV[k1*(k1+1)/2 + k2];;
    cout << endl;
  }
*/
}


/**********************************************************/
/*      Solve BCS Gap Equation                            */
/**********************************************************/
double BCSGapEquation(const bool delta, const int np, double lambda, Basis *basis, double *matV, SPEnergy *spe, double *d0)
{
  const double eps     = 1.0e-6;
  const double etrunc  = 8.0;
  const int itermax    = 500;

  /*** cut s.p. states above this energy */
  double truncate = etrunc + lambda;

  for(int k=0 ; k<basis->nv ; k++) d0[k] = 1.0;

  /*** iterate gap equation until chemical potential converges */
  for(int itr=0 ; itr<itermax ; itr++){

    bool pass1 = false, pass2 = false;

    /*** first, determine the highest single particle state */
    int lcut = 0;
    for(int k2=0 ; k2<basis->nv ; k2++){
      double f2 = (delta) ? BCSCutoffFactor1(spe->energy[k2],lambda)
                          : BCSCutoffFactor2(spe->energy[k2],lambda);
      if((spe->energy[k2] > truncate) && (fabs(f2) < 1.0e-5)){
        lcut = k2; break;
      }
    }

    /*** state dependent energy gap
         -1/2 sum V(ij) f(i) D(i) / sqrt{(e(j)-lambda)^2 + f(j) D(j)^2} */

    double deld = 0.0;
    for(int k1=0 ; k1<lcut ; k1++){
      double sum = 0.0;
      for(int k2=0 ; k2<lcut ; k2++){
        double f2 = (delta) ? BCSCutoffFactor1(spe->energy[k2],lambda)
                            : BCSCutoffFactor2(spe->energy[k2],lambda);
        double e  = BCSDenominator(spe->energy[k2],lambda,f2,d0[k2]);

        int k = (k1 >= k2) ? k1*(k1+1)/2 + k2 : k2*(k2+1)/2 + k1;
        sum += matV[k] * f2 * d0[k2] / e;
      }
      double x = -0.5*sum;
      deld += fabs(x - d0[k1]);
      d0[k1] = x;
    }
    if( (deld/lcut) < eps) pass1 = true;

    /*** chemical potential */
    double x1 = 0.0;
    double x2 = 0.0;

    for(int k=0 ; k<lcut ; k++){
      double f = (delta) ? BCSCutoffFactor1(spe->energy[k],lambda)
                         : BCSCutoffFactor2(spe->energy[k],lambda);
      double e = BCSDenominator(spe->energy[k],lambda,f,d0[k]);

      x1 += 1.0 - spe->energy[k]/e;
      x2 += 1.0/e;
    }

    double x = lambda;
    lambda = (np-x1)/x2;

    if( fabs(lambda - x) < eps ) pass2 = true;

    if(pass1 && pass2) break;
  }

  return(lambda);
}


/**********************************************************/
/*      Cut Off Factor                                    */
/**********************************************************/
static inline double BCSCutoffFactor1(const double e, const double lambda)
{
  const double cutenergy = 5.0;
  const double diffuse   = 0.5;

  double x = (fabs(e - lambda) - cutenergy)/diffuse;
  double f = 1.0/sqrt(1.0 + exp(x));

  return(f);
}

static inline double BCSCutoffFactor2(const double e, const double lambda)
{
  const double cutenergy = 6.0;
  const double diffuse   = 0.2;

  double x = (e - (lambda + cutenergy))/diffuse;
  double f = 1.0/sqrt(1.0 + exp(x));

  return(f);
}


/**********************************************************/
/*      sqrt{(e-lambda)^2 + Delta^2}                      */
/**********************************************************/
static inline double BCSDenominator(const double e, const double lambda, const double f, const double d)
{
  double x = sqrt((e-lambda)*(e-lambda) + f*d*d);
  return(x);
}

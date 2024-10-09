/******************************************************************************/
/*  HOpolynomial.cpp                                                          */
/*        Hermite and Laguerre polynomials on each Gauss quadrature mesh      */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "HObasis.h"

static void HOHermite (const int, HermitePolynomial *);
static void HOLaguerre (const int, const int, LaguerrePolynomial *);
static void HOHermiteRoot (const int, double *);
static void HOLaguerreRoot (const int, double *, double);

extern double *fact;

#undef DEBUG


/**********************************************************/
/*      Hermite and Laguerre Polynomials                  */
/**********************************************************/
void HOPolynomials(const int nr, const int nz, Basis *basis, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  /*** Laguerre Polynomials */
  fl->setGridSize(basis->nrmax,nr,basis->lambdamax);
  fl->memalloc();
  HOLaguerre(nr,basis->lambdamax,fl);

  /*** Hermite Polynomials */
  fh->setGridSize(basis->nzmax,nz);
  fh->memalloc();
  HOHermite(nz,fh);
}


/**********************************************************/
/*      Hermite Polynomials                               */
/*      --------                                          */
/*      x[i]     X-axis points from -Nz to Nz (2Nz + 1)   */
/*      P0[n][i] n-th order Hermite polynomials at x[i]   */
/*      P1[n][i] = -x P0[n+1][i] + dP0[n][i]/dx           */
/*      W[i]     Gauss-Hermite weight, sum W = sqrt(PI)   */
/*      G[i]     weight with exp(x^2)                     */
/**********************************************************/
void HOHermite(const int nz, HermitePolynomial *f)
{
  double *x;
  x = new double [f->getNsize()];

  /*** find zeros of the Hermite polynomials */
  HOHermiteRoot(nz, x);

  /*** store values at each Z-grid position, from -nz to nz */
  for(int i=1 ; i<=f->nz ; i++){
    f->x[ i] = x[f->nz-i];
    f->x[-i] = x[f->nz+i-1];

    double x0 = f->x[i];

    f->P0[0][ i] =  1.0;      // H(0,i) = 1
    f->P0[1][ i] =  2.0 * x0; // H(1,i) = 2x

    f->P1[0][ i] = -x0;
    f->P1[1][ i] =  2.0 - 2.0*x0*x0;

    /*** note: this recurrance form causes less of significance
         H(n+1) = 2{x H(n) - n H(n-1)} */
    for(int m=1 ; m<f->getMsize()-1 ; m++){
      f->P0[m+1][i] =  2.0*(x0 * f->P0[m][i] - m * f->P0[m-1][i]);
      f->P1[m+1][i] = -x0 * f->P0[m+1][i] + 2.0 * (m+1) * f->P0[m][i];
    }

    /*** weight: 2^(n-1) n! sqrt(pi) / n^2 H(n-1,i)^2 */
    f->W[i] = f->W[-i] = SQRTPI * pow(2.0,(double)(nz-1))
                       * exp(fact[nz-1])/(nz * f->P0[nz-1][i] * f->P0[nz-1][i]);
    f->G[i] = f->G[-i] = f->W[i] * exp(x0*x0);
  }

  /*** normalization by Hermite orthogonality
       int Hn(x) Hm(x) exp(-x*x) dx = sqrt(pi) 2^n n! delta(n,m) */
  for(int i=1 ; i<=f->nz ; i++){
    for(int m=0 ; m<f->getMsize() ; m++){

      double c = sqrt(SQRTPI * pow(2.0,(double)m) * exp(fact[m]));

      f->P0[m][ i] = f->P0[m][i] / c;
      f->P0[m][-i] = f->P0[m][i] * ((m%2 == 0) ? 1.0 : -1.0);
      f->P1[m][ i] = f->P1[m][i] / c;
      f->P1[m][-i] = f->P1[m][i] * (((m+1)%2 == 0) ? 1.0 : -1.0);
    }
  }


#ifdef DEBUG
  for(int i=-f->nz ; i<=f->nz ; i++){
    std::cout << std::setw(12) << f->x[i];
    for(int n=0 ; n<f->getMsize() ; n++)  std::cout << std::setw(12) << f->P0[n][i];
    std::cout << std::endl;
  }
  std::cout << std::endl;
  for(int i=1 ; i<=f->nz ; i++){
    std::cout << std::setw(12) << f->x[i];
    for(int n=0 ; n<f->getMsize() ; n++)  std::cout << std::setw(12) << f->P1[n][i];
    std::cout << std::endl;
  }
#endif

  delete [] x;
}


/**********************************************************/
/*      Laguerre Polynomials                              */
/*      --------                                          */
/*      x[i]     X-axis points from 0 to Nr (Nr + 1)      */
/*      P0[n][l][i] n-th and l-th order Laguerre at x[i]  */
/*      P1[n][l][i] = m P0 - x P0 + 2x dP0/dx             */
/*      W[i]     Gauss-Laguerre weight                    */
/*      G[i]     weight with exp(x)                       */
/**********************************************************/
void HOLaguerre(const int nr, const int lmax, LaguerrePolynomial *f)
{
  double lambda = 0.0;

  /*** find zeros of Laguerre polynomials */
  HOLaguerreRoot(nr, f->x, lambda);

  /*** store values at each R-grid position, from 0 to nr */
  for(int i=0 ; i<nr ; i++){
    double x0 = f->x[i];

    for(int l=0 ; l<=lmax ; l++){
      f->P0[0][l][i] = 1.0;
      f->P0[1][l][i] = l + 1.0 - x0;
      f->P1[0][l][i] = l - x0;
      f->P1[1][l][i] = x0*x0 - (2.0*l + 3.0)*x0 + l*(l+1.0);

      for(int m=1 ; m<f->getMsize()-1 ; m++){
        f->P0[m+1][l][i] = ((2.0*m + l + 1.0 - x0)*f->P0[m][l][i] - (m+l)*f->P0[m-1][l][i]) / (m+1.0);
        f->P1[m+1][l][i] = (2.0*(m+1.0) + l - x0)*f->P0[m+1][l][i] - 2.0*(m + 1.0 + l)*f->P0[m][l][i];
      }
    }
  }

  /*** normalization constant for Laguerre functions, sqrt{n! / (n+Lambda)!} */
  for(int i=0 ; i<nr ; i++){
    for(int l=0 ; l<=lmax ; l++){
      for(int m=0 ; m<f->getMsize() ; m++){
        double c = sqrt( exp(fact[m] - fact[m+l]) );
        f->P0[m][l][i] *= c;
        f->P1[m][l][i] *= c;
      }
    }

    double x0 = f->x[i];
    f->W[i] = -x0/(nr*nr * f->P0[nr-1][0][i]*(f->P0[nr][0][i] - f->P0[nr-1][0][i]));
    f->G[i] = f->W[i] * exp(x0);
  }

#ifdef DEBUG
  for(int l=0 ; l<=lmax ; l++){
    std::cout << "# Lambda  = " << l << " " << nr << " " << f->getMsize() << std::endl;
    for(int i=0 ; i<nr ; i++){
      std::cout << std::setw(12) << f->x[i];
      for(int n=0 ; n<f->getMsize() ; n++) std::cout << std::setw(12) << f->P0[n][l][i];
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
#endif
}



/**********************************************************/
/*      Zeros of Hermite Polynomials                      */
/**********************************************************/
void HOHermiteRoot(const int n, double *x)
{
  const double eps = 1.0e-4;
  const int maxit = 100;

  double p14 = pow(PI, -0.25);

  double z = 0.0;
  for(int i=0 ; i<(n+1)/2 ; i++){
    if(i == 0){
      z = sqrt(2.0*n+1.0) - 1.85575 * pow((2.0*n+1.0),-1.0/6.0);
    }
    else if(i == 1){
      z = z - 1.14 * pow((double)n, 0.426)/z;
    }
    else if(i == 2){
      z = 1.86 * z - 0.86 * x[0];
    }
    else if(i == 3){
      z = 1.91 * z - 0.91 * x[1];
    }
    else{
      z = 2.0 * z - x[i-2];
    }

    double zt = z + 1.0/eps;
    int niter = 0;
    while( fabs(z - zt) > eps ){
      if(niter++ > maxit) break;

      double p1 = p14, p2 = 0.0, p3 = 0.0;
      for(int j=0 ; j<n ; j++){
        p3 = p2;
        p2 = p1;
        double xj = j+1.0;
        p1 = z * sqrt(2.0/xj) * p2 - sqrt(j/xj) * p3;
      }
      p3 = sqrt(2.0*n) * p2;
      zt = z;
      z  = zt - p1/p3;
    }
    x[i] = z;
    x[n-1-i] = -z;
  }
}


/**********************************************************/
/*      Zeros of Laguerre Polynomials                     */
/**********************************************************/
void HOLaguerreRoot(const int n, double *x, double lambda)
{
  const double eps = 1.0e-8;
  const int maxit = 100;

  double z = 0.0;
  for(int i=0 ; i<n ; i++){
    if(i == 0){
      z = (1.0 + lambda)*(3.0 + 0.92*lambda)/(1.0 + 2.4*n + 1.8*lambda);
    }
    else if(i == 1){
      z = z + (15.0 + 6.25*lambda) / (1.0 + 0.9*lambda + 2.5*n);
    }
    else{
      double p = i - 1.0;
      z = z + ((1.0 + 2.55*p)/(1.9*p)
            + 1.26*p*lambda/(1.0 + 3.5*p)) * (z - x[i-2])/(1.0 + 0.3*lambda);
    }

    double zt = z + 1.0/eps;
    int niter = 0;
    while( fabs(z - zt) > eps ){
      if(niter++ > maxit) break;

      double p1 = 1.0, p2 = 0.0, p3 = 0.0;
      for(int j=1 ; j<=n ; j++){
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*j - 1.0 + lambda - z)*p2 - (j - 1.0 + lambda)*p3)/(double)j;
      }
      p3 = (n*p1 - (n+lambda)*p2)/z;
      zt = z;
      z = zt - p1/p3;
    }
    x[i] = z;
  }
}


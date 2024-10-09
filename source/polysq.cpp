/******************************************************************************/
/**     POLYSQ, Least-Squares Fitting by Polynomial or Legendre              **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "physicalconstant.h"
#include "polysq.h"
#include "terminate.h"

static void   polyDesignMatrixP (const int, const int, double *, double *);
static void   polyDesignMatrixL (const int, const int, double *, double *);
static double legendre (const int, const double);

static const double eps1 = 1.0e-10;
static const double eps2 = 1.0e-03;

/**********************************************************/
/*      Least-Squares Fitting of Polynomials to Data      */
/**********************************************************/
int LSQPolynomial(
  const int n,     // number of data points
  const int m,     // number of parameters (order)
  double *xdata,
  double *ydata,   // input data array (x,y) to be fitted
  double *a)       // output coefficients
{
  double  *f, *v, *x;
  
  v = new double [n*(n+1)/2];
  x = new double [m*(m+1)/2];
  f = new double [n*m];

  /*** design matrix */
  polyDesignMatrixP(n,m,xdata,f);

  /*** covariance matrix */
  for(int i=0 ; i<n ; i++){
    for(int j=0 ; j<=i ; j++)  v[i*(i+1)/2+j] = (i==j) ? 1.0 : 0.0;
  }

  if( (LSQCalc(n,m,ydata,a,v,x,f)) < 0.0 ){
    message << "least-squares equation not solved";
    cohTerminateCode("LSQPolynomial");
  }

  for(int j=0 ; j<m ; j++) if(fabs(a[j]) < eps1) a[j] = 0.0;

  delete [] v;
  delete [] x;
  delete [] f;

  return(0);
}


/**********************************************************/
/*      Least-Squares Fitting of Polynomials to Data      */
/**********************************************************/
int LSQLegendre(
  const bool opt,  // flag to adjust the order automatically
  const int n,     // number of data points
  const int m,     // number of parameters (order)
  double *xdata,
  double *ydata,   // input data array (x,y) to be fitted
  double *a)       // output coefficients
{
  double  *f, *v, *x;
  int     order = 0;
  
  v = new double [n*(n+1)/2];
  x = new double [m*(m+1)/2];
  f = new double [n*m];

  /*** optimize the highest order */
  if(opt){
    for(int mopt=3 ; mopt<=m ; mopt+=2){

      /*** design matrix */
      polyDesignMatrixL(n,mopt,xdata,f);

      /*** covariance matrix, just diagonal */
      for(int i=0 ; i<n ; i++){
        for(int j=0 ; j<=i ; j++)  v[i*(i+1)/2+j] = (i==j) ? 1.0 : 0.0;
      }

      if( (LSQCalc(n,mopt,ydata,a,v,x,f)) < 0.0 ){
        message << "least-squares equation not solved";
        cohTerminateCode("LSQLegendre");
      }

      double chi2 = 0.0;
      for(int i=0 ; i<n ; i++){
        double xx = 0.0;
        for(int j=0 ; j<mopt ; j++) xx += f[i*mopt+j]*a[j];
        if(xx > 0.0) xx = ydata[i]/xx - 1.0;
        chi2 += xx*xx;
      }

      if( chi2 <= eps2 ){
        order = mopt;
        break;
      }
    }
  }

  /*** fixed order */
  else{
    /*** design matrix */
    polyDesignMatrixL(n,m,xdata,f);

    /*** covariance matrix, just diagonal */
    for(int i=0 ; i<n ; i++){
      for(int j=0 ; j<=i ; j++)  v[i*(i+1)/2+j] = (i==j) ? 1.0 : 0.0;
    }

    if( (LSQCalc(n,m,ydata,a,v,x,f)) < 0.0 ){
      message << "least-squares equation not solved";
      cohTerminateCode("LSQLegendre");
    }

    for(int j=0 ; j<m ; j++) if(fabs(a[j]) <= eps1) a[j] = 0.0;
  }
/*
  for(int i=0 ; i<180 ; i++){
    double z = a[0];
    for(int j=1 ; j<m ; j++){
      z += a[j]*legendre(j,(i+1.0));
    }
    cout << setw(12) << i+1 << setw(12) << z <<endl;
  }
*/
  delete [] v;
  delete [] x;
  delete [] f;

  return(order);
}


/**********************************************************/
/*      Least-Squares Fitting to Data                     */
/**********************************************************/
void polyDesignMatrixP(const int n, const int m, double *x, double *f)
{
  for(int i=0 ; i<n ; i++){
    f[i*m] = 1.0;
    for(int j=1 ; j<m ; j++) f[i*m+j] = pow(x[i],j);
  }
}

void polyDesignMatrixL(const int n, const int m, double *x, double *f)
{
  for(int i=0 ; i<n ; i++){
    f[i*m] = 1.0;
    for(int j=1 ; j<m ; j++) f[i*m+j] = legendre(j,x[i]);
  }
}


/**********************************************************/
/*     Legendre Function                                  */
/**********************************************************/
double legendre(const int n, const double t)
{
  double x  = cos(PI*t/180.0);
  double p0 = 1.0;
  double p1 = x;
  double p2 = (3*x*x-1)/2;
  double p3 = 0.0;

  double p  = 0.0;
  if(n == 0)        p=p0;
  else if(n == 1)   p=p1;
  else if(n == 2)   p=p2;
  else{
    for(int i=2 ; i<n ; i++){
      double pn = (double)i;
      p3 = ((2*pn+1)*x*p2-pn*p1)/(pn+1);
      p1=p2;  p2=p3;
    }
    p=p3;
  }
  return(p);
}



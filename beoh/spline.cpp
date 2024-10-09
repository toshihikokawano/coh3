/******************************************************************************/
/*  spline.cpp                                                                */
/*        B-spline interpolation                                              */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "spline.h"


/*************************************************/
/*    Calculate B-Spline Support                 */
/*************************************************/
void splineSupport(const int m, Spline *spl)
{
  double *a = new double [m];
  double *b = new double [m];
  double *p = new double [m];
  double *q = new double [m];

  p[0] = 0.0;
  for(int i=1 ; i<m ; i++) p[i] = spl->GetXdat(i)-spl->GetXdat(i-1);

  a[0]   = 0.0;
  a[m-1] = 1.0;

  b[0]   = 3.0*(spl->GetYdat(1)   - spl->GetYdat(0  ))/p[1];
  b[m-1] = 3.0*(spl->GetYdat(m-1) - spl->GetYdat(m-2))/p[m-1];

  for(int i=1 ; i<m-1 ; i++){
    a[i] = p[i+1]/(p[i]+p[i+1]);
    b[i] = 3.0*((spl->GetYdat(i  )-spl->GetYdat(i-1))*a[i]/p[i])
         + 3.0*((spl->GetYdat(i+1)-spl->GetYdat(i  ))*(1.0-a[i])/p[i+1]);
  }

  p[0]  = 2.0;
  q[0]  = b[0];
  for(int i=1 ; i<m ; i++){
    double c = a[i]/p[i-1];
    p[i] = 2.0  - c*(1.0 - a[i-1]);
    q[i] = b[i] - c*q[i-1];
  }

  spl->SetSpline(m-1,q[m-1]/p[m-1]);
  for(int i=m-2 ; i>=0 ; i--)
    spl->SetSpline(i, (q[i] - (1.0-a[i]) * spl->GetSpline(i+1))/p[i]);

  delete [] a;
  delete [] b;
  delete [] p;
  delete [] q;

  return;
}


/*************************************************/
/*    Interpolation by B-Spline                  */
/*************************************************/
double splineInterpolation(const int n, const double u, Spline *spl)
{
  int i=0, j=0, k=0;

  for(i=0 ; i<n-1 ; i++){
    if( (spl->GetXdat(i) <= u) && (u < spl->GetXdat(i+1)) ){
      j = i;
      k = i+1;
      break;
    }
  }
  if(i == n-1){ k = i; }

  double a = spl->GetXdat(k) - u;
  double b = u - spl->GetXdat(j);
  double h = spl->GetXdat(k) - spl->GetXdat(j);
  double v = (spl->GetSpline(j)*a*a*b     - spl->GetSpline(k)*b*b*a    )/(h*h)
           + (spl->GetYdat(j)*a*a*(2*b+h) + spl->GetYdat(k)*b*b*(2*a+h))/(h*h*h);

  return(v);
}

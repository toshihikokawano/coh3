/******************************************************************************/
/*  etc.cpp                                                                   */
/*        miscellaneous functions                                             */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "physicalconstant.h"
#include "etc.h"


static double expint1(const double);
static double expint2(const double);
static double expint3(const double);

extern  double  *fact;


/**********************************************************/
/*     Laguerre Polynomial Function                       */
/**********************************************************/
double laguerre(const int n, const double a, const double x)
{
  double p2=0.0;

  double p0 = 1.0;       if(n == 0) return(p0);
  double p1 = 1.0+a-x;   if(n == 1) return(p1);

  for(int i=2 ; i<=n ; i++){
    p2=((-x+2*i+a-1.0)*p1-(i+a-1.0)*p0)/(double)i;
    p0=p1; p1=p2;
  }

  return(p2);
}

/**********************************************************/
/*     Gamma Function                                     */
/**********************************************************/
double loggamma(double x)
{
  double v, w;

  v = 1.0;
  while (x < 8) {  v *= x;  x++;  }
  w = 1.0 / (x * x);
  return (((((((-0.0295506535947712  * w + 0.0064102564102564) * w
                -0.0019175269175269) * w + 0.0008417508417508) * w
                +0.0005952380952381) * w + 0.0007936507936508) * w
                -0.0027777777777778) * w + 0.0833333333333333) / x
                + 0.5 * LOG_2PI - log(v) - x + (x - 0.5) * log(x);
}

double gam(const double x)
{
  if(x < 0) return( PI / (sin(PI * x) * exp(loggamma(1 - x))) );
  else      return( exp(loggamma(x)) );
}

/**********************************************************/
/*     Legendre Function                                  */
/**********************************************************/
double legendre(const int n, const double x)
{
  double p = 0.0;

  if(n <= 7){
    double x2 = x*x;
    switch(n){
    case 0: p = 1.0; break;
    case 1: p = x  ; break;
    case 2: p = (3.0*x2 - 1.0) * 0.5;    break;
    case 3: p = x*(5.0*x2 - 3.0) * 0.5;  break;
    case 4: p = (x2*(35.0*x2 - 30.0) + 3.0)/8.0;     break;
    case 5: p = x*(x2*(63.0*x2 - 70.0) + 15.0)/8.0;  break;
    case 6: p = (x2*(x2*(231.0*x2 - 315.0) + 105.0) - 5.0)/16.0;   break;
    case 7: p = x*(x2*(x2*(429.0*x2 - 693.0) + 315.0) - 35.0)/16.0; break;
    default: break;
    }
  }
  else{
    double p1 = x;
    double p2 = (3*x*x-1)/2;
    double p3 = 0.0;
    for(int i=2 ; i<n ; i++){
      double pn = (double)i;
      p3 = ((2*pn+1)*x*p2-pn*p1)/(pn+1);
      p1 = p2;  p2 = p3;
    }
    p = p3;
  }
  return(p);
}

double legendre1(const int n, const double x)
{
  double d = sqrt(1.0-x*x);
  if(n == 0 || d == 0.0) return(0.0);
  else return( n*(-x*legendre(n,x)+legendre(n-1,x))/d );
}

double dlegendre(const int n, const double x)
{
  if(n == 0) return(0.0);
  if(n == 1) return(1.0);
  if(x == 1.0) return(n*(n+1)/2.0);

  double d = n/(x*x-1.0);
  double dp = d*(x*legendre(n,x) - legendre(n-1,x));
  return(dp);
}

/**********************************************************/
/*   Associated Legendre Polynomials   P^m_L              */
/**********************************************************/
double assocLegendrePol(const int l, const int m, const double x)
{
  double p1 = 1.0;
  double p2 = sqrt((1.0-x)*(1.0+x));
  double p3 = 1.0;

  for(int i=1 ; i<=m ; i++){
    p1 *= -p3*p2;
    p3 += 2.0;
  }
  if(l == m) return(p1);

  p2 = x*(2*m+1)*p1;
  if(l == m+1) return(p2);

  for(int i=m+2 ; i<=l ; i++){
    p3 = (x*(2*i-1)*p2-(i+m-1)*p1)/(i-m);  p1 = p2;  p2 = p3;
  }
  return(p3);
}


/**********************************************************/
/*   Spherical Harmonics Real of Y^m_L(theta,phi)         */
/**********************************************************/
double SphericalHarmonics(const int l, const int m, const double theta, const double phi)
{
  double f = 0.0;
  if( (m < -l) || (m > l) ) return(f);

  int    n  = std::abs(m);
  double y1 = sqrt( exp(fact[l-n] - fact[l+n]) * (2.0*l+1.0)/PI4 );
  double y2 = (((m+n)/2)%2 == 0) ? 1.0 : -1.0;
  double y3 = assocLegendrePol(l,n,cos(theta));

  f = y1 * y2 * y3;

  if(m >= 0) f *= cos(m*phi);
  else       f *= sin(-m*phi);

  return(f);
}


/**********************************************************/
/*     Bessel Function I_{n+1/2}(x)                       */
/**********************************************************/
static const double EPSBESS = 1.0e-10;
double bessi2(const int n, const double x)
{
  double bi=0.0,x2,am,an,g1,g2,g3,s;
  int m,k;

  if(x <= 1.0){
    x2=0.25*x*x;
    s=g1=1.0/(n+0.5);
    for(k=1;;k++){
      g1 *= x2/( k * (n+1+k*0.5) );
      s  += g1;
      if(fabs(g1) < EPSBESS) break;
    }
    bi=pow(0.5*x,(double)n+0.5)*s/gam(n+0.5);
  }
  else{
    if(     x <=  10.0) { am = 1.6  *x+11;  an = 1.1  *x+ 4; }
    else if(x <=  40.0) { am = 0.65 *x+18;  an = 0.466*x+11; }
    else if(x <= 100.0) { am = 0.367*x+31;  an = 0.25 *x+21; }
    else                { am = 0.28 *x+39;  an = 0.2  *x+26; }
    m = (int)(am + ( ((double)n>an) ? n-an : 0.0 ));
    x2=2.0/x;
    s=g3=0.0; g2=1.e-36;
    for(k=m;k>=0;k--){
      g1 = (k+0.5)*x2*g2+g3;  if(k == n) bi=g2;
      s += (k+0.5)*g2;
      g3=g2;
      g2=g1;
    }
    bi *= sqrt(x/(2*PI))*exp(x)/s;
  }

  return(bi);
}


/**********************************************************/
/*     Hankel Function K_{n+1/2}(x)                       */
/**********************************************************/
double bessk2(const int n, const double x)
{
  double bk,x2,g,s;
  int k;

  x2=1.0/(2*x);
  s=g=1.0;
  for(k=1;k<=n;k++){
    g *= x2*(n+k)*(n-k+1)/(double)k;
    s += g;
  }
  bk=sqrt(PI*x2)*s*exp(-x);
  return(bk);
}



/**********************************************************/
/*     Exponential Integral Ei(x)                         */
/**********************************************************/
#include<cfloat>
double expint(const double x)
{
  double ei = 0.0;

  if(x == 0.0) return(-DBL_MAX);

  if( (x < -5.0) || (x > 50.0) ) ei = expint1(x);
  else if(x < 6.8)               ei = expint2(x);
  else                           ei = expint3(x);

  return(ei);
}


static double expint1(double x)
{
  const double eps = DBL_EPSILON;
  double am1 = 1.0, a0 = 0.0, bm1 = 0.0, b0 = 1.0;
  double a = exp(x);
  double b = -x + 1.0;
  double ap1 = b * a0 + a * am1;
  double bp1 = b * b0 + a * bm1;

  int j = 1;
  a = 1.0;
  while(fabs(ap1 * b0 - a0 * bp1) > eps*fabs(a0 * bp1)){
    if(fabs(bp1) > 1.0){
      am1 = a0 / bp1;  a0 = ap1 / bp1;  bm1 = b0 / bp1;  b0 = 1.0l;
    }
    else{
      am1 = a0;        a0 = ap1;        bm1 = b0;        b0 = bp1;
    }
    a = -j * j;
    b += 2.0;
    ap1 = b * a0 + a * am1;
    bp1 = b * b0 + a * bm1;
    j += 1;
  }
  return (-ap1 / bp1);
}


static double expint2(const double x)
{ 
  const double eps = DBL_EPSILON;
  double xn = -x;
  double sn = -x;
  double sm1 = 0.0;
  double hsum = 1.0;
  double g = 0.5772156649015328606065121;
  double y = 1.0;
  double f = 1.0;
  
 
  while(fabs(sn - sm1) > eps*fabs(sm1)){
    sm1 = sn;
    y += 1.0;
    xn *= (-x);
    f *= y;
    hsum += (1.0 / y);
    sn += hsum * xn / f;
  }
  return(g + log(fabs(x)) - exp(x) * sn);
}


static double expint3(const double x)
{
  static long double ei[] = {
      1.915047433355013959531e2,  4.403798995348382689974e2,  1.037878290717089587658e3,
      2.492228976241877759138e3,  6.071406374098611507965e3,  1.495953266639752885229e4,
      3.719768849068903560439e4,  9.319251363396537129882e4,  2.349558524907683035782e5,
      5.955609986708370018502e5,  1.516637894042516884433e6,  3.877904330597443502996e6,
      9.950907251046844760026e6,  2.561565266405658882048e7,  6.612718635548492136250e7,
      1.711446713003636684975e8,  4.439663698302712208698e8,  1.154115391849182948287e9,
      3.005950906525548689841e9,  7.842940991898186370453e9,  2.049649711988081236484e10,
      5.364511859231469415605e10, 1.405991957584069047340e11, 3.689732094072741970640e11,
      9.694555759683939661662e11, 2.550043566357786926147e12, 6.714640184076497558707e12,
      1.769803724411626854310e13, 4.669055014466159544500e13, 1.232852079912097685431e14,
      3.257988998672263996790e14, 8.616388199965786544948e14, 2.280446200301902595341e15,
      6.039718263611241578359e15, 1.600664914324504111070e16, 4.244796092136850759368e16,
      1.126348290166966760275e17, 2.990444718632336675058e17, 7.943916035704453771510e17,
      2.111342388647824195000e18, 5.614329680810343111535e18, 1.493630213112993142255e19,
      3.975442747903744836007e19, 1.058563689713169096306e20};

  const double eps = DBL_EPSILON;

  int  k = (int) (x + 0.5);
  int  j = 0;
  double xx = (double)k;
  double dx = x - xx;
  double xxj = xx;
  double edx = exp(dx);
  double sm = 1.0;
  double sn = (edx - 1.0) / xxj;
  double term = DBL_MAX;
  double f = 1.0;
  double dxj = 1.0;

  while(fabs(term) > eps*fabs(sn)){
    j++;
    f *= (double)j;
    xxj *= xx;
    dxj *= -dx;
    sm +=  dxj / f;
    term = (f * (edx * sm - 1.0))/xxj;
    sn += term;
  }
   
  return(ei[k-7] + sn * exp(xx));
}

/**********************************************************/
/*     Exponential Integral E1(x)                         */
/**********************************************************/
std::complex<double> expintE1(std::complex<double> z)
{
  const int nmax = 10;
  std::complex<double> p(1.0,0.0),q(0.0,0.0);
  double f = 1.0;

  for(int k=1 ; k<nmax ; k++){
    p *= -z;
    f *= (double)k;
    q += p/((double)k*f);
  }
  p = -EULER - log(z) - q;

  return(p);
}

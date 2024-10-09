/******************************************************************************/
/*  amplitud.cpp                                                              */
/*        scattering amplitudes                                               */
/******************************************************************************/

#include <cmath>

#include "optical.h"
#include "etc.h"

static inline double abs_square(std::complex<double> x)
{ return( x.real() * x.real() + x.imag() * x.imag() ); }
static inline std::complex<double> omCoulombScattering (const double, const double, std::complex<double>);

static const std::complex<double> ci(0.0,1.0);


/***********************************************************/
/*      Scattering Amplitude for Neutron                   */
/***********************************************************/
double omScatteringAmplitudeN(const int lmax, double *cx, const double angle, const double wave, std::complex <double> *s)
{
  std::complex<double> a(0.0,0.0), b(0.0,0.0);

  for(int l=0 ; l<=lmax ; l++){
    double xl = (double)l;
    double p0 = legendre(l,angle);
    double p1 = legendre1(l,angle);
    a += ((xl+1.0)*(1.0 - s[l*3]) + xl*(1.0 - s[l*3+1]))*p0;
    b += (               -s[l*3]  +           s[l*3+1] )*p1;
  }
  double x = 1.0/(2*wave);
  a *= x * ci;
  b *= x;

  cx[0] = abs_square(a)+abs_square(b);
  cx[1] = -2*(a.real()*b.real() + a.imag()*b.imag())/cx[0];
  cx[2] = 0.0;

  return(cx[0]);
}


/***********************************************************/
/*      Scattering Amplitude for Proton                    */
/***********************************************************/
double omScatteringAmplitudeP(const int lmax, double *cx, const double angle, const double wave, const CoulombData cdata, std::complex<double> *s)
{
  std::complex<double> a(0.0,0.0), b(0.0,0.0);
  std::complex<double> cou(cos(2*cdata.sigma0),sin(2*cdata.sigma0));

  for(int l=0 ; l<=lmax ; l++){
    double xl = (double)l;
    double p0 = legendre(l,angle);
    double p1 = legendre1(l,angle);

    a += cou*((xl+1.0)*(1.0 - s[l*3]) + xl*(1.0 - s[l*3+1]))*p0;
    b += cou*(               -s[l*3]  +           s[l*3+1] )*p1;

    cou = omCoulombScattering(xl+1,cdata.eta,cou);
  }
  double x = 1.0/(2*wave);
  a *= x * ci;
  b *= x;

  std::complex<double> f = omRutherfordAmplitude(angle,wave,cdata);
  a += f;

  cx[0] = abs_square(a)+abs_square(b);
  cx[1] = -2*(a.real()*b.real() + a.imag()*b.imag())/cx[0];
  cx[2] = abs_square(f);

  return(cx[0]);
}


/***********************************************************/
/*      Scattering Amplitude for Alpha Particle            */
/***********************************************************/
double omScatteringAmplitudeA(const int lmax, double *cx, const double angle, const double wave, const CoulombData cdata, std::complex<double> *s)
{
  std::complex<double> a(0.0,0.0);
  std::complex<double> cou(cos(2*cdata.sigma0),sin(2*cdata.sigma0));

  for(int l=0 ; l<=lmax ; l++){
    double xl = (double)l;
    double p0 = legendre(l,angle);

    a += cou*((2*xl+1.0)*(1.0 - s[l*3]))*p0;

    cou = omCoulombScattering(xl+1,cdata.eta,cou);
  }
  double x = 1.0/(2*wave);
  a *= x * ci;

  std::complex<double> f = omRutherfordAmplitude(angle,wave,cdata);
  a += f;

  cx[0] = abs_square(a);
  cx[1] = 0.0;
  cx[2] = abs_square(f);

  return(cx[0]);
}


/***********************************************************/
/*      Scattering Amplitude for Deuteron                  */
/***********************************************************/
double omScatteringAmplitudeD(const int lmax, double *cx, const double angle, const double wave, const CoulombData cdata, std::complex<double> *s)
{
  std::complex<double> a(0.0,0.0), b(0.0,0.0), c(0.0,0.0), d(0.0,0.0), e(0.0,0.0);
  std::complex<double> cou(cos(2*cdata.sigma0),sin(2*cdata.sigma0));

  for(int l=0 ; l<=lmax ; l++){
    double xl = (double)l;
    double xx = xl*(xl+1.0);
    double p0 = legendre(l,angle);
    double p1 = legendre1(l,angle);
    double p2 = (fabs(angle) == 1.0) ? 0.0
              : 2.0*angle/(sqrt(1.0-angle*angle))*p1-l*(l+1.0)*p0;

    a += cou*((xl+1.0)*(1.0 - s[l*3]) + xl*(1.0 - s[l*3+1]))*p0;
    b += cou*( (    xl+2.0)*(1.0 - s[l*3  ])
              +(2.0*xl+1.0)*(1.0 - s[l*3+1])
              +(    xl-1.0)*(1.0 - s[l*3+2]) )*p1;
    if(l > 0){
       c += cou*(s[l*3] - s[l*3+2])*p1;
       d += cou*(  xl*(xl+2.0) * (1.0 - s[l*3  ])
                 -(2.0*xl+1.0) * (1.0 - s[l*3+1])
                 -( xl*xl-1.0) * (1.0 - s[l*3+2]) ) * p1/xx;
    }
    if(l > 1){
       e += cou*(      xl      * (1.0 - s[l*3  ])
                 -(2.0*xl+1.0) * (1.0 - s[l*3+1])
                 +(    xl+1.0) * (1.0 - s[l*3+2]) ) * p2/xx;
    }

    cou = omCoulombScattering(xl+1,cdata.eta,cou);
  }

  double x1 = 1.0/(2*wave);
  double x2 = x1/2.0;
  double x3 = x1/sqrt(2.0);
  a *= x1 * ci;
  b *= x2 * ci;
  c *= x3 * ci;
  d *= x3 * ci;
  e *= x2 * ci;

  std::complex<double> f = omRutherfordAmplitude(angle,wave,cdata);
  a += f;
  b += f;

  cx[0] = (abs_square(a) +2*( abs_square(b)+abs_square(c)
                             +abs_square(d)+abs_square(e) ))/3.0;
  cx[1] = sqrt(8.0) *(  a.imag()*c.real() - a.real()*c.imag()
                      + b.imag()*d.real() - b.real()*d.imag()
                      + d.imag()*e.real() - d.real()*e.imag());
  cx[2] = abs_square(f);

  return(cx[0]);
}


/***********************************************************/
/*      Coulomb Scattering Amplitude                       */
/***********************************************************/
std::complex<double> omCoulombScattering(const double xl1, const double eta, std::complex<double> s)
{
  double g  = xl1*xl1 + eta*eta;
  double a1 = (xl1*xl1 - eta*eta)/g;
  double a2 = (2*eta*xl1)/g;

  std::complex<double> x(a1*s.real() - a2*s.imag(), a1*s.imag() + a2*s.real());

  return(x);
}


/***********************************************************/
/*      Ratherford Scattering                              */
/***********************************************************/
std::complex<double> omRutherfordAmplitude(const double angle, const double wave, const CoulombData cdata)
{
  double x1 = (1.0-angle) / 2.0;
  double x2 = -cdata.eta*log(x1) + 2*cdata.sigma0;
  x1 = -cdata.eta / (2*wave*x1);

  std::complex<double> f(x1*cos(x2),x1*sin(x2));

  return(f);
}

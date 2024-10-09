// miscellaneous functions, complex variable manipulation carried over from old C code

#ifndef __COMPLEX_H__
#define __COMPLEX_H__
#include <complex>
#endif


static inline double cfmin(double x, double y){ return( (x < y) ? x : y); }
static inline double cfmax(double x, double y){ return( (x > y) ? x : y); }


/**************************************/
/*      etc.cpp                       */
/**************************************/
double  laguerre (const int, const double, const double);
double  gam (const double);
double  loggamma (double);
double  legendre (const int, const double);
double  legendre1 (const int, const double);
double  dlegendre (const int, const double);
double  assocLegendrePol (const int, const int, const double);
double  SphericalHarmonics (const int, const int, const double, const double);
double  bessi2 (const int, const double);
double  bessk2 (const int, const double);
double  expint (const double);
std::complex<double> expintE1 (std::complex<double>);


/******************************************************************************/
/*  dsdoverlap.cpp                                                            */
/*        overlap integral <R|f|x>                                            */
/******************************************************************************/

#include <iostream>

#include "structur.h"
#include "optical.h"
#include "dsd.h"

static std::complex<double> dsdVibCouplingFactor (const int, const double, const double, const double);

/**********************************************************/
/*      Direct/Semidirect Formfactors                     */
/**********************************************************/
void    dsdFormfactor(const int c, Optical *omp, Potential *pot, const double v1, const double w1, std::complex<double> *f1, std::complex<double> *f2)
{
  for(int i=0 ; i<=pot->n_match ; i++){
    double r=(double)i*(pot->width);

    /*** direct formfactor,  <R|r|X> */
    f1[i] = std::complex<double>(r,0.0);

    /*** semidirect formfactor,  <R|h(r)|X> */
    std::complex<double> h = dsdVibCouplingFactor(c, r,omp->R0 ,omp->a0);
    f2[i] = std::complex<double>(v1*h.real(),w1*h.imag());
  }
}


/**********************************************************/
/*      Semidirect Coupling Factor                        */
/**********************************************************/
std::complex<double> dsdVibCouplingFactor(const int c, const double x, const double r, const double a)
{
  std::complex<double> f;

  double z = exp((x-r)/a);
  switch(c){
  case 0:    // J.P. Boisson and S.Jang,  df(r)/dr
             // Nucl. Phys. A189, 334 (1972)
    f = std::complex<double>(z/((1+z)*(1+z)),0.0);
    break;
  case 1:    // H. Kitazawa, rf(r)
             // Nucl. Phys. A307, 1 (1978)
    f = std::complex<double>(x/(1+z), 0.0);
    break;
  case 2:    // M. Potokar, rf(r)-i W 4a df(r)/dr
             // Phys. Lett., 46B, 346 (1973)
    f = std::complex<double>(x/(1+z), x*z/((1+z)*(1+z)) *4);
    break;
  default:
    f = std::complex<double>(0.0, 0.0);
  }

  return(f);
}


/**********************************************************/
/*      Overlap Integral                                  */
/*      Distorted Wave and Bound Wave                     */
/**********************************************************/
std::complex<double> dsdRadialIntegral(const int nint, const double width, std::complex<double> *form, double *bw, std::complex<double> *wfn)
{
  std::complex<double> sum;

  double f0=width    /3.0;
  double f1=width*2.0/3.0;
  double f2=width*4.0/3.0;

  sum = wfn[nint] * form[nint] * bw[nint] * f0;

  for(int i=1 ; i<nint ; i++){
    sum += wfn[i] * form[i] * bw[i] * ((i%2 == 0) ? f1 : f2);
  }

  return(sum);
}

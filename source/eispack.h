#ifndef __COMPLEX_H__
#define __COMPLEX_H__
#include <complex>
#endif

int EISPACKUnitaryDiag (const int, double *, double *, double **);
int EISPACKHermiteDiag (const int, std::complex<double> *, double *, std::complex<double> *);

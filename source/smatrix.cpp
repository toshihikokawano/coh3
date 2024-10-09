/******************************************************************************/
/*  smatrix.cpp                                                               */
/*        s matrix elements                                                   */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "optical.h"
#include "etc.h"

static inline std::complex<double>
lagrange(double h, std::complex<double> a,std::complex<double> b, std::complex<double> c,
                   std::complex<double> e,std::complex<double> f, std::complex<double> g)
{
  return( ((g-a)/60.0+0.15*(b-f)+0.75*(e-c))/h );
}


/***********************************************************/
/*      S Matrix Element for Spherical Calculation         */
/***********************************************************/
std::complex<double> omSmatrix(const int m, const double mesh, const double rm, std::complex<double> d1, std::complex<double> d2, Wavefunc *wfn)
{
  std::complex<double> win, dwin, f, u1, u2, s;

  win  = wfn->internal[m-3];
  dwin = lagrange(mesh,wfn->internal[m-6],wfn->internal[m-5],
                       wfn->internal[m-4],wfn->internal[m-2],
                       wfn->internal[m-1],wfn->internal[m  ]);
  f = rm * dwin/win;

  u1 = std::complex<double>(f.real() - d1.real(), f.imag() + d1.imag());
  u2 = std::complex<double>(f.real() - d1.real(), f.imag() - d1.imag());

  s = u1 / u2 * d2;

  return(s);
}


/***********************************************************/
/*      S and C Matrix Elements for CC Calculation         */
/***********************************************************/
int ccSmatrix(const int m, const double mesh, CCdata *cdt, Wavefunc *wfn, std::complex<double> *d, std::complex<double> *s, std::complex<double> *c)
{
  std::complex<double> win, dwin, h, dh;

  for(int i=0 ; i<m ; i++){
    int k  = cdt[i].level;
    int l  = cdt[i].chn.l;
    h  = wfn[k].external[l];
    dh = wfn[k].extderiv[l];
    double  x = sqrt(fabs(cdt[i].lev->wavesq));

    dh *= x;
    x *= -1;

    for(int j=0 ; j<m ; j++){
      int ij = i*m+j;
      int ji = j*m+i;

      win  = wfn[ij].internal[3];
      dwin = lagrange(mesh,wfn[ij].internal[0],wfn[ij].internal[1],
                           wfn[ij].internal[2],wfn[ij].internal[4],
                           wfn[ij].internal[5],wfn[ij].internal[6]);

      double ar = (win.real()*dh.real() - dwin.real()*h.real())/x;
      double ai = (win.imag()*dh.real() - dwin.imag()*h.real())/x;

      double sr = (win.real()*dh.imag() - dwin.real()*h.imag())/x;
      double si = (win.imag()*dh.imag() - dwin.imag()*h.imag())/x;

      s[ji] = std::complex<double>(sr,si);
      d[ji] = std::complex<double>(-ar+si,-ai-sr);
    }
  }

  for(int k=0 ; k<m-1 ; k++){
    int    kk = k*m+k;
    double x  = d[kk].real()*d[kk].real() + d[kk].imag()*d[kk].imag();
    if(x == 0.0) return(-1);

    for(int i=k+1 ; i<m ; i++){
      int ik = i*m+k;
      double er, ei, fr, fi;

      fr = (d[ik].real()*d[kk].real() + d[ik].imag()*d[kk].imag())/x;
      fi = (d[ik].imag()*d[kk].real() - d[ik].real()*d[kk].imag())/x;
       
      for(int j=k ; j<m ; j++){
        int ij = i*m+j;
        int kj = k*m+j;
        er = d[ij].real() - (fr*d[kj].real() - fi*d[kj].imag());
        ei = d[ij].imag() - (fr*d[kj].imag() + fi*d[kj].real());
        d[ij] = std::complex<double>(er,ei);
      }

      for(int j=0 ; j<m ; j++){
        int ij = i*m+j;
        int kj = k*m+j;
        er = s[ij].real() - (fr*s[kj].real() - fi*s[kj].imag());
        ei = s[ij].imag() - (fr*s[kj].imag() + fi*s[kj].real());
        s[ij] = std::complex<double>(er,ei);
      }
    }
  }

  int    kk =(m-1)*m+m-1;
  double x  = d[kk].real()*d[kk].real() + d[kk].imag()*d[kk].imag();
  if(x == 0.0) return(-1);

  for(int j=0 ; j<m ; j++){
    int kj =(m-1)*m+j;
    double cr = (s[kj].real()*d[kk].real() + s[kj].imag()*d[kk].imag())/x;
    double ci = (s[kj].imag()*d[kk].real() - s[kj].real()*d[kk].imag())/x;
    c[kj] = std::complex<double>(cr,ci);
  }

  for(int i=m-2 ; i>=0 ; i--){
    int    ii = i*m+i;
    double x  = d[ii].real()*d[ii].real() + d[ii].imag()*d[ii].imag();
    if(x == 0.0) return(-1);

    for(int j=0 ; j<m ; j++){ 
      int ij = i*m+j;
      double er = 0.0, ei = 0.0, fr = 0.0, fi = 0.0;
      for(int k=i+1 ; k<m ; k++){
        int ik = i*m+k;
        int kj = k*m+j;
        er += d[ik].real()*c[kj].real() - d[ik].imag()*c[kj].imag();
        ei += d[ik].real()*c[kj].imag() + d[ik].imag()*c[kj].real();
      }
      fr = s[ij].real() - er;
      fi = s[ij].imag() - ei;

      double cr = (fr*d[ii].real() + fi*d[ii].imag())/x;
      double ci = (fi*d[ii].real() - fr*d[ii].imag())/x;
      c[ij] = std::complex<double>(cr,ci);
    }
  }

  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++ ){
      int ij = i*m+j;
      double x = pow(fabs(cdt[j].lev->wavesq)/fabs(cdt[i].lev->wavesq),0.25);
      c[ij] *= x;
      s[ij] = std::complex<double>(((i==j) ? 1 : 0) - 2*c[ij].imag(),2*c[ij].real());
    }
  }

  return(0);
}

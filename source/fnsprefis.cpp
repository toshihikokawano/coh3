/******************************************************************************/
/*  fnprefis.cpp                                                              */
/*        Neutron Evaporation before Fission                                  */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "fns.h"

/**********************************************************/
/*      PreFission Neutron Spectum                        */
/*      ----                                              */
/*          assume neutron spectra are the same as        */
/*          the neutron emission from CN(A -1)            */
/*          times fission probability at each bin         */
/**********************************************************/
double fnsPreFissionSpectrum(const int ne, double *eout, double *pref, double *w,  Nucleus *n)
{
  if(n->ncont <= 1) return(0.0);

  /*** if only 2 points, connect zero and the first point linearly */
  if(n->ncont == 2){
    if(pref[1] == 0.0) return(0.0);

    for(int i=0 ; i<ne ; i++){
      w[i] = (eout[i] > n->de) ? 0.0 : pref[1] / n->de * eout[i];
    }
  }
  /*** general case */
  else{

    if( (pref[1] == 0.0) || (pref[2] == 0.0) ) return(0.0);

    /*** determine evaporation parameters at low emission energies */
    double tp = n->de / (log(pref[1]/pref[2]) - log(0.5));
    double cp = pref[1]/n->de * exp( n->de/tp );

    /*** energy point interpolation */
    for(int i=0 ; i<ne ; i++){

      /*** first bin, below 2 dE */
      if(eout[i] <= 2*n->de){
        /*** use evaporation form for spec */
        w[i] = cp*eout[i]*exp(-eout[i]/tp);
      }
      /*** set zero above Emax */
      else if(eout[i] >= n->max_energy){
        w[i] = 0.0;
      }
      /*** linear interpolation */
      else{
        int k = 0;
        for(k=1 ; k<n->ncont ; k++){
          if(n->de*k > eout[i]) break;
        }
        double p0 = pref[k-1];
        double p1 = pref[k];
        w[i] = (p1 - p0)/n->de * (eout[i] - n->de*(k-1)) + p0;
        if(w[i] < 0.0) w[i] = 0.0;
      }
    }
  }

  /*** save pre-fission spectrum in pref array  */
  for(int i=0 ; i<MAX_ENERGY_BIN ; i++){
    pref[i] = (i < ne) ? w[i] : 0.0;
  }

  fnsNormalizeSpectrum(ne,eout,pref);
  double ep = fnsAverageSpectrumEnergy(ne,eout,pref);

  return(ep);
}

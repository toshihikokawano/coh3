/******************************************************************************/
/*  fnssystematics.cpp                                                        */
/*        phenomenological model pameters for fission calculation             */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "fnssystematics.h"
#include "masstable.h"

static const int TKEmodel = 0; // Viola for TKE as default


/**********************************************************/
/*      <Nu> Calculation                                  */
/**********************************************************/
double fnsNubar(const double txe, const double ecms, const double esep, const double eg)
{
  double nu = (txe - eg) / (esep + ecms);
  if(nu < 0.0) nu = 0.0;

  return(nu);
}


/**********************************************************/
/*      Simple TKE Estimates                              */
/**********************************************************/
double fnsTKESystematics(const int zc, const int ac)
{
  double tke = 0.0;

  double a13 = pow((double)ac,1/3.0);
  double c   = zc * zc / a13;

  /*** Proc. 3rd IAEA Symp. Phys. and Chem. of Fission, Vol 2., p 19 */
  if(TKEmodel == 1){
    tke = 0.13323  * c - 11.64;
  }
  /*** V.E. Viola Jr, Nucl. Data Tables A1, 391 (1966), also given in Viola (1985) */
  else if (TKEmodel == 2){
    tke = 0.10710  * c + 22.2;
  }
  /*** M.L. Walker and V.E. Viola Jr, Indiana Nuclear Chemistry Report, No. INC-40007-6 (1981), ditto */
  else if(TKEmodel == 3){
    tke = 0.1166   * c + 9.0;
  }

  /*** V.E. Viola, K. Kwiatkowski, M. Walker,  Phys. Rev. C 31, 1550 (1985) */
  else{
    tke = 0.1189   * c + 7.3;
  }

  return(tke);
}


/**********************************************************/
/*      Rough Estimate of TKE Slope from 235U data        */
/**********************************************************/
double fnsTKEEdepSystematics(const int ac)
{
  double dtke = -0.18 - (236 - ac) * 0.02;
  return(dtke);
}


/**********************************************************/
/*      Simple TKE Energy Dependence                      */
/**********************************************************/
double fnsTKEEnergyDependence(double tke0, double tke1, double sn, const int n, double *ex, double *pf)
{
  double tke = tke0;

  double x0 = 0.0, x1 = 0.0;
  for(int i=0 ; i<n ; i++){
    tke = (ex[i] < sn) ? tke0 : tke0 + tke1 * (ex[i] - sn);
    if(tke < 0.0) tke = 0.0;

    x0 += pf[i];
    x1 += pf[i] * tke;
  }
  if(x0 > 0.0) tke = x1 / x0;
  else tke = tke0 + tke1 * ex[0];

  return(tke);
}


/*** parameters given in CGMF */
double fnsTKEEnergyDependenceCGMF(const int zc, const int ac, double sn, const int n, double *ex, double *pf)
{
  double tke = 0.0;
  double a=0.0, e=0.0, b=0.0, d=0.0;

  if(zc != 92) tke = fnsTKESystematics(zc,ac);
  else{
    if(ac == 233){      a = 171.23; e = 0.75; b = 0.5000; d = -0.4500;  }
    else if(ac == 234){ a = 170.65; e = 0.60; b = 0.6218; d = -0.0557;  }
    else if(ac == 235){ a = 170.09; e = 0.90; b = 0.7559; d = -0.2371;  }
    else if(ac == 236){ a = 171.74; e = 0.75; b = 0.7181; d = -0.0975;  }

    double c = a + (b - d)*e;

    tke = a;
    double x0 = 0.0, x1 = 0.0;
    for(int i=0 ; i<n ; i++){
      double en = ex[i] - sn;

      if(en < 0.0)    tke = a;
      else if(en < e) tke = a + b * en;
      else            tke = c + d * en;
      
      x0 += pf[i];
      x1 += pf[i] * tke;
    }
    if(x0 > 0.0) tke = x1 / x0;
  }

  return(tke);
}


/**********************************************************/
/*      Simple Egamma Estimate                            */
/**********************************************************/
double fnsEgammaSystematics(const int ac)
{
  /*** prompt gamma-energy systematics by E.Fort, 1990 IAEA Vienna Mtg. */
  double eg = 0.02772*ac + 0.0891;
  return(eg);
}



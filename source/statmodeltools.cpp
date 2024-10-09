/******************************************************************************/
/*  statmodeltools.cpp                                                        */
/*        some utilitiy functions to calcualte statistical model              */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "nucleus.h"
#include "statmodel.h"
#include "global.h"


/**********************************************************/
/*      Find a Bin Index for Given Excitation             */
/**********************************************************/
int specFindEnergyBin(double e, const double de)
{
  double eps = 1.0e-5;

  if(fabs(e) <= eps) return(0);
  if( (e < 0.0) || (e > (MAX_ENERGY_BIN-1)*de) ) return(-1);

  if(e < de*0.5) return(0);

  /*** Add first energy bin, which has a half width,
       so that zero-energy emission is mapped onto k=0 */ 
  e += 0.5*de;

  /*** To avoid round-off error, a small number is added */
  int k = (int)((e+eps)/de);
  if( (de*k <= (e+eps)) && (e < de*(k+1)) ) return (k);
  if( fabs(de*k - (e+eps)) < 1e-10 ) return (k);

  std::cout << " out of range " << k << " " << de*k <<  " " << e << " " << e+eps << " " << de*(k+1) << " " << de*k - (e+eps) << std::endl;
  return(-1);
}


/**********************************************************/
/*      Sum Up All Particle Emission Spectra              */
/**********************************************************/
void specCumulativeSpectra(const int csize, const int nsize, double **spc, Nucleus *n)
{
  int cs = crx.getNspectra();
  int ns = crx.getNspecbin();

  for(int id=0 ; id<csize ; id++){
    if(n->cdt[id].status){
      for(int k=0 ; k<nsize ; k++){
        if((id < cs) && (k < ns)) crx.spectra[id][k] += spc[id][k];
      }
    }
  }
}


/**********************************************************/
/*      Gaussian Broadening Results by Given Resolution   */
/**********************************************************/
void specGaussianBroadening(const int m, const double de, const double dm, double *y0)
{
  double *y1 = new double [m];

  int d = 3 * (int)(dm / de) + 1;

  for(int i0=0 ; i0<m ; i0++){

    double z = 0.0;
    y1[i0] = 0.0;
    for(int i1=i0-d ; i1<=i0+d ; i1++){
      if((i1 < 0) || (i1 >= m)) continue;
      double e0 = i0 * de;
      double e1 = i1 * de;
      double w = exp(-(e0-e1)*(e0-e1)/(2*dm*dm));
      z += w;
      y1[i0] += w * y0[i1];
    }
    y1[i0] = (z > 0.0) ? y1[i0]/z : 0.0;
  }

  for(int i=0 ; i<m ; i++) y0[i] = y1[i];

  delete [] y1;
}


/**********************************************************/
/*      Add Gaussian Broadened Line to Continuum          */
/**********************************************************/
void specGaussianBroadening(const int m, const double de, const double dm, double *y0, const double x1, const double y1)
{
  int d = 3 * (int)(dm / de) + 1;
  double *y2 = new double [2*d+1];

  int i0 = -1;
  for(int i1=0 ; i1<m-1 ; i1++){
    double e0 = i1 * de;
    double e1 = e0 + de;
    if((e0 <= x1) && (x1 < e1)){
      i0 = i1;
      break;
    }
  }
  if(i0 < 0) return;

  double z = 0.0;
  for(int k = 0 ; k<=2*d ; k++){
    y2[k] = 0.0;

    int i1 = i0 - d + k;
    if((i1 < 0) || (i1 >= m)) continue;

    double e1 = i1 * de;
    double w = exp(-(x1-e1)*(x1-e1)/(2*dm*dm));
    z += w;
    y2[k] = w * y1;
  }

  if(z > 0.0){
    z *= de;
    for(int k = 0 ; k<=2*d ; k++) y2[k] = y2[k] / z;
  }

  for(int k = 0 ; k<=2*d ; k++) y0[i0-d+k] += y2[k];

  delete [] y2;
}


/**********************************************************/
/*      Primary Gamma-Ray Spectra                         */
/*      -----                                             */
/*             The primary gamma-ray emission cross       */
/*             section should be the same as the level    */
/*             population for the C0=0, K0=0 case         */
/**********************************************************/
int specStorePrimaryGamma(Nucleus *n, double **pg, GammaProduction *gp)
{
  int ng = 0;

  /*** copy non-zero production cross section to PG array */
  for(int i=0 ; i<n->ndisc ; i++){
    double eg = n->excitation[0] - n->lev[i].energy;
    if(n->lpop[i] > 0.0){
      pg[0][ng] = eg;
      pg[1][ng] = n->lpop[i];
      ng++;
    }
  }

  /*** primary gamma lines distinguished by negative energy */
  for(int i=0 ; i<ng ; i++) gp->push(n->za,-pg[0][i],pg[1][i]);
  
  return gp->getN();
}


/**********************************************************/
/*      Store Cascading Gamma Lines                       */
/**********************************************************/
int specStoreDiscreteGamma(Nucleus *n, GammaProduction *gp)
{
  for(int i0=n->ndisc-1 ; i0>0 ; i0--){
    for(int j=0 ; j<n->lev[i0].ngamma ; j++){
      int i1 = n->lev[i0].fstate[j];

      /*** discrete gamma-rays including internal conversion factor */
      double p = n->lev[i0].branch[j] * n->lpop[i0];
      if(opt.internalconversion) p *= n->lev[i0].gratio[j];
      gp->push(n->za,n->lev[i0].energy-n->lev[i1].energy,p);
    }
  }

  return gp->getN();
}


/**********************************************************/
/*      Sort Discrete Gamma Lines                         */
/**********************************************************/
void specSortDiscreteGamma(GammaProduction *gp)
{
  GammaLine tmp;

  for(int j=0 ; j<gp->getN() ; j++){
    int l = j;
    for(int i=j ; i<gp->getN() ; i++){
      if(gp->line[i].energy < gp->line[l].energy) l = i;
    }
    tmp = gp->line[j];
    gp->line[j] = gp->line[l];
    gp->line[l] = tmp;
  }
}

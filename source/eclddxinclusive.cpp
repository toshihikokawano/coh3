/******************************************************************************/
/*  eclddxinclusive.cpp                                                       */
/*        calculate inclusive double differential cross sections              */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "physicalconstant.h"
#include "structur.h"
#include "nucleus.h"
#include "statmodel.h"
#include "eclipse.h"
#include "ddxkalbach.h"
#include "parameter.h"
#include "etc.h"
#include "terminate.h"

static void   ddxMainInclusiveSpectrum (const int, const int, const int, const double, double **);
static int    ddxSumDiscreteTransition (const int, const int, const double, double *);

static double    **cleg;      // Legendre expansion coefficients
static double    **fleg;      // pre-calculated Legendre polynomials
static double     *epar;      // secondary particle energies
static double     *xang;      // secondary particle angles

static const double DANG = 30.0; // calculation angle step
static int nang = 0;


/**********************************************************/
/*      Inclusive Double Differential Cross Sections      */
/**********************************************************/
void eclDDXInclusive(const int km, System *sys, double **fmsd)
{
  nang = (int)(180.0 / DANG);
  if(sys->incident.pid == neutron){
    nang ++; // include zero degree in the case of neutron
  }

  /*** bin number for elastic when neutron induced */
  int kel = -1;
  if(sys->incident.pid == neutron){
    kel = specFindEnergyBin(ncl[sys->target_id].max_energy,sys->energy_bin);
  }

  try{
    /*** cleg[LegOrder][Energy], fleg[LegOrder][Angle] */
    epar = new double [NCNT];
    xang = new double [nang];
    cleg = new double * [MAX_J];
    fleg = new double * [MAX_J];
    for(int l=0 ; l<MAX_J ; l++){
      cleg[l] = new double [km];
      fleg[l] = new double [nang];
    }
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("eclDDXinclusive");
  }

  /*** secondary particle energy */
  for(int k=0 ; k<NCNT ; k++) epar[k] = (k + 0.5) * sys->energy_bin;

  /*** secondary particle angle, and pre-calculated Legendre function */
  for(int n=0 ; n<nang ; n++){
    if(sys->incident.pid == neutron) xang[n] = n * DANG;
    else xang[n] = (n+1) * DANG;

    double th = cos( (double)xang[n] * PI / 180.0 );
    for(int l=0 ; l<MAX_J ; l++){
      fleg[l][n] = legendre(l,th);
    }
  }

  /*** Generate DDX from Kalbach systematics */
  ddxKalbachSetParm(sys->cms_energy,sys->compound,sys->incident.za);
  ddxMainInclusiveSpectrum(km,sys->max_channel,kel,sys->energy_bin,fmsd);


  for(int l=0 ; l<MAX_J ; l++){
    delete [] cleg[l];
    delete [] fleg[l];
  }
  delete [] cleg;
  delete [] fleg;
  delete [] epar;
  delete [] xang;
}


/**********************************************************/
/*      Main Part for Inclusive DDX                       */
/**********************************************************/
void ddxMainInclusiveSpectrum(const int km, const int cm, const int kel, const double de, double **fmsd)
{
  double *ctmp = new double [MAX_J];
  double *espc = new double [km];
  int    *lmax = new int [km];

  for(int c=0 ; c<cm ; c++){

    if(!ncl[0].cdt[c].status) continue;
    int id = ncl[0].cdt[c].next;

    /*** find non-zero element from the high-side */
    int kmax = km;
    for(int k=km-1 ; k>=0 ; k--){
      if(crx.spectra[c][k] > 0.0){
        kmax = k + 2;
        break;
      }
    }

    if(kmax <= 0) continue;

    for(int k=0 ; k<km ; k++){
      espc[k] = 0.0;
      for(int l=0 ; l<MAX_J ; l++) cleg[l][k] = 0.0;
    }

    /*** for each energy bin */
    for(int k=0 ; k<kmax-1 ; k++){

      if(c == 0) lmax[k] = 1;  // for photons, no angular distribution
      else{
        /*** continuum spectrum, use Kalback systematics */
        if(k < ncl[id].ncont){
          ddxKalbach(c,epar[k],fmsd[c][k],NLEG,ctmp);
          for(int l=0 ; l<NLEG ; l++) cleg[l][k] = 2 * ctmp[l];
          lmax[k] = NLEG;
        }
        /*** replace by discrete angular distributions */
        else{
          lmax[k] = ddxSumDiscreteTransition(c,k,de,ctmp);
          if(lmax[k] > 0){
            for(int l=0 ; l<lmax[k] ; l++) cleg[l][k] = ctmp[l];
          }
        }
      }
      cleg[0][k] = 1.0;
      espc[k] = crx.spectra[c][k] / PI4;

      /*** add shape elastic scattering ***/
      if(k == kel) espc[k] += crx.elastic / PI4;
    }

    if(c != 0){
      /*** Gaussian broadening */
      double gw = parmGetValue(parmBROD);
      if(gw > 0.0){
        int kmax2 = specGaussianBroadening(km,de,gw,espc);

        /*** copy legendre coefficients at the highest bin */
        for(int k=kmax-1 ; k<kmax2 ; k++){
          lmax[k] = lmax[kmax-2];
          for(int l=0 ; l<lmax[k] ; l++) cleg[l][k] = cleg[l][kmax-2];
        }
        kmax = kmax2 + 1;
      }
    }

    /*** output DDX */
    eclOutInclusiveSpectrum(c,kmax,nang,xang,lmax,de,cleg,fleg,espc);
  }

  delete [] espc;
  delete [] ctmp;
  delete [] lmax;
}


/**********************************************************/
/*      Legendre Coefficient for Discrete Transitions      */
/**********************************************************/
int ddxSumDiscreteTransition(const int c, int k, const double de, double *ctmp)
{
  int id = ncl[0].cdt[c].next;

  for(int j=0 ; j<MAX_J ; j++) ctmp[j] = 0.0;

  /*** look for discrete levels within the energy bin, and retrieve their Legendre coefficients */
  int nl = 0;
  for(int i=0 ; i<ncl[id].ndisc ; i++){

    double ep = ncl[id].max_energy - ncl[id].lev[i].energy;
    int    kp = specFindEnergyBin(ep,de);

    if(kp == k){
      for(int j=0 ; j<MAX_J ; j++) ctmp[j] += crx.legcoef[id][i][j];
      nl ++;
    }
  }

  int lmax = 1;

  if(nl > 0){
    for(int j=MAX_J-1 ; j>=0 ; j--){
      if(ctmp[j] != 0.0){
        lmax = j+1;
        break;
      }
    }

    for(int j=1 ; j<lmax ; j++){
      ctmp[j] = ctmp[j] / ctmp[0];
    }
  }
  ctmp[0] = 1.0;
  

  return lmax;
}



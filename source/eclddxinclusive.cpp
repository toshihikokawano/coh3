/******************************************************************************/
/*  eclddxinclusive.cpp                                                       */
/*        calculate inclusive double differential cross sections              */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "physicalconstant.h"
#include "structur.h"
#include "nucleus.h"
#include "global.h"
#include "statmodel.h"
#include "eclipse.h"
#include "ddxkalbach.h"
#include "parameter.h"
#include "etc.h"
#include "terminate.h"

static void   ddxMainInclusiveSpectrum (const int, System *, const int, double **);
static int    ddxSumDiscreteTransition (const int, const int, const double, double *);
static int    ddxEnergyArray (const int, const double, const double, double *);

static double    **cleg;      // Legendre expansion coefficients
static double    **fleg;      // pre-calculated Legendre polynomials
static double     *epar;      // secondary particle energies
static double     *spec;      // secondary particle normalized spectra
static double     *xang;      // secondary particle angles

static const double DANG = 15.0; // calculation angle step
static int nang = 0;
static double sigR = 0.0, dummy[7];

/**********************************************************/
/*      Inclusive Double Differential Cross Sections      */
/**********************************************************/
void eclDDXInclusive(System *sys, double **fmsd)
{
  nang = (int)(180.0 / DANG);
  if(sys->incident.pid == neutron){
    nang ++; // include zero degree in the case of neutron
  }

  int km = ncl[0].ntotal + 1; // include zero and highest points

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
    spec = new double [km];
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

  sigR = 0.0;
  for(int i=0 ; i<sys->max_compound ; i++) if(crx.prod[i].xsec > 0.0) sigR += crx.prod[i].xsec;
  if(sigR == 0.0) sigR = crx.reaction;

  if(ctl.endfspec){
    eclOutHead(sys->max_compound,sys->max_channel, sys->lab_energy);
    for(int c=0 ; c<7 ; c++) dummy[c] = 0.0;
    eclOutNucleusHead(0,sys->max_channel,dummy);
  }

  /*** Generate DDX from Kalbach systematics */
  ddxKalbachSetParm(sys->cms_energy,sys->compound,sys->incident.za);
  ddxMainInclusiveSpectrum(km,sys,kel,fmsd);


  for(int l=0 ; l<MAX_J ; l++){
    delete [] cleg[l];
    delete [] fleg[l];
  }
  delete [] cleg;
  delete [] fleg;
  delete [] epar;
  delete [] xang;
  delete [] spec;
}


/**********************************************************/
/*      Main Part for Inclusive DDX                       */
/**********************************************************/
void ddxMainInclusiveSpectrum(const int km, System *sys, const int kel, double **fmsd)
{
  double *ctmp = new double [MAX_J];
  double *espc = new double [km];
  int    *lmax = new int [km];

  for(int c=0 ; c<sys->max_channel ; c++){

    if(!ncl[0].cdt[c].status) continue;
    int id = ncl[0].cdt[c].next;

    /*** find non-zero element from the high-side */
    int kmax1 = km;
    for(int k=km-1 ; k>=0 ; k--){
      if(crx.spectra[c][k] > 0.0){
        kmax1 = k + 2;
        break;
      }
    }
    if(kmax1 <= 0) continue;

    for(int k=0 ; k<km ; k++){
      espc[k] = 0.0;
      for(int l=0 ; l<MAX_J ; l++) cleg[l][k] = 0.0;
    }

    if(ctl.endfspec){
      /*** particle multiplicity */
      double sum = 0.0;
      for(int k=0 ; k<km ; k++) sum += crx.spectra[c][k] * sys->energy_bin;
      for(int k=0 ; k<km ; k++){
        if(sum > 0.0) spec[k] = crx.spectra[c][k] / sum;
        else spec[k] = 0.0;
      }

      double np = sum/sigR;
      eclOutChannelHead(c,np,false);
      if(np == 0.0) continue;

      double epmax = ncl[c].max_energy;
      int kmax2 = ddxEnergyArray(kmax1,epmax,sys->energy_bin,epar);

      int nleg = (c == 0) ? 1 : NLEG;

      for(int k=1 ; k<kmax2 ; k++){
        if(c == 0) cleg[0][k] = spec[k];
        else{
          if(k < ncl[id].ncont) ddxKalbach(c,epar[k],fmsd[c][k-1],NLEG,ctmp);
          else{
            int lmax = ddxSumDiscreteTransition(c,k,sys->energy_bin,ctmp);
            for(int l=0 ; l<lmax ; l++) ctmp[l] *= 0.5;
          }

          /*** convert into ENDF-6 data, F(E') = 2 f/(2L+1) x (ds/dE) */
          /*** note that discrete contribution might have higher L component,
               but we cut the max L to NLEG */
          for(int l=0 ; l<nleg ; l++) cleg[l][k] = ctmp[l] / (l+0.5) * spec[k];
        }
      }

      if(c == 0) eclOutGammaLine(0,0,&dummy[0],&dummy[1]);

      if(kmax2 > 0) eclOutLegCoeff(kmax2,nleg,epar,cleg);
      else          eclOutLegCoeff(-1,0,epar,cleg);
    }
    else{
      /*** for each energy bin */
      for(int k=0 ; k<kmax1-1 ; k++){

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
            lmax[k] = ddxSumDiscreteTransition(c,k,sys->energy_bin,ctmp);
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

      /*** output DDX */
      if(c != 0){
        /*** Gaussian broadening */
        double gw = parmGetValue(parmBROD);
        if(gw > 0.0){
          int kmax2 = specGaussianBroadening(km,sys->energy_bin,gw,espc);

          /*** copy legendre coefficients at the highest bin */
          for(int k=kmax1-1 ; k<kmax2 ; k++){
            lmax[k] = lmax[kmax1-2];
            for(int l=0 ; l<lmax[k] ; l++) cleg[l][k] = cleg[l][kmax1-2];
          }
          kmax1 = kmax2 + 1;
        }
      }
      eclOutInclusiveSpectrum(c,kmax1,nang,xang,lmax,sys->energy_bin,cleg,fleg,espc);
    }
  }

  delete [] espc;
  delete [] ctmp;
  delete [] lmax;
}


/**********************************************************/
/*      Set Up Outgoing Secondary Particle Energies       */
/**********************************************************/
int ddxEnergyArray(const int km, const double epmax, const double de, double *ep)
{
  ep[0] = 0.0;
  ep[1] = de/4.0;

  /*** generate energy bins for spectrum */
  bool emx = true;
  int  km1 = km + 1;
  for(int k=1 ; k<km ; k++){
    ep[k+1] = de*k;
    if(ep[k+1] >= epmax){
      ep[k+1] = epmax;
      km1 = k+2;
      emx = false;
      break;
    }
  }
  /*** if km range does not cover the Epmax, add the max point */
  if(emx) ep[km1++] = epmax;

  return(km1);
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
    if(i >= MAX_ANGDISTLEVELS) continue;

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



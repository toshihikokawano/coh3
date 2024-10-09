/******************************************************************************/
/*  spectra.cpp                                                               */
/*        main particle emission spectra calculation                          */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "beoh.h"
#include "statmodel.h"
#include "beohstat.h"
#include "nucleus.h"
#include "global.h"

static void specCmsToLab(const int, const double, const double, double *, double *);

/**********************************************************/
/*      Decay of each compound nucleus                    */
/**********************************************************/
void    beohspectra(System *sys, Pdata *pdt, Transmission **tc, Transmission **td, double **tg, Spectra *spc, GammaProduction *gml)
{
  /*** Loop over each daughter compound nucleus */
  for(int c0=0 ; c0<sys->max_compound ; c0++){

    /*** Store particle transmission data into array */
    statStoreContinuumTransmission(c0,ncl[c0].excitation[0],pdt,tc);

    /*** Clear spectrum array */
    spc->memclear("cn");

    /*** Loop over parent excitation */
    for(int k0=0 ; k0<ncl[c0].ncont ; k0++){

      /*** Store gamma-ray transmission data into array */
      statStoreContinuumGammaTransmission(k0,tg,&ncl[c0]);

      /*** Transmission coefficients for discrete levels */
      statStoreDiscreteTransmission(c0,ncl[c0].excitation[k0],pdt,td);

      /*** All residual CN decay */
      specCompoundDecay(c0,k0,tc,td,tg,spc);
    }

    /*** Particle decay from levels if possible */
    specLevelDecay(c0,pdt,&ncl[c0],tc,td,spc);

    /*** Gamma-ray cascading */
    specGammaCascade(spc->cn[0],&ncl[c0]);

    /*** Ground state population + meta states when no gamma decay */
    crx.prod[c0].xsec += ncl[c0].lpop[0];
    for(int i=1 ; i<ncl[c0].ndisc ; i++){
      if(ncl[c0].lev[i].ngamma == 0) crx.prod[c0].xsec += ncl[c0].lpop[i];
    }

    /*** Add to total particle emission spectra */
    specCumulativeSpectra(spc->getCsize(),spc->getNsize(),spc->cn,&ncl[c0]);

    /*** Store discrete gamma lines */
    if(opt.finegammaspectrum) specStoreDiscreteGamma(&ncl[c0],gml);
  }
}


/**********************************************************/
/*      Convert emission spectrum into the LAB frame      */
/**********************************************************/
void    beohspectraLAB(System *sys, Pdata *pdt,
                       Transmission **tc, Transmission **td, double **tg, Spectra *spc, double exmean, double ekmean)
{
  /*** Use pe as a temporal LAB spectrum strage */
  spc->memclear("pe");

  /*** Loop over each daughter compound nucleus */
  for(int c0=0 ; c0<sys->max_compound ; c0++){

    /*** Store particle transmission data into array */
    statStoreContinuumTransmission(c0,ncl[c0].excitation[0],pdt,tc);

    /*** Loop over parent excitation */
    for(int k0=0 ; k0<ncl[c0].ncont ; k0++){

      /*** Kinetic energy is a difference between initial CN excitation and average kinetic energy */
      double ek = ekmean - (ncl[0].excitation[k0] - exmean); // ncl[0] is not a bug

      /*** Clear spectrum here in order to convert into the LAB */
      spc->memclear("cn");

      /*** Store gamma-ray transmission data into array */
      statStoreContinuumGammaTransmission(k0,tg,&ncl[c0]);

      /*** Transmission coefficients for discrete levels */
      statStoreDiscreteTransmission(c0,ncl[c0].excitation[k0],pdt,td);

      /*** All residual CN decay */
      specCompoundDecay(c0,k0,tc,td,tg,spc);

      if(ncl[c0].excitation[k0] >= ncl[c0].cdt[neutron].binding_energy){
        specCmsToLab(c0,ek,ncl[0].de,spc->cn[neutron],spc->pe[neutron]);
      }
      /*** Add to total particle emission spectrum */
      specCumulativeSpectra(spc->getCsize(),spc->getNsize(),spc->cn,&ncl[c0]);
    }

    /*** Clear cn-array one more time, since spc is already added in the k-loop */
    spc->memclear("cn");

    /*** Particle decay from levels if possible, ignore spec contribution */
    specLevelDecay(c0,pdt,&ncl[c0],tc,td,spc);

    /*** Gamma-ray cascading, ground state population */
    specGammaCascade(spc->cn[0],&ncl[c0]);

    crx.prod[c0].xsec += ncl[c0].lpop[0];

    /*** Add to total particle emission spectrum */
    specCumulativeSpectra(spc->getCsize(),spc->getNsize(),spc->cn,&ncl[c0]);
  }
}


/**********************************************************/
/*      Distribute p(eps) to pL(Elmin) - pL(Elmax)        */
/**********************************************************/
void specCmsToLab(const int c0, const double ek, const double de, double *spccms, double *spclab)
{
  double mass = (double)ncl[c0].za.getA();
  int n = beohZeroCut(spccms);
  if(n == 0) return;

  const int nt = 180;
  double dt = PI/nt;

  int nmax = n;
  for(int i=0 ; i<n ; i++){
    double ecms = (i == 0) ? 0.25*de : i*de;

    /*** Elab = Ecms + Ek/A + 2 sqrt(Ecms Ek / A) cos(th) = x1 + x2 cos(th) */
    double x1 = ecms + ek / mass;
    double x2 = 2.0 * sqrt(ecms * ek / mass);

    double w = 0.0;
    for(int j=0 ; j<nt ; j++)  w += cos(j*dt) - cos((j+1)*dt);
//  for(int j=0 ; j<nt ; j++)  w += sin((j+0.5)*dt);


    for(int j=0 ; j<nt ; j++){
      double t0  = j * dt;
      double t1  = t0 + dt;

      double el0 = x1 + x2*cos(t0);
      double el1 = x1 + x2*cos(t1);
      double el  = (el0 + el1) * 0.5;

      int kl = specFindEnergyBin(el,de);

      spclab[kl] += spccms[i] * (cos(t0) - cos(t1)) / w;
//    spclab[kl] += spccms[i] * sin((t0+t1)*0.5) / w;

      if(kl > nmax) nmax = kl;
    }
  }
}



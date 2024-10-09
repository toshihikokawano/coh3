/******************************************************************************/
/*  statlevdecay.cpp                                                          */
/*        decay of discrete level                                             */
/******************************************************************************/

#include <iostream>

#include "structur.h"
#include "statmodel.h"
#include "nucleus.h"


/**********************************************************/
/*      Particle Emission from Discrete Levels            */
/**********************************************************/
void specLevelDecay(const int c0, Pdata *pdt, Nucleus *n0, Transmission **tc, Transmission **td, Spectra *spc)
{
  bool   tstat[MAX_CHANNEL];

  /*** Loop over Discrete Levels, from top to bottom*/
  for(int k0=n0->ndisc-1 ; k0>=1 ; k0--){

    /*** check if the parent level population is not zero */
    double pop = n0->lpop[k0];
    if( pop==0.0 ) continue;

    /*** check if particle emission from the level is possible */
    tstat[gammaray] = true;
    for(int id=1 ; id<MAX_CHANNEL ; id++){
      tstat[id] = n0->cdt[id].status;
      if(tstat[id] && (n0->lev[k0].energy < n0->cdt[id].binding_energy)) tstat[id]=false;
    }

    bool ldecay = false;
    for(int id=1 ; id<MAX_CHANNEL ; id++) if(tstat[id] == true) ldecay = true;
    if(!ldecay) continue;

    /*** find bin number corresponding to that level */
    int kl = 0;
    for(int k=0 ; k<n0->ntotal ; k++){
      if(n0->excitation[k] <= n0->lev[k0].energy){
        kl = k;
        break;
      }
    }

    /*** transmission from the level to continuum/daughter levels */
    statStoreContinuumTransmission(c0,n0->lev[k0].energy,pdt,tc);
    statStoreDiscreteTransmission(c0,n0->lev[k0].energy,pdt,td);

    int j0 = (int)(2.0*n0->lev[k0].spin);
    int p0 = n0->lev[k0].parity;

    /*** Sum particle transmissions */
    double tp = 0.0;
    for(int id=1 ; id<MAX_CHANNEL ; id++){
      if(!tstat[id]) continue;
      Nucleus *n1 = &ncl[n0->cdt[id].next];
      tp += specLevelTransitionParticle(sumall,kl,id,p0,j0,tc[id],td[id],0.0,n0,n1,spc->cn[id],spc->dp[id]);
    }

    /*** calculate gamma-ray only if particle transition is possible
         the gamma transition between two levels are approximated
         by the strength function model */
    if(tp > 0.0){
      double tg = 0.0;
      /*** Sum runs over lower discrete levels by gamma */
      if(tstat[gammaray]) tg = specLevelTransitionGamma(sumall,k0,p0,j0,0.0,n0);

      /*** correct level population */
      double r = tg/(tp+tg);
      n0->lpop[k0] *= r;

      /*** increase populations in the daughter */
      double dp = pop/(tp+tg);
      for(int id=1 ; id<MAX_CHANNEL ; id++){
        if(!tstat[id]) continue;
        Nucleus *n1 = &ncl[n0->cdt[id].next];
        specLevelTransitionParticle(hauser,kl,id,p0,j0,tc[id],td[id],dp,n0,n1,spc->cn[id],spc->dp[id]);
      }
    }
  }
}


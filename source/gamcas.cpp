/******************************************************************************/
/*  gamcas.cpp                                                                */
/*        gamma-ray cascading                                                 */
/******************************************************************************/

#include <iostream>

#include "structur.h"
#include "statmodel.h"
#include "output.h"
#include "parameter.h"
#include "global.h"


/**********************************************************/
/*      Gamma Cascading From Each Level                   */
/**********************************************************/
int specGammaCascade(double *spc, Nucleus *n)
{
  double s = 0.0;
  int    k = 0, m = 0;

  double tisom = parmGetValue(parmISOM);
  bool   gcut  = (tisom > 0.0) ? true : false;

  /*** from a parent level */
  for(int i0=n->ndisc-1 ; i0>0 ; i0--){

    /*** if this level is long-lived, skip decay */
    if( gcut && (n->lev[i0].halflife > tisom) ) continue;

    /*** for each gamma transition to a daughter level */
    for(int j=0 ; j<n->lev[i0].ngamma ; j++){
      int i1 = n->lev[i0].fstate[j]; // final state index

      /*** increment the daughter level population */
      n->lpop[i1] += (s = n->lpop[i0]*n->lev[i0].branch[j]);

      if(s > 0.0){
        /*** add to spectrum */
        double e = n->lev[i0].energy - n->lev[i1].energy;
        if( (k = specFindEnergyBin(e,n->de)) >= 0 ){
          if(opt.internalconversion) s *= n->lev[i0].gratio[j];
          spc[k] += s/n->de;
          m++;
        }
      }
    }
  }
  return (m);
}


/**********************************************************/
/*      Gamma-rays From Levels Are Added to dPop          */
/**********************************************************/
int specGammaCascadeExclusive(const int k0, double *dlp, Nucleus *n)
{
  double s     = 0.0;
  double emin0 = n->excitation[0]-(k0+0.5)*n->de; if(emin0 < 0.0) emin0 = 0.0;
  double emax0 = n->excitation[0]-(k0-0.5)*n->de;

  for(int k1=k0 ; k1<n->ntotal ; k1++){

    double emin1 = n->excitation[0]-(k1+0.5)*n->de; if(emin1 < 0.0) emin1 = 0.0;
    double emax1 = n->excitation[0]-(k1-0.5)*n->de;

    for(int i0=n->ndisc-1 ; i0>0 ; i0--){

      if( (emin0 > n->lev[i0].energy) || (n->lev[i0].energy >= emax0) ) continue;

      for(int j=0 ; j<n->lev[i0].ngamma ; j++){
        int i1 = n->lev[i0].fstate[j];

        if( (emin1 > n->lev[i1].energy) || (n->lev[i1].energy >= emax1) ) continue;

        n->lpop[i1] += (s = n->lpop[i0]*n->lev[i0].branch[j]);

        if(opt.internalconversion) s *= n->lev[i0].gratio[j];
        if(s > 0.0) dlp[k1] += s/n->de;
      }
    }
  }
  return (0);
}

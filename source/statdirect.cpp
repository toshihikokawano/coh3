/******************************************************************************/
/*  statdirect.cpp                                                            */
/*        to incorporate direct reaction into statistical model               */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "nucleus.h"
#include "statmodel.h"
#include "global.h"
#include "terminate.h"


/**********************************************************/
/*     Add Direct Cross Section to Level Excite           */
/**********************************************************/
void statSetupDirectInelastic(Level *dir, Nucleus *n, CrossSection *crx)
{
  const double e1 = 0.01, e2 = 1e-6;
  double *ap[MAX_DIRECT];
  int     dlev[MAX_DIRECT];

  /*** count number of levels to which direct inelastic cross sections
       are given, elastic included */
  int    ndir = 1;
  for(int i=0 ; i<MAX_DIRECT ; i++){

    ap[i]    = crx->angdist[i];
    dlev[i]  = -1;
    if(dir[i].energy <= 0.0) continue;

    ndir++;

    /*** search for the same spin and energy levels */
    bool found = false;
    for(int j=0 ; j<n->ndisc ; j++){
      if( (fabs(     dir[i].energy -  n->lev[j].energy) <= e1) &&
          (fabs(fabs(dir[i].spin)  -  n->lev[j].spin  ) <= e2) ){
        found = true;
        dlev[i] = j;
        break;
      }
    }

    /*** if not found, check if DWBA weak coupled levels */
    if(!found){
      found = false;
      for(int j=0 ; j<n->ndisc ; j++){
        if(fabs(dir[i].energy -  n->lev[j].energy) <= e1){
          found = true;
          dlev[i] = j;
          break;
        }
      }
    }

    if(!found){
      message << i << "-th irect level energy " << dir[i].energy << " did not match";
      cohTerminateCode("statSetupDirectInelastic");
    }
  }
  dlev[0] = 0;

  /*** if shape elastic is only given, no more process */
  if(ndir == 1) return;

  int k = ndir; // lowest angdist pointer to be used for non-direct levels

  /*** for all discrete levels */
  for(int j=0 ; j<n->ndisc ; j++){

    bool found = false;
    int i = 0;
    for(i=0 ; i<ndir ; i++){
      if(i >= MAX_DIRECT) continue;
      if(dlev[i] == j){ found = true; break; }
    }
    if(found){
      /*** add direct excitation to the level population,
           exclude elastic cross section */

      if(i != 0) n->lpop[j] += fabs(crx->direct[i]);

      /*** point ang-dist arrays that contain DI cross sections */
      if(j<MAX_ANGDISTLEVELS && prn.angdist) crx->angdist[j] = ap[i];
    }
    else{
      /*** shift ang-dist pointer */
      if(k<MAX_ANGDISTLEVELS && prn.angdist) crx->angdist[j] = crx->angdist[k];
      k++;
    }
  }
}


/**********************************************************/
/*     Add Direct Cross Section to Total Spectrum         */
/**********************************************************/
void statDirectSpectra(const double ein, const double de, Level *dir, double *spc, double *cdir)
{
  for(int i=1 ; i<MAX_DIRECT ; i++){
    if(dir[i].energy <= 0.0) continue;

    double eout = ein - dir[i].energy;
    /*** for closed channels, the energies can be negative */
    if(eout < 0.0) continue;

    int kp = specFindEnergyBin(eout,de);
    if(kp >= 0) spc[kp] += fabs(cdir[i])/de;
  }
}

/******************************************************************************/
/*  statmaniplevel.cpp                                                        */
/*        manipulate discrete levels, esp. for including meta-stable states   */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "levden.h"
#include "kcksyst.h"
#include "parameter.h"
#include "mt19937.h"
#include "metastable.h"
#include "terminate.h"


/**********************************************************/
/*      Add Extra Discrete Level Data                     */
/*      --------                                          */
/*      number of extra levels given by parmLCUT          */
/**********************************************************/
bool statAddDiscreteLevels(NuclearStructure *d)
{
  /*** check the highest level stored in array */
  int nhigh = 1;
  for(int i=1 ; i<MAX_LEVELS ; i++){
    if(d->lev[i].energy == 0.0){ nhigh = i; break; }
  }

  int lcut = parmGetValue(parmLCUT,d->za);

  bool t = false;
  int  norg = d->nlevel;

  if(lcut <= nhigh){
    if(d->nlevel < lcut){
      norg = d->nlevel;
      d->nlevel = lcut;
      t = true;
    }
    statFixDiscreteLevels(d);
  }
  else{
    message << "level cut number " << lcut << " larger than stored data " << nhigh;
    cohWarningMessage("statAddDiscreteLevels");
  }
  
  if(t){
    message << "number of discrete levels for Z " << d->za.getZ() << " - A " << d->za.getA() << " extended from " << norg << " to " << lcut;
    cohNotice("statAddDiscreteLevels");
  }

  return(t);
}


/**********************************************************/
/*      Force Include Meta-Stable States                  */
/*      --------                                          */
/*      add meta-states if exist above the highest level  */
/**********************************************************/
bool statForceIncludeMeta(NuclearStructure *d, Nucleus *n)
{
  /*** search for a meta-stable state above the highest included level
       meta: index of the highest meta
       ntop: number of levels between the highest included level and meta */
  int meta = 0, ntop = 0;
  for(int i=d->nlevel ; i<MAX_LEVELS ; i++){
    if(d->lev[i].energy == 0.0){ ntop = meta - d->nlevel; break; }
    if(d->lev[i].halflife > thalfmin) meta = i; // remember the highest one
  }
  /*** when no meta-stable found */
  if(meta == 0) return false;

  /*** when meta is the first excited state and
       unknown if there is any other levels in between,
       we assume the meta-state is actually the first excited state */
  if((meta == 1) && (ntop == 0)){ d->nlevel = 2; return true; }

  double e0 = d->lev[d->nlevel-1].energy; // highest known level energy
  double e1 = d->lev[meta].energy;        // meta-state energy to be included
  double de = 0.01;                       // integration width, 10keV
  int    ediv  = (e1 - e0)/de;
  de = (e1 - e0)/(double)ediv;

  /*** general case, estimate missing levels from level density in [E0,E1] */
  double r0 = 0.0, r1 = 0.0, sr = 0.0;
  for(int k=1 ; k<=ediv ; k++){
    r1 = ldLevelDensity(k*de,(double)n->za.getA(),&n->ldp);
    sr += (r0 + r1) * 0.5 * de;
    r0 = r1;
  }

  /*** include extra artificial levels */
  int nx = (int)sr - ntop; // number of artificial levels

  if(nx > 0){
    de = (e1 - e0)/(nx + 1.0);
    int m = 0;
    for(int i=1 ; i<=nx ; i++){
      int idx = meta + i;
      if(idx >= MAX_LEVELS-1) break;

      d->lev[idx].energy    = e0 + de * i;
      d->lev[idx].spin      = -1.0;
      d->lev[idx].parity    = 0.0;
      /*** gamma branch, 100% decays to ground state */
      d->lev[idx].ngamma    = 1;
      d->lev[idx].fstate[0] = 0;
      d->lev[idx].branch[0] = 1.0;
      d->lev[idx].gratio[0] = 1.0;
      m++;
    }
    d->nlevel += ntop + m;
  }
  else d->nlevel = meta + 1;

  /*** make up spin and parity up to the highest meta */
  statFixDiscreteLevels(d);

  return true;
}


/**********************************************************/
/*      Fix Discrete Level Data If Not Assigned           */
/*      --------                                          */
/*      J-pi assignment based on level density            */
/**********************************************************/
int statFixDiscreteLevels(NuclearStructure *d)
{
  /*** quick scan if data are incomplete */
  bool needfix = false;
  for(int i=0 ; i<d->nlevel ; i++){
    if( (d->lev[i].spin < 0.0) || (d->lev[i].parity == 0) ) needfix = true;
  }
  if(!needfix) return(0);

  /*** reset random number seed by Z-A,
       so that the sequence of spin-parity will be the same all the time */
  unsigned long seed = d->za.getZ() * 1000 + d->za.getA();
  genrand_init(seed);

  double s0 = sqrt(kckSpinCutoff(d->za.getA()));
  double c  = sqrt(PI2)*s0*s0*s0;

  for(int i=0 ; i<d->nlevel ; i++){

    /*** re-assign spin, sampling from spin and parity distributions */
    if(d->lev[i].spin < 0.0){
      double r = genrand_real1();
      double s = 0.0;
      for(int j=0 ; j<MAX_J ; j++){
        double x = j+halfint(d->lev[0].spin);
        s += (x+0.5)/c * exp(-(x+0.5)*(x+0.5)/(2*s0*s0)) * (2*x+1.0);
        if(s > r){
          d->lev[i].spin = x;
          break;
        }
      }
    }

    /*** assume even distribution for parity */
    if(d->lev[i].parity == 0){
      double r = genrand_real1();
      d->lev[i].parity = (r > 0.5) ? 1 : -1;
    }
  }

  return(d->nlevel);
}


/**********************************************************/
/*      Population Trap in Meta-Stable States             */
/*      --------                                          */
/*      set number of gamma-ray branch zero,              */
/*      so that the level does not decay                  */
/**********************************************************/
void statPopTrapMeta(NuclearStructure *d)
{
  for(int i=1 ; i<MAX_LEVELS ; i++){
    if(d->lev[i].energy == 0.0) break;
    if(d->lev[i].halflife > thalfmin){
      d->lev[i].ngamma = 0;

      message << "discrete level of " << d->lev[i].energy << " marked as meta";
      cohNotice("statPopTrapMeta");
    }
  }
}


/**********************************************************/
/*      Remove Discrete Levels At Initial CN Energies     */
/**********************************************************/
void statCutDiscreteLevels(const double emax, Nucleus *n)
{
  int nhigh = n->ndisc;
  for(int i=1 ; i<n->ndisc ; i++){
    if(n->lev[i].energy >= emax){
      nhigh = i-1;

      message << "discrete level for Z " << n->za.getZ() << " - A " << n->za.getA() << " truncated at " << nhigh << " because Emax = " << emax;
      cohNotice("statCutDiscreteLevels");

      break;
    }
  }
  n->ndisc = nhigh;
}


/***********************************************************/
/*     Subtract Discrete Levels from Continuum             */
/*     --------                                            */
/*     modify bins that contain extra levels               */
/***********************************************************/
void statAdjustDensityByExtraLevels(Nucleus *n)
{
  /*** fix level density of each bin that contains extra levels*/
  for(int k=1 ; k<n->ncont ; k++){
    double ex = n->excitation[k];  if(ex == 0.0) continue;
    double e0 = ex - n->de*0.5;  if(e0 < 0.0) e0 = 0.0;
    double e1 = ex + n->de*0.5;

    /*** count number of extra levels inside an energy bin */
    int nl = 0;
    for(int i=1 ; i<n->ndisc ; i++){
      if(n->lev[i].energy < n->excitation[n->ncont]) continue;
      if( (e0 < n->lev[i].energy) && (n->lev[i].energy <= e1) ) nl++;
    }

    /*** subtract partial density corresponding to the levels,
         spin dependence ignored */
    if(nl > 0){
      double dr = (double)nl/n->de;
      double r  = ldLevelDensity(ex,(double)n->za.getA(),&n->ldp) - dr;
      if(r < 0.0) r = 0.0;

      for(int j=0 ; j<MAX_J ; j++){
        double sd = ldSpinDistribution(j+halfint(n->lev[0].spin),n->ldp.spin_cutoff,n->ldp.sigma0,1.0);
        n->density[k][j].even = r*sd*ldParityDistribution( 1);
        n->density[k][j].odd  = r*sd*ldParityDistribution(-1);
      }
    }
  }
}


/***********************************************************/
/*     Zero Level Density at Low Excitation Energies       */
/*     --------                                            */
/*     modify low energy bins regardless extra-levels      */
/***********************************************************/
void statClearDensityByExtraLevels(Nucleus *n)
{
  /*** reduce level densities from the bottom */
  /*** count total number of extra levels */
  double nl = 0.0;
  for(int i=1 ; i<n->ndisc ; i++){
    if(n->lev[i].energy >= n->excitation[n->ncont]) nl += 1.0;
  }

  if(nl > 0.0){
    /*** integrage level density until the sum becomes more than NL */
    double rt = 0.0, rx = 0.0;
    int kx = 0;
    for(int k=n->ntotal ; k>=0 ; k--){
      double ex = n->excitation[k]; // if(ex == 0.0) continue;
      rt += ldLevelDensity(ex,(double)n->za.getA(),&n->ldp) * n->de;

      if(rt >= nl){
        kx = k;                // below this bin, level density will be zero
        rx = (rt - nl)/n->de;  // extra density to be stored, (sum rho - NL)/dE
        break;
      }
    }

    if(kx > 0){
      /*** from the bottom to kx, store zero */
      for(int k=n->ntotal ; k>kx ; k--){
        for(int j=0 ; j<MAX_J ; j++){
          n->density[k][j].even = 0.0;
          n->density[k][j].odd  = 0.0;
        }
      }
      /*** at kx-th bin, store the rest density */
      for(int j=0 ; j<MAX_J ; j++){
        double sd = ldSpinDistribution(j+halfint(n->lev[0].spin),n->ldp.spin_cutoff,n->ldp.sigma0,1.0);
        n->density[kx][j].even = rx*sd;
        n->density[kx][j].odd  = rx*sd;
      }
    }
  }
}


/***********************************************************/
/*     Adjust Discrete/Continuum Boundary                  */
/***********************************************************/
void statAdjustContinuumBoundary(Nucleus *n)
{
  double meta = 0.0;
  for(int i=1 ; i<n->ndisc ; i++){
    if(n->lev[i].halflife > thalfmin) meta = n->lev[i].energy;
  }

  int nc = n->ncont;
  for(int k=0 ; k<n->ncont ; k++){
    if(n->excitation[k] < meta){ nc = k; break; }
  }
  n->ncont = nc;
}

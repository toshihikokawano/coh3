/******************************************************************************/
/*  setup.cpp                                                                 */
/*        setting up basic model parameters                                   */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "beoh.h"
#include "setupinit.h"
#include "statmodel.h"
#include "nucleus.h"
#include "masstable.h"
#include "ripl2levels.h"
#include "terminate.h"
#include "parameter.h"

static void setupEnergy (const bool, System *);


/**********************************************************/
/*     Clean Parameters                                   */
/**********************************************************/
void setupClearParameter(const int mc, System *sys, Pdata *pdt, GDR *gdr)
{
  /*** clear system parameters */
  sys->init();

  /*** clear channel parameter */
  for(int p=0 ; p<mc ; p++) pdt[p].omp = 0;

  /*** clear GDR parameters */
  for(int i=0 ; i<MAX_GDR ; i++) gdr[i].clear();

  /*** clear fission barriers */
  for(int i=0 ; i<MAX_FISS_CHANCE ; i++) fbr[i].init();

  crx.init();
  for(int i=0 ; i<MAX_COMPOUND ; i++) crx.prod[i].init();

  /*** clear gamma-ray production array */
  gml.init();
}


/**********************************************************/
/*     Setup System Parameters and Store                  */
/**********************************************************/
void setupGeneralParameter(const bool betacalc, System *sys)
{
  /*** fake incident particle, to be photon */
  sys->incident.za.setZA(0,0);
  sys->incident.pid = gammaray;
  setupEnergy(betacalc,sys);

  /*** precursor level data stored into nst[0] */
  nst[0].nlevel = riplReadDiscreteLevels(&sys->target, nst[0].lev, reassign);
  nst[0].za.setZA(sys->target.getZ(),sys->target.getA());

  if(nst[0].lev[0].spin < 0.0){
    message << "ground state spin not assigned";
    cohTerminateCode("setupGeneralParameter");
  }

  /*** look for isomeric states in the precursor */
  if(sys->target_level > 0){
    bool found = false;
    int  isomer = 0;
    for(int i=1 ; i<nst[0].nlevel ; i++){
      if(nst[0].lev[i].halflife > 0.0){
        isomer++;
        if(sys->target_level == isomer){
          if(nst[0].lev[i].spin < 0.0){
            message << "isomeric state spin not assigned";
            cohTerminateCode("setupGeneralParameter");
          }
          sys->excitation = nst[0].lev[i].energy;
          sys->target_level = i;
          found = true;
          break;
        }
      }
    }
    if(!found){
      message << "isomeric state " << sys->target_level << " not found";
      cohTerminateCode("setupGeneralParameter");
    }
  }

  /*** check adjustable parameters, if Gaussian widths are given */
  parmCheckFactor();
}


/**********************************************************/
/*     Setup Compound Reaction Parameters                 */
/**********************************************************/
int setupStatModel(const bool betacalc, const int jmax, System *sys, Pdata *pdt)
{
  /*** max spin in the compound (initial guess) */
  ncl[0].jmax = (jmax >= MAX_J) ? MAX_J-1 : jmax;

  /*** make reaction chain */
  sys->max_compound = setupChain(sys->max_particle,sys->ex_total,&sys->compound,pdt);
  sys->uniq_compound = setupCountCompound(sys->max_compound);

  /*** read in discrete levels */
  int nl = 0;

  /*** skip the first object that contains precursor in the beta-decay case */
  if(betacalc) nl++;

  for(int i=0 ; i<sys->max_compound ; i++){

    /*** check if the same nucleus exists */
    int k = -1;
    for(int j=0 ; j<i ; j++){
      if(ncl[i].za == ncl[j].za){ k = j; break; }
    }

    /*** if yes, point the same memory space */
    if(k >= 0){
      ncl[i].ndisc = ncl[k].ndisc;
      ncl[i].lev   = ncl[k].lev;
    }
    /*** if this is new, read database, and point the array */
    else{
      cohAllocateNuclearStructure();
      nst[nl].nlevel = riplReadDiscreteLevels(&ncl[i].za,nst[nl].lev,reassign);
      nst[nl].za.setZA(ncl[i].za.getZ(),ncl[i].za.getA());

      /*** if level cut-off number is given, and its less than
           nlevel (RIPL tells they are OK), move the nlevel to lcut */

      int lcut = parmGetValue(parmLCUT,nst[nl].za);
      if((1 <= lcut) && (lcut < nst[nl].nlevel)) nst[nl].nlevel = lcut;

      ncl[i].ndisc = nst[nl].nlevel;
      ncl[i].lev = nst[nl].lev;
      nl ++;
    }

    /*** clear population array */
    ncl[i].memclear();
  }

  /*** max J for each compound */
  ncl[0].jmax = jmax + (int)ncl[sys->target_id].lev[sys->target_level].spin;
  if(ncl[0].jmax >= MAX_J) ncl[0].jmax = MAX_J-1;
  for(int i=1 ; i<sys->max_compound ; i++) ncl[i].jmax = ncl[0].jmax;

  /*** count number of fission barriers given */
  for(int i=0 ; i<MAX_FISS_CHANCE ; i++){
    fbr[i].hump = 0;
    for(int k=0 ; k<MAX_HUMP ; k++) if(fbr[i].barrier[k].height > 0.0) fbr[i].hump ++;
  }

  /*** point fission to given fission barriers */
  for(int i=0 ; i<MAX_FISS_CHANCE ; i++){
    if(fbr[i].za.getA() == 0 || fbr[i].za.getZ() == 0) continue;
    for(int j=0 ; j<sys->max_compound ; j++){
      if(fbr[i].za == ncl[j].za){
        ncl[j].fission = &fbr[i];
        ncl[j].fissile = true;
      }
    }
  }

  /*** number of actual open channels */
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(pdt[j].omp == 0) continue;
    sys->max_channel = j;
  }
  sys->max_channel ++;

  return(sys->max_compound);
}


/**********************************************************/
/*     Set Energy                                         */
/**********************************************************/
void setupEnergy(const bool betacalc, System *sys)
{
  double qbeta = 0.0;

  /*** for beta-decay calculation, target is given as an input */
  if((sys->target.getZ() > 0) && (sys->target.getA() > 0) && betacalc){

    /*** determine beta- or beta+ */
    double qm = mass_excess(sys->target.getZ()  ,sys->target.getA())
              - mass_excess(sys->target.getZ()+1,sys->target.getA());

    double qp = mass_excess(sys->target.getZ()  ,sys->target.getA())
              - mass_excess(sys->target.getZ()-1,sys->target.getA());

    /*** when Qbeta(-) > 0, we assume this is beta minus, regardless of Qbeta(+) */
    if(qm > 0.0){
      sys->compound.setZA(sys->target.getZ()+1,sys->target.getA());
      qbeta = qm;
    }
    else if(qp > 0.0){
      sys->compound.setZA(sys->target.getZ()-1,sys->target.getA());
      qbeta = qp;
    }
    else{
      message << "precursor is stable";
      cohTerminateCode("setupEnergy");
    }

    ncl[0].za = sys->compound;
  }
  /*** for statistical decay calculation, compound is already set */
  else{
    sys->target.setZA(sys->compound.getZ(),sys->compound.getA());
  }

  /*** Qbeta, if not given */
  if(sys->ex_total == 0.0){
    if(betacalc){
      sys->ex_total = sys->excitation = qbeta;
    }
    else{
      message << "compound nucleus excitation energy zero";
      cohTerminateCode("setupEnergy");
    }
  }
  else{
    sys->excitation = sys->ex_total;
  }
}


/**********************************************************/
/*     Set Particle Data                                  */
/**********************************************************/
Pdata setupParticle(const Particle p)
{
  Pdata pdt;

  switch(p){
  case neutron:
    pdt.za.setZA(0,1);  pdt.mass = MNEUTRON;   pdt.mass_excess = ENEUTRON;
    break;
  case proton:
    pdt.za.setZA(1,1);  pdt.mass = MPROTON;    pdt.mass_excess = EPROTON;
    break;
  case alpha:
    pdt.za.setZA(2,4);  pdt.mass = MALPHA;     pdt.mass_excess = EALPHA;
    break;
  case deuteron:
    pdt.za.setZA(1,2);  pdt.mass = MDEUTERON;  pdt.mass_excess = EDEUTERON;
    break;
  case triton:
    pdt.za.setZA(1,3);  pdt.mass = MTRITON;    pdt.mass_excess = ETRITON;
    break;
  case helion:
    pdt.za.setZA(2,3);  pdt.mass = MHELIUM3;   pdt.mass_excess = EHELIUM3;
    break;
  case gammaray:
    // fall through
  default:
    pdt.za.setZA(0,0);  pdt.mass = 0.0;        pdt.mass_excess = 0.0;
    break;
  }

  return(pdt);
}


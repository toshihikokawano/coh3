/******************************************************************************/
/*  setup.cpp                                                                 */
/*        setting up basic model parameters                                   */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "coh.h"
#include "setupinit.h"
#include "nucleus.h"
#include "masstable.h"
#include "ripl2levels.h"
#include "global.h"
#include "terminate.h"
#include "parameter.h"

static void setupZA (System *);
static void setupEnergy (System *);

/**********************************************************/
/*     Clean Parameters                                   */
/**********************************************************/
void setupClearParameter(const int mc, System *sys, Pdata *pdt, Direct *dir, Dcapt *cap, FNSpec *fns)
{
  /*** clear system parameters */
  sys->init();

  /*** clear channel parameter */
  for(int p=0 ; p<mc ; p++) pdt[p].omp = 0;

  /*** clear direct reaction parameters */
  dir->init();

  /*** clear DSD parameters */
  cap->init();

  /*** clear fission barriers */
  for(int i=0 ; i<MAX_FISS_CHANCE ; i++) fbr[i].init();

  /*** clear fission neutron spectrum parameters */
  fns->init();

  /*** clear all cross section results */
  setupClearCrxArray();

  ctl.entrance       = false;
  ctl.deformed       = false;
  ctl.exciton        = false;
  ctl.dwba           = false;
  ctl.dsd            = false;
  ctl.fluctuation    = false;
  ctl.statmodel      = false;
  ctl.fission        = false;
  ctl.fns            = false;
}


/**********************************************************/
/*     Setup System Parameters and Store                  */
/**********************************************************/
void setupGeneralParameter(System *sys, Pdata *pdt, Direct *dir)
{
  /*** check incident particle */
  setupZA(sys);
  if(sys->incident.pid == unknown){
    message << "unknown incident particle " << sys->incident.pid << " detected";
    cohTerminateCode("setupGeneralParameter");
  }

  sys->incident = pdt[sys->incident.pid];

  setupEnergy(sys);
  if(sys->ex_total < 0.0){
    message << "target is unbound, the total energy " << sys->ex_total << " is negative";
    cohTerminateCode("setupGeneralParameter");
  }

  if((sys->incident.omp == 0) && (sys->incident.pid != gammaray)){
    message << "incident channel optical potential not defined";
    cohTerminateCode("setupGeneralParameter");
  }

  /*** target level data stored into ncl[0], to be changed later */
  nst[0].nlevel = riplReadDiscreteLevels(&sys->target, nst[0].lev, reassign);
  nst[0].za.setZA(sys->compound.getZ(),sys->compound.getA());
  ncl[0].ndisc = nst[0].nlevel;
  ncl[0].lev = nst[0].lev;

  /*** insert ground state for the direct reaction calculation */
  if(sys->excitation == 0.0){
    dir->lev[0].energy = 0.0;
    dir->lev[0].spin = ncl[0].lev[0].spin * ncl[0].lev[0].parity;
  }

  /*** if coupled levels are given, deformed calculation */
  double beta2r = dir->defstatic[0];
  double beta2v = dir->defdynamic[0];
  ctl.deformed = ( (beta2r != 0.0) || (beta2v != 0.0) ) ? true : false;
  if(ctl.deformed) sys->beta2 = fabs(beta2r);

  /*** target level, if excited nucleus */
  sys->target_level = 0;
  if(sys->excitation > 0.0){
    double eps = 0.0001; // if difference is less than 0.1 keV
    bool found = false;
    for(int i=1 ; i<ncl[0].ndisc ; i++){
      if( fabs(sys->excitation - ncl[sys->target_id].lev[i].energy) <= eps ){
        sys->target_level = i;
        found = true;
        break;
      }
    }
    if(!found){
      message << "excited target state " << sys->excitation << " not found in level scheme";
      cohTerminateCode("setupGeneralParameter");
    }
  }

  /*** this is for the entrance channel */
  ctl.entrance = true;

  /*** check adjustable parameters, if Gaussian widths are given */
  parmCheckFactor();
}


/**********************************************************/
/*     Setup Compound Reaction Parameters                 */
/**********************************************************/
int setupStatModel(const int jmax, System *sys, Pdata *pdt)
{
  /*** make reaction chain */
  sys->max_compound = setupChain(sys->max_particle,sys->ex_total,&sys->compound,pdt);
  sys->uniq_compound = setupCountCompound(sys->max_compound);

  /*** this is not an entrance channel */
  ctl.entrance = false;

  /*** look for target nucleus index */
  for(int id=0 ; id<MAX_CHANNEL ; id++){
    if(ncl[0].cdt[id].pid == sys->incident.pid){
      sys->target_id = ncl[0].cdt[id].next;
      sys->inc_id    = id;
      break;
    }
  }

  /*** read in discrete levels,
       data in the first object will be overwritten */
  int nl = 0;
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

  /*** eliminate levels in the initial compound nucleus,
       if their energies are larger than the separation energy (avoid resonance) */
  if(ncl[0].cdt[sys->inc_id].binding_energy <  ncl[0].lev[ncl[0].ndisc-1].energy){
    int ntop = ncl[0].ndisc;
    for(int j=ntop-1 ; j>=0 ; j--){
      if(ncl[0].lev[j].energy < ncl[0].cdt[sys->inc_id].binding_energy){
        ntop = j;
        break;
      }
    }
    ncl[0].ndisc = ntop;
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
        ctl.fission = true;
      }
    }
  }

  /*** number of actual open channels */
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(pdt[j].omp == 0) continue;
    sys->max_channel = j;
  }
  sys->max_channel ++;

  message << "max J for initial CN " << ncl[0].jmax;
  cohNotice("setupStatModel");

  return(sys->max_compound);
}


/**********************************************************/
/*     Fill Zero in Cross Section Data                    */
/**********************************************************/
void setupClearCrxArray()
{
  /*** clear cross section arrays */
  crx.init();

  for(int i=0 ; i<MAX_DIRECT ; i++){
    crx.direct[i] = 0.0;
    for(int j=0 ; j<MAX_J*3 ; j++) crx.transmission[i][j] = 0.0;
  }

  for(int i=0 ; i<MAX_ANGDISTLEVELS ; i++){
    for(int k=0 ; k<MAX_ANGDIST ; k++) crx.angdist[i][k] = 0.0;
  }

  for(int i=0 ; i<MAX_CHANNEL ; i++){
    for(int j=0 ; j<MAX_LEVELS   ; j++) crx.levexcite[i][j] = 0.0;
    for(int k=0 ; k<MAX_ANGDISTLEVELS ; k++){
      for(int j=0 ; j<MAX_J ; j++) crx.legcoef[i][k][j] = 0.0;
    }
  }

  for(int i=0 ; i<MAX_COMPOUND ; i++) crx.prod[i].init();

  /*** clear gamma-ray production array */
  gml.init();
}


/**********************************************************/
/*     Check Incident Particle                            */
/**********************************************************/
void setupZA(System *sys)
{
  switch(sys->incident.za.getZ()){
  case 0:
    if(     sys->incident.za.getA() == 0) sys->incident.pid = gammaray;
    else if(sys->incident.za.getA() == 1) sys->incident.pid = neutron;
    else                                  sys->incident.pid = unknown;
    break;
  case 1:
    if(     sys->incident.za.getA() == 1) sys->incident.pid = proton;
    else if(sys->incident.za.getA() == 2) sys->incident.pid = deuteron;
    else if(sys->incident.za.getA() == 3) sys->incident.pid = triton;
    else                                  sys->incident.pid = unknown;
    break;
  case 2:
    if(     sys->incident.za.getA() == 3) sys->incident.pid = helion;
    else if(sys->incident.za.getA() == 4) sys->incident.pid = alpha;
    else                                  sys->incident.pid = unknown;
    break;
  case 3: // this is temporary defined in dataread.cpp 
    sys->incident.pid = userdef;
    break;
  }
}


/**********************************************************/
/*     Set Energy                                         */
/**********************************************************/
void setupEnergy(System *sys)
{
  /*** compound Z and A */
  sys->compound = sys->target + sys->incident.za;

  double tmx = mass_excess(sys->target.getZ()  ,sys->target.getA()  );
  double cmx = mass_excess(sys->compound.getZ(),sys->compound.getA());

  double tm  = sys->target.getA() + tmx / AMUNIT;
  double im  = sys->incident.mass;

  /*** mass and energies for binary reaction */
  sys->reduced_mass = tm * im / (tm+im);

  if(sys->input_energy > 0.0){
    sys->lab_energy = sys->input_energy;
    sys->cms_energy = sys->lab_energy*tm/(tm+im);
  }else{
    sys->cms_energy = fabs(sys->input_energy);
    sys->lab_energy = sys->cms_energy*(tm+im)/tm;
  }

  sys->wave_number = (sys->incident.pid == gammaray) ?
    sys->cms_energy/HBAR/VLIGHT :
    sqrt(sys->cms_energy*sys->reduced_mass*2.0*AMUNIT)/VLIGHT/HBAR;

  sys->ex_total = sys->cms_energy + sys->incident.mass_excess + tmx - cmx;

  /*** in the case of excited target */
  if(sys->excitation>0.0) sys->ex_total += sys->excitation;

  message << "target mass " << tmx << " compound mass " << cmx;
  cohNotice("setupEnergy");
}


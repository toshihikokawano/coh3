/******************************************************************************/
/*  cohmain.cpp                                                               */
/*        CoH main calculation loop                                           */
/******************************************************************************/

#include <iostream>

#include "structur.h"
#include "nucleus.h"
#include "coh.h"
#include "setupinit.h"
#include "omcalc.h"
#include "global.h"
#include "parameter.h"
#include "output.h"
#include "fileout.h"
#include "fisbar.h"
#include "cohegrid.h"
#include "macs.h"
#include "terminate.h"

static int  cohReadInputFile    (MFTparm *);
static void cohSuperLazy        (const double, const int, const int);
static int  cohEntranceChannel  (void);
static void cohCheckRange       (ZAnumber *, const double);
static inline int cohEGrid      (const int, double *);
static inline int macsEGrid     (Particle);

static System      sys;                  // system parameters
static Dcapt       cap;                  // radiative capture parameters for both DSD and Stat
static FNSpec      fns(MAX_FISS_CHANCE); // fission neutron spectrum parameters
static Pdata       pdt[MAX_CHANNEL+1];   // emitting particle data 
static Direct      dir;                  // direct reaction parameters
static char        buf[256];             // input data buffer
static double      ein[MAX_EINCIDENT];   // multi-energy point calculation
static double      mkt[MAX_MACSTEMP];    // MACS temperatures

/**********************************************************/
/*      Main Loop                                         */
/*      --------                                          */
/*      Read input and repeat calculations for each block */
/*      All the data will be initialized at each time     */
/**********************************************************/
void cohMainLoop(const double elab, const int targz, const int targa, const unsigned long nsim)
{
  MFTparm     mfp;              // mean-field theory parameters
  HFInterface hfsp[2];          // mean-field single particle states
  int     sec = 0, jmax = 0;

  for(int i=0 ; i<MAX_EINCIDENT ; i++) ein[i] = 0.0;
  for(int i=0 ; i<MAX_MACSTEMP  ; i++) mkt[i] = 0.0;

  /*** store parameters which are constant throughtout the calculation */
  setupInitSystem(MAX_CHANNEL+1,pdt);

  /*** scattering angle, if printed */
  if(prn.angdist) setupScatteringAngle();

  /*** super lazy calculation */
  if(ctl.superlazy){
    cohSuperLazy(elab,targz,targa);
    return;
  }

  /*** read STDIO until EOF */
  do{
    if ( (sec = readHead(buf)) < 0 ) break;

    /*** start over */
    if(sec == BEGIN){
      /*** clear data array */
      setupClearParameter(MAX_CHANNEL,&sys,pdt,&dir,&cap,&fns);
      mfp.init();

      if(prn.system) outTitle(buf);

      /*** read input data */
      if(cohReadInputFile(&mfp) < 0) break;

      /*** if MACS temperature is given, set up default energy grid */
      if(ctl.macs){
        macsEGrid(sys.incident.pid);
        crx.astroalloc(N_REACTION_RATE,N_MACS_EGRID);
      }

      /*** check given parameter range,
           generate random parameter if Gaussian width given */
      parmCheckFactor();

      /*** mean-field theory calculation */
      if(ctl.meanfield){
        ctl.expansion = true;
        hfsp[0].memalloc();
        hfsp[1].memalloc();
        mfp.target = sys.target;
        mftCalc(&mfp,hfsp,dir.defstatic);

        /*** point MFT results for DSD calculation */
        cap.hfsp = (sys.incident.pid == neutron) ? &hfsp[0] : &hfsp[1];
        cap.hf = true;
      }

      /*** override incident energy if command line value is given */
      if(elab != 0.0){ ein[0] = elab;  ein[1] = 0.0; }

      /*** loop for each incident energy */
      for(int i=0 ; i<MAX_EINCIDENT ; i++){ if(ein[i] == 0.0) break;
        sys.input_energy = ein[i];

        /*** clear all cross section results */
        setupClearCrxArray();

        /*** set up energy dependent parameters */
        setupGeneralParameter(&sys,pdt,&dir);

        cohCheckRange(&sys.target,sys.lab_energy);
        outSystem(&sys,prn.system);

        /*** entrance channel optical model, ignore photo reaction */
        if(sys.incident.pid != gammaray){
          jmax = cohEntranceChannel();
          if(jmax <= 0) continue;
        }
        else{
          jmax = MAX_MULTIPOL / 2 + 1;
        }

        /*** DWBA */
        if(ctl.dwba){
          int iret = dwbaCalc(sys.cms_energy,sys.excitation,&sys.incident,&sys.target,sys.reduced_mass,&dir,&crx);
          if(iret < 0) break;
        }

        /*** DSD, Pre-equilibrium, and Hauser-Feshbach */
        if(ctl.dsd || ctl.exciton || ctl.statmodel){
          /*** create nucleus, generate reaction chain */
          setupStatModel(jmax,&sys,pdt);
          statModel(&sys,pdt,&cap,&dir,&fns,nsim);
        }

        if(prn.system) outParameter();

        /*** write production cross sections on file */
        if(ctl.fileout){
          cohFileOut(sys.max_compound,sys.target_id,sys.lab_energy);
        }

        /*** store cross sections for MACS calculation */
        if(ctl.macs) macsStoreCrossSection(i,sys.max_compound);

        /*** if command-line Ein is given, force single energy calc. */
        if(elab != 0) break;

        /*** prepare fpr the next energy point, reset ncl and dsc pointers */
        cohResetAllocationPointer();
      }

      /*** MACS at given temperature */
      if(ctl.macs) macsReactionRate(sys.target_id,sys.incident.pid,mkt);
    }
  }while(!std::cin.eof());
}


/**********************************************************/
/*      Scan Input Data from BEGIN to END                 */
/**********************************************************/
int cohReadInputFile(MFTparm *mfp)
{
  int sec = 0; 

  do{
    if ( (sec = readHead(buf)) < 0 ) return(-1);

    /*** DATA section, OM calc. */
    if(sec == DATA){
      dir.ncc = readSystem(buf,&sys,pdt,&dir,ein,mkt);
      if(dir.ncc == 0) dir.ncc++;

      int eset = (int)parmGetValue(parmESET);
      if(eset > 0) cohEGrid(eset,ein);
    }

    /*** DWBA section */
    else if(sec == DWBA){
      dir.nlev = readDwba(buf,&dir);
      ctl.dwba = true;
    }

    /*** DSD parameter section */
    else if(sec == DSD){
      readDsd(buf,&cap);
      ctl.dsd = true;
    }

    /*** PREEQ Pre-equilibrium parameter section */
    else if(sec == PREEQ){
      readExciton(buf);
      ctl.exciton = true;
    }

    /*** HFS Statistical model parameter section */
    else if(sec == HFS){
      readStatModel(buf,cap.gdr,fbr);
      ctl.statmodel = true;
    }

    /*** FNS Fission neutron spectrum model parameter section */
    else if(sec == FNS){
      if(!ctl.exclusive) ctl.fns = true;
      ctl.exclusive = true;
      readFns(buf,&fns);
    }

    /*** MFT Hartree-Fock and FRDM model parameter section */
    else if(sec == MFT){
      ctl.meanfield = true;
      readMeanfield(buf,mfp);
    }

    /*** Perform Hauser-Feshbach calculation */
    else if(sec == END){
      break;
    }

  }while(!std::cin.eof());

  return 0;
}


/**********************************************************/
/*      SuperLazy Mode (only for neutron incident)        */
/**********************************************************/
void cohSuperLazy(const double elab, const int targz, const int targa)
{
  for(int i=0 ; i<MAX_DIRECT ; i++) dir.lev[i].energy = 0.0;

  setupClearParameter(MAX_CHANNEL,&sys,pdt,&dir,&cap,&fns);
  if(prn.system) outTitle(buf);

  sys.input_energy = elab;
  sys.target.setZA(targz,targa);
  sys.incident.pid = neutron;
  sys.incident.za.setZA(0,1);

  /*** always include neutron */
  pdt[neutron].omp  =  6; // Koning-Delaroche
  /*** ignore charged particles if heavier than Nd */
  if(targz < 60){
    if(MAX_CHANNEL>=2)  pdt[proton].omp   =  6; // Koning-Delaroche
    if(MAX_CHANNEL>=3)  pdt[alpha].omp    = 15; // Avrigeanu2009
    if(MAX_CHANNEL>=4)  pdt[deuteron].omp = 18; // Bojowald
//  if(MAX_CHANNEL>=5)  pdt[triton].omp   =  2; // Becchetti
//  if(MAX_CHANNEL>=6)  pdt[helion].omp   =  2; // Becchetti
  }

  /*** default gamma-ray strength function normalization */
  adj.addparm(parmGSTR,-1.0,0.0);

  /*** modules to be performed */
  ctl.dwba        = false;
  ctl.exciton     = true;
  ctl.dsd         = true;
  ctl.fluctuation = true;
  ctl.statmodel   = true;

  /*** if heavier than Ac-225, include fission */
  if(targz >= 89 && targa >= 225){
    ctl.fission   = true;
    fisDataRead(&sys.target,fbr);
  }
  else{
    ctl.fission   = false;
  }

  setupGeneralParameter(&sys,pdt,&dir);
  cohCheckRange(&sys.target,sys.lab_energy);
  outSystem(&sys,prn.system);

  int jmax = cohEntranceChannel();
  if(jmax <= 0) return;

  setupStatModel(jmax,&sys,pdt);
  statModel(&sys,pdt,&cap,&dir,&fns,0L);
}


/**********************************************************/
/*      Entrance Channel Optical Model Calculation        */
/**********************************************************/
int cohEntranceChannel()
{
  int jmax = 0;

  /*** coupled channels calculation */
  if(ctl.deformed){
    jmax = ccCalc(sys.cms_energy,sys.excitation,&sys.incident,&sys.target,sys.reduced_mass,&dir,crx.transmission,&crx);
  }

  /*** spherical optical model */
  else
    jmax = omCalc(sys.cms_energy,&sys.incident,&sys.target,sys.reduced_mass,crx.transmission[0],&crx);

  if(crx.reaction == 0.0) jmax = 0;

  return(jmax);
}


/**********************************************************/
/*      Check Z/A/E Range                                 */
/*      -----                                             */
/*      Valid Z and A numbers, and Energy range           */
/*      are empirically defined as 5 < Z < 118,           */
/*      11< A < 300, and E < 100 MeV.                     */
/**********************************************************/
void cohCheckRange(ZAnumber *t, const double e)
{
  const double       emax = 100.0;
  const unsigned int zmin =     5;
  const unsigned int zmax =   118;
  const unsigned int amin =    11;
  const unsigned int amax =   300;

  std::string msg;
  bool rangechk = false;
  if((e <= 0.0) || (e > emax)){
    message << "incident energy " << e << " out of range";
    rangechk = true;
  }
  else if((t->getZ() < zmin) || (t->getZ() > zmax)){
    message << "target Z-number " << t->getZ() << " out of range";
    rangechk = true;
  }
  else if((t->getA() < amin) || (t->getA() > amax)){
    message << "target A-number " << t->getA() << " out of range";
    rangechk = true;
  }

  if(rangechk) cohTerminateCode("cohCheckRagne");
}


/**********************************************************/
/*      Copy Energy Grid for Cross Section calculation    */
/**********************************************************/
static inline int cohEGrid(const int k, double *ein)
{
  int n = 0;
  double *e;

  if((1 <= k) && (k <= MAX_GRID_TYPE)){
    switch(k){
    case  1: n = N_EGRID1; e = coh_egrid1; break;
    case  2: n = N_EGRID2; e = coh_egrid2; break;
    case  3: n = N_EGRID3; e = coh_egrid3; break;
    case  4: n = N_EGRID4; e = coh_egrid4; break;
    default: n = 0       ; e = coh_egrid1; break;
    }
    for(int i=0 ; i<MAX_EINCIDENT ; i++) ein[i] = (i < n) ? e[i] : 0.0;
  }
  else{
    message << "default energy grid option " << k << " not defined";
    cohTerminateCode("cohEgrid");
  }

  return(n);
}


/**********************************************************/
/*      Copy Energy Grid for MACS calculation             */
/**********************************************************/
static inline int macsEGrid(Particle p)
{
  /*** photon incident case */
  if( p == gammaray ){
    for(int i=0 ; i<MAX_EINCIDENT ; i++){
      ein[i] = (i < N_MACS_GGRID) ? macs_ggrid[i] : 0.0;
    }
  }
  /*** neutron (or other particles) incident case */
  else{
    /*** negative energies given stand for CMS */
    for(int i=0 ; i<MAX_EINCIDENT ; i++){
      ein[i] = (i < N_MACS_EGRID) ? -macs_egrid[i] : 0.0;
    }
  }
  ctl.macs = true;

  return(N_MACS_EGRID);
}

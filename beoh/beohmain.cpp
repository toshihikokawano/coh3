/******************************************************************************/
/*  beohmain.cpp                                                              */
/*        BeoH main calculation loop for beta or statistical decay mode       */
/******************************************************************************/

#include <iostream>

#include "beoh.h"
#include "beohoutput.h"
#include "global.h"
#include "statmodel.h"
#include "beohstat.h"
#include "nucleus.h"
#include "setupinit.h"
#include "terminate.h"

static int  beohReadInputFile (System *, Pdata *, GDR *, Beta *, BData *, FFragData *);
static void beohDecayMain (System *, Pdata *, GDR *, Beta *, BData *, const unsigned long);
static void beohBetaSetting (System *, Beta *, Beta *, BetaProfile *);
static void beohStoreBetaPopulation (System *, Beta *, Beta *, BetaProfile *);
static void beohStoreStatPopulation (BData *);
static void beohStatCalculation (System *, Pdata *, const unsigned long);
static void beohPrintResult (System *, BData *);


static char buf[256];
static bool fisdat = false;

static Spectra       spc;
static Transmission  **tc;
static Transmission  **td;
static double        **tg;

/**********************************************************/
/*      Loop Over Input File                              */
/**********************************************************/
void beohMainLoop(const unsigned long nsim)
{
  System      sys;                  // system parameters
  BData       bdt;                  // BeoH specific system parameters
  Pdata       pdt[MAX_CHANNEL];     // emitting particle data 
  GDR         gdr[MAX_GDR];         // photon strength function parameters
  Beta        gts;                  // beta strength distribution from inputs
  FFragData   fdt;                  // fission fragment decay data


//---------------------------------------
//      Initial Setting

  /*** store parameters which are constant throughtout the calculation */
  setupInitSystem(MAX_CHANNEL,pdt);


//---------------------------------------
//      Main Loop

  int sec = 0;
  do{
    if ((sec = readHead(buf)) < 0) break;

    if(sec == BEGIN){

      /*** initialize all data */
      cohResetNucleus();
      cohResetNuclearStructure();
      setupClearParameter(MAX_CHANNEL,&sys,pdt,gdr);

      bdt.init();
      fdt.init();

      /*** read input file */
      if(beohReadInputFile(&sys,pdt,gdr,&gts,&bdt,&fdt) < 0) break;

      /*** fission fragment decay mode */
      if(bdt.getMode() == fissiondecay || bdt.getMode() == fissionspec || bdt.getMode() == cumulativeyield){
        beohFissionFragmentMain(&sys,&bdt,&fdt,pdt[neutron].omp);
      }
      /*** beta-decay or statistical decay mode */
      else{
        beohDecayMain(&sys,pdt,gdr,&gts,&bdt,nsim);
      }
    }
  }while(!std::cin.eof());
}


/**********************************************************/
/*      Mail Calculation for Beta and Stat Decay Modes    */
/**********************************************************/
void beohDecayMain(System *sys, Pdata *pdt, GDR *gdr, Beta *gts, BData *bdt, const unsigned long nsim)
{
  Beta        ens;    // beta strength distribution from ENSDF
  BetaProfile bpf;    // mixed beta strength distribution

  int jmax = MAX_J-1;

  /*** allocate memory */
  beohAllocateMemory();

  if(bdt->exmean > 0.0) sys->ex_total = bdt->exmean + EXRANGE_FACTOR * bdt->exwidth;

  setupGeneralParameter(bdt->isBetaCalc(),sys);

  setupStatModel(bdt->isBetaCalc(),jmax,sys,pdt);

  beohParameterSetting(sys,gdr,false);

  /*** print system parameters */
  if(prn.system){
    outTitle(buf);
    outSystem(bdt->getMode(),sys);
    outLevelDensity(sys->max_compound,0.0);
    outGDR(false,sys->max_compound);
    if(fisdat) outFissionBarrier(sys->max_compound);
    outCompound(sys->max_compound,pdt);
  }
 
  if(bdt->isBetaCalc()){
    /*** set up beta-decay strength */
    beohStoreBetaPopulation(sys,gts,&ens,&bpf);

    /*** print beta-decay strength parameters */
    if(prn.system){
      outTargetState(&nst[0],sys->target_level);
      outBeta(&sys->compound,gts,&ens);
    }
  }
  else{
    /*** store initial population for statistical decay */
    beohStoreStatPopulation(bdt);
  }

  /*** statistical decay main calculation */
  beohStatCalculation(sys,pdt,nsim);

  /*** electron and neutrino spectra */
  if(bdt->isBetaCalc()){
    beohElectronSpectrum(ncl[0].max_energy,bpf.rcont,bpf.rdisc,&ncl[0],crx.spectra[MAX_CHANNEL],crx.spectra[MAX_CHANNEL+1]);
  }

  /*** output results */
  beohPrintResult(sys,bdt);

  /*** free allocated memory */
  beohDeleteAllocated();
}


/**********************************************************/
/*      BeoH for CGM Compatible Mode                      */
/**********************************************************/
void beohCGMCompatibleMode(const int targZ, const int targA, const double iniJ, const int iniP, const double energy, const double exwidth, const double spinfact, const double targE, const double beta2, const unsigned long nsim)
{
  System sys;
  BData  bdt;
  Pdata  pdt[MAX_CHANNEL];
  GDR    gdr[MAX_GDR];
  int    jmax = MAX_J-1;

  beohAllocateMemory();
  setupInitSystem(MAX_CHANNEL,pdt);
  setupClearParameter(MAX_CHANNEL,&sys,pdt,gdr);

  sys.compound.setZA(targZ,targA);
  sys.ex_total     = energy;
  sys.beta2        = beta2;
  pdt[neutron].omp = 6; // Koning-Delaroche potential as default

  /*** set initial condition given by command-line options */
  bdt.setMode(statdecay);
  bdt.initspin     = iniJ;
  bdt.initparity   = iniP;
  bdt.spinfactor   = spinfact;
  bdt.ekmean       = targE;

  /*** when excitation energy is spread, max energy is set to
       mean_energy + factor times spread_width
       where the factor is defined in beoh.h */
  bdt.exmean       = energy;
  bdt.exwidth      = exwidth;
  sys.ex_total     = bdt.exmean + EXRANGE_FACTOR * bdt.exwidth;

  setupGeneralParameter(bdt.isBetaCalc(),&sys);
  setupStatModel(false,jmax,&sys,pdt);
  beohParameterSetting(&sys,gdr,false);

  if(prn.system){
    outSystem(bdt.getMode(),&sys);
    outLevelDensity(sys.max_compound,0.0);
    outGDR(false,sys.max_compound);
    outCompound(sys.max_compound,pdt);
  }

  beohStoreStatPopulation(&bdt);
  beohStatCalculation(&sys,pdt,nsim);
  beohPrintResult(&sys,&bdt);

  beohDeleteAllocated();
}


/**********************************************************/
/*      Scan Input Data from BEGIN to END                 */
/**********************************************************/
int beohReadInputFile(System *sys, Pdata *pdt, GDR *gdr, Beta *gts, BData *bdt, FFragData *fdt)
{
  int sec = 0;

  /*** first, mode is determined by the keyword for (Z,A) in the DATA section.
       then this is re-determine by existing BETA or FFRAG section. */
  do{
    if((sec = readHead(buf)) < 0) return(-1);

    /*** DATA section */
    if(sec == DATA){
      readSystem(buf,sys,pdt,fdt,bdt);
    }
    /*** BETA section */
    else if(sec == BETA){
      bdt->setMode(betadecay);
      readBeta(buf,gts);
    }
    /*** HFS section */
    else if(sec == HFS){
      readStatModel(buf,gdr,fbr);
    }
    /*** FFRAG section */
    else if(sec == FFRAG){
      bdt->setMode(fissiondecay);
      readFissionFragment(buf,fdt,bdt);
    }
    /*** MCF section */
    else if(sec == MCF){
      /*** this section will be repeated for all the multi-chance fission blocks */
      int nc = fdt->getFissionChance();
      if((bdt->getMode() == fissiondecay || bdt->getMode() == fissionspec) && (nc < fdt->getMaxFissionChance())){
        readMultiChanceFission(buf,&fdt->mc[nc]);
        fdt->inclFissionChance();
      }
    }
    /*** CFY section */
    else if(sec == CFY){
      bdt->setMode(cumulativeyield);
      readFissionYield(buf,fdt);
    }
    /*** END */
    else if(sec == END){
      break;
    }
  }while(!std::cin.eof());

  return 0;
}


/**********************************************************/
/*      Store Initial Population for Beta-decay Case      */
/**********************************************************/
void beohStoreBetaPopulation(System *sys, Beta *gts, Beta *ens, BetaProfile *bpf)
{
  /*** store beta strengths in arrays */
  beohBetaSetting(sys,gts,ens,bpf);

  /*** beta strength profile in the daughter nucleus */
  beohStrengthProfile(gts,ens,bpf);
  if(prn.system) outBetaProfile(bpf);

  /*** precursor ground or isomeric state spin and parity */
  double pcs = nst[0].lev[sys->target_level].spin;
  int    pcp = nst[0].lev[sys->target_level].parity;

  /*** initial population in the continuum and discrete levels */
  beohClearInitialPopulation(&ncl[0]);
  beohInitialPopulation(pcs,pcp,bpf,&ncl[0]);
}


/**********************************************************/
/*      Store Initial Population for General Stat. Calc.  */
/**********************************************************/
void beohStoreStatPopulation(BData *bdt)
{
  beohClearInitialPopulation(&ncl[0]);

  if(ncl[0].ncont > 0){
    /*** store initial population in the top bin */
    beohInitialPopulationStat(1.0,bdt->initspin,bdt->initparity,bdt->spinfactor,bdt->exmean,bdt->exwidth,&ncl[0]);
  }
  else{
    /*** when no continuum, store all population at the highest level */
    int nstart = beohStoreLevelExcite(1.0,&ncl[0]);
    ncl[0].ndisc = nstart;
  }
}


/**********************************************************/
/*      Main Statistical Decay Calculation                */
/**********************************************************/
void beohStatCalculation(System *sys, Pdata *pdt, const unsigned long nsim)
{
  crx.specclear();
  spc.memclear();
  gml.init();

  if(nsim == 0L){
    /*** regular Hauser-Feshbach calculation */
    if(ncl[0].ncont > 0) beohspectra(sys,pdt,tc,td,tg,&spc,&gml);
    else                 specGammaCascade(spc.cn[0],&ncl[0]);
  }
  else{
    /*** Monte Carlo Hauser-Feshbach calculation */
    beohspectraMC(pdt,tc,td,tg,&spc,nsim);
  }
}


/**********************************************************/
/*      Output Calculated Results                         */
/**********************************************************/
void beohPrintResult(System *sys, BData *bdt)
{
  if(prn.spectra){
    outSpectrum(bdt->isBetaCalc(),ncl[0].de,crx.spectra,&ncl[0]);
    outSpectrumSum(bdt->isBetaCalc(),ncl[0].de,crx.spectra,&ncl[0]);

    if(bdt->ekmean > 0.0){
      double *spl = new double [MAX_ENERGY_BIN];
      double ef = bdt->ekmean / sys->compound.getA();
      beohLabSpectrum(MAX_ENERGY_BIN,ef,beohZeroCut(crx.spectra),ncl[0].de,crx.spectra[1],spl);
      outSpectrumLab(beohZeroCut(spl),ncl[0].de,spl);
      delete [] spl;
    }
  }

  if(prn.xsection){
    outGSProduction(sys->max_compound);
    outPrepTotalResidual(sys->max_compound,sys->uniq_compound);
    if(ncl[0].fissile) outFission(sys->max_compound);
  }

  if(prn.gcascade){
    for(int i=0 ; i<sys->max_compound ; i++) outGammaCascade(1.0,&ncl[i]);
  }
}


/***********************************************************/
/*      Eliminate Zeros at High Energies                   */
/***********************************************************/
int beohZeroCut(double *spc)
{
  int k0 = 0;

  /*** cut zeros at higher energy bins */
  for(int k=MAX_ENERGY_BIN-1 ; k>0 ; k--){
    if(spc[k] > 0.0){
      k0 = k+1;
      break;
    }
  }

  return(k0);
}


int beohZeroCut(double **spc)
{
  int k0 = 0;

  /*** cut zeros at higher energy bins */
  for(int k=MAX_ENERGY_BIN-1 ; k>0 ; k--){

    bool nonzero = false;
    for(int c=0 ; c<MAX_CHANNEL+2; c++) if(spc[c][k] > 0.0) nonzero = true;
    if(nonzero){
      k0 = k+1;
      break;
    }
  }

  return(k0);
}


/**********************************************************/
/*      Initialize Parameters for Compound Reactions      */
/**********************************************************/
void beohParameterSetting(System *sys, GDR *gdr, const bool meta)
{
  /*** bin width */
  if(sys->energy_bin_in == 0.0) sys->energy_bin = ENERGY_BIN;
  else                          sys->energy_bin = sys->energy_bin_in;

  /*** continuum binning and level density */
  for(int i=0 ; i<sys->max_compound ; i++){
    ncl[i].de = sys->energy_bin;
    statSetupEnergyBin(&ncl[i]); 
    statSetupLevelDensity(&ncl[i],&ncl[i].ldp);
  }

  /*** clean discrete levels */
  for(int i=0 ; i<sys->uniq_compound ; i++){
    /*** re-set ndisc, if extra levels are included */
    if(statAddDiscreteLevels(&nst[i])){
      for(int j=0 ; j<sys->max_compound ; j++){
        if(nst[i].za == ncl[j].za){
          ncl[j].ndisc = nst[i].nlevel;
          statAdjustDensityByExtraLevels(&ncl[j]);
        }
      }
    }
  }

  /*** include meta state at higher excitation energies */
  if(meta){
    for(int i=0 ; i<sys->uniq_compound ; i++){
      int k = 0;
      for(k=0 ; k<sys->max_compound ; k++){ if(nst[i].za == ncl[k].za) break; }

      if(statForceIncludeMeta(&nst[i],&ncl[k])){
        for(int j=0 ; j<sys->max_compound ; j++){
          if(nst[i].za == ncl[j].za){
            ncl[j].ndisc = nst[i].nlevel;
            statAdjustContinuumBoundary(&ncl[j]);
          }
        }
      }
      /*** population trap in meta-stable state */
      statPopTrapMeta(&nst[i]);
    }
  }

  /*** GDR parameters for the first CN */
  statSetupGdrParameter(&ncl[0],gdr,sys->beta2);
  /*** all others, we always use systematics, but assume the same deformation */
  for(int i=1 ; i<sys->max_compound ; i++) statSetupGdrSystematics(&ncl[i],sys->beta2);

  /*** fission level density */
  for(int i=0 ; i<sys->max_compound ; i++){
    if(ncl[i].fissile){
      statSetupFissionParameter(ncl[i].fission);
      statSetupFissionLevelDensity(&ncl[i],ncl[i].fission);
      fisdat = true;
    }
  }
}


/**********************************************************/
/*      Initialize Parameters for Beta Decay              */
/**********************************************************/
void beohBetaSetting(System *sys, Beta *gts, Beta *ens, BetaProfile *bpf)
{
  /*** read ENSDF */
  beohENSDFRead(&sys->compound, ens);

  /*** renormalize the probability */
  double sum = 0.0;
  for(int j=0 ; j<gts->nstate ; j++) sum += gts->br[j].getR();
  if(sum==0.0) gts->nstate = 0;

  for(int j=0 ; j<gts->nstate ; j++) gts->br[j].scaleR(1.0/sum);

  if( (gts->nstate == 0) && (ens->nstate == 0) ){
    message << "beta-decay strength = 0";
    cohTerminateCode("beohBetaSetting");
  }

  /*** beta strength profile setup */
  bpf->lev        = ncl[0].lev;
  bpf->ntotal     = ncl[0].ntotal;
  bpf->ncont      = ncl[0].ncont;
  bpf->ndisc      = ncl[0].ndisc;
  bpf->excitation = ncl[0].excitation;
}


/**********************************************************/
/*      Allocate T and Spec Array                         */
/**********************************************************/
void beohAllocateMemory(void)
{
  try{
    /*** allocate particle transmission arrays */
    tc = new Transmission * [MAX_CHANNEL];
    td = new Transmission * [MAX_CHANNEL];
    for(int j=1 ; j<MAX_CHANNEL ; j++){
      tc[j] = new Transmission [MAX_ENERGY_BIN];
      td[j] = new Transmission [MAX_LEVELS    ];

      for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[j][k].memalloc(MAX_J);
      for(int k=0 ; k<MAX_LEVELS     ; k++) td[j][k].memalloc(MAX_J);
    }

    /*** allocate photon transmission array */
    tg = new double * [MAX_ENERGY_BIN];
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tg[k] = new double [MAX_MULTIPOL];

    /*** allocate energy spectra arrays */
    spc.memalloc(MAX_CHANNEL,MAX_ENERGY_BIN,MAX_LEVELS);
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("beohAllocateMemory");
  }
}


/**********************************************************/
/*      Free Allocated T and Spec Array                   */
/**********************************************************/
void beohDeleteAllocated(void)
{
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[j][k].memfree();
    for(int k=0 ; k<MAX_LEVELS     ; k++) td[j][k].memfree();
    delete [] tc[j];
    delete [] td[j];
  }
  delete [] tc;
  delete [] td;

  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) delete [] tg[k];
  delete [] tg;

  spc.memfree();
}

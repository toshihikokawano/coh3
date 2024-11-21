/******************************************************************************/
/*  statmodel.cpp                                                             */
/*        main Hauser-Feshbach model calculation                              */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "nucleus.h"
#include "coh.h"
#include "statmodel.h"
#include "eclipse.h"
#include "fns.h"
#include "global.h"
#include "output.h"
#include "outext.h"
#include "terminate.h"
#include "parameter.h"

static void statPreequilibrium (System *, Pdata *);
static void statHauserFeshbach (System *, Pdata *, Direct *);
static void statParameterSetting (System *, Pdata *, GDR *);
static void statInitialPopulation (System *, Level *, Transmission *);
static void statReplaceTransmission (const int, Level *, Nucleus *, int *, Transmission *, double **);
static void statAllocateMemory (void);
static void statDeleteAllocated (void);

static Spectra       spc;          // particle and gamma emission spectra
static Transmission  tin;          // transmission coeff. for incident channel
static Transmission  **tc;         // transmission coeff. for continuum
static Transmission  **td;         // transmission coeff. for discrete levels
static double        **tg;         // transmission coeff. for gamma-ray


/**********************************************************/
/*      Hauser-Feshbach Calculation                       */
/*         Set-up Discrete Level Data, Level Densities,   */
/*         GDR, Initial Population.                       */
/**********************************************************/
int statModel
( System        *sys,              // system parameters
  Pdata         *pdt,              // data for all emitting partiles
  Dcapt         *cap,              // radiative capture parameters
  Direct        *dir,              // discrete levels for which direct c.s. given
  FNSpec        *fns,              // fission neutron spectrum data
  const unsigned long nsim)        // number of Monte Carlo simulations
{

//---------------------------------------
//      Initial Setting

  /*** memory allocation */
  statAllocateMemory();

  /*** parameter setting */
  statParameterSetting(sys,pdt,cap->gdr);

  if(ctl.exclusive) eclAllocateMemory(sys->max_compound,sys->max_channel);

  /*** initial population in the first CN */
  statInitialPopulation(sys,dir->lev,&tin);
  spc.memclear("cn");
  spc.memclear("pe");
  crx.specclear();


//---------------------------------------
//      Statistical Model Main Part

  /*** calculate transmission coefficients for all channels */
  statStoreContinuumTransmission(0,ncl[0].excitation[0],pdt,tc);
  statStoreDiscreteTransmission(0,ncl[0].excitation[0],pdt,td);

  /*** in the case of deformed nucleus, replace Tj by the CC calc. */
  if(ctl.deformed && sys->inc_id != gammaray){
    statReplaceTransmission(sys->target_level,dir->lev,&ncl[sys->target_id],sys->bandidx,td[sys->target_id],crx.transmission);
  }

  /***  direct/semidirect nucleon capture */
  if(ctl.dsd){
    dsdDirectCaptureModel(sys->inc_id,sys->cms_energy,&sys->incident,&sys->target,sys->reduced_mass,cap,sys->beta2,spc.cn[0]);
    specCumulativeSpectra(spc.getCsize(),spc.getNsize(),spc.cn,&ncl[0]);
  }

  /*** preequilibrium spectrum */
  if(ctl.exciton) statPreequilibrium(sys,pdt);

  /*** Hauser-Feshbach */
  if(ctl.statmodel) statHauserFeshbach(sys,pdt,dir);

  if(pex.gammaline){
    /*** if Gaussian width is given, broaden the photon spectrum */
    double gw = parmGetValue(parmBROD);
    extGammaLine(sys->lab_energy,crx.spectra[gammaray],&ncl[0],gw);
  }


//---------------------------------------
//     Exclusive Particle Energy Spectra

  if(ctl.exclusive) eclCalc(sys,spc.pe,nsim);


//---------------------------------------
//     Madland-Nix Fission Spectrum Calculation

  if(ctl.fns) fnsFissionNeutronModel(sys,&pdt[neutron],spc.cn[0],fns);


//---------------------------------------
//      Free Allocated Memory

  statDeleteAllocated();
  if(ctl.exclusive) eclDeleteAllocated(sys->max_compound,sys->max_channel);

  return 0;
}


/**********************************************************/
/*      Main Part for Preequilibrium Calculation          */
/**********************************************************/
void statPreequilibrium(System *sys, Pdata *pdt)
{
  /*** exciton model */
  preqExcitonModel(sys,pdt,tc,td,&spc);

  if(sys->incident.pid != gammaray){
    /*** transfer reaction */
    preqExcitonTransfer(sys,pdt,tc,td,&spc);

    /*** alpha particle knockout */
    preqAlphaKnockout(sys,pdt,tc,td,&spc);
  }

  /*** renormalize pre-equilibrium spectra, if more than 99.9% of reaction cs */
  double xpre = preqTotalPopulation();
  double xdir = 0.0;

  if(ctl.dsd) xdir += crx.dsd;
  for(int i=1 ; i<MAX_DIRECT ; i++){
    if(crx.direct[i] > 0.0) xdir += crx.direct[i];
  }
  if( (xdir+xpre) > 0.999*crx.reaction ){
    double f = (0.999*crx.reaction - xdir)/xpre;
    preqRenormalizePopulation(f);
  }

  /*** energy Spectrum from population */
  preqExcitonSpectra(&spc);

  /*** add to total particle energy spectra */
  specCumulativeSpectra(spc.getCsize(),spc.getNsize(),spc.pe,&ncl[0]);

  if(!ctl.fns){
    if(prn.spectra)  outSpectrum(1,spc.pe,&ncl[0]);
    if(prn.xsection) outSpectrumSum(1,spc.pe,&ncl[0]);
  }
}


/**********************************************************/
/*      Main Part for Hauser-Feshbach Calculation         */
/**********************************************************/
void statHauserFeshbach(System *sys, Pdata *pdt, Direct *dir)
{
  spc.memclear("cn");

  /*** add direct inelastic scattering */
  statSetupDirectInelastic(dir->lev,&ncl[sys->target_id],&crx);
  statDirectSpectra(sys->cms_energy,ncl[sys->inc_id].de,dir->lev,spc.cn[sys->inc_id],crx.direct);
  specCumulativeSpectra(spc.getCsize(),spc.getNsize(),spc.cn,&ncl[0]);

  /*** compound decay chain start */
  if(ctl.exclusive){
    /*** copy PE, DSD, Inel spectra to the population increment */
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      for(int k=0 ; k<MAX_ENERGY_BIN ; k++) spc.dp[j][k] = crx.spectra[j][k];
    }
    spectraExclusive(sys,pdt,&tin,tc,td,tg,dir,&spc);
  }
  else
    spectra(sys,pdt,&tin,tc,td,tg,dir,&spc);

  /*** before output, store re-calculated sigmaR in output area */
  outSetSigmaReaction(sys->max_compound);
  
  /*** output total emission spectra */
  if(prn.spectra & !ctl.fns){
    /*** if Gaussian width is given, broaden the photon spectrum */
    double gw = parmGetValue(parmBROD);
    if(gw > 0.0) specGaussianBroadening(MAX_ENERGY_BIN,ncl[0].de,gw,crx.spectra[0]);
    outSpectrum(2,crx.spectra,&ncl[0]);

    /*** if finegammaspectrum option, print gamma-ray spectrum at finer energy grid */
    if(opt.finegammaspectrum) outSpectrumFineGamma(crx.spectra[gammaray],&gml,&ncl[0]);
  }

  if(prn.xsection){
//  outSpectrumSum(2,crx.spectra,&ncl[0]);

    /*** output production cross section */
    outReaction(sys->max_compound,sys->target_id,sys->excitation);

    /*** calculate isotope production cross section */
    outPrepTotalResidual(sys->max_compound,sys->uniq_compound);

    /*** output fission cross section */
    if(ctl.fission) outFission(sys->max_compound);

    /*** output particle production cross section */
    if(!ctl.fns) outParticleProduction(sys->max_compound,ncl[0].cdt,crx.spectra);
  }
}


/**********************************************************/
/*      Initialize Parameters for Compound Reactions      */
/**********************************************************/
void statParameterSetting(System *sys, Pdata *pdt, GDR *gdr)
{
  /*** bin width */
  if(sys->energy_bin_in == 0.0) sys->energy_bin = statBinWidth(sys->lab_energy);
  else                          sys->energy_bin = sys->energy_bin_in;

  message << "default energy bin " << sys->energy_bin;
  cohNotice("statParameterSetting");

  /*** suppress resonances in the initial CN */
  statCutDiscreteLevels(sys->ex_total,&ncl[0]);

  /*** continuum binning and level density */
  for(int i=0 ; i<sys->max_compound ; i++){
    ncl[i].de = sys->energy_bin;
    statSetupEnergyBin(&ncl[i]);
    statSetupLevelDensity(&ncl[i],&ncl[i].ldp);
  }

  /*** adjust level density, if D0 is given */
  statAdjustD0(ncl[0].cdt[sys->inc_id].binding_energy,&ncl[0],&ncl[sys->target_id]);

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

    /*** if level delete flag is given, remove them from calculation */
    statSkipDiscreteLevel(&nst[i]);
  }
  
  /*** GDR parameters for the first CN */
  statSetupGdrParameter(&ncl[0],gdr,sys->beta2);
  /*** all others, we always use systematics, but assume the same deformation */
  for(int i=1 ; i<sys->max_compound ; i++) statSetupGdrSystematics(&ncl[i],sys->beta2);

  /*** read photo-absorption table when this option is activated */
  if(opt.readphotoabsorption){
    gdrAbsorptionDataRead();
    statSetupResetGdrParameter(gdrAbsorptionDataYsize(),&ncl[0]);
  }

  /*** photon strength function */
  double d0 = 0.0, gg = 0.0;
  if(sys->incident.pid == neutron){
    /*** calculate D0 for the target nucleus */
    double sn = ncl[0].cdt[sys->inc_id].binding_energy;
    d0 = statCalculateD0(sn,&ncl[0],&ncl[sys->target_id]);

    /*** E1 renormalization */
    int jt = (int)(2.0*ncl[sys->target_id].lev[sys->target_level].spin);
    int pt =           ncl[sys->target_id].lev[sys->target_level].parity;
    gg = gdrRenormE1StrengthFunction(pt,jt,sys->target.getA(),d0,sn,&ncl[0],tg);
  }

  /*** if E > neutron separation, ignore fluctuation */
  if(sys->incident.pid == neutron){
    double s2n = ncl[ncl[0].cdt[neutron].next].cdt[neutron].binding_energy;
    if(s2n > 0.0 && sys->cms_energy > s2n){
      ctl.fluctuation = false;
    }else{
      ctl.fluctuation = true;
    }
  }
  else ctl.fluctuation = false;

  message << "width fluctuation is ";
  if(ctl.fluctuation) message << "ON"; else message << "OFF";
  cohNotice("statParameterSetting");

  /*** check Engelbrecht-Weidenmueller transformation option */
  if(opt.ewtransformation){
    ctl.ewtransform = true;
    if(!ctl.deformed)    ctl.ewtransform = false;
    if(!ctl.fluctuation) ctl.ewtransform = false;
    /*** when angular distribution is on, we don't calculate EWT */
    if(prn.angdist)      ctl.ewtransform = false;
  }

  /*** fission level density */
  if(ctl.fission){
    for(int i=0 ; i<sys->max_compound ; i++){
      if(ncl[i].fissile){
        statSetupFissionParameter(ncl[i].fission);
        statSetupFissionLevelDensity(&ncl[i],ncl[i].fission);
      }
    }
  }

  /*** replace calculated level densities by those in files in current directory */
  if(opt.readdensity){
    for(int i=0 ; i<sys->max_compound ; i++)  statReadLevelDensity(&ncl[i]);

    if(ctl.fission){
      for(int i=0 ; i<sys->max_compound ; i++){
        if(ncl[i].fissile){
          statReadFissionLevelDensity(&ncl[i]);
          statSetupFissionLevelDensity(&ncl[i],ncl[i].fission);
        }
      }
    }
  }


//---------------------------------------
//      Output System Parameters

  /*** system parameters */
  if(prn.system){
    outTargetState(sys->target_id,sys->target_level,sys->excitation);
    outCompound(sys->max_compound,pdt);
    outLevelDensity(sys->max_compound,d0);
    outGDR(false,sys->max_compound,gg);
    if(ctl.fission) outFissionBarrier(sys->max_compound);
  }

  /*** discrete levels, level densities, and photon strength function (non-standard) */
  if(pex.gbranch)   extGammaBranch();
  if(pex.cumulevel) extCumulativeLevels(&ncl[0]);
  if(pex.density)   extContinuumDensity(&ncl[0]);
  if(ctl.fission & pex.fisdensity) extFissionDensity(sys->max_compound);
  if(pex.gammastf){
    statStoreContinuumGammaTransmission(0,tg,&ncl[0]);
    extGammaStrengthFunction(ncl[0].ntotal,ncl[0].de,tg);
  }
}


/**********************************************************/
/*      Initial Population for Particle/Photo-Incident    */
/**********************************************************/
void statInitialPopulation(System *sys, Level *dir, Transmission *tin)
{
  double coef = NORM_FACT*PI/(sys->wave_number*sys->wave_number);

  /*** for photo-induced reactions */
  if(sys->incident.pid == gammaray){
    ctl.fluctuation = false;
    ctl.dsd = false;
    tin->tran = crx.transmission[0];
    crx.reaction = photoAbsorption(sys->target_level,sys->target_id,sys->cms_energy,coef,&ncl[sys->target_id],tin);
  }
  /*** for particle induced reactions */
  else{
    tin->lmax = ncl[0].jmax;
    tin->ecms = sys->cms_energy;
    if(!ctl.deformed) tin->tran = crx.transmission[0];
    else{
      /*** look for target state, and obtain target index in CC calculation */
      int  kt = 0;
      for(int k=0 ; k<MAX_DIRECT ; k++){
        if(dir[k].energy == sys->excitation){
          kt = k;
          break;
        }
      }
      tin->tran = crx.transmission[kt];
    }
    statStoreInitPopulation(sys->incident.spin2,sys->target_level,sys->target_id,coef,tin);
  }


  if(sys->incident.pid == gammaray){
    statAdjustInitPopulation(crx.reaction);

    double sigQDB = crx.reaction * photoQDratio();
    double sigGDR = crx.reaction * (1.0 - photoQDratio());

    if(prn.xsection)  outCrossSection(2, sigGDR, sigQDB, crx.reaction);
  }
}


/**********************************************************/
/*      Replace Tlj by Entrance Channel Calculations      */
/*      -----                                             */
/*             This is only needed if the entrance Tlj's  */
/*             are calculated with the CC method.         */
/*             For the spherical case, the same values    */
/*             are re-calcuated in the subroutine         */
/*             statStoreDiscreteTransmission()            */
/**********************************************************/
void statReplaceTransmission(const int targlev, Level *dir, Nucleus *n, int *bidx, Transmission *td, double **transinc)
{
  const double e1 = 0.01, e2 = 1e-6;

  /*** first element should be the elastic channel */
  for(int j=0 ; j<3*MAX_J ; j++) td[targlev].tran[j] = transinc[0][j];
  bidx[targlev] = 0;

  int k = 1;
  for(int i=1 ; i<MAX_DIRECT ; i++){

    if(dir[i].energy <= 0.0) continue;

    /*** search for the levels included in the CC method */
    /*** CC direct cross sections are given by negative values */
    if(crx.direct[i] >= 0.0) continue;

    /*** search for the same spin and energy levels */
    bool found = false;
    for(int j=k ; j<n->ndisc ; j++){
      if( (fabs(     dir[i].energy - n->lev[j].energy)<=e1) &&
          (fabs(fabs(dir[i].spin)  - n->lev[j].spin  )<=e2) ){
        found = true;
        k = j;
        break;
      }
    }

    if(found){
      /*** set coupled-level number in the level array */
      bidx[k] = i;
      /*** copy transmission coefficients */
      for(int j=0 ; j<3*MAX_J ; j++) td[k].tran[j] = transinc[i][j];
    }
    else{
      message << "direct level energy " << dir[i].energy << " did not match";
      cohTerminateCode("statReplaceTransmission");
    }
  }
}


/**********************************************************/
/*      Memory Allocation                                 */
/**********************************************************/
void statAllocateMemory()
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
    crx.spectalloc(MAX_CHANNEL,MAX_ENERGY_BIN);
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("statAllocateMemory");
  }
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void statDeleteAllocated()
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
  crx.spectfree();
}

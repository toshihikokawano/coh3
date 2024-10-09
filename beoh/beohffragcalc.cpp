/******************************************************************************/
/*  beohffragcalc.cpp                                                         */
/*        Main calculation part for fission fragment decay                    */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "beoh.h"
#include "beohoutput.h"
#include "nucleus.h"
#include "masstable.h"
#include "beta2.h"
#include "global.h"
#include "setupinit.h"
#include "terminate.h"

#define CALC_MONITOR

static bool PRINT_INDIV_RESULT        = false;
static bool PRINT_INDIV_SPEC          = false;
static bool PRINT_INDIV_CHI           = false;
static bool PRINT_INDIV_GAMMACASCADE  = false;
static bool PRINT_INDIV_GAMMALINE     = false;

static bool FFPSelectZA (const unsigned int, const unsigned int, ZAPair, FFragData *);
static void FFPMultiChanceYield (const int, const int);
static int  FFPMultiChanceFindIndex (const int, const int, const int, char *);
static int  FFPMultiChanceFindIndexL (const int, const int);
static int  FFPMultiChanceFindIndexH (const int, const int);
static void FFPSetExcitationEnergy (const int, BData *, BData *, FissionFragmentPair *, double *);

static void FFPSplitCalc1 (const int, System *, System *, BData *, BData *, FFragData *, Pdata *);
static void FFPSplitCalc2 (const int, System *, System *, BData *, BData *, FFragData *, Pdata *);
static void FFPPost (const int, System *, BData *, double, FragmentObservable *);
static void FFPCombined (void);
static void FFPCombinedSpectra (void);

#ifdef CALC_MONITOR
static void   FFPCalcMonitor (const int, const int, const int, const int, BData *, BData *);
#endif

extern FissionFragmentPair *ffp;
extern FissionObservable   fob;
extern GammaProduction     gmlL, gmlH;

static double yieldL = 0.0;
static double yieldH = 0.0;


/**********************************************************/
/*      For Each Fission Fragment Pair                    */
/**********************************************************/
void FFPSplit(const CalcMode mode, System *sys, FFragData *fdt, const int ompid)
{
  BData       *bdtL, *bdtH;
  System      sysL, sysH;
  Pdata       pdt[MAX_CHANNEL];

  setupInitSystem(MAX_CHANNEL,pdt);
  pdt[neutron].omp = ompid;

  /*** Energy etc. data for both light and heavy fragments including multi-chance fission component */
  bdtL = new BData [fdt->getFissionChance()];
  bdtH = new BData [fdt->getFissionChance()];

  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    bdtL[nc].init();
    bdtH[nc].init();
    bdtL[nc].setMode(statdecay);
    bdtH[nc].setMode(statdecay);
    bdtL[nc].spinfactor = bdtH[nc].spinfactor = fdt->mc[nc].spinfactor;
    bdtL[nc].fraction   = bdtH[nc].fraction   = fdt->mc[nc].fraction;
  }

  /*** count number of total fragment pairs to be included */
  int ntotal = 0;
  for(int k=0 ; k<ffp[0].getN() ; k++){
    FFPMultiChanceYield(fdt->getFissionChance(),k);
    if( (yieldL + yieldH)*0.5 > fdt->ycutoff ) ntotal ++;
  }

  /*** loop over fission fragment pairs */
  int ncount = 1;
  for(int k=0 ; k<ffp[0].getN() ; k++){

    /*** if Z and/or A is selected */
    if(FFPSelectZA(sys->compound.getZ(),sys->compound.getA(),ffp[0].fragment[k],fdt)) continue;

    FFPMultiChanceYield(fdt->getFissionChance(),k);
    if( (yieldL + yieldH)*0.5 <= fdt->ycutoff ) continue;

    /*** Z,A for each fragment */
    sysL.init();
    sysH.init();

    sysL.compound.setZA(ffp[0].fragment[k].getZl(),ffp[0].fragment[k].getAl());
    sysH.compound.setZA(ffp[0].fragment[k].getZh(),ffp[0].fragment[k].getAh());

    /*** find beta2 deformation parameter */
    sysL.beta2 = beta2Read(&sysL.compound);
    sysH.beta2 = beta2Read(&sysH.compound);

    /*** copy energy bin width */
    sysL.energy_bin_in = sys->energy_bin_in;
    sysH.energy_bin_in = sys->energy_bin_in;

    /*** main calculation for each pair */
    if(mode != fissionspec){
      FFPSplitCalc1(k,&sysL,&sysH,bdtL,bdtH,fdt,pdt);
      /*** postprocess for L and H combined quantities */
      FFPCombined();
    }
    else{
      FFPSplitCalc2(k,&sysL,&sysH,&bdtL[0],&bdtH[0],fdt,pdt);
      /*** total spectra from L and H */
      FFPCombinedSpectra();
    }

#ifdef CALC_MONITOR
    /*** calculation monitor */
    FFPCalcMonitor(ntotal,ncount,k,fdt->getFissionChance(),bdtL,bdtH);
#endif

    /*** print individual result if needed */
    if(ncl[0].de > 0.0){
      if(PRINT_INDIV_SPEC){
        
        FFPOutputIndividualSpectrum(&ffp[0].fragment[k],fob.getNsize(),ncl[0].de,fob.L.spectrum[gammaray],fob.H.spectrum[gammaray],fob.spectrum[gammaray]);

        FFPOutputIndividualSpectrum(&ffp[0].fragment[k],fob.getNsize(),ncl[0].de,fob.L.spectrum[neutron],fob.H.spectrum[neutron],fob.spectrum[neutron]);
      }

      if(PRINT_INDIV_CHI) FFPOutputIndividualSpectrum(&ffp[0].fragment[k],fob.getNsize(),ncl[0].de,fob.L.speclab,fob.H.speclab,fob.chi);
    }

    ncount ++;
  }

  /*** add prefission neutrons */
  fob.eprefis = 0.0;
  fob.nprefis = 0.0;
  for(int nc=1 ; nc < fdt->getFissionChance() ; nc++){
    fob.nprefis += nc * fdt->mc[nc].fraction;
    fob.eprefis += nc * fdt->mc[nc].fraction * fdt->mc[nc].eprefis;
  }

  delete [] bdtL;
  delete [] bdtH;
}


/**********************************************************/
/*      Selected Fragment Calculation                     */
/**********************************************************/
bool FFPSelectZA(const unsigned int zcn, const unsigned int acn, ZAPair frag, FFragData *fdt)
{
  bool skip = false;
  if( (fdt->selectA == 0) && (fdt->selectZ == 0) ) return skip;

  if(fdt->selectA > 0){
    int asl = fdt->selectA;
    if(asl > (int)acn / 2) asl = acn - asl;

    /*** both A and Z selected  */
    if(fdt->selectZ > 0){
      int zsl = fdt->selectZ;
      if(zsl > (int)zcn / 2) zsl = zcn - zsl;
      if(frag.getZl() != (unsigned int)zsl || frag.getAl() != (unsigned int)asl) skip = true;
    }
    /*** A selected but Z not */
    else{
      if(frag.getAl() != (unsigned int)asl) skip = true;
    }
  }
  /*** only Z selected */
  else{
    if(fdt->selectZ > 0){
      int zsl = fdt->selectZ;
      if(zsl > (int)zcn / 2) zsl = zcn - zsl;
      if(frag.getZl() != (unsigned int)zsl) skip = true;
    }
  }

  return skip;
}


/**********************************************************/
/*      Fragment Yield Including MultiChance Components   */
/**********************************************************/
void FFPMultiChanceYield(const int nc, const int k0)
{
  if(nc == 1){
    yieldL = yieldH = ffp[0].fragment[k0].yield;
    return;
  }
  else{
    yieldL = yieldH = ffp[0].fragment[k0].yield * ffp[0].fraction;
  }

  char z;
  for(int i=1 ; i<nc ; i++){
    /*** for the light fragemnt */
    int kl = FFPMultiChanceFindIndex(i,ffp[0].fragment[k0].getZl(),ffp[0].fragment[k0].getAl(),&z);
    if(kl >= 0) yieldL += ffp[i].fragment[kl].yield * ffp[i].fraction;
    /*** for the heavy fragemnt */
    int kh = FFPMultiChanceFindIndex(i,ffp[0].fragment[k0].getZh(),ffp[0].fragment[k0].getAh(),&z);
    if(kh >= 0) yieldH += ffp[i].fragment[kh].yield * ffp[i].fraction;
  }
}


/**********************************************************/
/*      Index for Same Z,A in Multi-Chance List           */
/**********************************************************/
int FFPMultiChanceFindIndex(const int c, const int z0, const int a0, char *d)
{
  ZAnumber za0(z0,a0);

  int kx = -1;
  *d = 'X';

  for(int k1=0 ; k1<ffp[c].getN() ; k1++){
    ZAnumber zaL(ffp[c].fragment[k1].getZl(),ffp[c].fragment[k1].getAl());
    ZAnumber zaH(ffp[c].fragment[k1].getZh(),ffp[c].fragment[k1].getAh());

    if(za0 == zaL){ kx = k1; *d = 'L'; break;}
    if(za0 == zaH){ kx = k1; *d = 'H'; break;}
  }
  return kx;
}


int FFPMultiChanceFindIndexL(const int c, const int k0)
{
  int k = -1;
  for(int k1=0 ; k1<ffp[c].getN() ; k1++){
    if( (ffp[0].fragment[k0].getZl() == ffp[c].fragment[k1].getZl()) &&
        (ffp[0].fragment[k0].getAl() == ffp[c].fragment[k1].getAl()) ){ k = k1; break; }
  }
  return k;
}

int FFPMultiChanceFindIndexH(const int c, const int k0)
{
  int k = -1;
  for(int k1=0 ; k1<ffp[c].getN() ; k1++){
    if( (ffp[0].fragment[k0].getZh() == ffp[c].fragment[k1].getZh()) &&
        (ffp[0].fragment[k0].getAh() == ffp[c].fragment[k1].getAh()) ){ k = k1; break; }
  }
  return k;
}


/**********************************************************/
/*      Fragment Pair for Yield Calculation               */
/*      Energy integration performed by energy            */
/*      distribution of initial population                */
/**********************************************************/
void FFPSplitCalc1(const int k, System *sysL, System *sysH, BData *bdtL, BData *bdtH, FFragData *fdt, Pdata *pdt)
{
  /*** for the first chance fission */
  BData bL, bH;
  if((ffp[0].fragment[k].el > 0.0) && (ffp[0].fragment[k].eh > 0.0)){
    /*** when data are given in a file */
    bL.exmean  = ffp[0].fragment[k].el;
    bL.exwidth = ffp[0].fragment[k].wl;
    bH.exmean  = ffp[0].fragment[k].eh;
    bH.exwidth = ffp[0].fragment[k].wh;

    /*** average kinetic energies */
    double x = ffp[0].fragment[k].getAl() + ffp[0].fragment[k].getAh();
    bL.ekmean = ffp[0].fragment[k].ek * ffp[0].fragment[k].getAh() / x;
    bH.ekmean = ffp[0].fragment[k].ek * ffp[0].fragment[k].getAl() / x;
  }
  else{
    /*** use RT model to calculate energy sharing */
    FFPSetExcitationEnergy(k,&bL,&bH,&ffp[0],fdt->mc[0].rta);
  }
  bdtL[0].exmean  = sysL->ex_total = bL.exmean;
  bdtH[0].exmean  = sysH->ex_total = bH.exmean;
  bdtL[0].exwidth = bL.exwidth;
  bdtH[0].exwidth = bH.exwidth;
  bdtL[0].ekmean  = bL.ekmean;
  bdtH[0].ekmean  = bH.ekmean;

  /*** multi-chance fission case */
  char z;
  for(int nc=1 ; nc < fdt->getFissionChance() ; nc++){

    /*** find the same ZA in the multi-chance chain */
    int kl = FFPMultiChanceFindIndex(nc,ffp[0].fragment[k].getZl(),ffp[0].fragment[k].getAl(),&z);
    FFPSetExcitationEnergy(kl,&bL,&bH,&ffp[nc],fdt->mc[nc].rta);

    if(z == 'L'){
      bdtL[nc].exmean  = bL.exmean;  if(bL.exmean > sysL->ex_total) sysL->ex_total = bL.exmean;
      bdtL[nc].exwidth = bL.exwidth;
      bdtL[nc].ekmean  = bL.ekmean;
    }
    else{
      bdtL[nc].exmean  = bH.exmean;  if(bH.exmean > sysL->ex_total) sysL->ex_total = bH.exmean;
      bdtL[nc].exwidth = bH.exwidth;
      bdtL[nc].ekmean  = bH.ekmean;
    }

    int kh = FFPMultiChanceFindIndex(nc,ffp[0].fragment[k].getZh(),ffp[0].fragment[k].getAh(),&z);
    FFPSetExcitationEnergy(kh,&bL,&bH,&ffp[nc],fdt->mc[nc].rta);
    if(z == 'L'){
      bdtH[nc].exmean  = bL.exmean;  if(bL.exmean > sysH->ex_total) sysH->ex_total = bL.exmean;
      bdtH[nc].exwidth = bL.exwidth;
      bdtH[nc].ekmean  = bL.ekmean;
    }
    else{
      bdtH[nc].exmean  = bH.exmean;  if(bH.exmean > sysH->ex_total) sysH->ex_total = bH.exmean;
      bdtH[nc].exwidth = bH.exwidth;
      bdtH[nc].ekmean  = bH.ekmean;
    }
  }

  /*** exec Hauser-Feshbach for each fragment, post processing, collecting calculated results */
  fob.L.clear();
  gmlL.init();
  beohFFDecay(fdt->getFissionChance(),sysL,pdt,bdtL,&gmlL,PRINT_INDIV_RESULT);
  FFPPost(fdt->getFissionChance(),sysL,bdtL,yieldL,&fob.L);

  if(PRINT_INDIV_GAMMACASCADE){
    for(int i=0 ; i<sysL->max_compound ; i++) outGammaCascade(yieldL,&ncl[i]);
  }
  if(PRINT_INDIV_GAMMALINE && opt.finegammaspectrum) outDiscreteGamma(yieldL,&gmlL);

  fob.H.clear();
  gmlH.init();
  beohFFDecay(fdt->getFissionChance(),sysH,pdt,bdtH,&gmlH,PRINT_INDIV_RESULT);
  FFPPost(fdt->getFissionChance(),sysH,bdtH,yieldH,&fob.H);

  if(PRINT_INDIV_GAMMACASCADE){
    for(int i=0 ; i<sysH->max_compound ; i++) outGammaCascade(yieldH,&ncl[i]);
  }
  if(PRINT_INDIV_GAMMALINE && opt.finegammaspectrum) outDiscreteGamma(yieldH,&gmlH);
}


/**********************************************************/
/*      Store Average Excitation Energies                 */
/**********************************************************/
void FFPSetExcitationEnergy(const int k, BData *bL, BData *bH, FissionFragmentPair *f, double *rt)
{
  if(k < 0){
    bL->exmean = bL->exwidth = bL->ekmean = 0.0;
    bH->exmean = bH->exwidth = bH->ekmean = 0.0;
    return;
  }

  /*** for the first chance fission */
  ZAnumber cL, cH;
  cL.setZA(f->fragment[k].getZl(),f->fragment[k].getAl());
  cH.setZA(f->fragment[k].getZh(),f->fragment[k].getAh());

  /*** energy sorting by the RT model */
  double eratio =0.0;
  // pull the correct RT values, depending on whether it was written in the file or not
  if (rt[0] > 0){
    int ahmin = rt[0];
    int ahmax = rt[1];
    int ah = f->fragment[k].getAh();
    int aind = ah - ahmin;
    if (aind < 0) aind = ahmax - ahmin;

    eratio = beohExcitationEnergyDivide(&cL,&cH,rt[2+aind],f->fragment[k].ex);
  } 
  else {
    eratio = beohExcitationEnergyDivide(&cL,&cH,rt[2],f->fragment[k].ex);
  }

  if(eratio > 0.0 && f->fragment[k].ex > 0.0){
    bL->exmean = f->fragment[k].ex * eratio;
    bH->exmean = f->fragment[k].ex * (1.0 - eratio);
  }
  else{
    bL->exmean = 0.0;
    bH->exmean = 0.0;
  }

  /*** when dE is given for TXE, divide dE into two fragments */
  if(f->fragment[k].sigma > 0.0){
    double d = sqrt(bL->exmean * bL->exmean + bH->exmean * bH->exmean);
    double c = 0.0;
    if(d > 0.0) c = ffp[0].fragment[k].sigma / d;
    bL->exwidth = c * bL->exmean;
    bH->exwidth = c * bH->exmean;
  }
  else{
    bL->exwidth = 0.0;
    bH->exwidth = 0.0;
  }

  /*** average kinetic energies */
  double x = f->fragment[k].getAl() + f->fragment[k].getAh();
  bL->ekmean = f->fragment[k].ek * f->fragment[k].getAh() / x;
  bH->ekmean = f->fragment[k].ek * f->fragment[k].getAl() / x;
}


/**********************************************************/
/*      Fragment Pair for Fission Spectrum Calculation    */
/*      Perform exact energy intergration                 */
/**********************************************************/
void FFPSplitCalc2(const int k, System *sysL, System *sysH, BData *bdtL, BData *bdtH, FFragData *fdt, Pdata *pdt)
{
  double *speclab = new double [fob.getNsize()];
  double *tkedist = new double [fob.getKsize()];
  double *rt = fdt->mc[0].rta;
  double rtval = 0;

  for(int j=0 ; j<fob.getNsize() ; j++){
    fob.L.speclab[j] = fob.H.speclab[j] = speclab[j] = 0.0;
    for(int p=0 ; p<2 ; p++) fob.L.spectrum[p][j] = fob.H.spectrum[p][j] = 0.0;
  }

  for(int i=0 ; i<fob.getKsize() ; i++){
    fob.L.nTKE[i] = fob.L.yTKE[i] = fob.H.nTKE[i] = fob.H.yTKE[i] = tkedist[i] = 0.0;
  }

  /*** determine Jmax for each fragment at the mean energy */
  FFPSetExcitationEnergy(k,bdtL,bdtH,&ffp[0],fdt->mc[0].rta);
  sysL->ex_total = bdtL->exmean;
  sysH->ex_total = bdtH->exmean;

  // pull the correct RT values, depending on whether it was written in the file or not
  if (rt[0] > 0){
    int ahmin = rt[0];
    int ahmax = rt[1];
    int ah = ffp[0].fragment[k].getAh();
    int aind = ah - ahmin;
    if (aind < 0) aind = ahmax - ahmin;
    rtval = rt[2+aind];
  } 
  else {
    rtval = rt[2];
  }

  int jmaxL = beohFFDecayPrep(sysL,pdt);
  int jmaxH = beohFFDecayPrep(sysH,pdt);

  /*** kinetic energy increment */
  int kmax = ffp[0].fragment[k].ek + EXRANGE_FACTOR * ffp[0].fragment[k].sigma;
  int kmin = ffp[0].fragment[k].ek - EXRANGE_FACTOR * ffp[0].fragment[k].sigma; if(kmin < 0) kmin = 0;
  int nk   = (double)(kmax - kmin) + 1;

  /*** Gaussian weight */
  double sg = 0.0;
  for(int i=0 ; i<nk ; i++){
    double ek = kmin + i;
    if(ffp[0].fragment[k].ex + ffp[0].fragment[k].ek - ek < 0.0) {nk = i; break;}
    tkedist[i] = exp(-(ek - ffp[0].fragment[k].ek)*(ek - ffp[0].fragment[k].ek)
               /(2.0*ffp[0].fragment[k].sigma*ffp[0].fragment[k].sigma));
    sg += tkedist[i];
  }
  if(sg > 0.0){ for(int i=0 ; i<nk ; i++) tkedist[i] /= sg; }

  /*** for each TKE bin */
  double nubarL = 0.0, nubarH = 0.0;
  for(int i=0 ; i<nk ; i++){
    double tke = (double)(kmin + i);
    double txe = ffp[0].fragment[k].ex + ffp[0].fragment[k].ek - tke;

//  if(tke < 100.0 || tke > 120.0) continue; // adhoc

    /*** divide TKE and TXE into two fragments */
    double eratio = beohExcitationEnergyDivide(&sysL->compound,&sysH->compound,rtval,txe);
    bdtL->exmean = txe * eratio;
    bdtH->exmean = txe * (1.0 - eratio);

    bdtL->ekmean = tke * sysH->compound.getA() / (sysL->compound.getA() + sysH->compound.getA());
    bdtH->ekmean = tke * sysL->compound.getA() / (sysL->compound.getA() + sysH->compound.getA());

    /*** statistical decay calculation and postprocessing for the light fragment */
    beohFFDecaySpec(sysL,pdt,bdtL,&gmlL,jmaxL,fob.getNsize(),speclab);
    beohAverageEnergy(MAX_ENERGY_BIN,fob.L.multiplicity,fob.L.eaverage,fob.L.etotal,crx.spectra);
    for(int j=0 ; j<fob.getNsize() ; j++){
      for(int p=0 ; p<2 ; p++) fob.L.spectrum[p][j] += crx.spectra[p][j] * tkedist[i];
      fob.L.speclab[j] += speclab[j] * tkedist[i];
    }

    /*** statistical decay calculation and postprocessing for the heavy fragment */
    beohFFDecaySpec(sysH,pdt,bdtH,&gmlH,jmaxH,fob.getNsize(),speclab);
    beohAverageEnergy(MAX_ENERGY_BIN,fob.H.multiplicity,fob.H.eaverage,fob.H.etotal,crx.spectra);
    for(int j=0 ; j<fob.getNsize() ; j++){
      for(int p=0 ; p<2 ; p++) fob.H.spectrum[p][j] += crx.spectra[p][j] * tkedist[i];
      fob.H.speclab[j] += speclab[j] * tkedist[i];
    }

    fob.L.nTKE[kmin+i] = fob.L.multiplicity[neutron];
    fob.H.nTKE[kmin+i] = fob.H.multiplicity[neutron];

    nubarL += fob.L.multiplicity[neutron] * tkedist[i];
    nubarH += fob.H.multiplicity[neutron] * tkedist[i];
  }

  fob.L.multiplicity[neutron] = nubarL;
  fob.H.multiplicity[neutron] = nubarH;

  for(int i=0 ; i<nk ; i++){
    fob.yTKE[kmin+i] += ffp[0].fragment[k].yield;
    fob.nTKE[kmin+i] += (fob.L.nTKE[kmin+i] + fob.H.nTKE[kmin+i]) * ffp[0].fragment[k].yield;
  }

  /*** just in case if TKE look above was skipped */
  cohResetNucleus();
  cohResetNuclearStructure();

  delete [] speclab;
  delete [] tkedist;
}


/**********************************************************/
/*      Accumulate Decay Calculation Results              */
/**********************************************************/
void FFPPost(const int nc, System *sys, BData *bdt, const double yield, FragmentObservable *fo)
{
  /*** photon and neutron spectra */
  for(int p=0 ; p<2 ; p++){
    for(int i=0 ; i<fob.getNsize() ; i++) fo->spectrum[p][i] = crx.spectra[p][i];
  }

  /*** simple conversion of neutron cms-spectrum to lab */
  for(int i=0 ; i<fob.getNsize() ; i++) fo->speclab[i] = crx.spectra[neutron][i];
  beohIntegCmsSpectrum(fob.getNsize(),nc,sys->compound.getA(),bdt,ncl[0].de,fo->speclab);

  /*** average energy and multiplicity */
  beohAverageEnergy(fob.getNsize(),fo->multiplicity,fo->eaverage,fo->etotal,fo->spectrum);

  /*** neutron multiplicity distributions, P(nu) */
  for(int i=0 ; i<fob.getLsize() ; i++) fo->Pn[i] = 0.0;
  int a0 = ncl[0].za.getA();
  for(int j=0 ; j<sys->max_compound ; j++){
    int da = a0 - ncl[j].za.getA();
    if(da < fob.getLsize()) fo->Pn[da] = yield * crx.prod[j].xsec;
  }

  /*** multiplicity and average energy distributions as function of A */
  for(int p=0 ; p<2 ; p++){
    fob.multiplicitydist[p][sys->compound.getA()] += yield * fo->multiplicity[p];
    fob.eaveragedist[p][sys->compound.getA()]     += yield * fo->eaverage[p];
  }

  /*** pre neutron emission mass yield */
  fob.preMassYield[sys->compound.getA()] += yield;

  /*** post neutron emission mass yield */
  for(int j=0 ; j<sys->max_compound ; j++){
    fob.postMassYield[ncl[j].za.getA()] += yield * crx.prod[j].xsec;
  }

  /*** add long-lived isotope productions */
  outPrepCumulativeResidual(sys->max_compound,yield);

  /*** average spin of initial compound */
  fob.javeragedist[sys->compound.getA()] += yield * bdt[0].averagespin;
}


/**********************************************************/
/*      Postprocess for L and H Combined Quantities       */
/**********************************************************/
void FFPCombined(void)
{
  /*** average photon and neutron multiplicities */
  for(int p=0 ; p<2 ; p++){
    double mtot = fob.L.multiplicity[p] + fob.H.multiplicity[p];
    double etot = yieldL * fob.L.multiplicity[p] * fob.L.eaverage[p] + yieldH * fob.H.multiplicity[p] * fob.H.eaverage[p];
    double eave = (mtot > 0.0) ? etot / mtot : 0.0;

    fob.multiplicity[p] += yieldL * fob.L.multiplicity[p] + yieldH * fob.H.multiplicity[p];
    fob.eaverage[p]     += eave;
    fob.etotal[p]       += etot;
  }

  /*** calculate neutron emission probability distribution
       P(0) = Pl(0) * Ph(0), P(1) = Pl(0) * Ph(1) + Pl(1) * Ph(0),
       P(2) = Pl(0) * Ph(2) + Pl(1) * Ph(1) + Pl(0) * Ph(2), ... */
  for(int i=0 ; i<fob.getLsize() ; i++){
    for(int j=0 ; j<=i ; j++){
      fob.Pn[i] += fob.L.Pn[j] * fob.H.Pn[i-j];
    }
  }

  /*** average spectra */
  FFPCombinedSpectra();
}


void FFPCombinedSpectra(void)
{
  /*** average spectra */
  for(int i=0 ; i<fob.getNsize() ; i++){
    for(int p=0 ; p<2 ; p++){
      fob.L.spectrum[p][i] *= yieldL;
      fob.H.spectrum[p][i] *= yieldH;
      fob.spectrum[p][i]   += fob.L.spectrum[p][i] + fob.H.spectrum[p][i];
    }
    fob.L.speclab[i] *= yieldL;
    fob.H.speclab[i] *= yieldH;
    fob.chi[i]       += fob.L.speclab[i] + fob.H.speclab[i];
  }

  /*** fine energy grid gamma */
  if(opt.finegammaspectrum){
    for(int i=0 ; i<gmlL.getN() ; i++){
      gml.push(gmlL.line[i].za, gmlL.line[i].energy, gmlL.line[i].production * yieldL);
    }
    for(int i=0 ; i<gmlH.getN() ; i++){
      gml.push(gmlH.line[i].za, gmlH.line[i].energy, gmlH.line[i].production * yieldH);
    }

    if(gml.getN() == MAX_GAMMALINE){
      message << "total number of discrete gamma lines exceeded the maximum of " << MAX_GAMMALINE;
      cohNotice("FFPCombinedSpectra");
    }
  }
}


/**********************************************************/
/*      Calculation Monitor                               */
/**********************************************************/
#ifdef CALC_MONITOR
static bool first_call= true;
void FFPCalcMonitor(const int m, const int n, const int k, const int nc, BData *bl, BData *bh)
{
  static bool detailed = true;

  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4);

  if(first_call){
    std::cout << "# Ntotal = " << m << std::endl;
    std::cout << "#    n    k  Zl  Al  Zh  Ah";
    std::cout << "  YieldL      YieldH      TKE[MeV]    TXE[MeV]    El[MeV]     Wl[MeV]     Eh[MeV]     Wh[MeV]" << std::endl;
    first_call = false;
  }

  double ek = 0.0, ex = 0.0;
  double blm = 0.0, blw = 0.0, bhm = 0.0, bhw = 0.0;
  if(nc == 1){
    ex = ffp[0].fragment[k].ex;
    ek = ffp[0].fragment[k].ek;
    blm = bl[0].exmean;
    bhm = bh[0].exmean;
    blw = bl[0].exwidth;
    bhw = bh[0].exwidth;
  }
  else{
    for(int i=0 ; i<nc ; i++){
      ex  += ffp[i].fraction * ffp[i].fragment[k].ex;
      ek  += ffp[i].fraction * ffp[i].fragment[k].ek;
      blm += ffp[i].fraction * bl[i].exmean;
      bhm += ffp[i].fraction * bh[i].exmean;
      blw += ffp[i].fraction * bl[i].exwidth * bl[i].exwidth;
      bhw += ffp[i].fraction * bh[i].exwidth * bh[i].exwidth;
    }
    blw = sqrt(blw);
    bhw = sqrt(bhw);
  }

  std::cout << "# " << std::setw(4) << n << std::setw(5) << k;
  std::cout << std::setw(4) << ffp[0].fragment[k].getZl() << std::setw(4) << ffp[0].fragment[k].getAl();
  std::cout << std::setw(4) << ffp[0].fragment[k].getZh() << std::setw(4) << ffp[0].fragment[k].getAh();
  std::cout << std::setw(12) << yieldL;
  std::cout << std::setw(12) << yieldH;
  std::cout << std::setw(12) << ek;
  std::cout << std::setw(12) << ex;
  std::cout << std::setw(12) << blm;
  std::cout << std::setw(12) << blw;
  std::cout << std::setw(12) << bhm;
  std::cout << std::setw(12) << bhw << std::endl;

  if(nc > 1 && detailed){
    for(int i=0 ; i<nc ; i++){
      std::cout << "#              ";
      std::cout << std::setw(12) << ffp[i].fraction;

      int kl = FFPMultiChanceFindIndexL(i,k);
      int kh = FFPMultiChanceFindIndexH(i,k);
      std::cout << std::setw(12) << ffp[i].fragment[kl].yield;
      std::cout << std::setw(12) << ffp[i].fragment[kh].yield;
      std::cout << std::setw(12) << ffp[i].fragment[k].ek;
      std::cout << std::setw(12) << ffp[i].fragment[k].ex;
      std::cout << std::setw(12) << bl[i].exmean;
      std::cout << std::setw(12) << bl[i].exwidth;
      std::cout << std::setw(12) << bh[i].exmean;
      std::cout << std::setw(12) << bh[i].exwidth << std::endl;
    }
  }
}
#endif

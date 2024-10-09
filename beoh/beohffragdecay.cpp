/******************************************************************************/
/*  beohffragdecay.cpp                                                        */
/*        Statistical decay of fission fragment pair                          */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "beoh.h"
#include "beohoutput.h"
#include "statmodel.h"
#include "beohstat.h"
#include "nucleus.h"
#include "global.h"
#include "terminate.h"

static void beohStoreFFragPopulation(const int, BData *);


/* data for statistical model calculation */
static Spectra       spc;
static Transmission  **tc;
static Transmission  **td;
static double        **tg;
static GDR           gdr[MAX_GDR];


/**********************************************************/
/*      Decay of Each Fission Fragment                    */
/**********************************************************/
void beohFFDecay(const int nc, System *sys, Pdata *pdt, BData *bdt, GammaProduction *gml, const bool printindiv)
{
  int jmax = MAX_J-1;

  for(int i=0 ; i<MAX_GDR ; i++) gdr[i].clear();

  /*** initialization */
  cohResetNucleus();
  cohResetNuclearStructure();
  crx.specclear();

  setupGeneralParameter(false,sys);

  /*** the highest excitation energy = <Ex> + f x dE, where f is given in beoh.h */
  sys->ex_total = bdt[0].exmean + EXRANGE_FACTOR * bdt[0].exwidth;

  setupStatModel(false,jmax,sys,pdt);
  beohParameterSetting(sys,gdr,true);

  /*** prepare initial population in the first compound nucleus*/
  beohStoreFFragPopulation(nc,bdt);

  /*** statistical decay calculation */
  if(ncl[0].ncont > 0) beohspectra(sys,pdt,tc,td,tg,&spc,gml);
  else{
    specGammaCascade(spc.cn[0],&ncl[0]);
    if(opt.finegammaspectrum) specStoreDiscreteGamma(&ncl[0],gml);
  }

  /*** print stat calc result for individual CN */
  if(printindiv){
    outSystem(bdt->getMode(),sys);
    outLevelDensity(sys->max_compound,0.0);
    outGDR(false,sys->max_compound);
    outCompound(sys->max_compound,pdt);
    outSpectrumSum(false,ncl[0].de,crx.spectra,&ncl[0]);
    outGSProduction(sys->max_compound);
  }
}


/**********************************************************/
/*      Initial Population for Fission Frag. Decay Calc.  */
/**********************************************************/
void beohStoreFFragPopulation(const int nc, BData *bdt)
{
  beohClearInitialPopulation(&ncl[0]);

  if(ncl[0].ncont > 0){
    /*** add each chance fission contribution */
    int jmax = 0;
    for(int i=0 ; i<nc ; i++){
      beohInitialPopulationStat(bdt[i].fraction,bdt[i].initspin,bdt[i].initparity,bdt[i].spinfactor,bdt[i].exmean,bdt[i].exwidth,&ncl[0]);
      if(ncl[0].jmax > jmax) jmax = ncl[0].jmax;
    }
    ncl[0].jmax = jmax;
  }
  /*** when no continuum, start calculation from the highest level */
  else{
    int nstart = beohStoreLevelExcite(1.0,&ncl[0]);
    ncl[0].ndisc = nstart;
  }

  double j0 = halfint(ncl[0].lev[0].spin);
  double xp = 0.0, xj = 0.0;
  for(int k=0 ; k<ncl[0].ncont ; k++){
    for(int j=0 ; j<MAX_J ; j++){
      xp +=  ncl[0].pop[k][j].even + ncl[0].pop[k][j].odd;
      xj += (ncl[0].pop[k][j].even + ncl[0].pop[k][j].odd) * (j + j0);
    }
  }
  bdt[0].averagespin = (xp > 0.0) ? xj / xp : 0.0; // although xp should be unity
}


/**********************************************************/
/*      Decay Calculation for Neutron Spectrum            */
/**********************************************************/
int beohFFDecayPrep(System *sys, Pdata *pdt)
{
  for(int i=0 ; i<MAX_GDR ; i++) gdr[i].clear();

  /*** determin Jmax at mean energy */
  setupGeneralParameter(false,sys);
  setupStatModel(false,MAX_J-1,sys,pdt);
  beohParameterSetting(sys,gdr,false);
  int jmax = beohPopinitJmax(sys->ex_total,&ncl[0]);

  return jmax;
}


void beohFFDecaySpec(System *sys, Pdata *pdt, BData *bdt, GammaProduction *gml, const int jmax, const int n, double *speclab)
{
  for(int i=0 ; i<MAX_GDR ; i++) gdr[i].clear();

  cohResetNucleus();
  cohResetNuclearStructure();
  crx.specclear();

  bdt->exwidth = 1.0;
  sys->ex_total = bdt->exmean;

  setupStatModel(false,jmax,sys,pdt);
  beohParameterSetting(sys,gdr,false);

  if(ncl[0].ncont > 0){
    beohStoreFFragPopulation(1,bdt);
    beohspectra(sys,pdt,tc,td,tg,&spc,gml);
  }
  else specGammaCascade(spc.cn[0],&ncl[0]);

  /*** convert CMS spectrum into LAB */
  double ef = bdt->ekmean / sys->compound.getA();
  beohLabSpectrum(n,ef,beohZeroCut(crx.spectra),ncl[0].de,crx.spectra[1],speclab);
}


/**********************************************************/
/*      Calculate <E> and Multiplicity from Spectrum      */
/**********************************************************/
void beohAverageEnergy(const int n, double *m, double *e, double *t, double **spc)
{
  double de = ncl[0].de;
  for(int p=0 ; p<2 ; p++){
    m[p] = t[p] = 0.0;
    for(int k=0 ; k<n ; k++){
      double e0 = (k > 0)  ? (k - 0.5) * de : 0;
      double e1 = (k + 0.5) * de;
      double ex = (e0 + e1) * 0.5;

      m[p] += spc[p][k] * de;
      t[p] += spc[p][k] * de * ex;
    }
    e[p] = (m[p] > 0.0) ? t[p]/m[p] : 0.0;
  }
}


/**********************************************************/
/*      Allocate Arrays                                   */
/**********************************************************/
void beohFFDecayAllocateMemory(void)
{
  try{
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
    cohTerminateCode("beohFFDecayAllocateMemory");
  }
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void beohFFDecayDeleteAllocated(void)
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


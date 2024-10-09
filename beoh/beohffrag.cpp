/******************************************************************************/
/*  beohffrag.cpp                                                             */
/*        Calculation of statistical decay of fission fragments               */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "beoh.h"
#include "beohoutput.h"
#include "statmodel.h"
#include "ffpmultichance.h"
#include "nucleus.h"
#include "masstable.h"
#include "global.h"
#include "terminate.h"
#include "dir.h"

static bool PRINT_FISSION_PAIR        = false;
static bool PRINT_FISSION_PAIR_SUPPL  = false;

static double finegrid = 0.01; // fine grid for gamma-ray output, 10 keV

static void beohFissionFragmentInit (System *, BData *, FFragData *);
static void beohEnergyDependentParameters (FFragData *);
static void beohFFPAllocateMemory (const int);
static void beohFFPDeleteAllocated (const int);
static void beohReadInputFiles(System *, FFragData *);


/**********************************************************/
/*      Global Scope Data Arrays                          */
/**********************************************************/
/* data for fission yield Y(Z,A,TXE) */
FissionFragmentPair *ffp;

/* data for post-processing */
FissionObservable fob;        // post-processed quantities for data comparison

/* data for discrete gamma output */
GammaProduction gmlL, gmlH;

/**********************************************************/
/*      Fission Fragment Decay Calculation                */
/**********************************************************/
void beohFissionFragmentMain(System *sys, BData *bdt, FFragData *fdt, unsigned int ompid)
{
  if(fdt->getFissionChance() < 1) return;

  if( (bdt->getMode() == fissionspec) && (fdt->getFissionChance() > 1) ){
    message << "fission specturm calculation for multi-chance case not supported";
    cohTerminateCode("beohFissionFragmentMain");
  }

  if(fdt->ycutoff == 0.0) fdt->ycutoff = YIELD_CUTOFF;
  if(ompid == 0) ompid = 6; // Koning-Delaroche potential as default


  /*** allocate Fragment Pair Object for each chance-fission */
  ffp = new FissionFragmentPair [fdt->getMaxFissionChance()];

  /*** determine fissioning system for the first chance fission and set binding energies */
  beohFissionFragmentInit(sys,bdt,fdt);
  beohReadInputFiles(sys,fdt);
  beohEnergyDependentParameters(fdt);

  /*** prepare fission fragment pair object and allocate memory */
  beohFFPAllocateMemory(fdt->getMaxFissionChance());
  beohFFDecayAllocateMemory();

  /*** create fission fragment pairs and store them in the FissionFragmentPair object */
  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){

    /*** 'ffydata'==external - YAmodel, TKEmodel, RTmodel exist */
    if (! fdt->yafile.empty()){
   
      /*** construct Gaussian parameters from the ones that were read in */
      double en = ffp[nc].ex - ffp[nc].bn;
      FFPConstructYieldParameters(fdt->mc[nc].yap,en,fdt->mc[nc].GaussSigma,fdt->mc[nc].GaussDelta,fdt->mc[nc].GaussFract,&ffp[nc]);
      if (! fdt->ytkefile.empty()) {
        FFPGeneratePairs(fdt->spontaneous,fdt->mc[nc].tke,&ffp[nc],fdt->mc[nc].ZpFactor,fdt->mc[nc].tkea, fdt->mc[nc].stkea,fdt->mc[nc].getNFTune(),fdt->mc[nc].ftuneMass,fdt->mc[nc].ftuneFactor);
      }
      else {
        FFPGeneratePairs(fdt->spontaneous,fdt->mc[nc].tke,&ffp[nc],fdt->mc[nc].ZpFactor,fdt->mc[nc].getNFTune(),fdt->mc[nc].ftuneMass,fdt->mc[nc].ftuneFactor);
      }
    }

    /*** 'ffydata'== input, Wahl or internal */
    else if(fdt->mc[nc].ffydata == "input" || fdt->mc[nc].ffydata == "Wahl" || fdt->mc[nc].ffydata == "internal"){

      /*** retrieve Gaussian parameters*/
      FFPMassYieldParameters(fdt->spontaneous,fdt->mc[nc].ffydata,fdt->mc[nc].GaussSigma,fdt->mc[nc].GaussDelta,fdt->mc[nc].GaussFract,&ffp[nc]);

      /*** generate fission fragment pairs and re-calculate TKE */
      FFPGeneratePairs(fdt->spontaneous,fdt->mc[nc].tke,&ffp[nc],fdt->mc[nc].ZpFactor,fdt->mc[nc].getNFTune(),fdt->mc[nc].ftuneMass,fdt->mc[nc].ftuneFactor);
    }

    /*** when data are given by experimental Y(A,TKE) data */
    else if(fdt->mc[nc].ffydata == "YATKE"){
      FFPReadExpdata(fdt->mc[nc].ffydata,&ffp[nc],fdt->mc[nc].ZpFactor);
    }

    /*** from externa file */
    else{
      /*** read data in external file */
      FFPReadPairs(fdt->mc[nc].ffydata,&ffp[nc]);
    }

    /*** recalculate TKE from fission product pairs, only for printing */
    fdt->mc[nc].tke = FFPAverageTKE(&ffp[nc]);
  }

  /*** ffp[0] will have a union list of the all FF appearing at all multi-chances */
  if(fdt->getFissionChance() > 1) FFPExpandList(fdt->getFissionChance(),ffp);

  /*** calculate decay of fission fragments */
  FFPSplit(bdt->getMode(),sys,fdt,ompid);

  /*** print parameters */
  if(PRINT_FISSION_PAIR){
    for(int nc=0 ; nc < fdt->getFissionChance() ; nc++) FFPOutputPairs(fdt->ycutoff,&ffp[nc]);
  }
  outSystem(bdt->getMode(),sys);
  outFissionFragment(sys->compound.getA(),fdt);

  /*** main output */
  double de = (sys->energy_bin_in == 0.0) ? ENERGY_BIN : sys->energy_bin_in;

  if(bdt->getMode() == fissionspec) FFPOutputSpec(de,fdt->maxtemperature,&fob);
  else{
    FFPOutput(de,fdt->massresolution,fdt->maxtemperature,&fob);

    if(opt.finegammaspectrum){
      outSpectrumFineGamma(fob.spectrum[0],&gml,de,finegrid);
      specSortDiscreteGamma(&gml);
      outDiscreteGamma(1.0,&gml);
    }

    /*** print production cross sections of ground and isomeric states */
    if(bdt->getMode() == cumulativeyield){
      beohFragmentBetaDecay(fdt->maxhalflife,fdt->branchdatafile);
    }
    else{
      outPrepPrintResidual();
      outPrepPrintIsomericRatio();
    }
  }

  /*** supplemental output */
  if(PRINT_FISSION_PAIR_SUPPL){
    for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
      FFPOutputYield(ffp[0].nmass,ffp[nc].ac,ffp[nc].mass_first,fdt->mc[nc].GaussSigma,fdt->mc[nc].GaussDelta,fdt->mc[nc].GaussFract);
      FFPOutputTKEDistribution(3,&ffp[nc]);
      FFPOutputTXEDistribution(&ffp[nc]);
    }
  }

  /*** free allocated memory */
  beohFFPDeleteAllocated(fdt->getMaxFissionChance());
  beohFFDecayDeleteAllocated();

  delete [] ffp;
}


/**********************************************************/
/*      Determine Fissioning System                       */
/**********************************************************/
void beohFissionFragmentInit(System *sys, BData *bdt, FFragData *fdt)
{
  /*** spontaneous fission case, no incident given */
  if(fdt->spontaneous){
    sys->compound = sys->target;
    ffp[0].en = ffp[0].ex = 0.0;
    ffp[0].setZA(sys->compound.getZ(),sys->compound.getA());
    fdt->resetFissionChance(1);
  }

  /*** particle induced fission cases */
  else{
    /*** determine incident particle */
    sys->incident = setupParticle(sys->incident.pid);
     
    /*** first-chance fissioning nucleus = (Z,A) + incident ZA */
    sys->compound.setZA(sys->target.getZ() + sys->incident.za.getZ(), sys->target.getA() + sys->incident.za.getA());

    /*** total excitation energy Binding + En(cms) */
    double tmx = mass_excess(sys->target.getZ()  ,sys->target.getA()  );
    double cmx = mass_excess(sys->compound.getZ(),sys->compound.getA());

    double tm  = sys->target.getA() + tmx / AMUNIT;
    double im  = sys->incident.mass;

    ffp[0].en = sys->input_energy = bdt->energy;

    sys->reduced_mass = tm * im / (tm+im);
    sys->lab_energy   = sys->input_energy;
    sys->cms_energy   = sys->lab_energy * tm/(tm+im);
    sys->ex_total     = sys->cms_energy + sys->incident.mass_excess + tmx - cmx;

    for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
      /*** ZA of fissioning compound nucleus */
      ffp[nc].setZA(sys->compound.getZ(),sys->compound.getA() - nc);

      /*** neutron binding energy */
      ffp[nc].bn = mass_excess(ffp[nc].zf,ffp[nc].af-1) + ENEUTRON - mass_excess(ffp[nc].zf,ffp[nc].af);
    }

    /*** multi-chance fission data from file */
    if(fdt->mcffile != ""){
      double dat[16];
      int ncdat = FFPReadMultichanceFissionData(fdt->mcffile,ffp[0].en,dat);
      for(int nc=0 ; nc < ncdat ; nc++){
        if(fdt->mc[nc].fraction == 1.0) fdt->mc[nc].fraction = dat[nc    ];
        if(fdt->mc[nc].tke      == 0.0) fdt->mc[nc].tke      = dat[nc + 4];
        if(fdt->mc[nc].exfis    == 0.0) fdt->mc[nc].exfis    = dat[nc + 8];
        if(fdt->mc[nc].eprefis  == 0.0) fdt->mc[nc].eprefis  = dat[nc +12];
      }
      /*** when data file contains more fission-chance, copy the first chance data */
      if(fdt->getFissionChance() < ncdat){
        for(int nc=fdt->getFissionChance() ; nc < ncdat ; nc++){
          ffp[nc].setZA(sys->compound.getZ(),sys->compound.getA() - nc);
          ffp[nc].bn = mass_excess(ffp[nc].zf,ffp[nc].af-1) + ENEUTRON - mass_excess(ffp[nc].zf,ffp[nc].af);
          fdt->mc[nc].spinfactor = fdt->mc[0].spinfactor;
          fdt->mc[nc].rt         = fdt->mc[0].rt;
          fdt->mc[nc].ffydata    = fdt->mc[0].ffydata;
        }
      }

      /*** summation check, if 5th chance or more is possible */
      double sumf = 0.0;
      for(int nc=0 ; nc < ncdat ; nc++) sumf += fdt->mc[nc].fraction;
      if(sumf < 1.0) fdt->mc[ncdat-1].fraction += 1.0 - sumf;

      fdt->resetFissionChance(ncdat);
    }

    /*** enforce single-chance calculation when Pf=1.0 */
    if(fdt->mc[0].fraction == 1.0) fdt->resetFissionChance(1);

    /*** fission probabilities */
    for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
      ffp[nc].fraction = fdt->mc[nc].fraction;
      if(ffp[nc].fraction == 0.0){ fdt->resetFissionChance(nc); break; }
    }

    /*** check closed channel and average exication energy causing fission */
    double bnsum = 0.0;
    for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){

      /*** when not given, we assume fission happens at the highest excitation energy */
      if(fdt->mc[nc].exfis == 0.0){ fdt->mc[nc].exfis = sys->ex_total - bnsum; }
      bnsum += ffp[nc].bn;

      ffp[nc].ex = fdt->mc[nc].exfis;
      /*** check closed channel */
      if(ffp[nc].ex < 0.0){ fdt->resetFissionChance(nc); break; }
    }
  }

  /*** check TKE */
  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    if(fdt->mc[nc].tke == 0.0){
      fdt->mc[nc].tke = FFPSystematics_TKE_En(fdt->spontaneous,ffp[nc].zf,ffp[nc].af,ffp[nc].ex,ffp[nc].bn);
    }
  }

  /*** determine mass-range to be calculated */
  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    ffp[nc].setZARange(MASS_FIRST,CHARGE_RANGE);
  }
}


/**********************************************************/
/*      Calculate Energy-Dependent Parameters             */
/**********************************************************/
void beohEnergyDependentParameters(FFragData *fdt)
{
  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    double en = ffp[nc].ex - ffp[nc].bn;
    if(en < 0.0) en = 0.0;

    if (! fdt->ytkefile.empty()) {
      double a = fdt->mc[nc].tkep[0];
      double e0 = fdt->mc[nc].tkep[1];
      double b = fdt->mc[nc].tkep[2];
      double d = fdt->mc[nc].tkep[3];
      double c = a + (b-d)*e0;
      if (en<=e0) fdt->mc[nc].tke = a+b*en;
      if (en>e0) fdt->mc[nc].tke = c+d*en;
    }
    else{

    if(fdt->mc[nc].tke1 != 0.0){
      fdt->mc[nc].tke = fdt->mc[nc].tke + en * fdt->mc[nc].tke1;
      if(fdt->mc[nc].tke < 0.0) fdt->mc[nc].tke = 0.0;
    }
    }

    /*** when RT < 1, then RT is always 1.0 */
    if(fdt->mc[nc].rt1 != 0.0){
      fdt->mc[nc].rt  = fdt->mc[nc].rt + en * fdt->mc[nc].rt1;
      if(fdt->mc[nc].rt < 1.0) fdt->mc[nc].rt = 1.0;
    }

    /*** spin factor energy dependence */
    if(fdt->mc[nc].spinfactor1 != 0.0){
      fdt->mc[nc].spinfactor = fdt->mc[nc].spinfactor + en * fdt->mc[nc].spinfactor1;
      /*** set lower limit of fJ to be 0.1 */
      if(fdt->mc[nc].spinfactor < 0.1) fdt->mc[nc].spinfactor = 0.1;
    }

    /*** ZpFactor energy dependence */
    for(int i=0 ; i<2 ; i++){
      if(fdt->mc[nc].ZpFactor[i+2] != 0.0){
        fdt->mc[nc].ZpFactor[i] = fdt->mc[nc].ZpFactor[i] + en * fdt->mc[nc].ZpFactor[i+2];
        /*** set Zp parameters lower limit 1e-10, since zero does not change calculations */

        if(fdt->mc[nc].ZpFactor[i] <= 0.0) fdt->mc[nc].ZpFactor[i] = 1e-10;
      }
    }


    // Y(A) energy dependence
    // if we read from a file, construct the Gaussians this way
    if (! fdt->yafile.empty()) {
       // parametrization from YAmodel.dat
       
    }
    // Gaussian values from 
    else {
    for(int i=0 ; i<4 ; i++){
      if(fdt->mc[nc].GaussFract[i+4] != 0.0){
        fdt->mc[nc].GaussFract[i] = fdt->mc[nc].GaussFract[i] + en * fdt->mc[nc].GaussFract[i+4];
        if(fdt->mc[nc].GaussFract[i] < 0.0){
          message << i+1 << "-th primary gaussian fraction became negative";
//        cohTerminateCode("beohEnergyDependentParameters");
          fdt->mc[nc].GaussFract[i] = 0.0;
        }
      }
      if(fdt->mc[nc].GaussSigma[i+4] != 0.0){
        fdt->mc[nc].GaussSigma[i] = fdt->mc[nc].GaussSigma[i] + en * fdt->mc[nc].GaussSigma[i+4];
        if(fdt->mc[nc].GaussSigma[i] < 0.0){
          message << i+1 << "-th primary gaussian width became negative";
          cohTerminateCode("beohEnergyDependentParameters");
        }
      }
      if(i == 3) continue; // symmetric part does not have mass shift
      if(fdt->mc[nc].GaussDelta[i+4] != 0.0){
        fdt->mc[nc].GaussDelta[i] = fdt->mc[nc].GaussDelta[i] + en * fdt->mc[nc].GaussDelta[i+4];
        if(fdt->mc[nc].GaussDelta[i] < 0.0){
          message << i+1 << "-th primary gaussian mass shift became negative";
          cohTerminateCode("beohEnergyDependentParameters");
        }
      }
      }
    }
  }
}


/**********************************************************/
/*      Read the input files		                  */
/**********************************************************/
void beohReadInputFiles(System *sys, FFragData *fdt)
{
  static std::string datadir = DATADIR;
  int ZAID = 1000*sys->compound.getZ() + sys->compound.getA();
  if (fdt->spontaneous) ZAID = -1*ZAID;
  bool found = false;

  // read the Y(A|E) parameter input file
  // for each nth-chance fission, if it exists
  // otherwise, use the previous chance parameters - the means should be shifted automatically
  std::string ya = datadir + "/" + PARAMETERSDIRECTORY + fdt->yafile;
  for (int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    int ZAIDt = ZAID - nc;
    if (! fdt->yafile.empty()){
      found = false;
      std::ifstream fp(ya.c_str(),std::ios::in);
      std::string line;
      int ZAIDtemp;

      if (!fp.is_open()) std::cout << "Parameters are not what you think for YA!!" << std::endl;

      while (getline(fp,line)) {
        if (line.find("#")>0){
          std::istringstream ss(line);
          ss.str(line);
          ss.clear();
          ss >> ZAIDtemp;
          if (ZAIDtemp==ZAIDt){
            found = true;
            for (int i=0; i<18; i++){
              ss >> fdt->mc[nc].yap[i];
            }
          }
        }
      }
      fp.close();
      if (! found && nc >0) {
        for (int i=0; i<18; i++) {
          fdt->mc[nc].yap[i] = fdt->mc[nc-1].yap[i];
        }
      }
    }
  }

  // read the Y(TKE) and TKE(E) parameter input file
  // for each nth-chance fission, if it exists
  // otherwise, use the previous nth-chance parameters
  std::string tke = datadir + "/" + PARAMETERSDIRECTORY + fdt->ytkefile;
  for (int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    int ZAIDt = ZAID - nc;
    if (! fdt->ytkefile.empty()) {
      found = false;
      std::ifstream fp(tke.c_str(),std::ios::in);
      std::string line;
      int ZAIDtemp;

      if (!fp.is_open()) std::cout << "Parameters are not what you think for TKE!!" << std::endl;

      while (getline(fp,line)){
        if (line.find("#")>0){
          std::istringstream ss(line);
          ss.str(line);
          ss.clear();
          ss >> ZAIDtemp;
          if (ZAIDtemp==ZAIDt){
            found = true;
            for (int i=0; i<4; i++){
              ss >> fdt->mc[nc].tkep[i];
            }
            for (int i=0; i<4; i++){
              ss >> fdt->mc[nc].tkea[i];
            }
            for (int i=0; i<3; i++){
              ss >> fdt->mc[nc].stkea[i];
            }
          }
        }
      }
      fp.close();
      if (! found && nc > 0){
	 for (int i=0; i<4; i++) fdt->mc[nc].tkep[i] = fdt->mc[nc-1].tkep[i];
	 for (int i=0; i<4; i++) fdt->mc[nc].tkea[i] = fdt->mc[nc-1].tkea[i];
	 for (int i=0; i<3; i++) fdt->mc[nc].stkea[i] = fdt->mc[nc-1].stkea[i];
      }
    }
  }

  // read the RT(A) input file
  // for each nth-chance fission, if it exists
  // otherwise, fill in the value from the given rt value
  std::string rta = datadir + "/" + PARAMETERSDIRECTORY + fdt->rtafile;
  for(int nc=0 ; nc < fdt->getFissionChance() ; nc++){
    // read RTA file, if it is given
    int ZAIDt = ZAID - nc;
    if (! fdt->rtafile.empty()) {
       std::ifstream fp(rta.c_str(),std::ios::in);
       std::string line;
       int ZAIDtemp;
       int nmax;
       found = false;

       if (!fp.is_open()) std::cout << "Parameters are not what you think for RTA!!" << std::endl;

       while (getline(fp,line)){
         if (line.find("#")>0){
           std::istringstream ss(line);
           ss.str(line);
           ss.clear();
           ss >> ZAIDtemp;
           if (ZAIDtemp==ZAIDt){
	     found = true;
	     ss >> nmax;
             for (int i=0; i<nmax+2; i++){
		ss >> fdt->mc[nc].rta[i];
             }
           }
         }
       }
       fp.close();
       if (! found){
	  fdt->mc[nc].rta[0] = fdt->mc[nc-1].rta[0];
	  fdt->mc[nc].rta[1] = fdt->mc[nc-1].rta[1];
	  for (int i=2; i<80; i++) fdt->mc[nc].rta[i] = fdt->mc[nc].rt;
       }
    } else { 
       // if we don't read the file, add the rt constant value, so we can just use this array
       fdt->mc[nc].rta[0] = 0;
       for (int i=1; i<80; i++) fdt->mc[nc].rta[i] = fdt->mc[nc].rt;
    }
  }
   
}


/**********************************************************/
/*      Allocate Arrays                                   */
/**********************************************************/
void beohFFPAllocateMemory(const int nc)
{
  try{
    for(int i=0 ; i<nc ; i++) ffp[i].memalloc(MAX_FISSION_PAIR);
    fob.memalloc(MAX_MULTIPLICITY_DIST,ffp[0].mass_first+abs(ffp[0].nmass+1),MAX_ENERGY_BIN);

    /*** allocate long-lived states */
    outPrepAllocateResidual(MAX_FISSION_PAIR * 2);

    /*** allocate discrete gamma lines from light and heavy fragments */
    gmlL.memalloc(MAX_GAMMALINE,CUT_GAMMALINE);
    gmlH.memalloc(MAX_GAMMALINE,CUT_GAMMALINE);
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("beohFFPAllocateMemory");
  }
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void beohFFPDeleteAllocated(const int nc)
{
  for(int i=0 ; i<nc ; i++) ffp[i].memfree();
  fob.memfree();
  outPrepFreeResidual();
  gmlL.memfree();
  gmlH.memfree();
}


#ifndef __STRUCTUR_H__
#define __STRUCTUR_H__
#include "structur.h"
#endif

#ifndef __PHYSICALCONSTANT_H__
#define __PHYSICALCONSTANT_H__
#include "physicalconstant.h"
#endif

#ifndef __BETACHAIN_H__
#define __BETACHAIN_H__
#include "betachain.h"
#endif

#ifndef __FFP_H__
#define __FFFP_H__
#include "ffp.h"
#endif

#ifndef __BEOHFFRAG_H__
#define __BEOHFFRAG_H__
#include "beohffrag.h"
#endif

/****************************/
/*   Defaults               */
/****************************/
/*** constants for fission-decay calculation */
static const double EXRANGE_FACTOR =    2.5; // total excitation will be <E> + FACTOR x width
//static const double EXRANGE_FACTOR =    2.0; // total excitation will be <E> + FACTOR x width
static const double YIELD_CUTOFF   = 1.0e-6; // lower limit for fragment yield

typedef enum {undef = 0, betadecay = 1, statdecay = 2, fissiondecay = 3, fissionspec = 4, cumulativeyield = 5} CalcMode;


/****************************/
/*   BeoH Specific Data     */
/****************************/
class BData{
 private:
  CalcMode  mode;
 public:
  int    initparity;          // parity of the initial compound
  double initspin;            // spin of the initial compound
  double spinfactor;          // multiplication factor for initial spin distribution
  double ekmean;              // average kinetic energy of the moving compound
  double exmean;              // average excitation energy of the initial compound
  double ekwidth;             // width of kinetic energy
  double exwidth;             // width of the excitation energy
  double energy;              // incident neutron energy
  double fraction;            // fractional contribution
  std::string yafile;	      // user defined Y(A) model parameters
  std::string ytkefile;	      // user defined Y(TKE|A) and TKE(E) model parameters
  std::string rtafile;        // user defined R_T(A) distribution
  double averagespin;         // average J of initial compound

  BData(){
    init();
  }

  void init(){
    mode           = undef;
    initparity     = 0;
    initspin       = 0.0;
    spinfactor     = 0.0;
    ekmean         = 0.0;
    exmean         = 0.0;
    ekwidth        = 0.0;
    exwidth        = 0.0;
    energy         = 0.0;
    fraction       = 1.0;
    yafile	   = "";
    ytkefile       = "";
    rtafile	   = "";
    averagespin    = 0.0;
  }

  void setMode(CalcMode m){
    if(m == betadecay) mode = betadecay;
    else if(m == statdecay) mode = statdecay;
    else if(m == fissiondecay) mode = fissiondecay;
    else if(m == fissionspec) mode = fissionspec;
    else if(m == cumulativeyield) mode = cumulativeyield;
    else mode = undef;
  }

  CalcMode getMode(){ return(mode); }

  bool isBetaCalc(){
    if(mode == betadecay) return(true);
    else return(false);
  }
};


/****************************/
/*   Beta Strength Dist.    */
/****************************/
class BetaProfile{
 public:
    int     ntotal;           // nubmer of total continuum bins
    int     ncont;            // nubmer of continuum bins
    int     ndisc;            // nubmer of discrete states
    double *excitation;       // continuum bin excitation energies
    double *rcont;            // continuum distribution
    double *rdisc;            // distribution in the levels
    double  tcont;            // sum of continuum distribution
    double  tdisc;            // sum of discrete distribution
    std::string  source;      // data source (GT, ENSDF, or MIX)
    Level  *lev;              // discrete level data

    BetaProfile(){
      ntotal     = 0;
      ncont      = 0;
      ndisc      = 0;
      rcont      = new double [MAX_ENERGY_BIN];
      rdisc      = new double [MAX_LEVELS];
      source     = "";
      for(int k=0 ; k<MAX_LEVELS     ; k++) rdisc[k] = 0.0;
      for(int k=0 ; k<MAX_ENERGY_BIN ; k++) rcont[k] = 0.0;
    }

    ~BetaProfile(){
      delete [] rcont;
      delete [] rdisc;
    }
};


enum section{BEGIN,DATA,BETA,HFS,FFRAG,MCF,CFY,END};
typedef enum section Section;

/**************************************/
/*      beoh.cpp                      */
/**************************************/
void    cohAllocateNucleus (void);
void    cohAllocateNuclearStructure (void);
void    cohResetAllocationPointer (void);

void    cohResetNucleus (void);
void    cohResetNuclearStructure (void);

/**************************************/
/*      beohmain.cpp                  */
/**************************************/
void    beohMainLoop (const unsigned long);
void    beohCGMCompatibleMode (const int, const int, const double, const int, const double, const double, const double, const double, const double, const unsigned long);
void    beohStoreStatPopulation (const int, const int, BData *);
void    beohParameterSetting (System *, GDR *, const bool);
int     beohZeroCut (double *);
int     beohZeroCut (double **);
void    beohAllocateMemory (void);
void    beohDeleteAllocated (void);

/**************************************/
/*      beohffragenergy.cpp           */
/**************************************/
double  beohExcitationEnergyDivide (ZAnumber *, ZAnumber *, double, const double);
double  beohRTDamping (const double, const double);

/**************************************/
/*      beohffragbeta.cpp             */
/**************************************/
void    beohFragmentBetaDecay (double, std::string);

/**************************************/
/*      beohsetup.cpp                 */
/**************************************/
void    setupClearParameter (const int, System *, Pdata *, GDR *);
void    setupGeneralParameter (const bool, System *);
int     setupStatModel (const bool, const int, System *, Pdata *);
Pdata   setupParticle (const Particle);

/**************************************/
/*      beohdataread.cpp              */
/**************************************/
int     readHead (char *);
int     readSystem (char *, System *, Pdata *, FFragData *, BData *);
int     readBeta (char *, Beta *);
int     readStatModel (char *, GDR *, Fission *);
int     readFissionFragment (char *, FFragData *, BData *);
int     readMultiChanceFission (char *, MultiChanceFissionData *);
int     readFissionYield (char *, FFragData *);

/**************************************/
/*      beohlabspectrum.cpp           */
/**************************************/
void    beohIntegCmsSpectrum (const int, const int, const unsigned int, BData *, double , double *);
void    beohLabSpectrum (const int, const double, const int, const double, double *, double *);
double  beohLabConvert (const int, const double, const double, const double, const double, double *);

/**************************************/
/*      beohensdf.cpp                 */
/**************************************/
int     beohENSDFRead (ZAnumber *, Beta *);

/**************************************/
/*      beohffrag.cpp                 */
/**************************************/
void    beohFissionFragmentMain (System *, BData *, FFragData *, unsigned int);

/**************************************/
/*      beohffragdecay.cpp            */
/**************************************/
void    beohFFDecay (const int, System *, Pdata *, BData *, GammaProduction *, const bool);
int     beohFFDecayPrep (System *, Pdata *);
void    beohFFDecaySpec (System *, Pdata *, BData *, GammaProduction *, const int, const int, double *);
void    beohAverageEnergy (const int, double *, double *, double *, double **);
void    beohFFDecayAllocateMemory (void);
void    beohFFDecayDeleteAllocated (void);

/**************************************/
/*      beohffragcalc.cpp             */
/**************************************/
void    FFPSplit (const CalcMode, System *, FFragData *, const int);

/**************************************/
/*      beohdataread.cpp              */
/**************************************/
int     readFissionFragment (char *, FFragData *, BData *);
int     readMultiChanceFission (char *, MultiChanceFissionData *);


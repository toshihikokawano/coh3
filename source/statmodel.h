// definition of transmission data array
// functions for statistical model calculations

/*** spin int / half-int distinction */
static inline double halfint(double x)
{
  return( (((int)(x*2.0+0.001) % 2)==1) ? 0.5 : 0.0 );
}

/*** transmission coefficient index calculator */
static inline int tj_index(int a, int b, int c){
  int k = a*3/2;
  if(c != 0) k += (c-b)/2;
  return(k);
}

/*** parity calculator */
#ifndef __PARITY_FUNC__
#define __PARITY_FUNC__
static inline int parity(int x){ return((((x)/2)%2==0) ? 1 : -1); }
#endif


typedef enum {sumall=0, wfactor=1, hauser=2, fluctuation=3} Statcalcmode;


/**************************************/
/*      Transmission Coefficients     */
/**************************************/
class Transmission{
 private:
  int       jsize;            // size of Tj array
  bool      allocated;
 public:
  int       lmax;            // max L
  double    ecms;            // cms energy
  double    sigr;            // reaction cross section
  double   *tran;            // Tj [l+s][l][l-s]

  Transmission(){
    lmax = 0;
    ecms = 0.0;
    sigr = 0.0;
    allocated = false;
  }

  ~Transmission(){
    memfree();
  }

  void memalloc(int j){
    if(!allocated){
      jsize = j;
      tran = new double [3*jsize];
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] tran;
      allocated = false;
    }
  }

  void memclear(){
    lmax = 0;
    ecms = 0.0;
    sigr = 0.0;
    if(allocated){
      for(int j=0 ; j<3*jsize ; j++) tran[j] = 0.0;
    }
  }

  int getJsize(){return jsize;}
};


/**************************************/
/*      Spectrum Array                */
/**************************************/
class Spectra{
 private:
  int       csize;            // number of channels
  int       nsize;            // number of energy bins
  int       lsize;            // number of discrete levels
  bool      allocated;
 public:
  double   **pe;              // preequilibrium spectrum
  double   **cn;              // spectrum from compound decay
  double   **dp;              // population increment
  double   *pg[2];            // primary gamma-ray
  int      npg;               // number of PG lines

  Spectra(){
    allocated = false;
    npg   = 0;
    csize = 0;
    nsize = 0;
  }

  ~Spectra(){
    memfree();
  }

  void memalloc(int c, int n, int l){
    if(!allocated){
      csize = c;
      nsize = n;
      lsize = l;
      pe    = new double * [csize];
      cn    = new double * [csize];
      dp    = new double * [csize];
      pg[0] = new double [lsize];
      pg[1] = new double [lsize];
      for(int j=0 ; j<csize ; j++){
        pe[j] = new double [nsize];
        cn[j] = new double [nsize];
        dp[j] = new double [nsize];
      }

      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      for(int j=0 ; j<csize ; j++){
        delete [] pe[j];
        delete [] cn[j];
        delete [] dp[j];
      }
      delete [] pe;
      delete [] cn;
      delete [] dp;
      delete [] pg[0];
      delete [] pg[1];
      allocated = false;
    }
  }

  void memclear(){
    if(allocated){
      for(int j=0 ; j<csize ; j++){
        for(int n=0 ; n<nsize ; n++) pe[j][n] = cn[j][n] = dp[j][n] = 0.0;
      }
      for(int l=0 ; l<lsize ; l++) pg[0][l] = pg[1][l] = 0.0;
    }
  }

  void memclear(std::string name){
    if(allocated){
      if(name == "cn"){
        for(int j=0 ; j<csize ; j++){
          for(int n=0 ; n<nsize ; n++) cn[j][n] = 0.0;
        }
      }
      else if(name == "pe"){
        for(int j=0 ; j<csize ; j++){
          for(int n=0 ; n<nsize ; n++) pe[j][n] = 0.0;
        }
      }
      else if(name == "dp"){
        for(int j=0 ; j<csize ; j++){
          for(int n=0 ; n<nsize ; n++) dp[j][n] = 0.0;
        }
      }
      else if(name == "pg"){
        for(int l=0 ; l<lsize ; l++) pg[0][l] = pg[1][l] = 0.0;
      }
    }
  }

  int getCsize(){return csize;}
  int getNsize(){return nsize;}
  int getLsize(){return lsize;}
};


/**************************************/
/*      statparm.cpp                  */
/**************************************/
double  statBinWidth (const double);
int     statSetupEnergyBin (Nucleus *);
void    statSetupLevelDensity (Nucleus *, LevelDensity *);
void    statSetupGdrParameter (Nucleus *, GDR *, const double);
void    statSetupGdrSystematics (Nucleus *, const double);
void    statSetupResetGdrParameter (const int, Nucleus *);
void    statSetupFissionParameter (Fission *);


/**************************************/
/*      statdirect.cpp                */
/**************************************/
void    statSetupDirectInelastic (Level *, Nucleus *, CrossSection *);
void    statDirectSpectra (const double, const double, Level *, double *, double *);


/**************************************/
/*      stataveragespace.cpp          */
/**************************************/
double  statCalculateD0 (const double, Nucleus *, Nucleus *);
void    statAdjustD0 (const double, Nucleus *, Nucleus *);
void    statDecayWidth (System *, Pdata *, Transmission **, Transmission **, double **, Spectra *);


/**************************************/
/*      statmaniplevel.cpp            */
/**************************************/
bool    statAddDiscreteLevels (NuclearStructure *);
bool    statForceIncludeMeta (NuclearStructure *, Nucleus *);
int     statFixDiscreteLevels (NuclearStructure *);
void    statPopTrapMeta (NuclearStructure *);
void    statCutDiscreteLevels (const double, Nucleus *);
void    statAdjustDensityByExtraLevels (Nucleus *);
void    statClearDensityByExtraLevels (Nucleus *);
void    statAdjustContinuumBoundary (Nucleus *);


/**************************************/
/*      statskiplevel.cpp             */
/**************************************/
int     statSkipDiscreteLevel (NuclearStructure *);


/**************************************/
/*      fisdens.cpp                   */
/**************************************/
void    statSetupFissionLevelDensity (Nucleus *, Fission *);


/**************************************/
/*      stattrans.cpp                 */
/**************************************/
void    statStoreContinuumTransmission (const int, const double, Pdata *, Transmission **);
void    statStoreDiscreteTransmission (const int, const double, Pdata *, Transmission **);
void    statStoreContinuumGammaTransmission (const int, double **, Nucleus *);
void    statStoreDiscreteGammaTransmission (const double, double *, Nucleus *);


/**************************************/
/*      popinit.cpp                   */
/**************************************/
double  statStoreInitPopulation (const int, const int, const int, const double, Transmission *);
double  statStoreInitPopulationPhoton (const int, const int, const double, Transmission *);
void    statAdjustInitPopulation (const double);


/**************************************/
/*      gtrans.cpp                    */
/**************************************/
void    gdrM1norm (const double, const double, GDR *);
void    gdrParameterSave (GDR *, const unsigned int);
void    gdrParameterReset (void);
double  gdrGammaTransmission (const GammaMultipolarity, const double, const double, const double);
double  gdrRenormE1StrengthFunction (const int, const int, const int, const double, const double, Nucleus *, double **);


/**************************************/
/*      gtransread.cpp                */
/**************************************/
void    gdrAbsorptionDataRead (void);
double  gdrGammaTransmissionInterpolate(const GammaMultipolarity, const double);
int     gdrAbsorptionDataYsize (void);


/**************************************/
/*      spectra.cpp                   */
/**************************************/
void    spectra (System *, Pdata *, Transmission *, Transmission **, Transmission **, double **, Direct *, Spectra *);
void    spectraExclusive (System *, Pdata *, Transmission *, Transmission **, Transmission **, double **, Direct *, Spectra *);


/**************************************/
/*      statmodeltools.cpp            */
/**************************************/
int     specFindEnergyBin (double, const double);
void    specCumulativeSpectra (const int, const int, double  **, Nucleus *);
void    specGaussianBroadening (const int, const double, const double, double *);
void    specGaussianBroadening (const int, const double, const double, double *, const double, const double);
int     specStorePrimaryGamma (Nucleus *, double **, GammaProduction *);
int     specStoreDiscreteGamma (Nucleus *, GammaProduction *);
void    specSortDiscreteGamma (GammaProduction *);


/**************************************/
/*      statdecay.cpp                 */
/**************************************/
void    specCompoundElastic (System *, double, Transmission *, Transmission **, Transmission **, double **, Direct *, Spectra *);


/**************************************/
/*      statcompounddecay.cpp         */
/**************************************/
void    specCompoundDecay (const int, const int, Transmission **, Transmission **, double **, Spectra *);
double  specTransmissionSum (const Statcalcmode, bool *, const int, const int, const int, double **, Transmission **, Transmission **, const  double, Nucleus *, Spectra *);
void    specLostPopFraction (const int, const int, double, bool *, double **);


/**************************************/
/*      statlevdecay.cpp              */
/**************************************/
void    specLevelDecay (const int, Pdata *, Nucleus *, Transmission **, Transmission **, Spectra *);


/**************************************/
/*      statreaddensity.cpp           */
/**************************************/
void    statReadLevelDensity (Nucleus *);
void    statReadFissionLevelDensity (Nucleus *);


/**************************************/
/*      gampop.cpp                    */
/**************************************/
double  specTransitionGamma (const Statcalcmode, const int, const int, const int, double **, const double, Nucleus *, double *, double *);
double  specLevelTransitionGamma (const Statcalcmode, const int, const int, const int, const double, Nucleus *);


/**************************************/
/*      fispop.cpp                    */
/**************************************/
double  specTransitionFission (const Statcalcmode, const int, const int, const int, const double, Nucleus *);
double  specRotationalBand (const int, KBand *, double);


/**************************************/
/*      fpotential.cpp                */
/**************************************/
double  specFissionPotentialModel (const int, const int, const int, const double, double *, Nucleus *);


/**************************************/
/*      parpop.cpp                    */
/**************************************/
double  specTransitionParticle (const Statcalcmode, const int, const int, const int, const int, Transmission *,Transmission *, const double, Nucleus *, Nucleus *, double *, double *);
double  specLevelTransitionParticle (const Statcalcmode, const int, const int, const int, const int, Transmission *, Transmission *, const double, Nucleus *, Nucleus *, double *, double *);


/**************************************/
/*      parleg.cpp                    */
/**************************************/
void    specLegendreCoefficient (const int, const int, const int, const double, const int, const int, const int, const int, Transmission **);
void    specLegendreCoefficientContinuum (const double, const int, const int, const int, const int, Transmission **, double ***);
int     statInelAngularDistribution (const int, const int, Transmission *);
int     statAngularDistribution (const int, const int, Transmission *);


/**************************************/
/*      gamcas.cpp                    */
/**************************************/
int     specGammaCascade (double *, Nucleus *);
int     specGammaCascadeExclusive (const int, double *, Nucleus *);


/**************************************/
/*      moldauer.cpp                  */
/**************************************/
void    statWidthFluctuationSet (const int, const int, const int, const int, const double);
void    statWidthFluctuationSetT0 (const double);
void    statWidthFluctuationReset (const double);
void    statWidthFluctuation (const double, const double, const double);
void    wfcheck (void);
double  statMoldauer (const int, const int, const int, const int, const double, const double);
double  statElasticEnhancement (const double);
double  statMoldauerCC (const double, const double);
double  statMoldauerCorrelation (const double, const double);


/**************************************/
/*      exciton.cpp                   */
/**************************************/
void    preqExcitonModel (System *, Pdata *, Transmission **, Transmission **, Spectra *);
double  preqExcitonSpectra (Spectra *);
void    preqExcitonPopulation (const int, Transmission *, Transmission *, double *);
double  preqTotalPopulation (void);
void    preqRenormalizePopulation (const double);


/**************************************/
/*      excitransf.cpp                */
/**************************************/
void    preqExcitonTransfer (System *, Pdata *, Transmission **, Transmission **, Spectra *);

/**************************************/
/*      knockout.cpp                 */
/**************************************/
void    preqAlphaKnockout (System *, Pdata *, Transmission **, Transmission **, Spectra *);

/**************************************/
/*      dsd.cpp                       */
/**************************************/
int     dsdDirectCaptureModel (const int, const double, Pdata *, ZAnumber *, const double, Dcapt *, const double, double *);

/**************************************/
/*      photoabs.cpp                  */
/**************************************/
double  photoAbsorption (const int, const int, const double, const double, Nucleus *, Transmission *);
double  photoQDratio ();



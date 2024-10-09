#include <cmath>

#ifndef __PHYSICALCONSTANT_H__
#define __PHYSICALCONSTANT_H__
#include "physicalconstant.h"
#endif

#ifndef __FOBS_H__
#define __FOBS_H__
#include "fobs.h"
#endif


static const int  MAX_FISSION_PAIR      = 1000;  // FissionPair array size
static const int  MASS_FIRST            =   50;  // lowest A
static const int  CHARGE_RANGE          =    5;  // charge distribution range
static const int  MAX_MULTIPLICITY_DIST =   10;  // array size for P(n)


static inline double gaussian(double f, double x, double s)
{ return f / (s * sqrt(PI2)) * exp(-x*x/(2.0*s*s)); }


/****************************/
/*   Fission Pair           */
/****************************/
class ZAPair{
 private:
  unsigned int Zl;            // atomic number of light fragment
  unsigned int Al;            // mass number of light fragment
  unsigned int Zh;            // atomic number of heavy fragment
  unsigned int Ah;            // mass number of heavy fragment

 public:
  double yield;               // fraction (sum = 1.0)
  double qval;                // Q-value for this partition
  double ex;                  // average TXE
  double ek;                  // average TKE
  double sigma;               // TKE (TXE) distribution width
  double el;                  // excitation energy of light fragment
  double wl;                  // energy of width light fragment
  double eh;                  // excitation energy of heavy fragment
  double wh;                  // energy width of heavy fragment

  ZAPair(void){
    Zl    = 0;
    Al    = 0;
    Zh    = 0;
    Ah    = 0;
    yield = 0.0;
    qval  = 0.0;
    ex    = 0.0;
    ek    = 0.0;
    sigma = 0.0;
    el = wl = 0.0;
    eh = wh = 0.0;
  }

  void setPair(int zl, int al, int zh, int ah){
    Zl  = zl;
    Al  = al;
    Zh  = zh;
    Ah  = ah;
  }

  unsigned int getZl(){ return Zl; }
  unsigned int getAl(){ return Al; }
  unsigned int getZh(){ return Zh; }
  unsigned int getAh(){ return Ah; }
};


/****************************/
/*   Gaussian Parameters    */
/****************************/
class GaussParameter{
 public:
  double width;
  double center;
  double fraction;

  GaussParameter(){
    clear();
  }

  void clear(){
    width    = 0.0;
    center   = 0.0;
    fraction = 0.0;
  }

  void setVal(double x1, double x2, double x3){
    width    = x1;
    center   = x2;
    fraction = x3;
  }
};


/****************************/
/*   All Fragment Data      */
/****************************/
class FissionFragmentPair{
 private:
  bool allocated;
  int  max_pair;              // size of allocated fission product object
  int  ntotal;                // actual total number of nuclide
  int  gtotal;                // number of given gaussians
  static const int ng = 7;    // total number of Gaussian

 public:
  ZAPair         *fragment;
  GaussParameter  gauss[ng];  // Gaussian fitted parameters to experimental yield
  unsigned int    zf;         // compound Z = Zl + Zh
  unsigned int    af;         // compound A = Al + Ah
  double          ac;         // half of A, mid-point of FFrag distribution
  double          bn;         // neutron separation energy
  double          en;         // neutron incident energy
  double          ex;         // excitation energy
  double          fraction;   // fraction of this pair-set
  int             mass_first; // the lowest mass number in the data
  int             mass_last;  // the hightest mass number
  int             nmass;      // number of masses, mass_last - mass_first + 1
  int             ncharge;    // number of distributed Z for a given A


  FissionFragmentPair(){
    allocated = false;
    max_pair = 0;
    ntotal   = 0;
    gtotal   = 0;
    zf       = 0;
    af       = 0;
    ac       = 0.0;
    bn       = 0.0;
    en       = 0.0;
    ex       = 0.0;
    fraction = 1.0;

    mass_first = 0;
    mass_last  = 0;
    nmass    = 0;
    ncharge  = 0;

    for(int i=0 ; i<ng ; i++) gauss[i].clear();
  }

  ~FissionFragmentPair(){
    memfree();
  }

  void memalloc(const int n){
    if(!allocated){
      max_pair = n;
      fragment = new ZAPair [max_pair];
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] fragment;
      allocated = false;
    }
  }

  void setZA(int z, int a){
    zf = z;
    af = a;
    ac = af * 0.5;
  }

  void setGauss(double x1, double x2, double x3){
    if(gtotal < ng){
      gauss[gtotal].setVal(x1,x2,x3);
      gtotal++;
    }
  }

  double getGaussWidth(int k){
    double x = (k < gtotal) ? gauss[k].width : 0.0;
    return x;
  }

  double getGaussCenter(int k){
    double x = (k < gtotal) ? gauss[k].center : 0.0;
    return x;
  }

  double getGaussFraction(int k){
    double x = (k < gtotal) ? gauss[k].fraction : 0.0;
    return x;
  }

  void setN(int n){ ntotal = n; }
  void setG(int n){ gtotal = n; }
  void resetN(){ gtotal = 0; }
  int getN(){ return ntotal; }
  int getNg(){ return gtotal; }
  int getMaxPair(){ return max_pair; }

  void setZARange(const int a0, const int zrange){
    int am = af / 2.0;
    mass_first = a0;
    mass_last  = am + (am - a0);
    nmass      = mass_last - mass_first + 1;
    ncharge    = 2 * zrange + 1;
  }
};


/**************************************/
/*      ffpparameter.cpp              */
/**************************************/
void   FFPGaussianParameters (const bool, const int, const int, const double, const double, double *, double *, double *);
double FFPSystematics_TKE_A (const bool, const int, const int, const int, const double);
double FFPSystematics_TKE_A_Width (const bool, const int, const int, const int);
double FFPSystematics_TKE_En (const bool, const int, const int, const double, const double);
double FFPConstructTKEA(const int, const int, double *);
double FFPConstructSigmaTKEA(const int, const int, double *);

double FFPMultichance_FissionProbability (const int, const int, const double, const int);
double FFPMultichance_TotalKineticEnergy (const int, const int, const double, const int);
double FFPMultichance_ExcitationEnergy (const int, const int, const double, const int);


/**************************************/
/*      ffpmultichance.cpp            */
/**************************************/
int FFPReadMultichanceFissionData (const std::string, const double, double *);


/**************************************/
/*      ffpyield.cpp                  */
/**************************************/
//void   FFPMassChainYieldParameters (const bool, std::string, double *, double *, double *, FissionFragmentPair *);
void   FFPConstructYieldParameters(double *, double, double *, double *, double *, FissionFragmentPair *);
//void   FFPGeneratePairs (const bool, const double, FissionFragmentPair *, double *);
void   FFPGeneratePairs (const bool, const double, FissionFragmentPair *, double *, double *, double *, const int, int *, double *);
void   FFPMassYieldParameters (const bool, std::string, double *, double *, double *, FissionFragmentPair *);
void   FFPGeneratePairs (const bool, const double, FissionFragmentPair *, double *, const int, int *, double *);
void   FFPGeneratePairsExpdata (FissionFragmentPair *, double *, double **, const int, const int, const int, const int);
double FFPAverageTKE (FissionFragmentPair *);
void   FFPExpandList (const int, FissionFragmentPair *);


/**************************************/
/*      ffpyieldread.cpp              */
/**************************************/
void   FFPReadPairs (const std::string, FissionFragmentPair *);
void   FFPReadExpdata (const std::string, FissionFragmentPair *, double *);


/**************************************/
/*      ffpoutput.cpp                 */
/**************************************/
void   FFPOutput (const double, const double, const double, FissionObservable *);
void   FFPOutputSpec (const double, const double, FissionObservable *);
void   FFPOutputPairs (const double, FissionFragmentPair *);
void   FFPOutputTXEDistribution (FissionFragmentPair *);
void   FFPOutputTKEDistribution (const int, FissionFragmentPair *);
void   FFPOutputYield (const int, const int, const int, double *, double *, double *);
void   FFPOutputIndividualSpectrum (ZAPair *, const int, const double, double *, double *, double *);


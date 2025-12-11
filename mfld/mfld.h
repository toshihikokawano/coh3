#ifndef __PHYSICALCONSTANT_H__
#define __PHYSICALCONSTANT_H__
#include "physicalconstant.h"
#endif

#ifndef __HF_H__
#define __HF_H__
#include "HF.h"
#endif

#ifndef __COMPLEX_H__
#define __COMPLEX_H__
#include <complex>
#endif

/****************************/
/*      ARRAY SIZE          */
/****************************/

const int MAX_SPSTATE       =   300;     // max single particle states
const int MAX_ENERGY_BIN    =   501;
const int MAX_PARTICLE      =    10;
const int MAX_M             =    16;
const int MAX_J             =    40;

const double ENERGY_BIN     =   0.2;     // default energy bin in the continuum
const double INTEG_WIDTH    =   0.1;     // default wave func. integral width
const double ENERGY_AVERAGE =   2.0;     // default smoothing width for J-dist

typedef enum {neutron=0, proton=1} Particle;
typedef enum {even=0, odd=1} Parity;

/**********************************************************/
/*   Z and A numbers of Nucleus                           */
/**********************************************************/
class ZAnumber{ 
 private:
    unsigned int Z;
    unsigned int A;
 public:
    ZAnumber(){
      Z = 0;
      A = 0;
    }
    ZAnumber(int z, int a){
      Z = z;
      A = a;
    }
    void setZA(int z, int a){
      Z = z;
      A = a;
    }
    unsigned int getZ(){ return (Z); }
    unsigned int getA(){ return (A); }
    unsigned int getN(){ return (A-Z); }
    ZAnumber operator+(ZAnumber x){
      ZAnumber y;
      y.Z = Z + x.Z;
      y.A = A + x.A;
      return y;
    }
    ZAnumber operator-(ZAnumber x){
      ZAnumber y;
      y.Z = Z - x.Z;
      y.A = A - x.A;
      return y;
    }
    bool operator==(ZAnumber x){
      bool z = false;
      if( (Z == x.Z) && (A == x.A) ) z = true;
      return z;
    }
};


/**********************************************************/
/*    System Data                                         */
/**********************************************************/
class SystemData{
 public:
  ZAnumber   target;        // nucleus Z and A numbers
  int        lsize;         // max size for 
  int        nsize;         // max size for energy mesh
  int        psize;         // max size of p-h configuration
  int        msize;         // max size of m-state
  int        jsize;         // max size for spin
  int        phmax;         // actual max size of p-h configuration
  double     de;            // energy bin width
  double     emax;          // highest excitation energy
  double     eave;          // energy width for J-dist printing
  double     spmax;         // highest single particle energy
  bool       sharpcutoff ;  // no BCS flag

  SystemData(){
    sharpcutoff  = false;
    target.setZA(0,0);
    lsize = 0;
    nsize = 0;
    psize = 0;
    msize = 0;
    jsize = 0;
    phmax = 0;
    de    = 0.0;
    emax  = 0.0;
    spmax = 0.0;
    eave  = ENERGY_AVERAGE;
  }

  ~SystemData(){
  }
};


/**********************************************************/
/*    Output Data                                         */
/**********************************************************/
class PrintData{
 public:
  bool TotalDensity;      // print total level density
  bool PartialDensity;    // print partial level density
  bool JDependentDensity; // J-dependence of level density
  bool BothParity;        // add even and odd parities
  bool Plotting;          // for plotting the results
  bool SpinDistribution;  // calculate spin distribution
  bool SpinCutOff;        // extract spin cut-off parameter
  bool AverageConfig;     // average number of excitons

  PrintData(){
    TotalDensity      = false;
    PartialDensity    = false;
    JDependentDensity = false;
    BothParity        = false;
    Plotting          = false;
    SpinDistribution  = false;
    SpinCutOff        = false;
    AverageConfig     = false;
  }
};


/**********************************************************/
/*      Partial Level Density                             */
/**********************************************************/
class Density{
 private:
  int    psize;  // max number of particles
  int    nsize;  // max number of energy bins
  int    jsize;  // max number of J states
  bool   allocated;
 public:
  double    de;  // energy bin size
  double ***r0;  // even pariry
  double ***r1;  // odd parity

  Density(){
    psize = 0;
    nsize = 0;
    jsize = 0;
    de    = 0.0;
    allocated = false;
  }

  ~Density(){
    memfree();
  }

  void memalloc(int p, int n, int j){
    if(!allocated){
      psize = p;
      nsize = n;
      jsize = j;

      r0 = new double ** [psize];
      r1 = new double ** [psize];
      for(p=0 ; p<psize ; p++){
        r0[p] = new double * [nsize];
        r1[p] = new double * [nsize];
        for(n=0 ; n<nsize ; n++){
          r0[p][n] = new double [jsize];
          r1[p][n] = new double [jsize];
        }
      }
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      for(int p=0 ; p<psize ; p++){
        for(int n=0 ; n<nsize ; n++){
          delete [] r0[p][n];
          delete [] r1[p][n];
        }
        delete [] r0[p];
        delete [] r1[p];
      }
      delete [] r0;
      delete [] r1;
      
      allocated = false;
    }
  }

  void memclear(){
    if(allocated){
      for(int p=0 ; p<psize ; p++){
        for(int n=0 ; n<nsize ; n++){
          for(int j=0 ; j<jsize ; j++) r0[p][n][j] = r1[p][n][j] = 0.0;
        }
      }
    }
  }

  double binwidth(){ return de; }

  int getPsize(){return psize;}
  int getNsize(){return nsize;}
  int getJsize(){return jsize;}
};



/**********************************************************/
/*      Single Particle State for Combinatorial Calc.     */
/**********************************************************/
class State{
 private:
  int    m;  // spin projection
  int    p;  // parity
  int    t;  // isospin (doubled)
  double e;  // binding energy
  double v2; // ocupation probability
 public:
  State(){
    m = 0; p = 0; t = 0; e = 0.0; v2 = 0.0;
  }

  void set(int k1, int k2, int k3, double c1, double c2){
    m = k1; p = k2; t = k3; e = c1; v2 = c2; 
  }
  void   setE (double x) { e  = x; }
  void   setV2(double x) { v2 = x; }
  void   setM(int k) { m = k; }
  void   setP(int k) { p = k; }
  void   setT(int k) { t = k; }
  int    getM() { return(m);  }
  int    getP() { return(p);  }
  int    getT() { return(t);  }
  double getE() { return(e);  }
  double getV2(){ return(v2); }
};


class SingleParticle{
 private:
  int    nalloc;
  int    nlev;
  bool   allocated;
  int    ng;
  bool   gridallocated;
 public:
  int    p0;     // lowest particle state orbit number
  int    p1;     // highest particle state
  int    h0;     // lowest hole state orbit number
  int    h1;     // highest hole state
  double lambda; // chemical potential
  State *state;
  GridData *wf;

  SingleParticle(){
    p0 = p1 = h0 = h1 = 0;
    nalloc = 0;
    nlev = 0;
    lambda = 0.0;
    state = NULL;
    allocated = false;
    ng = 0;
    gridallocated = false;
  }

  ~SingleParticle(){
    memfree();
    gridfree();
  }

  void memalloc(int n){
    if(!allocated){
      nalloc = n;
      state = new State [nalloc];
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] state;
      allocated = false;
    }
  }

  void gridalloc(int n){
    if(!gridallocated){
      ng = n;
      wf = new GridData [ng];
      gridallocated = true;
    }
  }

  void gridfree(){
    if(gridallocated){
      delete [] wf;
      gridallocated = false;
    }
  }

  void set(int k1, int k2, int k3, double c1, double c2){
    if(nlev < nalloc) state[nlev++].set(k1,k2,k3,c1,c2);
  }
  void setE (int i, double x) { state[i].setE(x) ; }
  void setV2(int i, double x) { state[i].setV2(x); }
  void setM(int i, int k) { state[i].setM(k); }
  void setP(int i, int k) { state[i].setP(k); }
  void setT(int i, int k) { state[i].setT(k); }

  double getE (int i) { return state[i].getE() ; }
  double getV2(int i) { return state[i].getV2(); }
  int getM(int i) { return state[i].getM(); }
  int getP(int i) { return state[i].getP(); }
  int getT(int i) { return state[i].getT(); }

  int getNmax(){ return nlev; }
};



/**********************************************************/
/*      M-Counting                                        */
/*      ----------                                        */
/*  This object contains the total number of M states for */
/*  each different Np-Nh configuraion in N and Z shells.  */
/*  The structure is:                                     */
/*  m0[Energy][M(0:Max)] = number of M (even parity)      */
/*  m1[Energy][M(0:Max)] = number of M (odd parity)       */
/*  note that M should be symmetric, which is implicit.   */
/**********************************************************/
class MCount{
 private:
  int    msize;  // max number of m states
  int    nsize;  // max number of energy bins
  bool   allocated;
 public:
  unsigned long int **m0, **m1;

  MCount(){
    msize = 0;
    nsize = 0;
    allocated = false;
  }

  ~MCount(){
    memfree();
  }

  void memalloc(int m, int n){
    if(!allocated){
      msize = m;
      nsize = n;
      m0 = new unsigned long int * [nsize];
      m1 = new unsigned long int * [nsize];
      for(int j=0 ; j<nsize ; j++){
        m0[j] = new unsigned long [msize+1];
        m1[j] = new unsigned long [msize+1];
      }
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      for(int j=0 ; j<nsize ; j++){
        delete [] m0[j];
        delete [] m1[j];
      }
      delete [] m0;
      delete [] m1;
      allocated = false;
    }
  }

  void clear(){
    for(int j=0 ; j<nsize ; j++){
      for(int m=0 ; m<=msize ; m++) m0[j][m] = m1[j][m] = 0;
    }
  }

  int getMsize(){ return msize; }
  int getNsize(){ return nsize; }
};


/**********************************************************/
/*      M-State Configuration                             */
/*      ---------------------                             */
/*  This contains cumulative quantities for recursive     */
/*  call for the configuration: N-particle out of Nlev    */
/*  For example, mvalue array will have                   */
/*  mvalue[0] = M1 (one selected leved from Nlev)         */
/*        [1] = M1 + M2 (sum of two different levels)     */
/*        [2] = M1 + M2 + M3                              */
/*  therefore the last element contains the total M       */
/**********************************************************/
class MConfig{
 private:
  int    csize;  // max number of configrations
  bool   allocated;
 public:
  int    *id;       // running index for multi-fold loop
  int    *step;     // index increment = 1;
  int    *parity;   // porduct of parities
  int    *mvalue;   // cumulative M-values
  double *energy;   // cumulative excitation energies
  double *fraction; // product of fraction of states
  double  emax;
  double  de;

  MConfig(){
    csize = 0;
    emax  = 0.0;
    de    = 0.0;
    allocated =false;
  }

  ~MConfig(){
    memfree();
  }

  void memalloc(int k){
    if(!allocated){
      csize    = k + 1;
      id       = new int [csize];
      step     = new int [csize];
      parity   = new int [csize];
      mvalue   = new int [csize];
      energy   = new double [csize];
      fraction = new double [csize];
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] id;
      delete [] step;
      delete [] parity;
      delete [] energy;
      delete [] fraction;
      delete [] mvalue;
      allocated = false;
    }
  }

  void setstep(int k){
    for(int i=0 ; i<csize ; i++) step[i] = k;
  }

  void init(double x1, double x2){
    emax = x1;
    de   = x2;
    energy[0] = 0.0;
    id[0] = -1;
    fraction[0] = 1.0;
    parity[0] = 1;
    mvalue[0] = 0;
    setstep(1);
  }
};


/**********************************************************/
/*      Count Number of Configuration for Given M and Tz  */
/**********************************************************/
class ConfigCount{
 private:
  int    msize;  // max number of M states
  int    tsize;  // max Tz
  bool   allocated;
  double *bm0, **pm0;
  double *bm1, **pm1;
 public:
  double **m0;  // summed even parity configurations
  double **m1;  // summed odd parity configurations

  ConfigCount(){
    msize = 0;
    tsize = 0;
    allocated = false;
  }

  ~ConfigCount(){
    memfree();
  }

  void memalloc(int m, int t){
    if(!allocated){
      msize = m;
      tsize = t;
      int wsize = 2*tsize + 1;
      int xsize = wsize * (msize+1);
      bm0 = new double [xsize];
      bm1 = new double [xsize];

      pm0 = new double * [wsize];
      pm1 = new double * [wsize];

      for(int i=0 ; i<wsize ; i++){
        pm0[i] = &bm0[i*(msize+1)];
        pm1[i] = &bm1[i*(msize+1)];
      }

      m0 = &pm0[tsize];
      m1 = &pm1[tsize];

      memclear();
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] pm0;
      delete [] pm1;
      delete [] bm0;
      delete [] bm1;
      allocated = false;
    }
  }

  void memclear(){
    if(allocated){
      for(int t=-tsize ; t<=tsize ; t++){
        for(int m=0 ; m<=msize ; m++) m0[t][m] = m1[t][m] = 0.0;
      }
    }
  }

  void increment(int p, int t, int m, double x){
    if((0 <= m) && (m <= msize) && (abs(t) <= tsize)){
      if(p > 0){ m0[t][m] += x; }
      else     { m1[t][m] += x; }
    }
  }

  int getMsize(){return msize;}
  int getTsize(){return tsize;}
};



//   mflddataread.cpp
int  readSystem (char *, SystemData *, PrintData *, double *);

//   mfldFRDM.cpp
void FRDMCalc (MFTSystem *, HFInterface *);

//   mfldprint.cpp
void printLevelDensity (const bool, Density *);
void printPartialLevelDensity (const bool, const bool, Density *);
void printPlottingLevelDensity (const bool, Density *);
void printJDependentLevelDensity (const bool, const bool, Density *);
void printSpinDistribution (const bool, const double, Density *);
void printAverageConfiguration (Density *);


void printMState (std::string, const int, const int, MCount *);
void printCCount (const int, const int, ConfigCount *);
void printMCount (const int, const int, double **, double **);

//   mfldsd.cpp
void SDUnperturbed (SystemData *, HFInterface *, Density *);
void SDStoreState (const Particle, HFInterface *, SingleParticle *);
void SDCombination (const int, SystemData *, const int, State *, MCount *);
void SDRecursive (const int, const int, const int, State *, MConfig *, MCount *);

//   mfldphconfig.cpp
void PHCount (const bool, SingleParticle *);
double PHPairingGap (SingleParticle *);


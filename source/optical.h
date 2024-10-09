// structure definition used in optical model calculations
// define global constants
// functions used in optical model

#ifndef __COMPLEX_H__
#define __COMPLEX_H__
#include <complex>
#endif

#ifndef __PHYSICALCONSTANT_H__
#define __PHYSICALCONSTANT_H__
#include "physicalconstant.h"
#endif

#ifndef __CONSTANT_H__
#define __CONSTANT_H__
#include "constant.h"
#endif

#ifndef __PARITY_FUNC__
#define __PARITY_FUNC__
static inline int parity(int x){ return((((x)/2)%2==0) ? 1 : -1); }
#endif

/**************************************/
/*      Constants                     */
/**************************************/
const int MAX_ITER         =  100;     // maximum iteration number
const int MAX_L            =   60;     // maximum angular momentum + 2
const int MAX_L0           =  100;     // maximum iteration for F function
const int MAX_CCMATRIX     =  120;     // maximal number of coupled equations
const int MAX_LTRANS       =   10;     // maximum angular momentum transfer
const int MAX_POINTS       =  300;     // maximum radial points

const double INTEG_WIDTH   =  0.1    ; // default wave func. integral width
const double CRIT_MATCHING =  1.0e-06; // matching radius criterion
const double CRIT_LMAX     =  1.0e-10; // maximum angular momentum criterion
const double CRIT_LCUT     =  1.0e-06; // reaction cross section converge
const double CRIT_LCUTCC   =  1.0e-08; // reaction converge for CC
const double CRIT_SCUT     =  1.0e-06; // spin distribution cut off
const double CRIT_WRONSK   =  1.0e-06; // wronskian satisfaction criterion
const double CRIT_ITER     =  1.0e-16; // infinit iteration criterion
const double CRIT_LEVDEN   =  1.0e-06; // level density matching criterion
const double ECUT_CHARGED  =  1.0e-10; // cutoff energy for charged particles


/**************************************/
/*      Enum                          */
/**************************************/
enum nuclearmodel{
  spherical = 0, rotation = 1, vibration = 2};
typedef enum nuclearmodel NuclearModel;


/**************************************/
/*      Structure data                */
/**************************************/

/****************************/
/*   OM Potential Parameter */
/****************************/
class Optical{
 private:
  unsigned int index; // index and potential form for the potential
 public:
    double r0  ;double r0s ;double rv  ;double rs  ;double rvso;double rwso;
    double a0  ;double a0s ;double av  ;double as  ;double avso;double awso;
    double v1  ;double vs1 ;double wv1 ;double ws1 ;double vso1;double wso1;
    double v2  ;double vs2 ;double wv2 ;double ws2 ;double vso2;double wso2;
    double v3  ;double vs3 ;double wv3 ;double ws3 ;double vso3;double wso3;
    double v4  ;
    double rc  ;
    double R0  ;double R0s ;double Rv  ;double Rs  ;double Rvso;double Rwso;
    double Rc  ;
    std::complex<double> volume    ;
    std::complex<double> surface   ;
    std::complex<double> spin_orbit;

    Optical(){
      init();
    }

    void init(){
      index  = 0;
      r0  = 0.0; r0s = 0.0; rv  = 0.0; rs  = 0.0; rvso= 0.0; rwso= 0.0;
      a0  = 0.0; a0s = 0.0; av  = 0.0; as  = 0.0; avso= 0.0; awso= 0.0;
      v1  = 0.0; vs1 = 0.0; wv1 = 0.0; ws1 = 0.0; vso1= 0.0; wso1= 0.0;
      v2  = 0.0; vs2 = 0.0; wv2 = 0.0; ws2 = 0.0; vso2= 0.0; wso2= 0.0;
      v3  = 0.0; vs3 = 0.0; wv3 = 0.0; ws3 = 0.0; vso3= 0.0; wso3= 0.0;
      v4  = 0.0;
      rc  = 0.0;
      R0  = 0.0; R0s = 0.0; Rv  = 0.0; Rs  = 0.0; Rvso= 0.0; Rwso= 0.0;
      Rc  = 0.0;
    }

    void setIndex(unsigned int k){ index = k; }

    unsigned int getIndex()   { return(index); }
    unsigned int getPotForm() { return(index & 0x00ff); }
    unsigned int getRealForm(){ return( (index & 0x00ff) >> 4 ); }
    unsigned int getImagForm(){ return(index & 0x000f); }

    void setDepth(double e){
      double e2 = e*e;
      double e3 = e*e*e;
      double vr = v1   + v2   *e + v3   *e2 + v4 * e3;
      double vi = wv1  + wv2  *e + wv3  *e2;

      double dr = vs1  + vs2  *e + vs3  *e2;
      double di = ws1  + ws2  *e + ws3  *e2;

      double sr = vso1 + vso2 *e + vso3 *e2;
      double si = wso1 + wso2 *e + wso3 *e2;

      volume     = std::complex<double>(vr,vi);
      surface    = std::complex<double>(dr,di);
      spin_orbit = std::complex<double>(sr,si);
    }

    void setReducedRadius(double a3){
      R0   = r0  *a3;
      R0s  = r0s *a3;
      Rv   = rv  *a3;
      Rs   = rs  *a3;
      Rvso = rvso*a3;
      Rwso = rwso*a3;
      Rc   = rc  *a3;
    }
};


/****************************/
/*   Optical Potential      */
/****************************/
class Potential{
 private:
  bool alloc_potential;
  bool alloc_spinorbit;
  bool alloc_coulomb;
  bool alloc_radius;
 public:
  int       n_match;            // integration step number
  double    rad_match;          // matching radius
  double    rho_match;          // rad_match * wave_number
  double    width;              // integration step width
  double   *r2inv;              // 1.0/(r*r)
  double   *coulomb;            // Coulomb field
  std::complex<double>  *mean_field; // Optical potential V   + iW
  std::complex<double>  *spin_orbit; // Optical potenital Vso + iWso

  Potential(){
    init();
  }

  ~Potential(){
    memfree();
  }

  void init (){
    n_match   = 0;
    rad_match = 0.0;
    rho_match = 0.0;
    width     = INTEG_WIDTH;
    alloc_potential = false;
    alloc_coulomb   = false;
    alloc_spinorbit = false;
    alloc_radius    = false;
  }

  void memalloc(int n){
    memalloc_potential(n);
    memalloc_coulomb(n);
    memalloc_spinorbit(n);
    memalloc_radius(n);
  }
  
  void memalloc_potential(int n){
    if(!alloc_potential){
      mean_field = new std::complex<double> [n];
      alloc_potential = true;
    }
  }
  
  void memalloc_coulomb(int n){
    if(!alloc_coulomb){
      coulomb = new double [n];
      alloc_coulomb = true;
    }
  }
  
  void memalloc_spinorbit(int n){
    if(!alloc_spinorbit){
      spin_orbit = new std::complex<double> [n];
      alloc_spinorbit = true;
    }
  }
  
  void memalloc_radius(int n){
    if(!alloc_radius){
      r2inv = new double [n];
      alloc_radius = true;
    }
  }

  void memfree(){
   if(alloc_potential){
     delete [] mean_field;
     alloc_potential = false;
   }
   if(alloc_coulomb){
     delete [] coulomb;
     alloc_coulomb = false;
   }
   if(alloc_spinorbit){
    delete [] spin_orbit;
    alloc_spinorbit = false;
   }
   if(alloc_radius){
     delete [] r2inv;
     alloc_radius = false;
   }
  }

  bool memstat_potential(){ return alloc_potential; }
  bool memstat_spinorbit(){ return alloc_spinorbit; }
  bool memstat_coulomb()  { return alloc_coulomb; }
  bool memstat_radius()   { return alloc_radius; }
  
  void setMatching(int n){
    n_match   = n;
    rad_match = width * (n-3);
  }
  void setRho(double k){
    rho_match = rad_match * k;
  }
};


/****************************/
/*   Unbound State Data     */
/****************************/
class Wavefunc{
 private:
  bool alloc_ext;
  bool alloc_int;
 public:
  std::complex<double> *external;     // external wave func.  G, F at Rm
  std::complex<double> *extderiv;     // derivative ext.wave dG,dF at Rm
  std::complex<double> *internal;     // internal wave func. R=0 to Rm

  Wavefunc(){
    alloc_ext = false;
    alloc_int = false;
  }
  ~Wavefunc(){
    if(alloc_ext){
      delete [] external;
      delete [] extderiv;
      alloc_ext = false;
    }
    if(alloc_int){
      delete [] internal;
      alloc_int = false;
    }
  }

  void memalloc(int n, int m){
    memalloc_ext(n);
    memalloc_int(m);
  }
  
  void memalloc_ext(int n){
    external  = new std::complex<double> [n];
    extderiv  = new std::complex<double> [n];
    alloc_ext = true;
  }

  void memalloc_int(int m){
    internal  = new std::complex<double> [m];
    alloc_int = true;
  }

  void memfree(){
    if(alloc_ext){
      delete [] external;
      delete [] extderiv;
      alloc_ext = false; 
    }
    if(alloc_int){
      delete [] internal;
      alloc_int = false;
    }
  }
};


/****************************/
/*   Channel Spin           */
/****************************/
class Chspin{
 public:
    int      l;               // angular momentum for channel
    int      j2;              // 2 x j (j = l+s)
    int      pt;              // parity of target
    int      pi;              // parity of emitted particle
    double   si2;             // 2 x spin of emitted particle
    double   st2;             // 2 x spin of target

    Chspin(){
      l   = 0;
      j2  = 0;
      pt  = 0;
      pi  = 0;
      si2 = 0;
      st2 = 0;
    }
};


/****************************/
/*   Discrete State Data    */
/****************************/
struct coulombData{
    double  eta;              // Coulomb parameter
    double  sigma0;           // sigma0
};
typedef struct coulombData CoulombData;

class LevelData{
 public:
    int     phonon;           // number of phonons
    double  excitation;       // level excitation
    double  spin;             // state spin and parity
    double  energy;           // cms energy for this level
    double  reduced_mass;     // reduced mass
    double  wave_number;      // wave number for this level
    double  wavesq;           // square of wave number
    CoulombData coulomb;      // Coulomb data

    LevelData(){
      clear();
    }
    void clear(){
      phonon         = 0;
      excitation     = 0.0;
      spin           = 0.0;
      energy         = 0.0;
      reduced_mass   = 0.0;
      wave_number    = 0.0;
      wavesq         = 0.0;
      coulomb.eta    = 0.0;
      coulomb.sigma0 = 0.0;
    }
};


/****************************/
/*   Channel Data           */
/****************************/
class CCdata{
 public:
    bool     open;            // open / closed channel
    Chspin   chn;             // channel spin
    int      level;           // index for discrete level
    int      lmax;            // max angular momentum for channel
    int      jmax;            // max spin
    LevelData *lev;           // pointer to discrete level data
    double   xl_term;         // L(L+1)
    double   so_term;         // J(J+1)-L(L+1)-S(S+1)

    CCdata(){
      open    = true;
      level   = 0;
      lmax    = 0;
      jmax    = 0;
      xl_term = 0.0;
      so_term = 0.0;
    }
};


/****************************/
/*   Collective Excite      */
/****************************/
class Collective{
 private:
    bool    allocated;
    int     ndim;
    int     mdim;
 public:
    int     nlevel;           // number of collective levels
    int     max_lambda;       // maximal order of expansion
    int     target_index;     // which is the target state ?
    int     *pindex;          // index for radial part potential
    double  pspin;            // projectile spin
    double  *beta;            // deformation parameter pointer
    LevelData *lev;           // discrete level data
    double  **qmatrix;        // nuclear matrix element

    Collective(){
      nlevel        = 0;
      max_lambda    = 0;
      target_index  = 0;
      pspin         = 0.0;
      allocated     = false;
      ndim = mdim = 0;
    }

    int memalloc(int n, int m){
      try{
        lev = new LevelData [n];
        qmatrix = new double * [m];
        for(int i=0 ; i<m ; i++){
          qmatrix[i] = new double [n*(n+1)/2];
        }
        pindex = new int [n*(n+1)/2];
        ndim = n;
        mdim = m;
      }
      catch(std::bad_alloc &e){ return -1; }
      for(int i=0 ; i<n ; i++) lev[i].clear();

      allocated = true;
      return 0;
    }

    void memfree(){
      delete [] lev;
      for(int i=0 ; i<mdim ; i++){
        delete [] qmatrix[i];
      }
      delete [] qmatrix;
      delete [] pindex;
      allocated = false;
    }
};



/**************************************/
/*      amplitud.cpp                  */
/**************************************/
double  omScatteringAmplitudeN (const int, double *, const double, const double, std::complex<double> *);
double  omScatteringAmplitudeP (const int, double *, const double, const double, const CoulombData, std::complex<double> *);
double  omScatteringAmplitudeA (const int, double *, const double, const double, const CoulombData, std::complex<double> *);
double  omScatteringAmplitudeD (const int, double *, const double, const double, const CoulombData, std::complex<double> *);
std::complex<double> omRutherfordAmplitude (const double, const double, const CoulombData);


/**************************************/
/*      asympt.cpp                    */
/**************************************/
void    omAsymptotic (const double, const double, const double, double *, double *);
double  omAsymptoticClosed (const double, const double);


/**************************************/
/*      builtin.cpp                   */
/**************************************/
unsigned int  omp_library (const int, const int, const int, const int, const int, const double, Optical *);
unsigned int  find_omp (std::string);


/**************************************/
/*      ccamp.cpp                     */
/**************************************/
void    ccAllocateAmplitute (Collective *);
void    ccFreeAmplitude (Collective *);

void    ccScatteringAmplitude (const int, const int, double *, Collective *, CCdata *, std::complex<double> *);
double  ccAngularDistribution (const int, const int, double *, Collective *);


/**************************************/
/*      ccfunction.cpp                */
/**************************************/
int     ccNumberOfChannels (Collective *);
int     ccSetChannel (const int, const int, const int, Collective *, CCdata *);
void    ccMatrixElementRot (const int, Collective *);
void    ccMatrixElementVib (Collective *);
void    ccCouplingPotential (const NuclearModel, const int, const int, double **, std::complex<double> **, CCdata *, Collective *, Potential *, Potential **, Potential **);
int     ccExternalFunction (const double, Wavefunc *, Collective *);
int     ccFindTargetLevel (const double, Collective *);
int     ccBandheadSpin (Collective *);


/**************************************/
/*      ccsetform.cpp                 */
/**************************************/
void    ccPotentialFormRot (const int, const int, double *, const int, Optical *, double, Potential *, Potential **);
void    ccPotentialFormVib (const int, const int, Optical *, const double, Potential *, Potential **, Potential **);

  
/**************************************/
/*      extwave.cpp                   */
/**************************************/
int     omExternalFunction (int, const double, const double, Wavefunc *);
int     omExternalClosed (int, const double, const double, Wavefunc *);


/**************************************/
/*      intwave.cpp                   */
/**************************************/
void    omInternalFunction (const int, const double, const double, const double, const double, const double, Potential *, Wavefunc *);
void    omIrregularFunction (const int, const double, const double, const double, const double, const double, Potential *, Wavefunc *);
std::complex<double> omNormalizationFactor (const int, const int, const CoulombData, std::complex<double>, Wavefunc *);
std::complex<double> omCoulombPhaseFactor (const int, CoulombData);
void    omNormalization (const int, const double, std::complex<double>, Potential *, Wavefunc *);
void    ccInternalFunction (const int, const int, const double, CCdata *, Wavefunc *, std::complex<double> **, std::complex<double> *[]);


/**************************************/
/*      omsetform.cpp                 */
/**************************************/
int     omPotentialForm (const int, Optical *, LevelData *,Potential *);
void    omPotentialFixedLength (const int, Optical *, const double, Potential *);
double  omPotentialRadialCoulomb (const double, Optical *);
int     dwCollectiveForm (const double, Optical *, Potential *);
std::complex<double> omPotentialRadialMean (const double, Optical *);
std::complex<double> omPotentialRadialSpo (const double, Optical *);
std::complex<double> omPotentialRadialFirstDeriv (const double, Optical *);
std::complex<double> omPotentialRadialSecondDeriv (const double, Optical *);


/**************************************/
/*      omsetparm.cpp                 */
/**************************************/
void    omSetEnergy (const double, const int, const double, LevelData *);
int     ccSetEnergy (const double, const int, const double, Collective *);
unsigned int omSetOmp (const int, const double, const int, const int, const int, const int, Optical *);


/**************************************/
/*      smatrix.cpp                   */
/**************************************/
std::complex<double> omSmatrix (const int, const double, const double, const std::complex<double>, const std::complex<double>, Wavefunc *);
int     ccSmatrix (const int, const double, CCdata *, Wavefunc *, std::complex<double> *, std::complex<double> *, std::complex<double> *);


/**************************************/
/*      dwcalc.cpp                    */
/**************************************/
void    dwStoreDistortedWave (const int, const double, CCdata *, Potential *, Wavefunc *,std::complex<double> *);


/**************************************/
/*      strength.cpp                  */
/**************************************/
void    omStrengthFunction (const int, const int, const double, const double, Wavefunc *, std::complex<double> *);


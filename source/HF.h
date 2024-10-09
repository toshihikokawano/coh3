#ifndef __PHYSICALCONSTANT_H__
#define __PHYSICALCONSTANT_H__
#include "physicalconstant.h"
#endif

#ifndef __HFINTERFACE_H__
#define __HFINTERFACE_H__
#include "HFinterface.h"
#endif

#ifndef __MFT_H__
#define __MFT_H__
#include "MFT.h"
#endif


/**********************************************************/
/*      Skyrme Force Parameters                           */
/**********************************************************/
class Skyrme{
 public:
  double b1;  double b2;  double b3;  double b4;
  double b5;  double b6;  double b7;  double b8;
  double b9;
  double alpha;  double w;  double v0;

  Skyrme(){
    b1 = 0.0;    b2 = 0.0;    b3 = 0.0;    b4 = 0.0;
    b5 = 0.0;    b6 = 0.0;    b7 = 0.0;    b8 = 0.0;
    b9 = 0.0;
    alpha = 0.0;    w = 0.0;    v0 = 0.0;
  }
};


/**********************************************************/
/*      Multipole Moments and Constraints in Iteration    */
/**********************************************************/
class Moment{
 private:
  double value[5];
  double weight[5];
  int    nfree[5];
 public:
  double zcm;
  double q20;
  double q30;
  double q40;
  double qn;
  double zmin;

  Moment(){
    for(int i=0 ; i<5 ; i++){
      value[i]  = 0.0;
      weight[i] = 0.0;
      nfree[i]  = 0;
    }
    zcm    = 0.0;
    q20    = 0.0;
    q30    = 0.0;
    q40    = 0.0;
    qn     = 0.0;
    zmin   = 0.0;
  }

  void setConstraint(int i, double v, double w, int n){
    if( (0 <= i) && (i < 5) ){
      value[i]  = v;
      weight[i] = w;
      nfree[i]  = n;
    }
  }

  double getV(int i){
    double v = 0.0;
    if( (0 <= i) && (i < 5) ){ v = value[i]; }
    return(v);
  }

  double getW(int i){
    double w = 0.0;
    if( (0 <= i) && (i < 5) ){ w = weight[i]; }
    return(w);
  }

  int getN(int i){
    int n = 0.0;
    if( (0 <= i) && (i < 5) ){ n = nfree[i]; }
    return(n);
  }
};


/**********************************************************/
/*      Energies                                          */
/**********************************************************/
class HFEnergy{
 public:
  double kinetic;
  double volume;
  double surface;
  double spinorbit;
  double coulomb;
  double pairing;
  double total;
  double nucleon;

  HFEnergy(){
    kinetic   = 0.0;
    volume    = 0.0;
    surface   = 0.0;
    spinorbit = 0.0;
    coulomb   = 0.0;
    pairing   = 0.0;
    total     = 0.0;
    nucleon   = 0.0;
  }

  void set(double x1, double x2, double x3, double x4, double x5, double x6){
    kinetic   = x1;
    volume    = x2;
    surface   = x3;
    spinorbit = x4;
    coulomb   = x5;
    pairing   = x6;
    setTotal();
  }

  void copy(HFEnergy e){
    kinetic   = e.kinetic;
    volume    = e.volume;
    surface   = e.surface;
    spinorbit = e.spinorbit;
    coulomb   = e.coulomb;
    pairing   = e.pairing;
    nucleon   = e.nucleon;
    setTotal();
  }

  double setTotal(){
    total = kinetic + volume + surface + spinorbit + coulomb + pairing + nucleon;
    return(total);
  }

  void setNucleonMass(double mx){
    nucleon = mx;
  }
};


/**********************************************************/
/*      Grid Data - Quantities Defined on the Grid        */
/**********************************************************/
class GridData{
 private:
  double *buf;
  bool   allocated;
 public:
  int    Rdim;
  int    Zdim;
  double **p;

  GridData(){
    Rdim = 0;
    Zdim = 0;
    allocated = false;
  }
  ~GridData(){
    memfree();
  }

  int memalloc(int nr, int nz){
    if(allocated) return 0;
    try{
      buf = new double [nr*(nz+1)];
      p   = new double * [nr];

      for(int i=0 ; i<nr ; i++){
        unsigned int adr = i*(nz+1) + nz/2;
        p[i] = &buf[adr];
      }
    }
    catch(std::bad_alloc &e){ return -1; }

    Zdim = nz+1;
    Rdim = nr;

    clear();

    allocated = true;
    return 0;
  }
  void memfree(){
    if(allocated){
      delete [] buf;
      delete [] p;
      allocated = false;
    }
  }

  void clear(){
    for(int i=0 ; i<Rdim*Zdim ; i++){
      buf[i] = 0.0;
    }
  }
};


/**********************************************************/
/*      BCS Parameters                                    */
/**********************************************************/
class BCSParameter{
 public:
  double strength;
  double pairing_energy;
  double chemical_potential;
  double average_gap;
  std::string pairing;

  BCSParameter(){
    strength = -1250.0;
    pairing_energy = 0.0;
    chemical_potential = 0.0;
    average_gap = 0.0;
    pairing       = "";
  }
};


/**********************************************************/
/*      Function Definitions                              */
/**********************************************************/

#ifndef __HOBASIS_H__
#define __HOBASIS_H__
#include "HObasis.h"
#endif


#ifndef COH_TOPLEVEL
extern double *fact;
#endif

/*** HFmain.cpp */
void HartreeFock (MFTSystem *, Basis *, Skyrme *, BCSParameter *, Moment *, HFInterface *);

/*** skyrmeparameter.cpp */
void SkyrmeParameter (std::string, Skyrme *);

/*** HFpotential.cpp */
void HFInitialPotential (const bool, MFTSystem *, Basis *, OneBodyPotential *, LaguerrePolynomial *, HermitePolynomial *);
void HFPotential (double, Basis *, OneBodyPotential *, GridData *, GridData *, GridData *, GridData *, GridData *, Skyrme *, Moment *, LaguerrePolynomial *, HermitePolynomial *);

/*** HFsolvebcs.cpp */
void HFSolveBCS (const int, Basis *, BCSParameter *, SPEnergy *, double *, double *, Operator *, GridData *, GridData *, LaguerrePolynomial *, HermitePolynomial *);

/*** HFdensity.cpp */
void HFDensityMatrix (Basis *, double *, double *, Operator *);
void HFLocalDensity (Basis *, double *, GridData *, GridData *, GridData *, GridData *, LaguerrePolynomial *, HermitePolynomial *);

/*** HFmultipole.cpp */
void HFMultipoleMoment (Basis *, Moment *, GridData *, LaguerrePolynomial *, HermitePolynomial *);
void HFRenormalizeLocalDensity (const int, const int, GridData *, GridData *, GridData *, GridData *, LaguerrePolynomial *, HermitePolynomial *);

/*** HFcoulomb.cpp */
void HFCoulomb (Basis *, GridData *, GridData *, LaguerrePolynomial *, HermitePolynomial *);

/*** HFenergy.cpp */
void HFCalcEnergy (const int, Basis *, GridData *, GridData *, GridData *, GridData *, GridData *, Skyrme *, BCSParameter *, HFEnergy *, LaguerrePolynomial *, HermitePolynomial *);

/*** HFprint.cpp */
void HFMonitor (const int, const double, HFEnergy *, Moment *);
void HFPrintBCS (BCSParameter *);
void HFPrintEnergy (HFEnergy *);
void HFPrintMoment (Moment *);
void HFPrintPotential (Basis *, OneBodyPotential *, LaguerrePolynomial *, HermitePolynomial *);
void HFPrintSPState (Basis *, SPEnergy *);
void HFPrintDensity (int, int, double, double, GridData *, GridData *);
void HFPrintDensityCalc (const int, const int, const double, const double, Basis *, double *, GridData *);


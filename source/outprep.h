// perform summation of all emerging ground states and isomeric states
// and print the production cross sections of these states

#ifndef __METASTABLE_H__
#define __METASTABEL_H__
#include "metastable.h"
#endif

/**********************************************************/
/*      Cumulating Production of Each Final State         */
/**********************************************************/
class ResidualProduct{
 public:
  ZAnumber     za;            // Z and A number of the nucleus
  int          metaflag;      // index of meta state, 0 for ground state
  double       energy;        // excitation energy of the state
  double       spin;          // spin of the state (doubled)
  int          parity;        // parity of the state
  double       halflife;      // half-life
  double       production;    // production cross section

  ResidualProduct(){
    clear();
  }

  void clear(){
    za.setZA(0,0);
    metaflag   = 0;
    energy     = 0.0;
    spin       = 0.0;
    parity     = 0;
    halflife   = 0.0;
    production = 0.0;
  }
};


/**********************************************************/
/*      Production Data Array                             */
/**********************************************************/
class CumulativeResidualProduct{
 private:
  bool         allocated;
  int          ntotal;
  int          ncurrent;
 public:
  unsigned int amin;
  unsigned int amax;
  unsigned int zmin;
  unsigned int zmax;
  ResidualProduct *rp;

  CumulativeResidualProduct(){
    ntotal    = 0;
    ncurrent  = 0;
    amin      = 0;
    amax      = 0;
    zmin      = 0;
    zmax      = 0;
    allocated = false;
  }

  ~CumulativeResidualProduct(){
    memfree();
  }

  void memalloc(int n){
    if(!allocated){
      ntotal = n;
      rp = new ResidualProduct [ntotal];
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] rp;
      allocated = false;
    }
  }

  void clear(){
    if(allocated){
      for(int i=0 ; i<ntotal; i++) rp[i].clear();
      ncurrent = 0;
    }
  }

  void setData(ZAnumber za, int mf, double e, double j, int p, double hl, double sig){
    if(allocated && (ncurrent < ntotal)){
      rp[ncurrent].za         = za;
      rp[ncurrent].metaflag   = mf;
      rp[ncurrent].energy     = e;
      rp[ncurrent].spin       = j;
      rp[ncurrent].parity     = p;
      rp[ncurrent].halflife   = hl;
      rp[ncurrent].production = sig;
      ncurrent ++;
    }
  }

  void addProduction(int i, double sig){
    if(allocated && (i < ncurrent)){ rp[i].production += sig; }
  }

  bool exist(ZAnumber za){
    bool f = false;
    if(allocated){
      for(int i=0 ; i<ncurrent ; i++){
        if(rp[i].za == za){ f = true; break; }
      }
    }
    return f;
  }

  int getIndex(ZAnumber za, double e){
    int id = -1;
    if(allocated){
      for(int i=0 ; i<ncurrent ; i++){
        if((rp[i].za == za) && (rp[i].energy == e)){ id = i; break; }
      }
    }
    return id;
  }

  void setRange(){
    if(allocated){
      amin = rp[0].za.getA();
      zmin = rp[0].za.getZ();
      amax = amin;
      zmax = zmin;
      for(int i=0 ; i<ncurrent ; i++){
        if(rp[i].za.getA() < amin) amin = rp[i].za.getA();
        if(rp[i].za.getA() > amax) amax = rp[i].za.getA();
        if(rp[i].za.getZ() < zmin) zmin = rp[i].za.getZ();
        if(rp[i].za.getZ() > zmax) zmax = rp[i].za.getZ();
      }
    }
  }

  int getNtotal(){ return ntotal; }
  int getNcurrent(){ return ncurrent; }
  bool isAlloc(){ return allocated; }
};


/**************************************/
/*      outprep.cpp                   */
/**************************************/
void    outPrepTotalResidual (const int, const int);
void    outPrepCumulativeResidual (const int, const double);
void    outPrepPrintResidual (void);
void    outPrepPrintIsomericRatio (void);
void    outPrepAllocateResidual (const int);
void    outPrepFreeResidual (void);

CumulativeResidualProduct *outPrepObjectPointer (void);

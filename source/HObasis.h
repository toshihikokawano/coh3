#ifndef __PHYSICALCONSTANT_H__
#define __PHYSICALCONSTANT_H__
#include "physicalconstant.h"
#endif

#ifndef __HFINTERFACE_H__
#define __HFINTERFACE_H__
#include "HFinterface.h"
#endif

/**********************************************************/
/*      Operator, Eigenvalues, and Eigenvectors           */
/**********************************************************/
class Operator{
 private:
  int ndim;
  double *buf;
  bool   allocated;
 public:
  double  *matrix;
  double  *value;
  double **vector;

  Operator(){
    ndim = 0;
    allocated = false;
  }

  ~Operator(){
    memfree();
  }

  int memalloc(int n){
    if(allocated) return 0;
    try{
      matrix = new double [n*(n+1)/2];
      value  = new double [n];
      buf    = new double [n*n];
      vector = new double * [n];
      for(int i=0 ; i<n ; i++){
        vector[i] = &buf[i*n];
      }
    }
    catch(std::bad_alloc &e){ return -1; }

    ndim = n;
    allocated = true;
    clear();

    return 0;
  }

  void memfree(){
    if(allocated){
      delete [] matrix;
      delete [] value;
      delete [] buf;
      delete [] vector;
      allocated = false;
    }
  }

  void clear(){
    for(int i=0 ; i<ndim ; i++){
      value[i] = 0.0;
      for(int j=0 ; j<ndim ; j++) vector[i][j] = 0.0;
      for(int j=0 ; j<=i ; j++)   matrix[i*(i+1)/2+j] = 0.0;
    }
  }

  int getN(){
    return ndim;
  }

  void setVal(int i, int j, double x){
    int k = (i>=j) ? i*(i+1)/2 + j : j*(j+1)/2 + i;
    matrix[k] = x;
  }

  double getVal(int i, int j){
    int k = (i>=j) ? i*(i+1)/2 + j : j*(j+1)/2 + i;
    return(matrix[k]);
  }
};


/**********************************************************/
/*      Harmonic Oscillator Basis Parameters              */
/**********************************************************/
class QuantumDeformed{
 private:
  int  nr;
  int  nz;
  int  l;
  int  s;
  int  j;
  int  p;
 public:
  QuantumDeformed(){
    nr = 0 ; nz = 0 ; l = 0 ; s = 0 ;  j = 0 ;  p = 0 ;
  }
  void set(int k1, int k2, int k3, int k4, int k5, int k6){
    nr = k1; nz = k2; l = k3; s = k4;  j = k5;  p = k6;
  }
  int getNr(){ return(nr); }
  int getNz(){ return(nz); }
  int getL(){ return(l); }
  int getS(){ return(s); }
  int getJ(){ return(j); }
  int getP(){ return(p); }
};


class Basis{
 private:
  int n0;        // number of Harmonic Oscillator basis
  bool allocated;
 public:
  int nv;        // total number of eigenvalues
  int nb;        // number of blocks having the same omega
  int nzmax;
  int nrmax;
  int lambdamax;
  double b;      // b = sqrt(hbar/m w)
  double q;
  double bp;
  double bp2;
  double bz;
  double bz2;

  QuantumDeformed *state;
  int  *omega;
  int  *index;
  int  *ndata;

  Basis(){
    n0    = 14;
    nv    = 0;
    nb    = 0;
    nzmax     = 0;
    nrmax     = 0;
    lambdamax = 0;
    b     = 0.5;
    q     = 1.0;
    bp    = 0.0;
    bp2   = 0.0;
    bz    = 0.0;
    bz2   = 0.0;

    allocated = false;
  }

  ~Basis(){
    memfree();
  }

  void memalloc(){
    if(!allocated){
      if(nv > 0){
        state = new QuantumDeformed [nv];
      }
      if(nb > 0){
        index = new int [nb+1];
        omega = new int [nb];
        ndata = new int [nb];
      }
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] state;
      delete [] omega;
      delete [] index;
      delete [] ndata;
      allocated = false;
    }
  }

  void setN0(int n){ n0 = n; }
  int getN0(){return n0;}
};


/**********************************************************/
/*      Spherical Harmonic Oscillator Expansion           */
/**********************************************************/
class QuantumSpherical{
 private:
  int  n;
  int  l;
  int  s;
  int  j;
  int  m;
 public:
  QuantumSpherical(){
    n = 0 ; l = 0 ; s = 0 ; j=0 ; m = 0 ;
  }
  void set(int k1, int k2, int k3, int k4, int k5){
    n = k1; l = k2;  s = k3, j = k4;  m = k5;
  }
  int getN(){ return(n); }
  int getL(){ return(l); }
  int getS(){ return(s); }
  int getJ(){ return(j); }
  int getM(){ return(m); }
};


class SphericalBasis{
 private:
  int nm;
  int lm;
  bool allocated;
 public:
  int nv;
  QuantumSpherical *state;

  SphericalBasis(){
    nm = 0;
    lm = 0;
    nv = 0;
    state = nullptr;
    allocated = false;
  }

  ~SphericalBasis(){
    memfree();
  }

  void setLimit(int nmax, int lmax){
    nm = nmax;
    lm = lmax;
  }

  void memalloc(){
    if(!allocated && nv > 0){
      state = new QuantumSpherical [nv];
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] state;
      allocated = false;
    }
  }

  int getN(int i){ return state[i].getN(); }
  int getL(int i){ return state[i].getL(); }
  int getS(int i){ return state[i].getS(); }
  int getJ(int i){ return state[i].getJ(); }
  int getM(int i){ return state[i].getM(); }
  int getNmax(){ return nm; }
  int getLmax(){ return lm; }
};


/**********************************************************/
/*      Hartree Fock / FRDM Potential                     */
/**********************************************************/
class OneBodyPotential{
 private:
  double *bufc;
  double *bufm;
  double *bufs;
  bool   allocated;
 public:
  int    Rdim;   // NGr
  int    Zdim;   // NGz
  double **central;
  double **effmass;
  double **spinorb;
  double **coulomb;

  OneBodyPotential(){
    Rdim      = 0;
    Zdim      = 0;
    allocated = false;
  }

  ~OneBodyPotential(){
    memfree();
  }

  int memalloc(int nr, int nz){
    if(allocated) return 0;
    try{
      bufc = new double [nr*(nz+1)];
      bufm = new double [nr*(nz+1)];
      bufs = new double [nr*(nz+1)];

      central = new double * [nr];
      effmass = new double * [nr];
      spinorb = new double * [nr];
      coulomb = effmass;

      for(int i=0 ; i<nr ; i++){
        unsigned int p = i*(nz+1) + nz/2;
        central[i] = &bufc[p];
        effmass[i] = &bufm[p];
        spinorb[i] = &bufs[p];
      }
    }
    catch(std::bad_alloc &e){ return -1; }

    Zdim = nz+1;
    Rdim = nr;

    for(int i=0 ; i<Rdim*Zdim ; i++){
      bufc[i] = 0.0;
      bufm[i] = 0.0;
      bufs[i] = 0.0;
    }

    allocated = true;
    return 0;
  }

  void memfree(){
    if(allocated){
      delete [] bufc;
      delete [] bufm;
      delete [] bufs;
      delete [] central;
      delete [] effmass;
      delete [] spinorb;
      allocated = false;
    }
  }
};


/**********************************************************/
/*      Single Particle Energy and Occupation Probability */
/**********************************************************/
class SPEnergy{
 private:
  bool    allocated;
 public:
  int     ne;
  double  lambda;
  double  delta;
  double *energy;
  double *v2;
  double *parity;
  int    *index;
  int    *block;
  int    *subindex;
  char   *truncate;

  SPEnergy(){
    ne = 0;
    lambda = 0.0;
    delta = 0.0;
    allocated = false;
  }

  ~SPEnergy(){
    memfree();
  }

  int memalloc(int n){
    if(allocated) return 0;
    try{
      energy   = new double [n];
      v2       = new double [n];
      parity   = new double [n];
      index    = new int [n];
      block    = new int [n];
      subindex = new int [n];
      truncate = new char [n];
    }
    catch(std::bad_alloc &e){ return -1; }

    ne = n;

    for(int i=0 ; i<ne ; i++){
      energy[i]   = 0.0;
      v2[i]       = 0.0;
      parity[i]   = 0.0;
      index[i]    = 0;
      block[i]    = 0;
      subindex[i] = 0;
      truncate[i] = 0;

      allocated = true;
    }
    return(0);
  }

  void memfree(){
    if(allocated){
      delete [] energy;
      delete [] v2;
      delete [] parity;
      delete [] index;
      delete [] block;
      delete [] subindex;
      delete [] truncate;
      allocated = false;
    }
  }
};


/**********************************************************/
/*      Hermite Polynomials                               */
/*      P[i][j], i=0:n, j=-nz:nz                          */
/**********************************************************/
class HermitePolynomial{
 private:
  double *x_array;
  double **P0_array;
  double **P1_array;
  double *W_array;
  double *G_array;
  int    ndim; // array size [-nz,nz] = 2*nz + 1 
  int    mdim; // array size [0, mz] = mz + 1
  bool   allocated;
 public:
  int    nz;   // number of Gauss-Hermite integration points
  int    mz;   // Max m of Hermite polynomials
  double *x;
  double **P0;
  double **P1;
  double *W;
  double *G;

  HermitePolynomial(){
    ndim = 0;
    mdim = 0;
    nz   = 0;
    mz   = 0;
    allocated = false;
  }

  ~HermitePolynomial(){
    memfree();
  }

  void setGridSize(int i, int j){
    mdim = ((i > j) ? i : j) + 1;
    ndim = j + 1;
  }

  int memalloc(){
    if(allocated) return 0;
    if(ndim == 0 && mdim == 0) return -1;
    try{
      x_array = new double [ndim];
      W_array = new double [ndim];
      G_array = new double [ndim];

      int nc = (ndim-1)/2;  // point the center
      x = &x_array[nc];
      W = &W_array[nc];
      G = &G_array[nc];

      P0_array = new double * [mdim];
      P0       = new double * [mdim];
      P1_array = new double * [mdim];
      P1       = new double * [mdim];

      for(int m=0 ; m<mdim ; m++){
        P0_array[m] = new double [ndim];
        P1_array[m] = new double [ndim];

        P0[m] = &P0_array[m][nc];
        P1[m] = &P1_array[m][nc];
      }
    }
    catch(std::bad_alloc &E){ return -1; }

    /*** initialize arrays */
    for(int n=0 ; n<ndim ; n++){
      x_array[n] = 0.0;
      W_array[n] = 0.0;
      G_array[n] = 0.0;
      for(int m=0 ; m<mdim ; m++){
        P0_array[m][n] = 0.0;
        P1_array[m][n] = 0.0;
      }
    }

    nz = (ndim-1)/2;
    mz = mdim-1;

    allocated = true;
    return 0;
  }

  void memfree(){
    if(allocated){
      delete [] x_array;
      delete [] W_array;
      delete [] G_array;
      for(int m=0 ; m<mdim ; m++){
        delete [] P0_array[m];
        delete [] P1_array[m];
      }
      delete [] P0_array;
      delete [] P1_array;
      delete [] P0;
      delete [] P1;
      allocated = false;
    }
  }

  int getNsize(){ return ndim; }
  int getMsize(){ return ndim; }
};


/**********************************************************/
/*      Laguerre Polynomials                              */
/*      P[i][l][j], i=0:n, l=0:lambdamax, j=0:nr          */
/**********************************************************/
class LaguerrePolynomial{
 private:
  int    ndim;   // array size [0, nr] = nr + 1
  int    mdim;   // array size [0, mr] = mr + 1
  int    ldim;   // max of lambda, [0,lambdamax]
  bool   allocated;
 public:
  int    nr;     // number of Gauss-Laguerre integration points
  int    mr;     // max n of Laguerre polynomials
  int    lmax;   // max lambda
  double *x;
  double ***P0;
  double ***P1;
  double *W;
  double *G;

  LaguerrePolynomial(){
    ndim = 0;
    mdim = 0;
    ldim = 0;
    nr   = 0;
    mr   = 0;
    allocated = false;
  }

  ~LaguerrePolynomial(){
    memfree();
  }

  void setGridSize(int i, int j, int k){
    mdim = ((i > j) ? i : j) + 1;
    ndim = j + 1;
    ldim = k + 1;
  }

  int memalloc(){
    if(allocated) return 0;
    if(ndim == 0 && mdim == 0) return -1;
    try{
      x = new double [ndim];
      W = new double [ndim];
      G = new double [ndim];

      P0 = new double ** [mdim];
      P1 = new double ** [mdim];
      for(int i=0 ; i<mdim ; i++){
        P0[i] = new double * [ldim];
        P1[i] = new double * [ldim];
        for(int j=0 ; j<ldim ; j++){
          P0[i][j] = new double [ndim];
          P1[i][j] = new double [ndim];
        }
      }
    }
    catch(std::bad_alloc &e){ return -1; }

    for(int i=0 ; i<ndim ; i++){
      x[i] = 0.0;
      W[i] = 0.0;
      G[i] = 0.0;
    }

    for(int m=0 ; m<mdim ; m++){
      for(int l=0 ; l<ldim ; l++){
        for(int n=0 ; n<ndim ; n++){
          P0[m][l][n] = 0.0;
          P1[m][l][n] = 0.0;
        }
      }
    }

    nr   = ndim-1;
    mr   = mdim-1;
    lmax = ldim-1;

    allocated = true;
    return 0;
  }

  void memfree(){
    if(allocated){
      delete [] x;
      delete [] W;
      delete [] G;
      for(int m=0 ; m<mdim ; m++){
        for(int l=0 ; l<ldim ; l++){
          delete [] P0[m][l];
          delete [] P1[m][l];
        }
        delete [] P0[m];
        delete [] P1[m];
      }
      delete [] P0;
      delete [] P1;
      allocated = false;
    }
  }

  int getNsize(){ return ndim; }
  int getMsize(){ return mdim; }
  int getLsize(){ return ldim; }
};



/*** HObasis.cpp */
int  HOCylindricalBasis (Basis *);
void HOSphericalExpansion (Basis *, SPEnergy *, Operator *, HFInterface *);

/*** HOpolynomial.cpp */
void HOPolynomials (const int, const int, Basis *, LaguerrePolynomial *, HermitePolynomial *);

/*** HOhamiltonian.cpp */
void HOHamiltonian (const int, Basis *, OneBodyPotential *, Operator *, LaguerrePolynomial *, HermitePolynomial *);

/*** HOsorteigenvalue.cpp */
void HOSortEigenvalue (Basis *, double *, int *, int *, int *);


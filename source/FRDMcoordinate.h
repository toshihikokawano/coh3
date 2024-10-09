static const int NTHETA    =   19;
static const int NPHI      =   37;

static const int NZETA     =   51;
static const int NRHO      =   21;

#define GAUSS_LEGENDRE40
#include "gauss.h"


/**************************************/
/*   3D Polar Coordinate              */
/**************************************/
class PolarCoordinate{
 private:
  bool allocated;
 public:
  double **r;       // 2dim array for radial position (theta, phi)
  double  *t;       // array for theta
  double  *p;       // array for phi
  double  dt;       // theta increment (equi-angle)
  double  dp;       // phi increment
  int     nt;       // number of points in theta-direction
  int     np;       // number of points in phi-direction

  PolarCoordinate(){
    dt = 180.0/(NTHETA-1);
    dp = 360.0/(NPHI-1);
    nt = NTHETA;
    np = NPHI;
    allocated = false;

    memalloc(nt,np);
  }

  PolarCoordinate(int c1, int c2){
    nt = c1+1;  if(nt > NTHETA) nt = NTHETA;
    np = c2+1;  if(np > NPHI  ) np = NPHI;
    dt = 180.0 / nt;
    dp = 360.0 / np;
    allocated = false;

    memalloc(nt,np);
  }

  void memalloc(int n, int m){
    if(!allocated){
      nt = n;
      np = m;
      r = new double * [nt];
      for(int i=0 ; i<nt ; i++) r[i] = new double [np];
      t = new double [nt];
      p = new double [np];
      for(int i=0 ; i<nt ; i++){ t[i] = dt*i/180.0*PI;  }
      for(int j=0 ; j<np ; j++){ p[j] = dp*j/180.0*PI;  }
      for(int i=0 ; i<nt ; i++){
        for(int j=0 ; j<np ; j++) r[i][j] = 0.0;
      }
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      for(int i=0 ; i<nt ; i++) delete [] r[i];
      delete [] r;
      delete [] t;
      delete [] p;
      allocated = false;
    }
  }

  double getX(int i, int j){ return( sin(t[i]) * cos(p[j]) ); }
  double getY(int i, int j){ return( sin(t[i]) * sin(p[j]) ); }
  double getZ(int i       ){ return( cos(t[i]) ); }
  double getR(int i, int j){ return( r[i][j] ); }
};



/**************************************/
/*   Gauss-Legendre Integration Point */
/**************************************/
class PolarGLPoint{
 private:
  bool allocated;
  double  *gx;      // Gauss-Legendre points in [0,1]
  double  *ga;      // Gauss-Legendre weight
  double    f;      // scaling factor
 public:
  double   *t;      // theta
  double   *p;      // phi
  double   *w;      // Gauss weight at each point
  double  **r;      // 2-dim array for the radial points
  double  **c;      // curvature at r(t,p)
  double ***v;      // normal vector at r(t,p)
  double   st;      // scaling factor for Gauss-Legendre integ in theta diection
  double   sp;      // scaling facror in the phi direction
  int      ng;      // Gauss-Legendre order

  PolarGLPoint(){
    allocated = false;
    memalloc(MAX_GAUSSLEG);
    f  = 1.0;
    st = GLFactor(0.0,PI,t);
    sp = GLFactor(0.0,PI,p);
  }

  ~PolarGLPoint(){
    memfree();
  }

  void memalloc(int mgauss){
    if(!allocated){
      ng = 2*mgauss;
      gx = gaussleg_x;
      ga = gaussleg_a;

      w = new double [ng];
      t = new double [ng];
      p = new double [ng];
      for(int i=0 ; i<ng/2 ; i++) w[i] = w[ng-i-1] = ga[i];

      r = new double * [ng];
      c = new double * [ng];
      v = new double ** [ng];
      for(int i=0 ; i<ng ; i++){
        r[i] = new double [ng];
        c[i] = new double [ng];
        v[i] = new double * [ng];
        for(int j=0 ; j<ng ; j++){
          r[i][j] = 0.0;
          c[i][j] = 0.0;
          v[i][j] = new double [3];
          for(int k=0 ; k<3 ; k++) v[i][j][k] = 0.0;
        }
      }
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      for(int i=0 ; i<ng ; i++){
        for(int j=0 ; j<ng ; j++){
          delete [] v[i][j];
        }
        delete [] v[i];
        delete [] c[i];
        delete [] r[i];
      }

      delete [] v;
      delete [] c;
      delete [] r;
      delete [] t;
      delete [] p;
      delete [] w;
      allocated = false;
    }
  }

  void setNormalization(double x){ f = x; }
  double getNormalization(){ return(f); }

  double getX(int i, int j){ return( sin(t[i]) * cos(p[j]) ); }
  double getY(int i, int j){ return( sin(t[i]) * sin(p[j]) ); }
  double getZ(int i       ){ return( cos(t[i]) ); }
  double getR(int i, int j){ return( r[i][j] ); }

  double normalVector(int i, int j, double *vec){
    double y = 1.0;
    if(j < 0){
      j = -j;
      y = -1.0;
    }
    double x = 0.0;
    int k = 0;
    vec[k] = v[i][j][k]    ;  x += vec[k]*vec[k]; k++;
    vec[k] = v[i][j][k] * y;  x += vec[k]*vec[k]; k++;
    vec[k] = v[i][j][k]    ;  x += vec[k]*vec[k];
    x = (x > 0.0) ? sqrt(x) : 0.0;
    return(x);
  }

  double cartesianVector(int i, int j, double *vec){
    double y = 1.0;
    if(j < 0){
      j = -j;
      y = -1.0;
    }
    double x  = getR(i,j);
    vec[0] = x * getX(i,j);
    vec[1] = x * getY(i,j) * y;
    vec[2] = x * getZ(i);
    return(x);
  }

  double cartesianUnitVector(int i, int j, double *vec){
    double x = cartesianVector(i,j,vec);
    if(x > 0.0){
      for(int i=0 ; i<3 ; i++) vec[i] /= x;
    }
    return(x);
  }

  double inline GLFactor(double xmin, double xmax, double *x){
    double x0 = (xmax + xmin)*0.5;
    double d  = (xmax - xmin)*0.5;
    for(int i=0 ; i<ng/2 ; i++){
      x[i]      = x0 - gx[i]*d;
      x[ng-i-1] = x0 + gx[i]*d;
    }
    return(d);
  }

  double getScaleFactor(){  return (2.0*st*sp);  }
};


/**************************************/
/*   2D Cylindrical Coordinate        */
/**************************************/
class CylindricalCoordinate{
 private:
  bool allocated;
 public:
  double  *r;       // array for R
  double  *t;       // array for theta
  double  *z;       // array for Z
  double  dz;       // z increment
  double  dt;       // theta increment
  int     nz;       // number of z position
  int     nt;       // number of theta

  CylindricalCoordinate(){
    dt = 360.0/(NTHETA-1);
    dz = 2.0/(NZETA-1);
    nt = NTHETA;
    nz = NZETA;
    allocated = false;

    memalloc(nz,nt);
  }

  CylindricalCoordinate(int c1, int c2){
    nt = c1+1;  if(nt > NTHETA) nt = NTHETA;
    nz = c2+1;  if(nz > NZETA ) nz = NZETA;
    dt = 360.0 / (nt-1);
    dz =   2.0 / (nz-1);
    allocated = false;

    memalloc(nt,nz);
  }

  ~CylindricalCoordinate(){
    memfree();
  }

  void memalloc(int n, int m){
    if(!allocated){
      nt = n;
      nz = m;
      t = new double [nt];
      r = new double [nz];
      z = new double [nz];
      for(int i=0 ; i<nt ; i++){ t[i] = dt*i/180.0*PI; }
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] r;
      delete [] z;
      allocated = false;
    }
  }

  void setZ(double z0, double z1){
    dz = (z1 - z0)/(nz - 1);
    for(int i=0 ; i<nz ; i++){ z[i] = z0 + dz*i; }
  }

  double getX(int i){ return( sin(t[i]) ); }
  double getY(int i){ return( cos(t[i]) ); }
  double getZ(       int j){ return( z[j] ); }
  double getR(       int j){ return( r[j] ); }
};



/**************************************/
/*   2D Cylindrical Gauss-Legendre    */
/**************************************/
class CylindricalGLPoint{
 private:
  bool allocated;
  double  *gx;      // Gauss-Legendre points in [0,1]
  double  *ga;      // Gauss-Legendre weight
  double    f;      // scaling factor
 public:
  double   *w;
  double  *z0;      // z-axis on the left-hand side
  double  *z1;      // z-axis on the right-hand side
  double  *r0;      // r-coordinate
  double  *r1;      // 
  double *dr0;      // derivative
  double *dr1;      //
  double  *th;      // theta
  double  *r;       // radius
  double  sz0;      // Gauss-Legendre scaling factor
  double  sz1;      // 
  double  sr;
  double  sth;      
  double zmin;      // left side end
  double zmax;      // right side end
  double zmid;      // mid point
  int      ng;      // order of Gauss-Legendre

  CylindricalGLPoint(){
    allocated = false;
    memalloc(MAX_GAUSSLEG);
    f = 1.0;
    setTArray();
  }

  ~CylindricalGLPoint(){
    memfree();
  }

  void memalloc(int mgauss){
    if(!allocated){
      ng = 2*mgauss;
      gx = gaussleg_x;
      ga = gaussleg_a;

      w  = new double [ng];
      z0 = new double [ng];
      z1 = new double [ng];
      r0 = new double [ng];
      r1 = new double [ng];
      dr0 = new double [ng];
      dr1 = new double [ng];
      r  = new double [ng];
      th = r;
      for(int i=0 ; i<ng/2 ; i++) w[i] = w[ng-i-1] = ga[i];

      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] w;
      delete [] z0;
      delete [] z1;
      delete [] r0;
      delete [] r1;
      delete [] dr0;
      delete [] dr1;
      delete [] r;
      allocated = false;
    }
  }

  void setNormalization(double x){ f = x; }
  double getNormalization(){ return(f); }

  void setZSize(double x0, double x1, double x2){
    zmin = x0;
    zmid = x1;
    zmax = x2;
    sz0 = GLFactor(zmin,zmid,z0);
    sz1 = GLFactor(zmid,zmax,z1);
  }

  void setRArray(double x0){
    sr = GLFactor(0.0,x0,r);
  }

  void setTArray(){
    sth = GLFactor(0.0,PI,th);
  }

  double getZ0(int i){ return( z0[i] ); }
  double getZ1(int i){ return( z1[i] ); }
  double getR0(int i){ return( r0[i] ); }
  double getR1(int i){ return( r1[i] ); }

  double inline GLFactor(double xmin, double xmax, double *x){
    double x0 = (xmax + xmin)*0.5;
    double d  = (xmax - xmin)*0.5;
    for(int i=0 ; i<ng/2 ; i++){
      x[i]      = x0 - gx[i]*d;
      x[ng-i-1] = x0 + gx[i]*d;
    }
    return(d);
  }
};


/**************************************/
/*   Quadratic Surface Parameters     */
/**************************************/
class QSurface{
 private:
  double unit;      // unit conversion factor
  double sigma[3];  // sigma parmaeters
  double alpha[3];  // alpha parameters
  double sigmat;    // max sigma
  double R0;        // size scaling radius
  
 public:
  double center[3];     // center of parabolae
  double junction[2];   // connection points
  double q2[3];         // a^2 / c^2
  double a2[3];         // r-direction size
  double c2[3];         // z-direction size
  double xalpha;
  double xbeta;
  double ztop;          // right end position
  double zbot;          // left end position
  double neck;          // neck position
  double dmin;          // possible minimum neck size
  double dmax;          // possible maximum neck size
  double s2max;         // maximum sigma2
  double dpos;          // mid point
  double wigner;        // parameter for Wigner term

  QSurface(){
    R0 = 0.0;
    init();
  }

  void init(){
    sigmat = 0.0;
    dmin   = 0.0;
    dmax   = 0.0;
    s2max  = 0.0;
    ztop   = 0.0;
    zbot   = 0.0;
    dpos   = 0.0;
    wigner = 1.0;

    for(int i=0 ; i<3 ; i++){
      sigma[i]  = 0.0;
      alpha[i]  = 0.0;
    }

    reset();
  }

  void reset(){
    unit = 1.0;
    alpha[0] = 0.0;

    neck   = 0.0;
    xalpha = 0.0;
    xbeta  = 0.0;

    for(int i=0 ; i<2 ; i++){
      junction[i] = 0.0;
    }
    for(int i=0 ; i<3 ; i++){
      center[i] = 0.0;
      q2[i] = 0.0;
      a2[i] = 0.0;
      c2[i] = 0.0;
    }
  }

  void setUnit(double x){ unit = x; }
  void setR0  (double x){ R0   = x; }  

  void setParm(double x1, double x2, double x3, double y1, double y2, double y3){
    sigma[0] = x1;    sigma[1] = x2;    sigma[2] = x3;
    alpha[0] = y1;    alpha[1] = y2;    alpha[2] = y3;
  }

  void setSigma(int i, double x){ if((1 <= i) && (i <= 3)) sigma[i-1] = x; }
  void setSigmaT(double x){ sigmat = x; }
  void setAlpha(int i, double x){ if((1 <= i) && (i <= 3)) alpha[i-1] = x; }

  double getUnit(){ return unit; }
  double getR0(){ return R0; }

  double getSigma(int i){
    if( (1 <= i) && (i <= 3) ) return(sigma[i-1]);
    return(0.0);
  }
  double getSigmaT(){ return sigmat; }
  double getAlpha(int i){
    if( (1 <= i) && (i <= 3) ) return(alpha[i-1]);
    return(0.0);
  }

  void neckmin(){ dmin = sqrt(a2[2]); }
  void neckmax(){ dmax = sqrt(a2[2]); }

};

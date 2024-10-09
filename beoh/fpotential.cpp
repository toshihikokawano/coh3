/******************************************************************************/
/*  fpotential.cpp                                                            */
/*        fission penetration factor by directly solving Schroedinger Eq.     */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "levden.h"
#include "nucleus.h"
#include "terminate.h"


static int    fpotPotentialEnergy (const int, const double, const double, std::complex<double> *, Fission *);
static double fpotInternalFunction (const int, const double, const double, const double, const double, std::complex<double> *);
static void   fpotPrintNormalized (const int, const double, const double, const double, const double, std::complex<double> *);
static double fpotLevelCompression(double, double, const double);
static int    fpotPotentialExcited (const int, const double, const double, const double, std::complex<double > *, std::complex<double> *);

static inline std::complex<double>
lagrange(double h, std::complex<double> a,std::complex<double> b, std::complex<double> c,
                   std::complex<double> e,std::complex<double> f, std::complex<double> g)
{ return( ((g-a)/60.0+0.15*(b-f)+0.75*(e-c))/h ); }

static bool complexpotential = false;
static std::complex<double> cnrm(0.0,0.0);

static const bool printwfunc = false;
static const bool printtrans = false;
static const bool printpot   = false;

/**********************************************************/
/*      One-Dimensional Potential Energy                  */
/**********************************************************/
double  specFissionPotentialModel(const int k0, const int p0, const int j0, const double ex, double *ch, Nucleus *n)
{
  Fission *fb = n->fission;

  /*** mass parameter */
  double mu = 0.054 * pow((double)n->za.getA(),5.0/3.0);

  /*** wave number */
  double w0 = sqrt(2.0 * mu * ex);

  /*** level compression factor */
  double compfact1p = fb->fisenhance.compfact1p;
  double compfact2p = fb->fisenhance.compfact2p;
  double compfact1n = fb->fisenhance.compfact1n;
  double compfact2n = fb->fisenhance.compfact2n;

  if(compfact1n == 0.0){
    compfact1n = compfact1p;
    compfact2n = compfact2p;
  }


  /*** integration over beta up to 3.0, at each 0.2/k */
  double dx = 0.2 / w0;
  int    nx = 3.0 / dx;

  std::complex<double> *vpot0 = new std::complex<double> [nx]; // potential for ground state
  std::complex<double> *vpot1 = new std::complex<double> [nx]; // potential for excited state
  for(int i=0 ; i<nx ; i++) vpot0[i] = std::complex<double>(0.0,0.0);

  /*** connected palaboras */
  int nm = fpotPotentialEnergy(nx,mu,dx,vpot0,fb);

  if(printpot && (k0 == 0) && (p0 > 0) && (j0 <= 1)){
    for(int i=0 ; i<nm+7 ; i++){

    std::complex<double>v = (i >= nm) ? std::complex<double>(0.0,0.0) : vpot0[i];

    std::cout << std::setw(12) << i*dx;
    std::cout << std::setw(12) << v.real();
    std::cout << std::setw(12) << v.imag() << std::endl;
    }
  }

  double tj = 0.0;
  if(nm > 0){
    /*** for discrete transitions */
    ch[0] = 0.0;
    for(int k=0 ; k<n->ndisc ; k++){
      if( ((int)(2.0*n->lev[k].spin) != j0) || (n->lev[k].parity != p0) ) continue;

      /*** solve Schroedinger equation */
      if(p0 > 0) nm = fpotPotentialExcited(nx,compfact1p,compfact2p,n->lev[k].energy,vpot0,vpot1);
      else       nm = fpotPotentialExcited(nx,compfact1n,compfact2n,n->lev[k].energy,vpot0,vpot1);

      double trn = fpotInternalFunction(nm,ex,mu,dx,w0,vpot1);

      tj += trn;
      ch[0] += 1.0;
    }

    /*** for transitions through continuum */
    int idx = (j0-(int)(2.0*halfint(n->lev[0].spin)))/2; // index for spin J
    for(int k=k0 ; k<n->ncont ; k++){

      /*** number of states */
      double rho = ((p0 > 0) ? n->density[k][idx].even 
                             : n->density[k][idx].odd ) * n->de;

      /*** solve Schroedinger equation */ 
      if(p0 > 0) nm = fpotPotentialExcited(nx,compfact1p,compfact2p,n->excitation[k],vpot0,vpot1);
      else       nm = fpotPotentialExcited(nx,compfact1n,compfact2n,n->excitation[k],vpot0,vpot1);
      double trn = fpotInternalFunction(nm,ex,mu,dx,w0,vpot1);

      tj += trn * rho;
      ch[0] += rho;
    }
  }

  /*** print wave function for a few cases */
  if((printwfunc || printtrans) && (k0 == 0) && (p0 > 0) && (j0 <= 1)){

    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout << std::setprecision(4);

    const int ndiv = 0;
    double ediv = ex / (ndiv + 1);
    for(int n=0 ; n<=ndiv ; n++){
      double elev = ediv * n;

      nm = fpotPotentialExcited(nx,0.0,0.0,elev,vpot0,vpot1);
      double tr = fpotInternalFunction(nm,ex,mu,dx,w0,vpot1);

      double t1 = 1/(1+exp(PI2*(fb->barrier[0].height-ex-elev)/fb->barrier[0].curvature));
      double t2 = 1/(1+exp(PI2*(fb->barrier[1].height-ex-elev)/fb->barrier[1].curvature));
      if(printtrans){
        std::cout << std::setw(12) << ex;
        std::cout << std::setw(12) << tr;
        std::cout << std::setw(12) << norm(cnrm);
        std::cout << std::setw(12) << t1 << std::setw(12) << t2 << std::setw(12) << t1*t2/(t1+t2) << std::endl;
      }

      if(printwfunc) fpotPrintNormalized(nm,ex,mu,dx,w0,vpot1);
    }
  }

  delete [] vpot0;
  delete [] vpot1;

  return(tj);
}


/**********************************************************/
/*      One-Dimensional Potential Energy                  */
/**********************************************************/
int     fpotPotentialEnergy(const int nx, const double mu, const double dx, std::complex<double> *vpot, Fission *fb)
{
  const double beta0 = 0.05;
  double xbeta[5], xjunc[4], c[5], vhigh[5], homeg[5];

  /*** if triple-humped, there will be five regions, otherwise three */
  int nregion = 5;
  if(fb->barrier[2].height == 0.0) nregion = 3;

  vhigh[0] = fb->barrier[0].height;
  homeg[0] = fb->barrier[0].curvature;

  vhigh[1] = fb->potwell[0].height;
  homeg[1] = fb->potwell[0].curvature;

  vhigh[2] = fb->barrier[1].height;
  homeg[2] = fb->barrier[1].curvature;

  if(nregion == 5){
    vhigh[3] = fb->potwell[1].height;
    homeg[3] = fb->potwell[1].curvature;

    vhigh[4] = fb->barrier[2].height;
    homeg[4] = fb->barrier[2].curvature;
  }
  else{
    vhigh[3] = vhigh[4] = homeg[3] = homeg[4] = 0.0;
  }
  
  for(int i=0 ; i<5 ; i++) c[i] = mu * homeg[i] * homeg[i];

  /*** check if potential is std::complex */
  if( (fb->potwell[0].absorb != 0.0) || (fb->potwell[1].absorb != 0.0) ) complexpotential = true;


  /*** xbeta is the center of each curves
       first point, set zero at beta = 0, and shift by beta0 */
  xbeta[0] = sqrt(2.0 * vhigh[0] / c[0]) + beta0;
  for(int i=1 ; i<nregion ; i++){
    double dv = fabs(vhigh[i-1] - vhigh[i]);
    xbeta[i] = xbeta[i-1] + sqrt(2.0 * dv * (c[i-1] + c[i])/(c[i-1] * c[i]));
  }

  /*** calculate conjunction points */
  for(int i=0 ; i<4 ; i++) xjunc[i] = 0.0;
  for(int i=0 ; i<nregion-1 ; i++){
    xjunc[i] = (c[i]*xbeta[i] + c[i+1]*xbeta[i+1])/(c[i] + c[i+1]);
  }


  /*** check the connection point order */
  bool connect = false;
  if(   (xbeta[0] < xjunc[0]) && (xjunc[0] < xbeta[1])
     && (xbeta[1] < xjunc[1]) && (xjunc[1] < xbeta[2]) ) connect = true;

  if(connect && (nregion == 5)){
    if(   (xbeta[2] < xjunc[2]) && (xjunc[2] < xbeta[3])
       && (xbeta[3] < xjunc[3]) && (xjunc[3] < xbeta[4]) ) connect = true;
  }
  if(!connect){
    message << "parabolas cannot be connected : ";
    message << " " << xbeta[0] <<" "<< xjunc[0] << " " << xbeta[1] << " " << xjunc[1];
    message << " " << xbeta[2] <<" "<< xjunc[2] << " " << xbeta[3] << " " << xjunc[3];
    message << " " << xbeta[4];
    cohTerminateCode("fpotPotentialEnergy");
  }

  /*** initialize */
  for(int i=0 ; i<nx ; i++) vpot[i] = std::complex<double>(0.0,0.0);


  /*** shift potential, and set zero below i0 */
  int i0 = 0;
  for(int i=i0 ; i<nx ; i++){
    double vx = i * dx;
    double yr = vhigh[0] - 0.5*c[0]*(vx - xbeta[0])*(vx - xbeta[0]);
    double yi = 0.0;
    vpot[i] = std::complex<double>(yr,yi);
    if(vpot[i].real() < 0.0){
      vpot[i] = std::complex<double>(0.0,0.0);
    }else{
      i0 = i;
      break;
    }
  }


  /*** calculate potential energy */
  for(int i=i0 ; i<nx ; i++){
    double vx = i * dx;

    int j = 0;
    for(j=0 ; j<nregion-1 ; j++) if(vx < xjunc[j]) break;

    double yr = vhigh[j] + ((j%2==0) ? -1.0 : 1.0)*0.5*c[j]*(vx - xbeta[j])*(vx - xbeta[j]);
    double yi = 0.0;

    /*** class-II and class-III well, if imaginary part is given */
    if( (j == 1) || (j == 3) ){
      double w = (j == 1) ? fb->potwell[0].absorb : fb->potwell[1].absorb;
      yi = yr - vhigh[j] - w;
      if(yi > 0.0) yi = 0.0;
    }

    vpot[i] = std::complex<double>(yr,yi);
  }

  /*** distortion test */
//  for(int i=i0 ; i<nx ; i++){
//    double vx = i * dx;
//    std::cout << vx <<" " << vpot[i].real() << " ";
//    vpot[i]  -= std::complex<double>(vx/2.0,0.0);
//    std::cout << vpot[i].real() << std::endl;
//  }


  /*** find zero */
  int nm = 0;
  for(int i=i0 ; i<nx ; i++){
    if(vpot[i].real() < 0.0){ nm = i; break; }
  }

  return(nm);
}


/**********************************************************/
/*      Solve One-Dimensional Potential Penetration       */
/**********************************************************/
double  fpotInternalFunction(const int nm, const double ex, const double mu, const double dx, const double wave, std::complex<double> *vpot)
{
  if(nm == 0) return 0.0;

  double h2 = dx*dx/12.0;
  double c  = 2.0*mu;

  /*** wave function integration using Fox-Goodwin method */
  std::complex<double> p0, p1, p2, q0(0.0,0.0), q1(0.0,0.0), q2;
  p0 = std::complex<double>(cos(-dx*wave),-sin(-dx*wave));
  p1 = std::complex<double>(cos(0.0),sin(0.0));

  q0 = std::complex<double>(-ex*c,0.0);
  q1 = std::complex<double>(-ex*c,0.0);

  for(int i=0 ; i<nm ; i++){
    q2 = (vpot[i] - ex)*c;
    p2 = ((2.0 + 10.0*h2*q1)*p1 - (1.0 - h2*q0)*p0) / (1.0 - h2*q2);

    q0 = q1;  q1 = q2;
    p0 = p1;  p1 = p2;
  }

  /*** extra points outside the potential for wave function matching
       using 7-point Lagrange derivative */
  std::complex<double> px[7];
  for(int i=nm ; i<nm+7 ; i++){
    q2 = - ex * c;
    p2 = ((2.0 + 10.0*h2*q1)*p1 - (1.0 - h2*q0)*p0) / (1.0 - h2*q2);

    q0 = q1;  q1 = q2;
    p0 = p1;  p1 = p2;
    px[i-nm] = p2;
  }
  std::complex<double> f = lagrange(dx,px[0],px[1],px[2],px[4],px[5],px[6]) / px[3];

  /*** matching point, and free particle wave functions */
  double xm = (nm+3) * dx * wave;

  std::complex<double> u1(cos(xm), sin(xm));
  std::complex<double> u2(cos(xm),-sin(xm));
  std::complex<double> du1(-wave*sin(xm), wave*cos(xm));
  std::complex<double> du2(-wave*sin(xm),-wave*cos(xm));

  /*** S-matrix */
  std::complex<double> smat = (f*u2 - du2)/(f*u1 - du1);

  /*** normalization factor */
  cnrm = (u2 - smat * u1) / px[3];

  double s  = norm(smat);
  double t1 = (fabs(1.0 - s) < 1e-12) ? 0.0 : 1.0 - s;
  double t2 = norm(cnrm);

  double t = (complexpotential) ? t2 : t1;

  return(t);
}


/**********************************************************/
/*      Discrete Level Compression Parameter              */
/**********************************************************/
double  fpotLevelCompression(double f0, double f1, const double el)
{
  double c = 1.0;

  if((f0 != 0.0) && (f1 == 0.0)) c = f0;
  else{
    if( (f0 == 0.0) && (f1 == 0.0) ){
      f0 = 0.8;
      f1 = 0.2;
    }
    c = f0 + (1 - exp(-f1 * el)) * (1 - f0);
  }

  return c;
}


/**********************************************************/
/*      Potential for Excited State                       */
/**********************************************************/
int     fpotPotentialExcited(const int nx, const double f0, const double f1, const double el, std::complex<double > *v0, std::complex<double> *v1)
{
  /*** clear array */
  for(int i=0 ; i<nx ; i++) v1[i] = std::complex<double>(0.0, 0.0);

  /*** level compression factor */
  double c = fpotLevelCompression(f0,f1,el);

  /*** initial ofset */
  int i0 = 0;
  for(int i=i0 ; i<nx ; i++){
    if(v0[i].real() > 0.0){ i0 = i; break; }
  }

  /*** calculate g.s. potential + excitation energy */
  int nm = 0;
  for(int i=i0 ; i<nx ; i++){
    v1[i] = v0[i] + std::complex<double>(c * el, 0.0);

    if(v1[i].real() < 0.0){
      v1[i] = std::complex<double>(0.0,0.0);
      nm = i;
      break;
    }
  }

  if(nm == 0) return(nm);

  /*** add extra points for matching calculation */
  nm += 20;
  if(nm > nx-1) nm = nx-1;

  return(nm);
}


/**********************************************************/
/*      Normalized Wave Function                          */
/**********************************************************/
void    fpotPrintNormalized(const int nm,const  double ex, const double mu, const double dx, const double wave, std::complex<double> *vpot)
{
  double h2 = dx*dx/12.0;
  double c  = 2.0*mu;

  /*** wave function integration using Fox-Goodwin method */
  std::complex<double> p0, p1, p2, q0(0.0,0.0), q1(0.0,0.0), q2, px, v;
  p0 = std::complex<double>(cos(-dx*wave),-sin(-dx*wave));
  p1 = std::complex<double>(cos(0.0),sin(0.0));

  for(int i=0 ; i<nm+7 ; i++){

    v = (i >= nm) ? std::complex<double>(0.0,0.0) : vpot[i];

    q2 = (v - ex)*c;
    p2 = ((2.0 + 10.0*h2*q1)*p1 - (1.0 - h2*q0)*p0) / (1.0 - h2*q2);

    q0 = q1;  q1 = q2;
    p0 = p1;  p1 = p2;

    px = p2 * cnrm;

    std::cout << std::setw(12) << i*dx;
    std::cout << std::setw(12) << v.real()  << std::setw(12) << v.imag();
    std::cout << std::setw(12) << px.real() << std::setw(12) << px.imag();
    std::cout << std::setw(12) << norm(px)  << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;
}


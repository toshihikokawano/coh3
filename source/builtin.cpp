/******************************************************************************/
/*  builtin.cpp                                                               */
/*        built in global optical model potentials                            */
/******************************************************************************/

#include <iostream>
#include <cstring>

#include "omplib/omplib.h"
#include "etc.h"

#define GAUSS_LEGENDRE20
#include "gauss.h"

static int findname (std::string *, char *, const int);
static unsigned int add_library (const int, const double, const int, const int, const int, const int, Optical *);

static const int N_OMPLIB = 23;
static std::string lib_name[N_OMPLIB] =
  {"Wilmore", "Becchetti", "Rapaport", "Walter", "CH89", 
   "Koning", "Soukhovitskii", "Kunieda",
   "Perey", "Schwandt", "Madland",
   "Lemos", "Nolte", "Avrigeanu", "Avrigeanu2009", "Avrigeanu2014", "TALYS-A",
   "Bojowald", "An", "Han",
   "bound", "test", "spline"};

static const int ADD_CUSTOM_OMP = 18;
static std::string add_name[ADD_CUSTOM_OMP] =
  {"scratch", "userdef",
   "SmithA120","SmithA136","SmithA415",
   "Flap22","Young-Am","Young-Pu","Young-Re","modSoukhovitskii","Soukhovitskii2005",
   "Dave1p","WLH1","WLH2","WLH3","WLH4","WLH5"};


/***********************************************************/
/*      Search for Given OMP Names in the Library          */
/***********************************************************/
unsigned int find_omp(std::string str)
{
   unsigned int k = findname(lib_name,&str[0],N_OMPLIB);
   if(k == 0){
        if( (k = findname(add_name,&str[0],ADD_CUSTOM_OMP)) != 0 ) k += N_OMPLIB;
   }
   if(k == 0) std::cerr << "ERROR     : OMP name [ " << str << " ] not found" << std::endl;
   return(k);
}


/***********************************************************/
/*      Compare OMP Names Provided With Library            */
/***********************************************************/
int findname(std::string *list, char *name, const int n)
{
  int  i;

  for(i=0 ; i<n ; i++){
    if(!strcmp(name,list[i].c_str())) break;
  }
  if(i == n)  return(0);
  else        return(i+1);
}


/***********************************************************/
/* Potential form bit field code description               */
/*   bit: 76543210  :  0 Imag  Volume  (WS)                */
/*                     1       Surface (deriv. WS)         */
/*                     2       Surface (gaussian)          */
/*                     3       undef.                      */
/*                     4 Real  Volume  (WS)                */
/*                     5       Surface (deriv. WS)         */
/*                     6       Surface (gaussian)          */
/*                     7       undef.                      */
/*   example                                               */
/*        00010001 = 0x0011 : Real WS + Imag WS            */
/*        00010010 = 0x0012 : Real WS + Imag dWS           */
/*        00010011 = 0x0013 : Real WS + Imag WS +dWS       */
/***********************************************************/
/***********************************************************/
/*      Select Optical Potential Parameter                 */
/***********************************************************/
unsigned int omp_library(const int id, const int zt, const int at, const int zi, const int ai, const double e, Optical *omp)
{
  unsigned int potfm;
  switch(id){
  case  1: potfm = OMPWilmoreHodgson(      at,   ai,zi,omp); break;
  case  2: potfm = OMPBecchettiGreenlees(  at,zt,ai,zi,omp); break;
  case  3: potfm = OMPRapaport(          e,at,zt,ai,zi,omp); break;
  case  4: potfm = OMPWalterGuss(        e,at,zt,ai,zi,omp); break;
  case  5: potfm = OMPChapelHill89(      e,at,zt,ai,zi,omp); break;
  case  6: potfm = OMPKoningDelaroche(   e,at,zt,ai,zi,omp); break;
  case  7: potfm = OMPSoukhovitskii(     e,at,zt,   zi,omp); break;
  case  8: potfm = OMPKunieda(           e,at,zt,ai,zi,omp); break;
  case  9: potfm = OMPPerey(               at,zt,ai,zi,omp); break;
  case 10: potfm = OMPMadlandSchwandt(   e,at,zt,ai,zi,omp); break;
  case 11: potfm = OMPMadlandYoung(      e,at,zt,      omp); break;
  case 12: potfm = OMPLemos(                     ai,zi,omp); break;
  case 13: potfm = OMPNolte(               at,zt,ai,zi,omp); break;
  case 14: potfm = OMPAvrigeanu(         e,at,zt,ai,zi,omp); break;
  case 15: potfm = OMPAvrigeanu2009(     e,at,zt,ai,zi,omp); break;
  case 16: potfm = OMPAvrigeanu2014(     e,at,zt,ai,zi,omp); break;
  case 17: potfm = OMPTALYS_alpha(       e,at,zt,ai,zi,omp); break;
  case 18: potfm = OMPBojowald(          e,at,zt,ai,zi,omp); break;
  case 19: potfm = OMPAnHaixia(            at,zt,ai,zi,omp); break;
  case 20: potfm = OMPHanYinlu(            at,zt,ai,zi,omp); break;
  case 21: potfm = OMPBoundState(                   zi,omp); break;
  case 22: potfm = OMPtest(                         zi,omp); break;
  case 23: potfm = OMPspline(            e,at,zt,ai,zi,omp); break;
  default: potfm = 0                                       ; break;
  }

  if(id>N_OMPLIB) potfm = add_library(id-N_OMPLIB,e,zt,at,zi,ai,omp);

  return(potfm);
}


unsigned int add_library(const int id, const double e, const int zt, const int at, const int zi, const int ai, Optical *omp)
{
  unsigned int potfm;
  switch(id){
  case  1: potfm = OMPscratch(           e,at,zt,ai,zi,omp); break;
  case  2: potfm = OMPuserdef(           e,at,zt,ai,zi,omp); break;
  case  3: potfm = OMPSmithA120(         e,at,zt,      omp); break;
  case  4: potfm = OMPSmithA136(         e,at,zt,      omp); break;
  case  5: potfm = OMPSmithA415(           at,zt,      omp); break;
  case  6: potfm = OMPFlap22(            e,at,zt,      omp); break;
  case  7: potfm = OMPYoung_Am(          e,at,zt,      omp); break;
  case  8: potfm = OMPYoung_Pu(          e,at,zt,      omp); break;
  case  9: potfm = OMPYoung_Re(          e,at,zt,      omp); break;
  case 10: potfm = OMPModSoukhovitskii(  e,at,zt,   zi,omp); break;
  case 11: potfm = OMPSoukhovitskii2005( e,at,zt,   zi,omp); break;
  case 12: potfm = OMPDave1p(              at,zt,      omp); break;
//  case 13: potfm = OMPWLH1(              e,at,zt,ai,zi,omp); break;
//  case 14: potfm = OMPWLH2(              e,at,zt,ai,zi,omp); break;
//  case 15: potfm = OMPWLH3(              e,at,zt,ai,zi,omp); break;
//  case 16: potfm = OMPWLH4(              e,at,zt,ai,zi,omp); break;
//  case 17: potfm = OMPWLH5(              e,at,zt,ai,zi,omp); break;
  default: potfm = 0; break;
  }

  return(potfm);
}


/*************************************************/
/*    Dispersion Contribution to Surface Term    */
/*    J.M. Quesada et al. PRC 67, 067601 (2003)  */
/*************************************************/
double OMPdeltaSurface(const double ex, const double e0, const double as, const double bs, const double cs)
{
  const int m = 2;

  double ep = ex + e0;
  double em = ex - e0;

  double resp = - ep*ep/(ep*ep + bs*bs);
  double resm =   em*em/(em*em + bs*bs);

  std::complex<double> vs(0.0,0.0), p, z;
  for(int j=1 ; j<=m ; j++){
    p = (j==1) ? std::complex<double>(0.0,bs) : std::complex<double>(0.0,-bs);
    z = ex/m * p*(2.0*p + ep - em)/((p + e0)*(p + ep)*(p - em));
    std::complex<double> z1 = -p*cs;
    vs += z * exp(z1) * expintE1(z1);
  }
  vs = as*(vs - resp * exp( cs*ep) * expint(-cs*ep)
              - resm * exp(-cs*em) * expint( cs*em))/PI;

  return(vs.real());
}


/*************************************************/
/*    Dispersion Contribution to Volume Term     */
/*************************************************/
double OMPdeltaVolume(const double ex, const double e0, const double av, const double bv)
{
  const int m = 2;

  double ep = ex + e0;
  double em = ex - e0;

  double resp = - ep*ep/(ep*ep + bv*bv);
  double resm =   em*em/(em*em + bv*bv);

  std::complex<double> vv(0.0,0.0), p, z;
  for(int j=1 ; j<=m ; j++){
    p = (j==1) ? std::complex<double>(0.0,bv) : std::complex<double>(0.0,-bv);
    z = ex/m * p*(2.0*p + ep - em)/((p + e0)*(p + ep)*(p - em));
    vv += z * log(-p);
  }
  vv = -av*(vv + resp * log(fabs(ep)) + resm * log(fabs(em)))/PI;
  
  return(vv.real());
}


double OMPasymmetricVolume(const double e, const double ef, const double ea, const double av, const double bv)
{
  const int maxrange = 20;
  const double eps = 1.0e-4;

  const double alpha = 1.65;

  double erange = 1.0;
  double estart = ef - ea;

  /*** integration below Ef - Ea */
  double y1 = 0.0;
  for(int n=0 ; n<maxrange ; n++){

    double e1 = estart;
    double e2 = estart - erange;

    double x0 = (e1 + e2) * 0.5;
    double xw = (e1 - e2) * 0.5;

    double dy = 0.0;
    for(int j=0 ; j<MAX_GAUSSLEG ; j++){

      double p1  = x0 - gaussleg_x[j] * xw;
      double p2  = x0 + gaussleg_x[j] * xw;

      double ex1 = (p1 - ef) * (p1 - ef);
      double ex2 = (p2 - ef) * (p2 - ef);

      double ef1 = (ef - p1 - ea) * (ef - p1 - ea);
      double ef2 = (ef - p2 - ea) * (ef - p2 - ea);

      double wv1 = av * ex1 / (ex1 + bv*bv) * (1.0 - ef1 / (ef1 + ea*ea));
      double wv2 = av * ex2 / (ex2 + bv*bv) * (1.0 - ef2 / (ef2 + ea*ea));

      double x1 = 1.0 / (p1 - e) - 1.0 / (p1 - ef);
      double x2 = 1.0 / (p2 - e) - 1.0 / (p2 - ef);

      dy += (x1*wv1 + x2*wv2) * gaussleg_a[j] * xw;
    }

    y1 += dy;
    if(fabs(dy / y1) < eps) break;

    estart = e2;
    erange *= 10.0;
  }

  /*** integration above Ef + Ea */
  erange = 1.0;
  estart = ef + ea;

  double a = pow(ef + ea, 1.5) * 0.5;
  double b = sqrt(ef + ea) * 1.5;

  double y2 = 0.0;
  for(int n=0 ; n<maxrange ; n++){

    double e1 = estart;
    double e2 = estart + erange;

    double x0 = (e1 + e2) * 0.5;
    double xw = (e2 - e1) * 0.5;

    double dy = 0.0;
    for(int j=0 ; j<MAX_GAUSSLEG ; j++){

      double p1  = x0 - gaussleg_x[j] * xw;
      double p2  = x0 + gaussleg_x[j] * xw;

      double ex1 = (p1 - ef) * (p1 - ef);
      double ex2 = (p2 - ef) * (p2 - ef);

      double wv1 = av * ex1 / (ex1 + bv*bv) + alpha * (sqrt(p1) + a/p1 - b);
      double wv2 = av * ex2 / (ex2 + bv*bv) + alpha * (sqrt(p2) + a/p2 - b);

      double x1 = 1.0 / (p1 - e) - 1.0 / (p1 - ef);
      double x2 = 1.0 / (p2 - e) - 1.0 / (p2 - ef);

      dy += (x1*wv1 + x2*wv2) * gaussleg_a[j] * xw;
    }

    y2 += dy;
    if(fabs(dy / y2) < eps) break;

    estart = e2;
    erange *= 10.0;
  }

  return( (y1 + y2)/ PI );
}


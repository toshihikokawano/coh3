/******************************************************************************/
/*  moldauer.cpp                                                              */
/*        width fluctuation correction factor                                 */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "structur.h"
#include "statmodel.h"

#define GAUSS_LEGENDRE20
#include "gauss.h"

//#define MOLDAUER_SYSTEMATICS
//#define  ERNEBJERG_SYSTEMATICS
#define  GOE_SYSTEMATICS
//#define  FIXED_DEGREE_OF_FREEDOM

static double  statMoldauerSystematics  (double, double);

static double transum = 0.0, transin = 0.0, w[2*MAX_GAUSSLEG];
static int    lin=0, j2in=0, idin=0, levin=0;

/**********************************************************/
/*      Moldauer width fluctuation factor                 */
/*          P.A.Moldauer                                  */
/*            Phys. Rev., 123, 986(1961)                  */
/*            Rev. Mod. Phys., 36, 1079(1964)             */
/**********************************************************/
/**********************************************************/
/*      Set width fluctuation parameters                  */
/**********************************************************/
void statWidthFluctuationSet(const int l, const int j2, const int id, const int lev, const double tj)
{
  statWidthFluctuationSetT0(tj);
  lin     = l;                // L for incident channel
  j2in    = j2;               // 2*J for incident channel
  idin    = id;               // index of incident channel
  levin   = lev;              // index of target level
}

void statWidthFluctuationSetT0(const double tj)
{
  transin = tj;               // Tj for incident channel
}


/**********************************************************/
/*      Reset width fluctuation parameters                */
/**********************************************************/
void statWidthFluctuationReset(const double tsum)
{
  transum = tsum;             // sum of transmissions
  transin = 0.0;              // incident channel transmission
  lin     = 0;
  j2in    = 0;
  idin    = 0;
  levin   = 0;
  for(int i=0 ; i<2*MAX_GAUSSLEG ; i++) w[i] = 1.0;
}


/**********************************************************/
/*      Calculate W factor at each Gauss-Legendre point   */
/**********************************************************/
void statWidthFluctuation(const double t, const double den, const double nu)
{
  double nu2 = (nu == 0.0) ? statMoldauerSystematics(t,transum)/2.0 : nu/2.0;
  double x   = 1.0-t/(transum*nu2);

  for(int ig=0 ; ig<MAX_GAUSSLEG ; ig++){
    double x1 = 0.5*(1+gaussleg_x[ig]);
    double x2 = 0.5*(1-gaussleg_x[ig]);
    /*** W at Gauss-Legendre inegration */
    w[MAX_GAUSSLEG+ig  ] *= pow( ((1-x*x1)/x2),-nu2*den );
    w[MAX_GAUSSLEG-ig-1] *= pow( ((1-x*x2)/x1),-nu2*den );
  }
}


void wfcheck(void)
{
  std::cout << std::setw(7) << transum << std::endl;
  for(int i=0 ; i<MAX_GAUSSLEG ; i++) std::cout << " "<< std::setw(7) << w[i];
  std::cout << std::endl;
  for(int i=MAX_GAUSSLEG ; i<2*MAX_GAUSSLEG ; i++) std::cout << " "<< std::setw(7) << w[i];
  std::cout << std::endl;
}

/**********************************************************/
/*      Moldauer's width fluctuation correction factor    */
/**********************************************************/
double  statMoldauer(const int l, const int j2, const int id, const int lev, const double t1, const double df)
{
  double sum, wfc = 1.0;

  if(t1 == 0.0) return(wfc);

  double t0  = transin;
  double df0 = statMoldauerSystematics(t0,transum);
  double df1 = (df == 0.0) ? statMoldauerSystematics(t1,transum) : df;
  double x3  = 2.0/(df0*transum);
  double x4  = 2.0/(df1*transum);

  sum = 0.0;
  for(int ig=0 ; ig<MAX_GAUSSLEG ; ig++){
    double x1 = 0.5*(1+gaussleg_x[ig]);
    double x2 = 0.5*(1-gaussleg_x[ig]);
    double g  = 1.0 - t0*x3;
    double f  = 1.0 - t1*x4;
    sum +=      gaussleg_a[ig]* ( w[MAX_GAUSSLEG+ig  ]/(1-f*x1)/(1-g*x1)
                                 +w[MAX_GAUSSLEG-ig-1]/(1-f*x2)/(1-g*x2));
  }
  wfc = sum/2.0;

  if(l==lin && j2==j2in && id==idin && lev==levin){
    wfc *= 1.0 + 2.0/df0; 
  }

  return(wfc);
}


/**********************************************************/
/*      Elastic Enhancement Factor                        */
/**********************************************************/
double statElasticEnhancement(const double t)
{
  double df = statMoldauerSystematics(t,transum);
  double ef = 1.0 + 2.0/df; 
  return(ef);
}


/**********************************************************/
/*      Width fluctuation for Coupled-Channels            */
/**********************************************************/
double  statMoldauerCC(const double t0, const double t1)
{
  double sum, wfc = 1.0;

  if(t1 == 0.0) return(wfc);

  double df0 = statMoldauerSystematics(t0,transum);
  double df1 = statMoldauerSystematics(t1,transum);
  double x3  = 2.0/(df0*transum);
  double x4  = 2.0/(df1*transum);

  sum = 0.0;
  for(int ig=0 ; ig<MAX_GAUSSLEG ; ig++){
    double x1 = 0.5*(1+gaussleg_x[ig]);
    double x2 = 0.5*(1-gaussleg_x[ig]);
    double g  = 1.0 - t0*x3;
    double f  = 1.0 - t1*x4;
    sum +=      gaussleg_a[ig]* ( w[MAX_GAUSSLEG+ig  ]/(1-f*x1)/(1-g*x1)
                                 +w[MAX_GAUSSLEG-ig-1]/(1-f*x2)/(1-g*x2));
  }
  wfc = sum/2.0;

  return(wfc);
}


/**********************************************************/
/*      Correlation Term                                  */
/**********************************************************/
double  statMoldauerCorrelation(const double t0, const double t1)
{
  double df0 = statMoldauerSystematics(t0,transum);
  double df1 = statMoldauerSystematics(t1,transum);
  double cor = (2.0/df0 - 1.0) * (2.0/df1 - 1.0);

  if(cor > 0.0) cor = sqrt(cor);

  return(cor);
}


/**********************************************************/
/*      P.A.Moldauer                                      */
/*      Nucl. Phys., A344, 185(1980)                      */
/**********************************************************/
#ifdef MOLDAUER_SYSTEMATICS
double statMoldauerSystematics(const double t, const double st)
{
  if(t < 0.0) return(1.0);

  double nu = 1.78+(pow(t,1.212)-0.78)*exp(-0.228*st);
  if(nu>2.0) nu = 2.0;
  return(nu);
}
#endif

/**********************************************************/
/*      M.Ernebjerg and M.Herman                          */
/*      Proc. conf. Nuclear Data, 2004, Santa Fe          */
/**********************************************************/
#ifdef ERNEBJERG_SYSTEMATICS
double statMoldauerSystematics(const double t, const double st)
{
  double f = 0.177/(1-pow(t,20.337));
  double g = 1.0 + 3.148*t*(1-t);
  double nu = 2.0 - 1.0/(1+f*pow(st,g));
  if(nu>2.0) nu = 2.0;
  else if(nu<1.0) nu = 1.0;

  return(nu);
}
#endif


/**********************************************************/
/*      T. Kawano and P. Talou                            */
/*      Proc. conf. Nuclear Data, 2013, New York          */
/**********************************************************/
#ifdef GOE_SYSTEMATICS
double statMoldauerSystematics(const double t, const double st)
{
  double a = 0.0287892 * t + 0.245856;
  double b = 1+t*(1-t)*exp(-2*st)*2.5;
  double c = t*t - (st-2*t)*(st-2*t);
  double d = 1.0;
  if( (c > 0.0) && (st < 2*t) ) d = sqrt(c)/t;
  double nu = 2.0 - 1.0/(1.0 + a*(st+t)/(1-t)*b*d);

  return(nu);
}
#endif


/**********************************************************/
/*      Ad Hoc Fixed Nu Calculation                       */
/**********************************************************/
#ifdef FIXED_DEGREE_OF_FREEDOM
double statMoldauerSystematics(const double t, const double st)
{
  double a = 0.0287892 * t + 0.245856;
  double b = 1+t*(1-t)*exp(-2*st)*2.5;
  double c = t*t - (st-2*t)*(st-2*t);
  double d = 1.0;
  if( (c > 0.0) && (st < 2*t) ) d = sqrt(c)/t;
  double nu = 2.0 - 1.0/(1.0 + a*(st+t)/(1-t)*b*d);

  if(t > 0.1){
    nu = 1.04;
    //    nu = 2.67;
  }

  return(nu);
}
#endif


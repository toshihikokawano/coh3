/******************************************************************************/
/*  beohffragenergy.cpp                                                       */
/*        Divide TXE into excitation energies of light and heavy fragments    */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "structur.h"
#include "kcksyst.h"
#include "levden.h"

static double RTmodel(double, double, LevelDensity *, LevelDensity *, const double, const double);
static double LDmodel(double, double, LevelDensity *, LevelDensity *, const double);
static double RTmodelParameter1(const int);
static double RTmodelParameter2(const int);
static double RTmodelParameter3(const int);
static double RTmodelParameter4(LevelDensity *, LevelDensity *);
static double RTmodelParameter9(const int);

/**********************************************************/
/*      Excitation Energy Sharing Model                   */
/**********************************************************/
double beohExcitationEnergyDivide(ZAnumber *cL, ZAnumber *cH, double rt, const double txe)
{
  if(txe == 0.0) return 0.0;

  LevelDensity ldpL, ldpH;

  double mL = (double)cL->getA();
  double mH = (double)cH->getA();

  /*** read shell correction and pairing energies */
  kckDataRead(cL,&ldpL);
  kckDataRead(cH,&ldpH);

  ldpL.a = kckAsymptoticLevelDensity(mL);
  ldpH.a = kckAsymptoticLevelDensity(mH);

  /*** some built-in RT models */
  if(rt < 0.0){
    int k = (int)fabs(rt);
    switch(k){
    case  1: rt = RTmodelParameter1(cH->getA());  break;
    case  2: rt = RTmodelParameter2(cH->getA());  break;
    case  3: rt = RTmodelParameter3(cH->getA());  break;
    case  4: rt = RTmodelParameter4(&ldpL,&ldpH); break;
    case  9: rt = RTmodelParameter9(cH->getA());  break;
    default: rt = 1.0;                            break;
    }
  }

  double r = 0.0;
  if(rt > 0.0) r = RTmodel(mL,mH,&ldpL,&ldpH,rt,txe);
  else         r = LDmodel(mL,mH,&ldpL,&ldpH,txe);

  return r;
}


/**********************************************************/
/*      RT Damping by Excitation Energy                   */
/**********************************************************/
double beohRTDamping(const double rt, const double en)
{
  const double e = 5.0;
  const double w = 1.0;

  double rt0 = rt  / (1.0 + exp( (en - e)/w ));
  double rt1 = 1.0 / (1.0 + exp(-(en - e)/w ));

  return(rt0 + rt1);
}


/**********************************************************/
/*      Standard RT Model                                 */
/**********************************************************/
double RTmodel(double mL, double mH, LevelDensity *ldpL, LevelDensity *ldpH, const double rt, const double txe)
{
  const double eps = 1.0e-6;
  const int    max_itr = 100;

  /*** level density parameters, a, for each fragment at TXE/2 */
  double eL = 0.5 * txe;
  double eH = 0.5 * txe;

  double aL = ldShellCorrection(eL - ldpL->pairing_energy,ldpL->a,ldpL->shell_correct,mL);
  double aH = ldShellCorrection(eH - ldpH->pairing_energy,ldpH->a,ldpH->shell_correct,mH);

  double p0 = aL/aH * rt*rt;

  /*** iterate till energies converge, this mush be very quick */
  for(int itr=0 ; itr<max_itr ; itr++){

    eH = (txe - (ldpL->pairing_energy + ldpH->pairing_energy))/(1.0 + p0) + ldpH->pairing_energy;
    eL = txe - eH;

    aL = ldShellCorrection(eL - ldpL->pairing_energy,ldpL->a,ldpL->shell_correct,mL);
    aH = ldShellCorrection(eH - ldpH->pairing_energy,ldpH->a,ldpH->shell_correct,mH);

    double p1 = aL/aH * rt*rt;
    double x  = fabs(1.0 - p1/p0); //  cout << itr <<" " << p0 <<" " << p1 << " " << x << endl;

    if(x < eps) break;
    p0 = p1;
  }

  return(eL / (eL + eH));
}


/**********************************************************/
/*      Level Density Sorting Model                       */
/**********************************************************/
double LDmodel(double mL, double mH, LevelDensity *ldpL, LevelDensity *ldpH, const double txe)
{
  const double eps = 1.0e-6;
  const int    max_itr = 100;

  double de    = 0.01;
  double remin = 0.001;
  double remax = 0.999;

  double re, eL, eH;
  for(int itr=0 ; itr < max_itr; itr++){
    re = (remax + remin) * 0.5;
    eL = re         * txe;
    eH = (1.0 - re) * txe;

    double rL = 0.0;
    for(int k=0 ; ; k++){
      double ex = de * k;
      rL += ldLevelDensity(ex,mL,ldpL) * de;
      if(ex > eL) break;
    }

    double rH = 0.0;
    for(int k=0 ; ; k++){
      double ex = de * k;
      rH += ldLevelDensity(ex,mH,ldpH) * de;
      if(ex > eH) break;
    }

    double r = log(rL / rH);

    if(fabs(remin/remax - 1.0) < eps) break;

    if(r < 0.0) remin = re;
    else if(r > 0.0) remax = re;
  }

  return(eL / (eL + eH));
}


/**********************************************************/
/*      RT Parameterization of Manailescu et al.          */
/*      NPA 867, 12 (2011)                                */
/**********************************************************/
double RTmodelParameter1(const int mass)
{
  double rt = 1.0;

  if(mass <= 132){
    rt = (1.6-1.0)/(132.0-118.0) * (mass - 118) + 1.0;
  }
  else if(mass <= 135){
    rt = (1.0-1.6)/(135.0-132.0) * (mass - 132) + 1.6;
  }
  else if(mass <= 142){
    rt = 1.0;
  }
  else{
    rt = (0.6-1.0)/(163.0-142.0) * (mass - 142) + 1.0;
  }

  return rt;
}


/**********************************************************/
/*      RT Parameterization by Litaize for Cf252          */
/*      PRC 82, 054616 (2010)                             */
/**********************************************************/
double RTmodelParameter2(const int mass)
{
  double rt = 1.0;
  if(mass <= 132){
    rt = (1.6-1.0)/(132.0-118.0) * (mass - 118) + 1.0;
  }
  else{
    rt = (0.4-1.6)/(165.0-132.0) * (mass - 132) + 1.6;
  }

  /*** to avoid too different temperature */
  if(rt < 0.8) rt = 0.8;

  return rt;
}


/**********************************************************/
/*      A-dependent RT Parameterization                   */
/**********************************************************/
double RTmodelParameter3(const int mass)
{
  double rt = 1.0;

  static double Pu239_RTAh [51] = { // Adjusted to fit Apalin, 1965
        1.000, 1.130, 1.120, 1.160, 1.140, 1.190, 1.225, 1.240, 1.330, 1.305,
        1.260, 1.110, 1.105, 1.075, 1.020, 1.010, 1.000, 1.080, 1.100, 1.075,
        1.065, 1.045, 1.030, 1.015, 1.005, 0.985, 0.955, 0.935, 0.910, 0.895,
        0.875, 0.850, 0.825, 0.850, 0.830, 0.830, 0.830, 0.830, 0.830, 0.830,
        0.830, 0.830, 0.830, 0.830, 0.830, 0.830, 0.830, 0.830, 0.830, 0.830,
        0.830 };
  static int ah0 = 120;

  if(mass < ah0 || mass >= ah0+51) rt = 1.0;
  else rt = Pu239_RTAh[mass - ah0];

  return rt;
}



/**********************************************************/
/*      Shell Energy-dependent RT Parameterization        */
/**********************************************************/
double RTmodelParameter4(LevelDensity *ldpL, LevelDensity *ldpH)
{
  const double scale = 20.0;

  double sc = ldpL->shell_correct + ldpH->shell_correct;
  double rt = exp(- sc / scale);
//cout <<"RT " << ldpL->shell_correct <<" " << ldpH->shell_correct<< " " << rt << endl;
  return rt;
}


/**********************************************************/
/*      Generic RT Parameterization                       */
/**********************************************************/
double RTmodelParameter9(const int mass)
{
  double rt = 1.0;
  if(mass <= 132){
    rt = 1.6;
  }
  else{
    rt = (0.4-1.6)/(165.0-132.0) * (mass - 132) + 1.6;
  }

  return rt;
}



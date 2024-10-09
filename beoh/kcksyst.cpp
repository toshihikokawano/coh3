/******************************************************************************/
/*  kcksyst.cpp                                                               */
/*        level density parameters based on KTUY mass model                   */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "structur.h"
#include "kcksyst.h"
#include "terminate.h"

static std::string datadir = DATADIR;

/***********************************************************/
/*      Read Level Density Paraters from a File            */
/***********************************************************/
int kckDataRead(ZAnumber *za, LevelDensity *ldp)
{
  std::ifstream  fp;
  std::string    d;
  int            z,a;

  std::string str = KCKDATAFILE;

  /*** try current directry */
  fp.open(&str[0]);

  /*** then system data area */
  if(!fp){
    str = datadir + "/" + str;
    fp.open(&str[0]);
  }

  if(!fp){
    message << "level density parameter file " << str << " not found";
    cohTerminateCode("kckDataRead");
  }

  bool found=false;
  while(getline(fp,str)){
    if(str[0] == '#') continue;

    d=str.substr( 0, 5);    z = atoi(&d[0]);
    d=str.substr( 5, 6);    a = atoi(&d[0]);
    ZAnumber za1(z,a);

    if((za1.getZ()==za->getZ()) && (za1.getA()==za->getA())){

      d=str.substr(11,13);  ldp->pairing_energy = atof(&d[0]);
      d=str.substr(24,13);  ldp->shell_correct  = atof(&d[0]);
      d=str.substr(37,10);  ldp->match_energy   = atof(&d[0]);
      d=str.substr(47,10);  ldp->a              = atof(&d[0]);

      if(ldp->match_energy > 0.0){
        d=str.substr(57,10);  ldp->temperature  = atof(&d[0]);
        d=str.substr(67,10);  ldp->E0           = atof(&d[0]);
      }else{
        d=str.substr(77,10);  ldp->temperature  = atof(&d[0]);
        d=str.substr(87,10);  ldp->E0           = atof(&d[0]);
      }
      found = true;
      break;
    }
  }
  fp.close();


  if(!found){
    message << "level density parameter for Z " << za->getZ() << " - A " << za->getA() << " not found";
    cohTerminateCode("kckDataRead");
  }

  return(0);
}


/***********************************************************/
/*      Systematics for Asymptotic Level Density Parameter */
/***********************************************************/
double kckAsymptoticLevelDensity(const double a)
{
  double astar = 0.126181*a + 7.52191e-05*a*a;
  return(astar);
}


/***********************************************************/
/*      Systematics for Spin-Cutoff Parameter              */
/***********************************************************/
double kckSpinCutoff(const double a)
{
//double sigma2 = 0.0465 *a + 4.296;
  double sigma2 = -0.000303626 * a*a + 0.141231 * a - 0.264688;
//double sigma2 = -0.000210084 * a*a + 0.0968439 * a + 0.938983; //from ENSDF


  if(sigma2 < 3.70) sigma2 = 3.70; // A < 30

  return(sigma2);
}


/***********************************************************/
/*      Systematics for T Parameter                        */
/***********************************************************/
double kckTemperature(const double a, const double dw)
{
  double gamma = 0.1;
  double tsyst = 47.0582 * pow(a,-0.893714) * sqrt(1.0 - gamma*dw);
  return(tsyst);
}


/***********************************************************/
/*      Systematics for Spin-Cutoff Parameter              */
/***********************************************************/
double kckE0(const double a, const double p, const double dw)
{
  double gamma = 0.16;
  double tsyst = kckTemperature(a,dw);
  double esyst = p - gamma*dw + tsyst*(-0.01380*a - 1.379);
  return(esyst);
}

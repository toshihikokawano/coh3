/******************************************************************************/
/*  FRDMreadshape.cpp                                                         */
/*        read nuclear shape parameters (epsilons)                            */
/******************************************************************************/

#include <sstream>
#include <fstream>
#include <iostream>

#include "FRDM.h"
#include "terminate.h"

static std::string datadir = DATADIR;

/***********************************************************/
/*      Read Nuclear Shape Parameters from File            */
/***********************************************************/
void FRDMReadShape(MFTSystem *sys)
{
  std::ifstream  fp;
  std::string    d;
  int       zt,at;
  double    eps[6];

  std::string str = FRDMDATAFILE;

  /*** try current directry */
  fp.open(&str[0]);

  /*** then system data area */
  if(!fp){
    str = datadir + "/" + str;
    fp.open(&str[0]);
  }

  if(!fp){
    message << "FRDM shape parameter file " << str << " not found";
    cohTerminateCode("FRDMReadShape");
  }

  for(int i=0 ; i<6 ; i++) eps[i] = 0.0;

  bool found=false;
  while(getline(fp,str)){
    if(str[0] == '#') continue;

    d=str.substr( 0, 5);  zt = atoi(&d[0]);
    d=str.substr( 5, 5);  at = atoi(&d[0]);

    if(zt == sys->getZ() && at == sys->getA()){
      d=str.substr(10,10);  eps[1] = atof(&d[0]);
      d=str.substr(20,10);  eps[2] = atof(&d[0]);
      d=str.substr(30,10);  eps[3] = atof(&d[0]);
      d=str.substr(40,10);  eps[4] = atof(&d[0]);
      d=str.substr(50,10);  eps[5] = atof(&d[0]);
      found = true;
      break;
    }
  }
  fp.close();


  if(!found){
    message << "FRDM shape parameter for Z " << sys->getZ() << " - A " << sys->getA() << " not found";
    cohTerminateCode("FRDMReadShape");
  }

  for(int i=0 ; i<6 ; i++) sys->setEps(i+1,eps[i]);
}

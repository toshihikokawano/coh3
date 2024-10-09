/******************************************************************************/
/*  beta2.cpp                                                                 */
/*        read beta2 deformation parameters                                   */
/******************************************************************************/

#include <string>
#include <sstream>
#include <ostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "structur.h"
#include "beta2.h"
#include "terminate.h"

static std::string datadir = DATADIR;

/***********************************************************/
/*      Read Fission Barrier Parameters from a File        */
/***********************************************************/
double beta2Read(ZAnumber *za)
{
  std::ifstream fp;
  std::string   d;
  unsigned int  zt,at;
  double        beta2 = 0.0;

  std::string str = FRDMDATAFILE;

  /*** try current directry */
  fp.open(&str[0]);

  /*** then system data area */
  if(!fp){
    str = datadir + "/" + str;
    fp.open(&str[0]);
  }

  if(!fp){
    message << "FRDM shape parameter file " << FRDMDATAFILE << " not found";
    cohTerminateCode("beta2Read");
  }

  bool found=false;
  while(getline(fp,str)){
    if(str[0] == '#') continue;

    d=str.substr( 0, 5);  zt = atoi(&d[0]);
    d=str.substr( 5, 5);  at = atoi(&d[0]);

    if(zt == za->getZ() && at == za->getA()){
      d=str.substr(60,12);  beta2  = atof(&d[0]);
//    d=str.substr(72,12);  beta4  = atof(&d[0]);
//    d=str.substr(84,12);  beta6  = atof(&d[0]);
      found = true;
      break;
    }
  }
  fp.close();


  if(!found){
    message << "FRDM shape parameter for Z " << za->getZ() << " - A " << za->getA() << " not found";
    cohTerminateCode("beta2Read");
  }

  return(beta2);
}

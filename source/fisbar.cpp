/******************************************************************************/
/*  fisbar.cpp                                                                */
/*        read fission barrier parameters                                     */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>

#include "structur.h"
#include "fisbar.h"
#include "terminate.h"

static std::string datadir = DATADIR;

/***********************************************************/
/*      Read Fission Barrier Parameters from a File        */
/***********************************************************/
int fisDataRead(ZAnumber *za, Fission *fbr)
{
  std::ifstream  fp;
  std::string    d;
  int zt,at,ac;

  std::string str = FISBARDATAFILE;

  /*** try current directry */
  fp.open(&str[0]);

  /*** then system data area */
  if(!fp){
    str = datadir + "/" + str;
    fp.open(&str[0]);
  }

  if(!fp) cohTerminateCode("fission barrier parameter file not found");

  int ic = 0;
  while(getline(fp,str)){
    if(str[0] == '#') continue;

    d=str.substr( 4, 4);  zt = atoi(&d[0]);
    d=str.substr( 8, 4);  at = atoi(&d[0]);
    d=str.substr(12, 4);  ac = atoi(&d[0]);

    ZAnumber za1(zt,at);

    if(za1 == *za){
      
      fbr[ic].za.setZA(zt,ac);

      d=str.substr(16, 8);  fbr[ic].barrier[0].height    = atof(&d[0]);
      d=str.substr(24, 8);  fbr[ic].barrier[0].curvature = atof(&d[0]);
      d=str.substr(32, 8);  fbr[ic].barrier[1].height    = atof(&d[0]);
      d=str.substr(40, 8);  fbr[ic].barrier[1].curvature = atof(&d[0]);

      ic++;
      if(ic >= MAX_FISS_CHANCE) break;
    }
  }
  fp.close();


  if(ic == 0){
    std::ostringstream os;
    os << "fission barrier parameter for Z " << za->getZ() << " - A " << za->getA() << " not found";
    cohTerminateCode(os.str());
  }

  return(ic);
}

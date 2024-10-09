/******************************************************************************/
/*  ffpmultichance.cpp                                                        */
/*        read multi-chance fission data from file                            */
/******************************************************************************/

#include <string>
#include <sstream>
#include <ostream>
#include <fstream>
#include <iostream>

#include "terminate.h"

static const int MaxChance = 4;
static const int Nitem     = 4;

/**********************************************************/
/*      Multi-Chance Fission Data from Data File          */
/**********************************************************/
int FFPReadMultichanceFissionData(const std::string fname, const double p, double *q)
{
  std::ifstream fp;
  std::string d, str, file = fname;

  /*** current directry */
  fp.open(&file[0]);
  if(!fp){
    message << "multi-chance fission data file " << fname << " not found";
    cohTerminateCode("FFPReadMultichanceFissionData");
  }

  /*** check the number of energy-points first */
  int nx = 0;
  while(getline(fp,str)){
    if(str[0] == '#') continue;
    nx ++;
  }
  fp.close();

  fp.open(&file[0]);

  /*** read data and store in a temporal buffer */
  double *x, **y;
  x = new double [nx];
  y = new double * [nx];
  for(int i=0 ; i<nx ; i++) y[i] = new double [Nitem * MaxChance];

  int i = 0;
  while(getline(fp,str)){
    if(str[0] == '#') continue;
    std::istringstream ss(str);
    ss >> x[i];
    for(int j=0 ; j<Nitem * MaxChance ; j++) ss >> y[i][j];
    i++;
  }
  fp.close();


  /*** data interpolation */
  if(p < x[0]){
    for(int j=0 ; j<Nitem * MaxChance ; j++) q[j] = y[0][j];
  }
  else if(p >= x[nx-1]){
    for(int j=0 ; j<Nitem * MaxChance ; j++) q[j] = y[nx-1][j];
  }
  else{
    for(int i=0 ; i<nx-1 ; i++){
      if( (x[i] <= p) && (p < x[i+1]) ){
        for(int j=0 ; j<Nitem * MaxChance ; j++){
          q[j] = (y[i+1][j] - y[i][j]) / (x[i+1] - x[i]) * (p - x[i]) + y[i][j];
        }
        break;
      }
    }
  }

  for(int i=0 ; i<nx ; i++) delete [] y[i];
  delete [] y;
  delete [] x;


  const double pcut = 1e-6;
  int m = MaxChance;
  for(int j=0 ; j<MaxChance ; j++){
    if(q[j] < pcut){ m = j; break; }
  }

  return m;
}


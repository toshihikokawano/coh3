/******************************************************************************/
/*  statreaddensity.cpp                                                       */
/*        read level density, and replce currently stored in the arrays       */
/******************************************************************************/

#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "structur.h"
#include "statmodel.h"
#include "terminate.h"

static const int negrid = 100; // Energy points
static const int njgrid =  50; // maximum J
static       int nefile =   0; // actual number of given energy points
static       int njfile =   0; // actual number of J

#define LINEARINTERPOLATE
#undef  DEBUG_INTERPOLATE
static int  statReadDensityFile (std::ifstream *, const unsigned int, double *, double **, double **);
static void statInterpolateDensity (const int, double *, double **, double **, double *, Parity **);

static inline double readData(std::string s, const int i0, const int i1)
{
  if(s.length() <= (unsigned)(i0+i1)) return(0);
  std::string d = s.substr(i0,i1);
  return(atof(&d[0]));
}

static inline double interpol(const double x, const double x0, const double x1, const double y0, const double y1)
{
  return( (y1-y0)/(x1-x0)*(x-x0) + y0 );
}



/**********************************************************/
/*      Read Level Density File and Replace Density       */
/**********************************************************/
void statReadLevelDensity(Nucleus *n)
{
  std::ostringstream os;
  std::ifstream      fp;
  std::string        str,file;

  double *x, **y0, **y1;
  x  = new double [negrid];
  y0 = new double * [negrid];
  y1 = new double * [negrid];
  for(int k=0 ; k< negrid ; k++){
    y0[k] = new double [negrid];
    y1[k] = new double [negrid];
  }

  os << std::setw(4) << std::setfill('0') << n->za.getZ();
  file = os.str();
  file[0] = 'z';
  str = file + ".tab";

  /*** try current directry first */
  fp.open(&str[0]);
  if(!fp){
    message << "level density data file " << str << " not found";
    cohTerminateCode("statReadLevelDensity");
  }

  std::cout << "# Reading " << str << " for Z = " << n->za.getZ() << "  A = " << n->za.getA() << std::endl;

  if(statReadDensityFile(&fp,n->za.getA(),x,y0,y1) != 0) std::cout << "# A = " << n->za.getA() << " not found" << std::endl;
  else statInterpolateDensity(n->ntotal,x,y0,y1,n->excitation,n->density);

  fp.close();

  for(int k=0 ; k< negrid ; k++){
    delete [] y0[k];
    delete [] y1[k];
  }
  delete [] y0;
  delete [] y1;
  delete [] x;
}


/**********************************************************/
/*      Read Fission Level Density File                   */
/**********************************************************/
void statReadFissionLevelDensity(Nucleus *n)
{
  std::ostringstream os;
  std::ifstream      fp;
  std::string        str1,str2,file;

  double *x, **y0, **y1;
  x  = new double [negrid];
  y0 = new double * [negrid];
  y1 = new double * [negrid];
  for(int k=0 ; k< negrid ; k++){
    y0[k] = new double [negrid];
    y1[k] = new double [negrid];
  }

  os << std::setw(4) << std::setfill('0') << n->za.getZ();
  file = os.str();
  file[0] = 'z';
  str1 = file + "max1.tab";
  str2 = file + "max2.tab";

  /*** first barrier */
  fp.open(&str1[0]);
  if(!fp){
    message << "fission level density data file for the first barrier " << str1 << " not found";
    cohTerminateCode("statReadFissionLevelDensity");
  }

  message << "# Reading " << str1 << " for Z = " << n->za.getZ() << "  A = " << n->za.getA();
  cohNotice("statReadFissionLevelDensity");

  if(statReadDensityFile(&fp,n->za.getA(),x,y0,y1) != 0){
    message << "# A = " << n->za.getA() << " not found";
    cohNotice("statReadFissionLevelDensity");
  }
  else statInterpolateDensity(n->ntotal,x,y0,y1,n->excitation,n->fission->barrier[0].density);

  fp.close();


  /*** second barrier */
  fp.open(&str2[0]);
  if(!fp){
    message << "fission level density data file for the second barrier " << str2 << " not found";
    cohTerminateCode("statReadFissionLevelDensity");
  }

  message << "# Reading " << str2 << " for Z = " << n->za.getZ() << "  A = " << n->za.getA();
  cohNotice("statReadFissionLevelDensity");

  if(statReadDensityFile(&fp,n->za.getA(),x,y0,y1) != 0){
    message << "# A = " << n->za.getA() << " not found";
    cohNotice("statReadFissionLevelDensity");
  }
  else statInterpolateDensity(n->ntotal,x,y0,y1,n->excitation,n->fission->barrier[1].density);

  fp.close();

  for(int k=0 ; k< negrid ; k++){
    delete [] y0[k];
    delete [] y1[k];
  }
  delete [] y0;
  delete [] y1;
  delete [] x;
}


/**********************************************************/
/*      Read File and Copy Data into Arrays               */
/**********************************************************/
int statReadDensityFile(std::ifstream *fp, const unsigned int anum, double *x, double **y0, double **y1)
{
  const int d0 = 42;
  const int dw =  9;

  std::string str;

  bool found = false;
  do{
    getline(*fp,str);
    if((unsigned int)readData(str,31,3) == anum){
      found = true; break;
    }
  }while(*fp);
  if(!found) return(-1);

  /*** even parity */
  getline(*fp,str);
  getline(*fp,str);

  njfile = (str.length() - d0)/dw;
  if(njfile > njgrid) njfile = njgrid;

  nefile = 0;
  for(int k=0 ; ; k++){
    getline(*fp,str);
    if(str == "") break;
    if(k < negrid){
      x[k] = readData(str,0,7);
      for(int j=0 ; j<njgrid ; j++) y0[k][j] = readData(str,d0 + j*dw, dw);
      nefile = k + 1;
    }
  }

  /*** odd parity */
  do{
    getline(*fp,str);
    if((unsigned int)readData(str,31,3) == anum) break;
  }while(*fp);
  getline(*fp,str);
  getline(*fp,str);

  for(int k=0 ; ; k++){
    getline(*fp,str); 
    if(str == "") break;
    if(k < negrid){
      x[k] = readData(str,0,7);
      for(int j=0 ; j<njgrid ; j++) y1[k][j] = readData(str,d0 + j*dw, dw);
    }
  }

  return(0);
}


/**********************************************************/
/*      Interpolate Level Density                         */
/**********************************************************/
void statInterpolateDensity(const int ntotal, double *x, double **y0, double **y1, double *ex, Parity **density)
{
  for(int k0=0 ; k0<ntotal ; k0++){

    /*** when Ex is lower than x[0], linear interpolation */
    if(ex[k0] < x[0]){
      for(int j=0 ; j<MAX_J ; j++){
        if(j < njgrid){

#ifdef LINEARINTERPOLATE
          density[k0][j].even = y0[0][j] / x[0] * ex[k0];
          density[k0][j].odd  = y1[0][j] / x[0] * ex[k0];
#else
          density[k0][j].even = y0[0][j];
          density[k0][j].odd  = y1[0][j];
#endif
        }
      }
    }

    /*** when Ex is inside the table in file */
    else if((x[0] <= ex[k0]) && (ex[k0] <= x[nefile-1])){
      /*** find two points when Ex sits inside */
      int k2 = 0;
      for(int k1=0 ; k1<nefile-1 ; k1++){
        if((x[k1] <= ex[k0]) && (ex[k0] < x[k1+1])){ k2 = k1; break; }
      }

      for(int j=0 ; j<MAX_J ; j++){
        if(j < njgrid){
          if((y0[k2][j] == 0.0) || (y0[k2+1][j] == 0.0)){
            density[k0][j].even = interpol(ex[k0],x[k2],x[k2+1],y0[k2][j],y0[k2+1][j]);
          }
          else{
            double r0 = interpol(ex[k0],x[k2],x[k2+1],log(y0[k2][j]),log(y0[k2+1][j]));
            density[k0][j].even = exp(r0);
          }

          if((y1[k2][j] == 0.0) || (y1[k2+1][j] == 0.0)){
            density[k0][j].odd  = interpol(ex[k0],x[k2],x[k2+1],y1[k2][j],y1[k2+1][j]);
          }
          else{
            double r1 = interpol(ex[k0],x[k2],x[k2+1],log(y1[k2][j]),log(y1[k2+1][j]));
            density[k0][j].odd  = exp(r1);
          }
        }
      }
    }
    /*** outside the range, copy the highest value (ad hoc) */
    else{
      for(int j=0 ; j<MAX_J ; j++){
        if(j < njgrid){
          density[k0][j].even = y0[nefile-1][j];
          density[k0][j].odd  = y1[nefile-1][j];
        }
      }
    }
  }

#ifdef DEBUG_INTERPOLATE
  for(int k0=0 ; k0<ntotal ; k0++){
    double r = 0.0;
    for(int j=0 ; j<MAX_J ; j++) r += density[k0][j].even + density[k0][j].odd;

    std::cout << std::setw(11) << ex[k0];
    std::cout << std::setw(12) << r;
    for(int j=0 ; j<=3 ; j++){
      std::cout << std::setw(12) << density[k0][j].even;
      std::cout << std::setw(12) << density[k0][j].odd;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;

  for(int k1=0 ; k1<nefile ; k1++){
    double r = 0.0;
    for(int j=0 ; j<njgrid ; j++) r += y0[k1][j] + y1[k1][j];

    std::cout << std::setw(11) << x[k1];
    std::cout << std::setw(12) << r;
    for(int j=0 ; j<=3 ; j++){
      std::cout << std::setw(12) << y0[k1][j];
      std::cout << std::setw(12) << y1[k1][j];
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;
#endif
}

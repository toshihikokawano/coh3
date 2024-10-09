/******************************************************************************/
/*  ffpyieldread.cpp                                                          */
/*        Read external data file and construct Y(Z,A,TKE) distribution       */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "masstable.h"
#include "wahl.h"
#include "ffp.h"
#include "terminate.h"


/**********************************************************/
/*      Fragment Yield and Energy from External File      */
/**********************************************************/
void FFPReadPairs(const std::string fname, FissionFragmentPair *ffp)
{
  std::ifstream fp;
  std::string str, file = fname;

  /*** pre-scan file contents */
  fp.open(&file[0]);
  if(!fp){
    message << "fission fragment pair data file " << fname << " not found";
    cohTerminateCode("FFPReadPairs");
  }

  /*** check the number of pairs */
  int nm = 0;
  while(getline(fp,str)){
    if(str[0] == '#') continue;
    nm ++;
  }
  fp.close();

  if(nm > MAX_FISSION_PAIR){
    message << " number of fission pairs " << nm << " in data file " << fname << " too large";
    cohTerminateCode("FFPReadPairs");
  }

  /*** read data and store in the FFP object */
  int x[6];
  double y[8];
  fp.open(&file[0]);

  unsigned int massmin = 1000, massmax = 0;
  nm = 0;
  while(getline(fp,str)){
    if(str[0] == '#') continue;
    std::istringstream ss(str);

    for(int i=0 ; i<6 ; i++) ss >> x[i];
    for(int i=0 ; i<8 ; i++) ss >> y[i];

    if(x[3] <= x[5]){
      ffp->fragment[nm].setPair(x[2],x[3],x[4],x[5]);
      ffp->fragment[nm].el    = y[4];
      ffp->fragment[nm].wl    = y[5];
      ffp->fragment[nm].eh    = y[6];
      ffp->fragment[nm].wh    = y[7];
    }
    else{
      ffp->fragment[nm].setPair(x[4],x[5],x[2],x[3]);
      ffp->fragment[nm].el    = y[6];
      ffp->fragment[nm].wl    = y[7];
      ffp->fragment[nm].eh    = y[4];
      ffp->fragment[nm].wh    = y[5];
    }

    ffp->fragment[nm].yield = y[0];
    ffp->fragment[nm].ek    = y[2];
    ffp->fragment[nm].ex    = y[3];

    double c  = sqrt(ffp->fragment[nm].el*ffp->fragment[nm].el + ffp->fragment[nm].eh*ffp->fragment[nm].eh);
    double s0 = c * ffp->fragment[nm].wl / ffp->fragment[nm].el;
    double s1 = c * ffp->fragment[nm].wh / ffp->fragment[nm].eh;
    ffp->fragment[nm].sigma = (s0 + s1) * 0.5; // s0 = s1 expected

    if(ffp->fragment[nm].getAl() < massmin) massmin = ffp->fragment[nm].getAl();
    if(ffp->fragment[nm].getAh() > massmax) massmax = ffp->fragment[nm].getAh();

    nm ++;
  }
  fp.close();

  /*** set mass range */
  ffp->setN(nm);
  ffp->mass_first = massmin;
  ffp->mass_last  = massmax;

  double mcn = mass_excess(ffp->zf,ffp->af) + ffp->ex;

  for(int i=0 ; i<ffp->getN() ; i++){
    int zl = ffp->fragment[i].getZl();
    int al = ffp->fragment[i].getAl();
    int zh = ffp->fragment[i].getZh();
    int ah = ffp->fragment[i].getAh();

    /*** calculate Q-value for this split, probably we don't use it */
    ffp->fragment[i].qval = mcn - (mass_excess(zl,al) + mass_excess(zh,ah));
  }
}


/**********************************************************/
/*      Read Y(A,TKE) Data from External File             */
/**********************************************************/
void FFPReadExpdata(const std::string fname, FissionFragmentPair *ffp, double *fzn)
{
  std::ifstream fp;
  std::string str, file = fname;

  /*** pre-scan file contents */
  fp.open(&file[0]);
  if(!fp){
    message << "Y(A,TKE) data file " << fname << " not found";
    cohTerminateCode("FFPReadYATKE");
  }

  /*** check the data size, min A, min TKE */
  int mmin = 1000, mmax = 0, tmin = 1000, tmax = 0, m;
  double t, y;
  while(getline(fp,str)){
    if(str[0] == '#') continue;
    std::istringstream ss(str);
    ss >> m >> t >> y;

    if(m > mmax) mmax = m;
    if(m < mmin) mmin = m;

    if((int)t > tmax) tmax = (int)t;
    if((int)t < tmin) tmin = (int)t;
  }
  fp.close();

  int nt = tmax - tmin + 1;
  int nm = mmax - mmin + 1;

  double **buf = new double * [nm];
  for(int i=0 ; i<nm ; i++) buf[i] = new double [nt];

  fp.open(&file[0]);
  while(getline(fp,str)){
    if(str[0] == '#') continue;
    std::istringstream ss(str);
    ss >> m >> t >> y;
    buf[m - mmin][(int)t - tmin] = 0.01 * y; // data given by percent
  }
  fp.close();

  /*** generate fragment pairs from Y(A,TKE) data */
  FFPGeneratePairsExpdata(ffp,fzn,buf,mmin,mmax,tmin,tmax);

  for(int i=0 ; i<nm ; i++) delete [] buf[i];
  delete [] buf;
}

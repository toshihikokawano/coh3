/******************************************************************************/
/*  ffpyield.cpp                                                              */
/*        Calculate Y(Z,A,TKE) distribution                                   */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "masstable.h"
#include "wahl.h"
#include "ffp.h"
#include "terminate.h"

static const unsigned char DEBUG_MASS = 0x01;
static const unsigned char DEBUG_LIST = 0x02;
static const unsigned char DEBUG_NORM = 0x04;
static const unsigned char DEBUG_CHRG = 0x08;
static const unsigned char DEBUG_ENRG = 0x10;
static const unsigned char DEBUG_ATKE = 0x20;
static const unsigned char DEBUG_LUNG = 0x40;

//static unsigned char DEBUG_ID = DEBUG_MASS | DEBUG_LIST;
static unsigned char DEBUG_ID = 0x00;

static void FFPDebug (unsigned char, FissionFragmentPair *);

/* data for fission yield */
static double *tke, *dtke;                    // TKE(A) and its width
static double *mass_yield, **charge_dist;    // Y(A), Y(Z,A)
static int    *charge_first;                  // smallest Z number for Y(Z,A) data

static void FFPMassYieldTune (FissionFragmentPair *, const int, int *, double *);
static void FFPZAList (FissionFragmentPair *);
static void FFPYield (FissionFragmentPair *);
static void FFPTXEDistribution (FissionFragmentPair *);
static void FFPAllocateMemory (FissionFragmentPair *);
static void FFPDeleteAllocated (FissionFragmentPair *);


/**********************************************************/
/*      Mass Yield by Gaussians                           */
/**********************************************************/
void FFPMassYieldParameters(const bool sf, std::string data, double *sigma, double *delta, double *fract, FissionFragmentPair *ffp)
{
  /*** internal Gaussian mass distribution based on experimental data */
  if(data == "internal"){
    FFPGaussianParameters(sf,ffp->zf,ffp->af,ffp->ex,ffp->bn,sigma,delta,fract);
  }
  /*** mass distribution by Wahl, set nu= 0 */
  else if(data == "Wahl"){
    WahlChainYield(ffp->zf,ffp->af,ffp->ex,sigma,delta,fract);
  }

  /*** re-normalize fractions, just in case */
  double s = 2.0*(fract[0] + fract[1] + fract[2]) + fract[3];
  double d = 2.0 - s;
  for(int i=0 ; i<4 ; i++) fract[i] += d * fract[i]/s;

  /*** copy parameters to FissionFragmentProduct object */
  ffp->resetN();
  for(int i=0 ; i<4 ; i++){
    ffp->setGauss(sigma[i], delta[i],fract[i]);
    if(delta[i] != 0.0) ffp->setGauss(sigma[i],-delta[i],fract[i]);
  }
}


/**********************************************************/
/*      Gaussian parameters from external file            */
/**********************************************************/
void FFPConstructYieldParameters(double *ya, double en, double *sigma, double *delta, double *fract, FissionFragmentPair *ffp)
{

  /*** Construct the Gaussian parametrization from the input YA file */
  fract[0] = 1.0/(1.0 + exp((en - ya[0]) / ya[1]));
  fract[1] = 1.0/(1.0 + exp((en - ya[6]) / ya[7]));
  fract[2] = 0.0;
  fract[3] = 2.0 - 2.0*(fract[0] + fract[1]); // keeps the correct normalization

  // TO DO: set fract[3] and reconstruct fract[2]

  delta[0] = ya[2] + ya[3] * en;
  delta[1] = ya[8] + ya[9] * en;
  delta[2] = 0.0;
  delta[3] = ya[12] + ya[13] * en;

  sigma[0] = ya[4] + ya[5] * en;
  sigma[1] = ya[10] + ya[11] * en;
  sigma[2] = 0.0;
  sigma[3] = ya[14] + ya[15] *en;

  /*** copy parameters to FissionFragmentProduct object */
  ffp->resetN();
  for (int i=0 ; i<4 ; i++){
   ffp->setGauss(sigma[i],delta[i],fract[i]);
   if (delta[i] != 0.0) ffp->setGauss(sigma[i],-delta[i],fract[i]);
  }
}


/**********************************************************/
/*      Fission Fragment Pair and Yield Calculation       */
/**********************************************************/
void FFPGeneratePairs(const bool sf, const double tke_input, FissionFragmentPair *ffp, double *fzn, int nft, int *ftmass, double *ftfactor)
{
  /*** allocate data arrays */
  FFPAllocateMemory(ffp);

  /*** mass chain yield by Gaussians */
  for(int i=0 ; i<ffp->nmass ; i++){
    mass_yield[i] = 0.0;
    int a = i + ffp->mass_first;
    for(int k=0 ; k<ffp->getNg() ; k++){
      if(ffp->getGaussWidth(k) > 0.0){
        double da = a - ffp->ac + ffp->getGaussCenter(k);
        mass_yield[i] += gaussian(ffp->getGaussFraction(k),da,ffp->getGaussWidth(k));
      }
    }
  }

  /*** fine tune mass yield */
  if(nft > 0) FFPMassYieldTune(ffp,nft,ftmass,ftfactor);

  if(DEBUG_ID & DEBUG_MASS) FFPDebug(DEBUG_MASS,ffp);

  /*** TKE(A) and its width from the systematics */
  for(int i=0 ; i<ffp->nmass ; i++){
    tke[i] =  FFPSystematics_TKE_A(sf,ffp->zf,ffp->af,i+ffp->mass_first,tke_input);
    dtke[i] = FFPSystematics_TKE_A_Width(sf,ffp->zf,ffp->af,i+ffp->mass_first);
  }

  /*** Z-distribution by Wahl's Zp model */
  WahlZpModel(ffp->zf,ffp->af,ffp->mass_first,ffp->nmass,ffp->ncharge,ffp->ex,charge_dist,charge_first,fzn);

  /*** prepare Z and A list and reaction Q-value in FissionPair object */
  FFPZAList(ffp);
  if(DEBUG_ID & DEBUG_LIST) FFPDebug(DEBUG_LIST,ffp);

  /*** independent yield for a given set of Z and A */
  FFPYield(ffp);
  if(DEBUG_ID & DEBUG_NORM) FFPDebug(DEBUG_NORM,ffp);
  if(DEBUG_ID & DEBUG_CHRG) FFPDebug(DEBUG_CHRG,ffp);

  /*** TXE distributions */
  FFPTXEDistribution(ffp);

  /*** check average TKE and distribution */
  double t0 = tke_input;
  double t1 = FFPAverageTKE(ffp);

  /*** adjust TKE(A) to give the fixed value */
  // if(t0 != 0.0 && t1 != 0.0){
  //   for(int i=0 ; i<ffp->nmass ; i++) tke[i] = FFPSystematics_TKE_A(sf,ffp->zf,ffp->af,i+ffp->mass_first,tke_input) * t0 / t1;
  //   FFPZAList(ffp);
  //   FFPYield(ffp);
  //   FFPTXEDistribution(ffp);
  // }

  if(t0 != 0.0){
    double delta = t0 - t1;
    for(int i=0 ; i<ffp->nmass ; i++) tke[i] = FFPSystematics_TKE_A(sf,ffp->zf,ffp->af,i+ffp->mass_first,tke_input) + delta;
    FFPZAList(ffp);
    FFPYield(ffp);
    FFPTXEDistribution(ffp);
    t1 = FFPAverageTKE(ffp);
  }

  if(DEBUG_ID & DEBUG_ENRG) FFPDebug(DEBUG_ENRG,ffp);
  if(DEBUG_ID & DEBUG_ATKE) FFPDebug(DEBUG_ATKE,ffp);
  if(DEBUG_ID & DEBUG_LUNG) FFPDebug(DEBUG_LUNG,ffp);

  /*** free allocated memory */
  FFPDeleteAllocated(ffp);
}


/**********************************************************/
/*      Fission Fragment Pair When Y(A,TKE) Given         */
/**********************************************************/
void FFPGeneratePairsExpdata(FissionFragmentPair *ffp, double *fzn, double **dat, const int mmin, const int mmax, const int tmin, const int tmax)
{
  /*** allocate data arrays */
  FFPAllocateMemory(ffp);

  /*** min values of A and TKE given in data */
  int nt = tmax - tmin + 1;
  int nm = mmax - mmin + 1;

  /*** Y(A), TKE, dTKE from Y(A,TKE) */
  double s = 0.0;
  for(int i=0 ; i<ffp->nmass ; i++){
    mass_yield[i] = tke[i] = dtke[i] = 0.0;

    int al = i + ffp->mass_first;
    int ah = ffp->af - al;

    if((double)al > ffp->ac) break;

    /*** average Y and TKE */
    double c = 0.0, t = 0.0;
    for(int j=0 ; j<nm ; j++){
      for(int k=0 ; k<nt ; k++){
        if(j + mmin == ah){
          c += dat[j][k];
          t += dat[j][k] * (k + tmin);
        }
      }
    }
    mass_yield[i] = c;
    if(c > 0.0) tke[i] = t / c;
    s += c;

    /*** variance of TKE */
    double d = 0.0;
    for(int j=0 ; j<nm ; j++){
      for(int k=0 ; k<nt ; k++){
        if(j + mmin == ah){
          d +=  dat[j][k] * (k + tmin - tke[i]) * (k + tmin - tke[i]);
        }
      }
    }
    if(c > 0.0) dtke[i] = sqrt(d / c);
  }
  if(s > 0.0){
    for(int i=0 ; i<ffp->nmass ; i++) mass_yield[i] /= s;
  }
  if(DEBUG_ID & DEBUG_MASS) FFPDebug(DEBUG_MASS,ffp);

  WahlZpModel(ffp->zf,ffp->af,ffp->mass_first,ffp->nmass,ffp->ncharge,ffp->ex,charge_dist,charge_first,fzn);
  FFPZAList(ffp);
  if(DEBUG_ID & DEBUG_LIST) FFPDebug(DEBUG_LIST,ffp);

  FFPYield(ffp);
  if(DEBUG_ID & DEBUG_NORM) FFPDebug(DEBUG_NORM,ffp);
  if(DEBUG_ID & DEBUG_CHRG) FFPDebug(DEBUG_CHRG,ffp);

  FFPTXEDistribution(ffp);

  if(DEBUG_ID & DEBUG_ENRG) FFPDebug(DEBUG_ENRG,ffp);
  if(DEBUG_ID & DEBUG_ATKE) FFPDebug(DEBUG_ATKE,ffp);
  if(DEBUG_ID & DEBUG_LUNG) FFPDebug(DEBUG_LUNG,ffp);

  /*** free allocated memory */
  FFPDeleteAllocated(ffp);
}


/**********************************************************/
/*      Fission Fragment Pair and Yield Calculation       */
/*      When TKE(A) and sTKE(A) are read in from file     */
/**********************************************************/
void FFPGeneratePairs(const bool sf, const double tke_input, FissionFragmentPair *ffp, double *fzn, double *tkea, double *stkea, int nft, int *ftmass, double *ftfactor)
{
  /*** allocate data arrays */
  FFPAllocateMemory(ffp);

  /*** mass chain yield by Gaussians */
  for(int i=0 ; i<ffp->nmass ; i++){
    mass_yield[i] = 0.0;
    int a = i + ffp->mass_first;
    for(int k=0 ; k<ffp->getNg() ; k++){
      if(ffp->getGaussWidth(k) > 0.0){
        double da = a - ffp->ac + ffp->getGaussCenter(k);
        mass_yield[i] += gaussian(ffp->getGaussFraction(k),da,ffp->getGaussWidth(k));
      }
    }
  }

  /*** fine tune mass yield */
  if(nft > 0) FFPMassYieldTune(ffp,nft,ftmass,ftfactor);


  if(DEBUG_ID & DEBUG_MASS) FFPDebug(DEBUG_MASS,ffp);

  /*** TKE(A) and its width from the systematics */
  for(int i=0 ; i<ffp->nmass ; i++){
    if(sf){
      tke[i] =  FFPSystematics_TKE_A(sf,ffp->zf,ffp->af,i+ffp->mass_first,tke_input);
      dtke[i] = FFPSystematics_TKE_A_Width(sf,ffp->zf,ffp->af,i+ffp->mass_first);
    }
    else{
      tke[i] = FFPConstructTKEA(ffp->af,i+ffp->mass_first,tkea);
      dtke[i] = FFPConstructSigmaTKEA(ffp->af,i+ffp->mass_first,stkea);
    }
  }

  /*** Z-distribution by Wahl's Zp model */
  WahlZpModel(ffp->zf,ffp->af,ffp->mass_first,ffp->nmass,ffp->ncharge,ffp->ex,charge_dist,charge_first,fzn);

  /*** prepare Z and A list and reaction Q-value in FissionPair object */
  FFPZAList(ffp);
  if(DEBUG_ID & DEBUG_LIST) FFPDebug(DEBUG_LIST,ffp);

  /*** independent yield for a given set of Z and A */
  FFPYield(ffp);
  if(DEBUG_ID & DEBUG_NORM) FFPDebug(DEBUG_NORM,ffp);
  if(DEBUG_ID & DEBUG_CHRG) FFPDebug(DEBUG_CHRG,ffp);

  /*** TXE distributions */
  FFPTXEDistribution(ffp);

  /*** check average TKE and distribution */
  double t0 = tke_input;
  double t1 = FFPAverageTKE(ffp);

  /*** adjust TKE(A) to give the fixed value */
  // if(t0 != 0.0 && t1 != 0.0){
  //   for(int i=0 ; i<ffp->nmass ; i++) tke[i] = FFPSystematics_TKE_A(sf,ffp->zf,ffp->af,i+ffp->mass_first,tke_input) * t0 / t1;
  //   FFPZAList(ffp);
  //   FFPYield(ffp);
  //   FFPTXEDistribution(ffp);
  // }

  if(t0 != 0.0){
    double delta = t0 - t1;
    //for(int i=0 ; i<ffp->nmass ; i++) tke[i] = FFPSystematics_TKE_A(sf,ffp->zf,ffp->af,i+ffp->mass_first,tke_input) + delta;
    for(int i=0 ; i<ffp->nmass ; i++) tke[i] = FFPConstructTKEA(ffp->af,i+ffp->mass_first,tkea) + delta;
    FFPZAList(ffp);
    FFPYield(ffp);
    FFPTXEDistribution(ffp);
    t1 = FFPAverageTKE(ffp);
  }

  if(DEBUG_ID & DEBUG_ENRG) FFPDebug(DEBUG_ENRG,ffp);
  if(DEBUG_ID & DEBUG_ATKE) FFPDebug(DEBUG_ATKE,ffp);


  /*** free allocated memory */
  FFPDeleteAllocated(ffp);
}


/**********************************************************/
/*      Fine Tune Mass Yield at Given Mass Points         */
/**********************************************************/
void FFPMassYieldTune(FissionFragmentPair *ffp, const int nft, int *ftmass, double *ftfactor)
{
  const double width = 2.0;
  double *dy = new double [ffp->nmass];
  for(int i=0 ; i<ffp->nmass ; i++) dy[i] = 0.0;

  /*** for all given A numbers */
  int al = 0, ah = 0;
  for(int j=0 ; j<nft ; j++){
    if(ftmass[j] < ffp->ac){
      al = ftmass[j];
      ah = ffp->af - al;
    }
    else{
      ah = ftmass[j];
      al = ffp->af - ah;
    }

    /*** Gaussian centered at Al and Ah with the default width */
    for(int i=0 ; i<ffp->nmass ; i++){
      int a = i + ffp->mass_first;
      dy[i] += gaussian(ftfactor[j]-1.0, a-al, width);
      dy[i] += gaussian(ftfactor[j]-1.0, a-ah, width);
    }
  }

  /*** add to mass yield, and re-normalize it again */
  double s = 0.0;
  for(int i=0 ; i<ffp->nmass ; i++){
    mass_yield[i] += dy[i];
    s += mass_yield[i];
  }
  s = 2.0 / s;
  for(int i=0 ; i<ffp->nmass ; i++) mass_yield[i] *= s;

  delete [] dy;
}


/**********************************************************/
/*      Store Z and A Data                                */
/**********************************************************/
void FFPZAList(FissionFragmentPair *ffp)
{
  double mcn = mass_excess(ffp->zf,ffp->af) + ffp->ex;

  int k = 0;
  for(int i=0 ; i<ffp->nmass ; i++){
    int al = i + ffp->mass_first;
    int ah = ffp->af - al;

    if((double)al > ffp->ac) break;

    for(int j=0 ; j< ffp->ncharge ; j++){
      int zl = j + charge_first[i];
      int zh = ffp->zf - zl;

      double q =  mcn - (mass_excess(zl,al) + mass_excess(zh,ah));

      ffp->fragment[k].setPair(zl,al,zh,ah);
      ffp->fragment[k].qval  = q;
      ffp->fragment[k].yield = 0.0;

      k++;
      if(k > ffp->getMaxPair()){
        std::cerr << "too many fission pairs " << k << std::endl;
        return;
      }

      if((al == ah) && (zl == zh)) break;
    }
  }

  ffp->setN(k);
}


/**********************************************************/
/*      Joint Mass and Charge Distributions               */
/**********************************************************/
void FFPYield(FissionFragmentPair *ffp)
{
  /* find min/max Al in FissionPair object */
  int almin = 1000, almax = 0;
  for(int i=0 ; i<ffp->getN() ; i++){
    int al = ffp->fragment[i].getAl();
    if(al < almin) almin = al;
    if(al > almax) almax = al;
  }

  for(int al=almin ; al<=almax ; al++){
    int i = al - ffp->mass_first; // index for Al in the charge_dist array

    /* for each A, determine the elements those are within the Z-range
       and renormalize them */
  
    double sum = 0.0;
    for(int k=0 ; k<ffp->getN() ; k++){
      if(ffp->fragment[k].getAl() != (unsigned int)al) continue;
      int j = ffp->fragment[k].getZl() - charge_first[i]; // index for Zl for a given Al
      if((0 <= j) && (j < ffp->ncharge)) sum += charge_dist[i][j];
    }

    /* mass chain yield times F(Z) */
    if(sum > 0.0) sum = 1.0/sum;
    else sum = 0.0;

    /* when symmetric, divide by 2 to make a smooth mass chain yield */
    if(ffp->af%2 == 0 && al == almax) sum *= 0.5;

    for(int k=0 ; k<ffp->getN() ; k++){
      if(ffp->fragment[k].getAl() != (unsigned int)al) continue;
      int j = ffp->fragment[k].getZl() - charge_first[i];
      if((0 <= j) && (j < ffp->ncharge)) ffp->fragment[k].yield = mass_yield[i] * charge_dist[i][j] * sum;
    }
  }
}


/**********************************************************/
/*      Calculate average TXE                             */
/**********************************************************/
void FFPTXEDistribution(FissionFragmentPair *ffp)
{
  const double TXEcut = 0.5; // eliminate pairs if excitation energy is too low

  for(int k=0 ; k<ffp->getN() ; k++){

    /* average TKE */
    double tke0 = tke[ffp->fragment[k].getAl() - ffp->mass_first];

    /* find the same A and calculate charge-product */
    double s = 0.0;
    int n = 0;
    for(int l=0 ; l<ffp->getN() ; l++){
      if(ffp->fragment[k].getAl() == ffp->fragment[l].getAl()){
        s += (double)ffp->fragment[l].getZl() * (double)ffp->fragment[l].getZh();
        n++;
      }
    }
    /* we assume TKE for each separation is proportional to the Coulomb repulsion */
    s = tke0 * (double)n / s;
    ffp->fragment[k].ek = (double)ffp->fragment[k].getZl() * (double)ffp->fragment[k].getZh() * s;

    /* average TXE */
    ffp->fragment[k].ex = ffp->fragment[k].qval - ffp->fragment[k].ek;

    /* width of TXE is proportional to TKE */
    ffp->fragment[k].sigma =  ffp->fragment[k].ek * dtke[ffp->fragment[k].getAl() - ffp->mass_first] / tke0;

    /* if average TXE is too small, exclude this pair */
    if(ffp->fragment[k].ex < TXEcut) ffp->fragment[k].yield = 0.0;
  }

  /* re-adjust yields */
  double s = 0.0;
  for(int k=0 ; k<ffp->getN() ; k++) s += ffp->fragment[k].yield;
  for(int k=0 ; k<ffp->getN() ; k++) ffp->fragment[k].yield /= s;
}


/**********************************************************/
/*      Calculate Average TKE                             */
/**********************************************************/
double FFPAverageTKE(FissionFragmentPair *ffp)
{
  /* Y sum and average TKE */
  double y = 0.0, e = 0.0;
  for(int k=0 ; k<ffp->getN() ; k++){
    y += ffp->fragment[k].yield;
    e += ffp->fragment[k].yield * ffp->fragment[k].ek;
  }
  if(y > 0.0) e /= y;

  return e;
}


/**********************************************************/
/*      Expand List to Include All FFs                    */
/**********************************************************/
void FFPExpandList(const int nc, FissionFragmentPair *ffp)
{
  if(nc == 1) return;

  double mcn = mass_excess(ffp[0].zf,ffp[0].af) + ffp[0].ex;

  /*** ffp[0] will include all the FF pairs */
  
  for(int i=1 ; i<nc ; i++){
    for(int k1=0 ; k1<ffp[i].getN() ; k1++){
      int zl1 = ffp[i].fragment[k1].getZl();
      int al1 = ffp[i].fragment[k1].getAl();
      int zh1 = ffp[i].fragment[k1].getZh();
      int ah1 = ffp[i].fragment[k1].getAh();


      /*** fisrt, search for (Zl,Al) in ffp[0] */
      bool found = false;
      int km = ffp[0].getN();
      for(int k0=0 ; k0<km ; k0++){
        int zl0 = ffp[0].fragment[k0].getZl();
        int al0 = ffp[0].fragment[k0].getAl();
        if( (zl0 == zl1) && (al0 == al1) ){ found = true; break; }
      }

      if( (!found) && (km < ffp[0].getMaxPair()) ){
        int zh0 = ffp[0].zf - zl1;  // this is for the first CN
        int ah0 = ffp[0].af - al1;

        double q =  mcn - (mass_excess(zl1,al1) + mass_excess(zh0,ah0));

        if( (al1 == ah0) && (zl1 > zh0) ) ffp[0].fragment[km].setPair(zh0,al1,zl1,ah0); // special case for Zl > Zh
        else  ffp[0].fragment[km].setPair(zl1,al1,zh0,ah0);
        ffp[0].fragment[km].qval  = q;
        ffp[0].fragment[km].ek    = 0.0;
        ffp[0].fragment[km].ex    = 0.0;
        ffp[0].fragment[km].sigma = 0.0;
        ffp[0].fragment[km].yield = 0.0;
        km ++;
        ffp[0].setN(km);
      }

      /*** second, search for (Zh,Ah) in ffp[0] */
      found = false;
      for(int k0=0 ; k0<km ; k0++){
        int zh0 = ffp[0].fragment[k0].getZh();
        int ah0 = ffp[0].fragment[k0].getAh();
        if( (zh0 == zh1) && (ah0 == ah1) ){ found = true; break; }
      }

      if( (!found) && (km < ffp[0].getMaxPair()) ){
        int zl0 = ffp[0].zf - zh1;
        int al0 = ffp[0].af - ah1;

        double q =  mcn - (mass_excess(zl0,al0) + mass_excess(zh1,ah1));

        if( (al0 == ah1) && (zl0 > zh1) ) ffp[0].fragment[km].setPair(zh1,al0,zl0,ah1);
        else ffp[0].fragment[km].setPair(zl0,al0,zh1,ah1);
        ffp[0].fragment[km].qval  = q;
        ffp[0].fragment[km].ek    = 0.0;
        ffp[0].fragment[km].ex    = 0.0;
        ffp[0].fragment[km].sigma = 0.0;
        ffp[0].fragment[km].yield = 0.0;
        km ++;
        ffp[0].setN(km);
      }
    }
  }
//FFPDebug(DEBUG_NORM,ffp);

  /*** duplication test */
/*
  for(int k0=0 ; k0<ffp[0].getN() ; k0++){
    int zl0 = ffp[0].fragment[k0].getZl();
    int al0 = ffp[0].fragment[k0].getAl();
    int zh0 = ffp[0].fragment[k0].getZh();
    int ah0 = ffp[0].fragment[k0].getAh();

    for(int k1=k0+1 ; k1<ffp[0].getN() ; k1++){
      int zl1 = ffp[0].fragment[k1].getZl();
      int al1 = ffp[0].fragment[k1].getAl();
      int zh1 = ffp[0].fragment[k1].getZh();
      int ah1 = ffp[0].fragment[k1].getAh();

      if( (zl0 == zl1) && (al0 == al1) && (zh0 == zh1) && (ah0 == ah1) ){
        std::cout << "duplicated " << zl0 << " " << al0 << " " << zh0 << " " << ah0 << std::endl;
      }
    }
  }
*/
}


/**********************************************************/
/*      Allocate Temporal Arrays                          */
/**********************************************************/
void FFPAllocateMemory(FissionFragmentPair *ffp)
{
  try{
    tke          = new double [ffp->nmass];   // TKE(A)
    dtke         = new double [ffp->nmass];   // width of TKE(A)
    mass_yield   = new double [ffp->nmass];   // Mass Yield Y(A)
    charge_first = new int    [ffp->nmass];   // lowest Z number for each Z-dist
    charge_dist  = new double * [ffp->nmass]; // probability of Z for a given A
    for(int i=0 ; i<ffp->nmass; i++) charge_dist[i] = new double [ffp->ncharge];
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("FFPAllocateMemory");
  }

  for(int a=0 ; a<ffp->nmass; a++){
    tke[a] = dtke[a] = mass_yield[a] = 0.0;
    charge_first[a] = 0;
    for(int z=0 ; z<ffp->ncharge; z++) charge_dist[a][z] = 0.0;
  }
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void FFPDeleteAllocated(FissionFragmentPair *ffp)
{
  delete [] tke;
  delete [] dtke;
  delete [] mass_yield;

  for(int a=0 ; a<ffp->nmass; a++) delete [] charge_dist[a];
  delete [] charge_dist;
  delete [] charge_first;
}


/**********************************************************/
/*      for Debugging                                     */
/**********************************************************/
void FFPDebug(unsigned char id, FissionFragmentPair *ffp)
{
  int zlmin = 1000, zlmax = 0;
  int almin = 1000, almax = 0;
  for(int i=0 ; i<ffp->getN() ; i++){
    int zl = ffp->fragment[i].getZl();
    if(zl < zlmin) zlmin = zl;
    if(zl > zlmax) zlmax = zl;
    int al = ffp->fragment[i].getAl();
    if(al < almin) almin = al;
    if(al > almax) almax = al;
  }


  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4);

  /* print mass chain yield */
  if(id & DEBUG_MASS){
    double s = 0.0;
    for(int i=0 ; i<ffp->nmass ; i++){
      s += mass_yield[i];
      std::cout << std::setw(5) << i + ffp->mass_first;
      std::cout << std::setw(14) << mass_yield[i];
      std::cout << std::setw(14) << s << std::endl;
    }
  }

  /* print all list */
  if(id & DEBUG_LIST){
    for(int k=0 ; k<ffp->getN() ; k++){
      std::cout << std::setw(5) << k;
      std::cout << std::setw(5) << ffp->fragment[k].getZl();
      std::cout << std::setw(5) << ffp->fragment[k].getAl();
      std::cout << std::setw(5) << ffp->fragment[k].getZh();
      std::cout << std::setw(5) << ffp->fragment[k].getAh();
      std::cout << std::setw(14) << ffp->fragment[k].qval << std::endl;
    }
  }

  /* check normalization */
  if(id & DEBUG_NORM){
    double s = 0.0;
    for(int k=0 ; k<ffp->getN() ; k++){
      s += ffp->fragment[k].yield;
      std::cout << std::setw(5) << k;
      std::cout << std::setw(5) << ffp->fragment[k].getZl();
      std::cout << std::setw(5) << ffp->fragment[k].getAl();
      std::cout << std::setw(5) << ffp->fragment[k].getZh();
      std::cout << std::setw(5) << ffp->fragment[k].getAh();
      std::cout << std::setw(14) << ffp->fragment[k].yield;
      std::cout << std::setw(14) << s << std::endl;
    }
  }

  /* check charge distribution */
  if(id & DEBUG_CHRG){
    double x = 0.0;
    for(int zl=zlmin ; zl<=zlmax ; zl++){
      double y = 0.0;
      for(int k=0 ; k<ffp->getN() ; k++){
        if(ffp->fragment[k].getZl() == (unsigned int)zl) y += ffp->fragment[k].yield;
      }
      x += y;
      std::cout << std::setw(5) << zl;
      std::cout << std::setw(5) << ffp->zf - zl;
      std::cout << std::setw(14) << y;
      std::cout << std::setw(14) << x << std::endl;
    }
  }

  /* check TKE and TXE */
  if(id & DEBUG_ENRG){
    double s = 0.0;
    for(int k=0 ; k<ffp->getN() ; k++){
      if(ffp->fragment[k].yield == 0.0) continue;
      s += ffp->fragment[k].yield;
      std::cout << std::setw(5)  << k;
      std::cout << std::setw(5)  << ffp->fragment[k].getZl();
      std::cout << std::setw(5)  << ffp->fragment[k].getAl();
      std::cout << std::setw(5)  << ffp->fragment[k].getZh();
      std::cout << std::setw(5)  << ffp->fragment[k].getAh();
      std::cout << std::setw(14) << ffp->fragment[k].qval;
      std::cout << std::setw(14) << ffp->fragment[k].ek;
      std::cout << std::setw(14) << ffp->fragment[k].ex;
      std::cout << std::setw(14) << ffp->fragment[k].sigma;
      std::cout << std::setw(14) << ffp->fragment[k].yield;
      std::cout << std::setw(14) << s << std::endl;
    }
  }

  /* check TKE and its width */
  if(id & DEBUG_ATKE){
    for(int al=almin ; al<=almax ; al++){
      double tke0 = tke[al - ffp->mass_first];
      if(tke0 == 0.0) continue;

      double t = 0.0, w = 0.0, y = 0.0;
      for(int k=0 ; k<ffp->getN() ; k++){
        if(ffp->fragment[k].getAl() == (unsigned int)al){
          t += ffp->fragment[k].ek * ffp->fragment[k].yield;
          w += ffp->fragment[k].sigma * ffp->fragment[k].sigma * ffp->fragment[k].yield;
          y += ffp->fragment[k].yield;
        }
      }
      if(y > 0.0){
        t = t/y;
        w = sqrt(w/y);
      }
      std::cout << std::setw(5)  << al;
      std::cout << std::setw(5)  << ffp->af - al;
      std::cout << std::setw(14) << tke0;
      std::cout << std::setw(14) << t;
      std::cout << std::setw(14) << t/tke0;
      std::cout << std::setw(14) << w << std::endl;
    }
  }

  /* lung plot */
  if(id & DEBUG_LUNG){

    int tmax = 250;
    int tmin = 100;
    double *b0 = new double [tmax - tmin + 1];
    double *b1 = new double [tmax - tmin + 1];

    for(int al=almin ; al<=almax ; al++){
      for(int t=0 ; t<=tmax-tmin ; t++) b0[t] = 0.0;

      for(int k=0 ; k<ffp->getN() ; k++){
        if(ffp->fragment[k].getAl() != (unsigned int)al) continue;
        if(ffp->fragment[k].yield == 0.0) continue;

        for(int t=0 ; t<=tmax-tmin; t++) b1[t] = 0.0;
        double s = 0.0;
        for(int t=0 ; t<=tmax-tmin ; t++){
          double tx = t + tmin + 0.5;
          double x1 = ffp->fragment[k].ek - tx;
          double x2 = ffp->fragment[k].sigma * ffp->fragment[k].sigma * 2.0;
          b1[t] = exp( - x1 * x1 / x2 );
          s += b1[t];
        }
        if(s > 0.0){
          s = ffp->fragment[k].yield / s;
          for(int t=0 ; t<=tmax-tmin ; t++) b0[t] += s * b1[t];
        }
      }

      for(int t=0 ; t<=tmax-tmin ; t++){
        std::cout << std::setw(5)  << al;
        std::cout << std::setw(5)  << ffp->af - al;
        std::cout << std::setw(5)  << t + tmin;
        std::cout << std::setw(14) << b0[t] << std::endl;;
      }
      std::cout << std::endl;
    }

    delete [] b0;
    delete [] b1;
  }
}

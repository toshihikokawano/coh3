/******************************************************************************/
/*  beohffragbeta.cpp                                                         */
/*        Beta-decay of independent fission products                          */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "beoh.h"
#include "beohoutput.h"
#include "output.h"
#include "outformat.h"
#include "terminate.h"
#include "dir.h"
#include "ripl2levels.h"

/*** default decay data file name */
static std::string datadir = DATADIR;
static std::string DEFAULT_DECAY_FILE = "decaybranchENDF.dat";


/*** define the longest half-live for the long lived isotope */
static double DEFAULT_LONG_LIVED = 1000.0; // 1000 years

static int MAX_NUCLIDES =   2000;
static int MAX_NUKINDEX = 150000;

static Isotope      *nuk = nullptr;  // properties of each nuclide
static unsigned int *idx = nullptr;  // idx[MZA] = sequential number in the data file

static int  cfyCopyPromptYield  (void);
static int  cfyReadDecayData    (std::string, int, double);
static void cfyCumulativeYield  (const int);
static void cfyCheckData        (const int);

static bool PRN_CHECKDATA  = false;

/**********************************************************/
/*      Decay                                             */
/**********************************************************/
void beohFragmentBetaDecay(double long_lived, std::string file)
{
  if(file == "") file = DEFAULT_DECAY_FILE;
  if(long_lived == 0.0) long_lived = DEFAULT_LONG_LIVED;

  /*** allocate data for beta-decay calculation */
  nuk  = new Isotope [MAX_NUCLIDES];
  idx  = new unsigned int [MAX_NUKINDEX];

  /*** initialize nuk array */
  for(int i=0 ; i<MAX_NUCLIDES ; i++) nuk[i].init();

  /*** clear idx array, -1 means no data given */
  for(int i=0 ; i<MAX_NUKINDEX ; i++) idx[i] = -1;

  /*** copy yield data from BeoH result */
  int nfpy = cfyCopyPromptYield();

  /*** read decay data file, half life and decay mode */
  nfpy = cfyReadDecayData(file,nfpy,long_lived);

  if(PRN_CHECKDATA) cfyCheckData(nfpy);

  /*** time-independent cumulative yield calculation */
  cfyCumulativeYield(nfpy);

  /*** print results */
  outFissionProductYield(nfpy,nuk);
  outFissionProductChainYield(nfpy,nuk);
  outDelayedNeutronYield(nfpy,nuk);

  delete [] nuk;
  delete [] idx;
}


/**********************************************************/
/*      Copy Prompt Fission Yield into Nuk Object         */
/**********************************************************/
int cfyCopyPromptYield()
{
  /*** point local object */
  CumulativeResidualProduct *res = outPrepObjectPointer();

  /*** copy data */
  int ntotal = res->getNcurrent();
  if(ntotal >= MAX_NUCLIDES){
    message << "number of nuclides larger than the maximum; " << ntotal << " > " << MAX_NUCLIDES;
      cohTerminateCode("cfyCopyPromptYield");
  }

  for(int i=0 ; i<ntotal ; i++){
    int    z = res->rp[i].za.getZ();
    int    a = res->rp[i].za.getA();
    int    m = res->rp[i].metaflag;
    double y = res->rp[i].production;

    /*** MZA contains M, Z, and A in a field */
    unsigned int mza = nuk[i].setZA(m,z,a);
    nuk[i].setJP(res->rp[i].spin,res->rp[i].parity);
    nuk[i].initial = y;   // initial yield
    nuk[i].yield   = y;   // final yield
    nuk[i].std     = 0.0; // uncertainty

    if(mza > (unsigned int)MAX_NUKINDEX){
      message << "nuclide index exceeded the maximum " << mza;
      cohTerminateCode("cfyCopyPromptYield");
    }

    idx[mza] = i;

//  cout << i <<" " << nuk[i].getZ() << " " << nuk[i].getA() << " " << nuk[i].getM() << " " << nuk[i].initial << " " << mza << endl;
  }

  return(ntotal);
}


/**********************************************************/
/*      Read Decay Data                                   */
/*          decay mode, decay constant, branching ratios  */
/*          CINDER decay_infor format                     */
/**********************************************************/
int cfyReadDecayData(std::string decaydata, int nuktotal, double long_lived)
{
  std::ifstream fp;
  std::string   line, s;
  int      z0, z1, a0, a1, m0, m1, nbranch, next[MAX_DECAYMODE];
  unsigned long za0, za1[MAX_DECAYMODE], mza0, mza1;
  double  lambda, br[MAX_DECAYMODE];

  /*** terminate chain for long-lived FPs */
  if(long_lived > 0.0) long_lived = log(2.0)/ (long_lived * 31556926.0);


  /*** try current directry first */
  fp.open(&decaydata[0]);

  /*** then system data area */
  if(!fp){
    decaydata = datadir + "/" + decaydata;
    fp.open(&decaydata[0]);
  }

  if(!fp){
    message << "beta-decay branching ratio data file " << decaydata << " not found";
    cohTerminateCode("cfyReadDecayData");
  }

  while( getline(fp,line) ){
    if(fp.eof() != 0) break;

    /*** each decay data entry */
    s    = line.substr(0, 8);
    za0  = atoi(s.c_str());
    a0   = (int)(za0/10000.0);
    z0   = (int)((za0 - a0*10000)/10.0);
    m0   = (int)(za0 - a0*10000 - z0*10);
    mza0 = tomza(m0,z0,a0);

    /*** check if stable, if not calculate decay constant */
    s = line.substr(18,11);
    if(s == "   STABLE  ") lambda = 0.0;
    else                   lambda = log(2.0)/atof(s.c_str());

    /*** number of decay modes (number daughter nuclei) */
    s = line.substr(29, 6);
    nbranch = atoi(s.c_str());

    /*** for each decay branch */
    for(int i=0 ; i<nbranch ; i++){
      getline(fp,line);
      s = line.substr( 7,11);  br[i]  = atof(s.c_str()); // branching ratio
      s = line.substr(19, 8);  za1[i] = atoi(s.c_str()); // final nuclide ZAM
    }

    if(mza0 >= (unsigned int)MAX_NUKINDEX) continue;

    /*** if this nuclide appears in decaying nuclide list */
    int k0 = idx[mza0];

    if(k0 >= 0){

      for(int i=0 ; i<MAX_DECAYMODE ; i++) next[i] = 0;

      /*** eliminate if ID is outside the range */
      double s = 0.0; // sum of branching ratios
      int    n = 0;   // actual number of decays

      for(int i=0 ; i<nbranch ; i++){
        a1   = (int)(za1[i]/10000.0);
        z1   = (int)((za1[i] - a1*10000)/10.0);
        m1   = (int)(za1[i] - a1*10000 - z1*10);
        mza1 = tomza(m1,z1,a1);
        next[i] = (mza1 >= (unsigned int)MAX_NUKINDEX) ? -2 : idx[mza1];

        if(next[i] == -2){ s += br[i]; n++; }
        /*** if the daughter nuclide cannot be found, add it */
        else if(next[i] == -1){

          /*** read RIPL */
          NuclearStructure ns;
          ns.memalloc(MAX_LEVELS);
          ns.za = ZAnumber(z1,a1);

          int nlev = riplReadDiscreteLevels(&ns.za,ns.lev,reassign);

          /*** find the m1-th metastable state in RIPL */
          int km = 0;
          if(m1 > 0){
            for(int k=1 ; k<nlev ; k++){
              if(ns.lev[k].halflife > thalfmin) km ++;
              if(km == m1) break;
            }
          }
          
          nuk[nuktotal].setZA(m1,z1,a1);
          nuk[nuktotal].setJP(ns.lev[km].spin,ns.lev[km].parity);
          idx[mza1] = next[i] = nuktotal;
          nuktotal ++;
        }
      }

      if(n > 0){
        for(int i=0 ; i<nbranch ; i++){
          if(s == 1.0) br[i] = 0.0;
          else         br[i] = (next[i] == -2) ? 0.0 : br[i]/(1.0-s);
        }
        nbranch -= n;
      }

      /*** if half-life is longer than given years, change this into stable */
      if(lambda < long_lived){
        lambda = 0.0;
        nbranch = 0;
      }

      nuk[k0].setLambda(lambda);
      nuk[k0].mode.setNDM(nbranch);
      for(int i=0 ; i<nbranch ; i++) nuk[k0].mode.setData(i,br[i],next[i]);
//    cout << k0 << " " << z0 << " " << a0 << " " << m0 << " " << nbranch << endl;
    }
  }

  fp.close();

  return(nuktotal);
}


/**********************************************************/
/*      Cumulative Yield                                  */
/*          Starting with independent yields,             */
/*          cumulative yields are given by folllowing     */
/*          decay chains of all FPs until it converges.   */
/**********************************************************/
void cfyCumulativeYield(const int ntotal)
{
  int    nt  = 100;   // ad hoc number of total iteration, converge very quickly
  double eps = 1e-10; // convergence criteria

  double *y1 = new double [ntotal];
  double *y2 = new double [ntotal];

  for(int i=0 ; i<ntotal ; i++){
    y1[i] = nuk[i].yield;
    y2[i] = 0.0;
  }

  /*** repeat decay chain */
  for(int t=0 ; t<nt ; t++){

    /*** for all daughter nuclides, Y2 */
    for(int i=0 ; i<ntotal ; i++){
      y2[i] = 0.0;

      /*** for all given nuclides, which produce I-th nuclide */
      for(int j=0 ; j<ntotal ; j++){
        if(i == j) continue;

        /*** check if the daughter nucleus is I-th */
        for(int m=0 ; m<nuk[j].mode.getNDM() ; m++){
          if(nuk[j].mode.getNext(m) == i){
            /*** add Y1[j] times branching ratio to Y2[i] */
            y2[i] += y1[j] * nuk[j].mode.getBranch(m);
          }
        }
      }
    }

    /*** convergence test */
    double d = 0.0;
    for(int i=0 ; i<ntotal ; i++){
      double z = y2[i] + nuk[i].yield;
      if(z > 0.0) d += abs(y1[i]/z - 1.0);
      y1[i] = y2[i] + nuk[i].yield;
    }
    if(d < eps) break;
  }

  /*** overwrite the yield data by calculated cumulative yield */
  for(int i=0 ; i<ntotal ; i++) nuk[i].yield = y1[i];

  delete [] y1;
  delete [] y2;
}


/**********************************************************/
/*      Data Monitor                                      */
/**********************************************************/
void cfyCheckData(const int nuktotal)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);

  std::cout << "#  K0  Z0  A0  M0  Initial     Final       T1/2     ";
  std::cout << "  NBR  K1  Z1  A1  M1  BRatio" << std::endl;
  for(int k=0 ; k<nuktotal ; k++){

    double t2 = (nuk[k].getLambda() == 0.0) ? 0.0 : log(2.0)/nuk[k].getLambda();

    std::cout << std::setw( 5) << k;
    std::cout << std::setw( 4) << nuk[k].getZ() << std::setw(4) << nuk[k].getA() << std::setw(4) << nuk[k].getM();
    std::cout << std::setw(12) << std::setprecision(4) << nuk[k].initial;
    std::cout << std::setw(12) << std::setprecision(4) << nuk[k].yield;
    std::cout << std::setw(12) << std::setprecision(4) << t2;
    std::cout << std::setw( 4) << nuk[k].mode.getNDM() << std::endl;

    for(int i=0 ; i<nuk[k].mode.getNDM(); i++){
      std::cout << "                                                         ";
      int k1 = nuk[k].mode.getNext(i);
      std::cout << std::setw(4) << k1;
      if(k1 < 0) std::cout << " decay daughter nuclide not found" << std::endl;
      else{
        std::cout << std::setw(4) << nuk[k1].getZ() <<  std::setw(4) << nuk[k1].getA() << std::setw(4) << nuk[k1].getM();
        std::cout << std::setw(12) << std::setprecision(4) << nuk[k].mode.getBranch(i) << std::endl;
      }
    }
  }
}

/******************************************************************************/
/*  eclmc.cpp                                                                 */
/*        Monte Carlo simulation for patricle and gamma-ray emissions         */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "structur.h"
#include "eclipse.h"
#include "global.h"
#include "mt19937.h"

static const int MAX_CASCADE = 100;
//#define CALC_MONITOR


static int   eclMCSimulation (const int, int *, int **, const double, double ****, unsigned long ***, unsigned long **, unsigned long *, EXSpectra *);
static void  eclMCHistory (const int, int, int *, int *, unsigned long ***, unsigned long **, unsigned long *, const double);
static void  eclMCTally (const int, const int, int *, unsigned int, double **, unsigned long ***, unsigned long **, unsigned long *, EXSpectra *);

static const bool TAG_ON_FISSION = false;
static bool fission_event = false;


/**********************************************************/
/*      Start Monte Carlo Simulation                      */
/**********************************************************/
void eclMC(const int nmax, const int kmax, const int cmax, int *ktot, int **cidx, double ****ptbl, double **mtpl, const unsigned long nsim, double de, EXSpectra *dat)
{
  unsigned long ***hist;    // Monte Carlo history for particle and gamma emissions
  unsigned long  **npar;    // particle and gamma multiplicity
  unsigned long nevt[nmax]; // total events

  /*** allocate memory */
  hist = new unsigned long ** [nmax];
  npar = new unsigned long *  [nmax];
  for(int n=0 ; n<nmax ; n++){
    hist[n] = new unsigned long * [cmax];
    npar[n] = new unsigned long   [cmax];
    for(int c=0 ; c<cmax ; c++){
      hist[n][c] = new unsigned long [kmax];
    }
  }

  /*** clear arrays */
  for(int n=0 ; n<nmax ; n++){
    for(int c=0 ; c<cmax ; c++){
      npar[n][c] = 0;
      for(int k=0 ; k<kmax ; k++) hist[n][c][k] = 0;
    }
    nevt[n] = 0;
  }

  std::cout << std::setprecision(4) << std::setiosflags(std::ios::scientific);

  /*** perform simulation for Nsim times */
  unsigned long tmis = 0; // number of faild simulation
  for(unsigned long t=0 ; t<nsim ; t++){
#ifdef CALC_MONITOR
    if(t%10000000L==0) std::cerr << "+";
    else if(t%1000000L==0) std::cerr << "*";
#endif

    if( eclMCSimulation(cmax,ktot,cidx,de,ptbl,hist,npar,nevt,dat) < 0 ) tmis++;
  }

  /*** calculate multiplicity and spectra */
  eclMCTally(nmax,cmax,ktot,nsim-tmis,mtpl,hist,npar,nevt,dat);

  for(int n=0 ; n<nmax ; n++){
    for(int c=0 ; c<cmax ; c++) delete [] hist[n][c];
    delete [] hist[n];
    delete [] npar[n];
  }
  delete [] hist;
  delete [] npar;
}


/*************************************************/
/*  Simulation Main Part                         */
/*************************************************/
int eclMCSimulation(const int cmax, int *ktot, int **cidx, const double de, double ****ptbl, unsigned long ***hist, unsigned long **npar, unsigned long *nevt, EXSpectra *dat)
{
  int iout[MAX_CASCADE];    // particle and gamma indices at each cascade step
  int kout[MAX_CASCADE];    // emission energies at each cascade step
  for(int i=0 ; i<MAX_CASCADE ; i++) kout[i] = iout[i] = 0;

  /*** indices of current location */
  int  c0     = 0;    // nucleus
  int  i0     = 0;    // decay channel
  int  k0     = 0;    // excitation

  /*** start MC simulation */
  int    mc     = 0;    // Monte Carlo steps
  bool   repeat = true; // false, if no more decay
  fission_event = false;
  do{
    /*** fission check, if PE skip checking for c0=0 and k0=0 case */
    if(genrand_real3() < dat[c0].pfis[k0] ){
      fission_event = true;
      iout[mc] = 7;
      kout[mc] = 0;
      mc++;
      break;
    }

    /*** find the next location (c1,k1) by integrating the probability */
    int    c1 = 0;
    int    k1 = 0;
    double s  = 0.0;
    double r  = genrand_real3();
    for(i0=0 ; i0<cmax ; i0++){
      c1 = cidx[c0][i0];  if(c1 < 0) continue;

      /*** sum all probabilities until sum > random number */
      bool   loopout = false;
      for(k1=k0 ; k1<=ktot[c1] ; k1++){
        s += ptbl[c0][k0][i0][k1];
        if(s >= r){
          loopout = true;
          break;
        }
      }
      if(loopout) break;
    }

    /*** due to round-off, sum P<1.0 happens; discard this event */
    if(i0 == cmax) return(-1);

    /*** store decay path and emitted energy */
    iout[mc] = i0;
    kout[mc] = k1-k0;

    /*** if no further probability table found, or reached at the lowest */
    if(ptbl[c1][k1][0][k1] < 0.0) repeat = false;
    if(k1 == ktot[c1]) repeat = false;

    /*** for the next step */
    c0 = c1;
    k0 = k1;

    if((mc++) >= MAX_CASCADE){
      std::cerr << "max cascading exceeded at "<<c0<<" "<<k0<<" "<<c1<<" "<<k1<< std::endl;
      break;
    }

  }while(repeat);

  eclMCHistory(c0,mc,iout,kout,hist,npar,nevt,de);

  return(c0);
}


/*************************************************/
/*  History Array                                */
/*************************************************/
void eclMCHistory(const int c0, int mc, int *iout, int *kout, unsigned long ***hist, unsigned long **npar, unsigned long *nevt, const double de)
{
  static char pname[] = {'g','n','p','a','d','t','h','f',' '};

  /*** eliminate the last gamma if energy is zero */
  if( (iout[mc-1] == 0) && (kout[mc-1] == 0) ) mc--;
  
  /*** if terminated by fission */
  if(TAG_ON_FISSION and !fission_event) return;
  if(fission_event) mc--;

  /*** count number of gamma-rays,
       increment total history to the final residual nucleus */
  // int ng = 0;
  for(int i=0 ; i<mc ; i++){
    // if(iout[i] == 0) ng ++;

    hist[c0][iout[i]][kout[i]] ++;
    npar[c0][iout[i]] ++;
  }
  if(mc > 0) nevt[c0]++;

  if(prn.mchistory){
    std::cout << std::setw(5) << mc <<" : ";
    for(int i=0 ; i<mc ; i++){
      double ep = kout[i] * de;
      std::cout << std::setw(2) << pname[ iout[i] ] << std::setw(12) << ep;
    }
    if(fission_event) std::cout << std::setw(2) << pname[7] << std::setw(12) << 0.0;

    std::cout << std::endl;
  }
}


/*************************************************/
/*  Tally Events                                 */
/*************************************************/
void eclMCTally(const int nmax, const int cmax, int *ktot, unsigned int tsim, double **mtpl, unsigned long ***hist, unsigned long **npar, unsigned long *nevt, EXSpectra *dat)
{
  for(int n=0 ; n<nmax ; n++){

    double ns = (double)nevt[n];
    for(int c=0 ; c<cmax ; c++){
      for(int k=0 ; k<=ktot[n] ; k++) dat[n].spec[c][k] = (double)hist[n][c][k]/tsim;
      mtpl[n][c] = (ns>0.0) ? npar[n][c]/ns : 0.0;
    }
  }
}

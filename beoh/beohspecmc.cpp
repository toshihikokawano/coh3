/******************************************************************************/
/*  specmc.cpp                                                                */
/*        Monte Carlo for particle emission and gamma-ray cascading           */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "beoh.h"
#include "statmodel.h"
#include "beohstat.h"
#include "nucleus.h"
#include "global.h"
#include "mchf.h"
#include "mt19937.h"
#include "terminate.h"


static MCPoint mcStartPoint (Nucleus *);
static void    mcCompoundDecay (Pdata *, Transmission **, Transmission **, double **, Spectra *, MCPoint);
static void    mcGammaCascade(double *, MCPoint);
static void    mcWriteHistory ();

static inline void   mcSetContinuumLevel (const int, const int, const int, const int, MCPoint *);
static inline void   mcSetDiscreteLevel (const int, const int, MCPoint *);
static inline void   mcRecordGammaHistory (const bool, MCPoint, MCPoint, double *);
static inline void   mcRecordParticleHistory (const int, MCPoint, MCPoint, double *);
static inline void   mcAddSpectra (const double, double *);

static bool PerturbExcitationEnergy = true;
static bool PrintParent = false;
static bool CalcMonitor = true;

static int       np = 0;           // number of emitted gamma-rays and particles
static MCHistory hs[MAX_CASCADE];  // emitted energy and particle identifier
static Spectra   spcmc;            // particle emission spectra

#undef DEBUG_POPCHECK
#ifdef DEBUG_POPCHECK
static void specPopCheck (Nucleus *);
#endif

/**********************************************************/
/*      Monte Carlo Decay of Compound Nucleus             */
/**********************************************************/
void    beohspectraMC(Pdata *pdt, Transmission **tc, Transmission **td, double **tg, Spectra *spc, const unsigned long nsim)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(5);

  /*** copy initial population distribution */
  Nucleus nc0;
  nc0.memalloc(MAX_LEVELS,MAX_ENERGY_BIN);
  nc0.memclear();

  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    nc0.excitation[k] = ncl[0].excitation[k];
    for(int j=0 ; j<MAX_J ; j++) nc0.pop[k][j] = ncl[0].pop[k][j];
  }
  for(int i=0 ; i<MAX_LEVELS ; i++) nc0.lpop[i] = ncl[0].lpop[i];
  for(int id=0 ; id<MAX_CHANNEL ; id++) nc0.cdt[id] = nc0.cdt[id];

  nc0.ncont = ncl[0].ncont;
  nc0.ndisc = ncl[0].ndisc;
  nc0.jmax  = ncl[0].jmax;
  nc0.de    = ncl[0].de;

  /*** spectrum by Monte Carlo, we don't use spc becuase it will be changed by compounddecay */
  spcmc.memalloc(MAX_CHANNEL,MAX_ENERGY_BIN,MAX_LEVELS);
  spcmc.memclear();

  MCPoint minit;

  /*** set MC starting point, must be done before NSIM look, as pop will be changed */
  if(nc0.ncont == 0) mcSetDiscreteLevel(0,nc0.ndisc-1,&minit);

  /*** main MC event loop */
  for(unsigned long t=0 ; t<nsim ; t++){
    np = 0; // reset history

    /*** MC starts from discrete level */
    if(nc0.ncont == 0){
      mcGammaCascade(spcmc.cn[gammaray],minit);
    }
    /*** general case */
    else{
      minit = mcStartPoint(&nc0);
      mcCompoundDecay(pdt,tc,td,tg,spc,minit);
    }

    /*** print event history */
    if(prn.mchistory) mcWriteHistory();
    if(CalcMonitor){
      if(t%1000L == 0L) std::cerr << "+";
      else if(t%100L == 0L) std::cerr << "*";
    }
  }

  /*** normalize total particle emission spectra */
  double f = 1.0 / (double)nsim / nc0.de;
  for(int id=0 ; id<MAX_CHANNEL ; id++){
    if(nc0.cdt[id].status && id < spcmc.getCsize()){
      for(int k=0 ; k<spc->getNsize() ; k++) spcmc.cn[id][k] *= f;
    }
  }
  specCumulativeSpectra(spcmc.getCsize(),spcmc.getNsize(),spcmc.cn,&nc0);


  spcmc.memfree();
}


/**********************************************************/
/*      Look for Monte Carlo Starting Point               */
/**********************************************************/
MCPoint mcStartPoint(Nucleus *n)
{
  double r = genrand_real3();
  double s = 0.0;

  /*** look for continuum bin */
  int  p = 0, k = 0, j = 0;
  bool loopout = false;
  for(k=0 ; k<n->ncont ; k++){
    for(j=0 ; j<=n->jmax ; j++){
      s += n->pop[k][j].even; if(s >= r){ p =  1; loopout = true; break; }
      s += n->pop[k][j].odd;  if(s >= r){ p = -1; loopout = true; break; }
    }
    if(loopout) break;
  }

  MCPoint minit;
  if(loopout){
    minit.set(0,p,j,k,n->excitation[k]);
    minit.setContinuum();
  }
  else{
    /*** this part is for future extension */
    for(k=0 ; k<n->ndisc ; k++){
      s += n->lpop[k];
      if(s >= r) break;
    }
    mcSetDiscreteLevel(0,k,&minit);
  }

  return minit;
}


/**********************************************************/
/*      Monte Carlo Particle Emission From Continuum      */
/**********************************************************/
void mcCompoundDecay(Pdata *pdt, Transmission **tc, Transmission **td, double **tg, Spectra *spc, MCPoint minit)
{
  bool tstat[MAX_CHANNEL];
  MCPoint mcurr,mnext;
  Nucleus *n0, *n1;

  mcurr = minit;

  while(true){

    /*** parent compund and MC start point */
    int c0 = mcurr.getNucleus();
    int k0 = mcurr.getEk();
    int j0 = mcurr.getSpin();
    int p0 = mcurr.getParity();
    n0 = &ncl[c0];

    /*** Check if channel is closed */
    tstat[gammaray] = true;
    for(int id=1 ; id<MAX_CHANNEL ; id++){
      tstat[id] = n0->cdt[id].status;
      if(tstat[id] && (n0->excitation[k0] < n0->cdt[id].binding_energy)) tstat[id] = false;
    }

    /*** clear all pop arrays */
    beohClearInitialPopulation(n0);

    for(int id=0 ; id<MAX_CHANNEL ; id++){
      if(!tstat[id]) continue;
      n1 = &ncl[n0->cdt[id].next];
      beohClearInitialPopulation(n1);
    }
    /*** and set 1.0 at the starting point */
    if(p0 > 0) n0->pop[k0][j0].even = 1.0;
    else       n0->pop[k0][j0].odd  = 1.0;

    /*** Store gamma-ray transmission data into array */
    statStoreContinuumGammaTransmission(k0,tg,n0);

    for(int id=1 ; id<MAX_CHANNEL ; id++){
      if(!tstat[id]) continue;
      /*** Store particle transmission data into array */
      statStoreContinuumTransmission(c0,n0->excitation[k0],pdt,tc);
      /*** Transmission coefficients for discrete levels */
      statStoreDiscreteTransmission(c0,n0->excitation[k0],pdt,td);
    }

    /*** residual population gives the probability to move to the next point */
    specCompoundDecay(c0,k0,tc,td,tg,spc);
#ifdef DEBUG_POPCHECK
    specPopCheck(&ncl[c0]);
#endif

    /*** move to the next point */
    double r = genrand_real3();
    double s = 0.0;

    int  c1 = 0, k1 = 0, k2 = 0, j1 = 0, id = 0;
    int loopout = 0; // 0: direct to g.s., 1: cont to cont, 2: cont to disc
    for(id=0 ; id<MAX_CHANNEL ; id++){
      if(!tstat[id]) continue;
      n1 = &ncl[ c1 = n0->cdt[id].next ];

      /*** continuum to continuum */
      for(k1=k0+1 ; k1<n1->ncont ; k1++){
        for(j1=0 ; j1<=n1->jmax ; j1++){
          s += n1->pop[k1][j1].even;
          if(s >= r){ mcSetContinuumLevel(c1,k1,j1, 1,&mnext); loopout = 1; break; }
          s += n1->pop[k1][j1].odd;
          if(s >= r){ mcSetContinuumLevel(c1,k1,j1,-1,&mnext); loopout = 1; break; }
        }
        if(loopout) break;
      }
      if(loopout) break;

      /*** continuum to discrete */
      for(k2=0 ; k2<n1->ndisc ; k2++){
        s += n1->lpop[k2];
        if(s >= r){ mcSetDiscreteLevel(c1,k2,&mnext); loopout = 2; break; }
      }
      if(loopout) break;
    }

//    std::cout << loopout << " " << k1<< " "<<j1<<" " <<k2<< " " << mnext.getEk() <<" " <<mnext.getEnergy() << std::endl;

    /*** in the case no transition happened, means went to the ground state */
    if(loopout == 0){
      id = gammaray;
      k2 = 0;
    }

    /*** save history */
    if(id == 0) mcRecordGammaHistory(false,mcurr,mnext,spcmc.cn[gammaray]);
    else        mcRecordParticleHistory(id,mcurr,mnext,spcmc.cn[id]);

    /*** proceed to the next step */
    mcurr = mnext;

    if(loopout == 2){   // the final state is in a discrete level
      mcGammaCascade(spcmc.cn[gammaray],mcurr);
      break;
    }
  }
}


/**********************************************************/
/*      Monte Carlo Cascading From Discrete Level         */
/**********************************************************/
void mcGammaCascade(double *spc, MCPoint minit)
{
  MCPoint mcurr, mnext;

  /*** set current discrete level */
  mcurr = minit;
  int c = mcurr.getNucleus();

  while(1){
    /*** current level number */
    int i0 = mcurr.getEk();

    /*** move to the next discrete level by branching ratios */
    double r = genrand_real3();
    double s = 0.0;

    for(int j=0 ; j<ncl[c].lev[i0].ngamma ; j++){
      int i1 = ncl[c].lev[i0].fstate[j];
      s += ncl[c].lev[i0].branch[j];

      if(s >= r){
        /*** if Internal Concversion happens, store no gamma infor */
        bool ic = false;
        if( opt.internalconversion && (genrand_real3() > ncl[c].lev[i0].gratio[j]) ) ic = true;

        /*** point next state */
        mcSetDiscreteLevel(c,i1,&mnext);

        /*** save history */
        mcRecordGammaHistory(ic,mcurr,mnext,spc);
        mcurr = mnext;
        break;
      }
    }
    if(i0 == 0) break;
  }
}


/**********************************************************/
/*      Save Level Information in Object                  */
/**********************************************************/
static inline void  mcSetContinuumLevel(const int c, const int k, const int j, const int p, MCPoint *mcp){
  double ex = ncl[c].excitation[k];

  /*** perturbation to the final state energy within its energy bin */
  if(PerturbExcitationEnergy){
    double emax = (k == 0) ? ex : (ncl[c].excitation[k-1] +  ncl[c].excitation[k])*0.5;
    double emin = (ncl[c].excitation[k+1] +  ncl[c].excitation[k])*0.5;
    if(emin < ncl[c].lev[ncl[c].ndisc].energy) emin = ncl[c].lev[ncl[c].ndisc].energy;
    ex = emin + (emax-emin) * genrand_real2();
  }

  mcp->set(c,p,j,k,ex);
  mcp->setContinuum();
}

static inline void  mcSetDiscreteLevel(const int c, const int i, MCPoint *mcp){
  mcp->set(c,ncl[c].lev[i].parity,(int)ncl[c].lev[i].spin,i,ncl[c].lev[i].energy);
  mcp->setDiscrete();
}


/**********************************************************/
/*      Record Decay Gamma-Ray Energies                   */
/**********************************************************/
static inline void mcRecordGammaHistory(const bool ic, MCPoint m0, MCPoint m1, double *spc)
{
  double ecms = m0.getEnergy() - m1.getEnergy();
  if(ecms < 0.0) return;

  if(ic) hs[np++].set(unknown , ecms, m0);
  else   hs[np++].set(gammaray, ecms, m0);
//std::cout << m0.getEnergy() <<" " <<  m1.getEnergy() <<  " " << ecms << " " << np << " " << ic <<  std::endl;

  /*** terminate if too many MC loops */
  if(np >= MAX_CASCADE){
    message << "too many particles in Monte Carlo statistical decay";
    cohTerminateCode("mcGammaStoreHistory");
  }

  if(!ic) mcAddSpectra(ecms,spc);
}


/**********************************************************/
/*      Record Decay Particle Energies                    */
/**********************************************************/
static inline void mcRecordParticleHistory(const int id, MCPoint m0, MCPoint m1, double *spc)
{
  /*** subtract separation energy */
  double ecms = m0.getEnergy() - m1.getEnergy() - ncl[m0.getNucleus()].cdt[id].binding_energy;
  if(ecms < 0.0) return;

  hs[np++].set(id, ecms, m0);

  /*** terminate if too many MC loops */
  if(np >= MAX_CASCADE){
    message << "too many particles in Monte Carlo statistical decay";
    cohTerminateCode("mcRecordParticleHistory");
  }

  mcAddSpectra(ecms,spc);
}


/**********************************************************/
/*      Add Event to Spectrum                             */
/**********************************************************/
void mcAddSpectra(const double ecms, double *spc)
{
  if(ecms <= 0.0) return;
  int k = specFindEnergyBin(ecms,ncl[0].de);
  if(k > 0) spc[k] += 1.0;
}


/**********************************************************/
/*      Write Decay History                               */
/**********************************************************/
void mcWriteHistory()
{
//double  etot = 0.0;
  char    name[9] = {'g','n','p','a','d','t','h','f','e'};

  /*** print out particle energies */
  std::cout << std::setw(8) << np << " : ";
  for(int i=0 ; i<np ; i++){
    std::cout << std::setw(4)  << name[hs[i].getIndex()];
    std::cout << std::setw(13) << hs[i].getEnergy();
//  etot += hs[i].getEnergy();
  }
  std::cout << std::endl;

  if(PrintParent){
    /*** print out parent state */
    std::cout << "           ";
    for(int i=0 ; i<np ; i++){
      MCPoint m = hs[i].getMCPoint();
      std::cout << std::setw(3)  << m.getSpin();
      std::cout << std::setw(1)  << ((m.getParity() == -1) ? '-' : '+');
      std::cout << std::setw(13) << m.getEnergy();
    }
    std::cout << std::endl;
  }
}


/**********************************************************/
/*      For Monitoring Only                               */
/**********************************************************/
#ifdef DEBUG_POPCHECK
void specPopCheck(Nucleus *n)
{
  double s1 = 0.0, s2 = 0.0;

  std::cout << "#  ";
  std::cout << std::setw(3) << std::setfill('0') << n->za.getZ() << '-';
  std::cout << std::setw(3) << std::setfill('0') << n->za.getA() << std::endl;
  std::cout << std::setfill(' ');

  for(int k=0 ; k<n->ncont ; k++){
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout << std::setw(7) << std::setprecision(3) <<  n->excitation[k];
    std::cout.setf(std::ios::scientific, std::ios::floatfield);

    double s3 = 0.0;
    for(int j=0 ; j<10 ; j++) {
      s1 += n->pop[k][j].even+n->pop[k][j].odd;
      s3 += n->pop[k][j].even+n->pop[k][j].odd;
      std::cout << std::setprecision(2) << std::setw(9) << n->pop[k][j].even;
    }
    std::cout << std::endl;

    std::cout << "       ";
    for(int j=0 ; j<10 ; j++) {
      std::cout << std::setprecision(2) << std::setw(9) << n->pop[k][j].odd;
    }
    std::cout << std::setprecision(4) << std::setw(11) << s3 << std::endl;
  }

  s2 = 0.0;
  for(int k=0;k<n->ndisc;k++){
    s2 += n->lpop[k];
    std::cout << std::setw(5) << k
         << std::setprecision(4) << std::setw(13) << n->lev[k].energy
         << std::setw(13) << n->lpop[k] << std::endl;
  }

  std::cout << "# SUMcont " << s1 << std::endl;
  std::cout << "# SUMdisc " << s2 << std::endl;
  std::cout << "# SUMall  " << s1+s2 << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}
#endif

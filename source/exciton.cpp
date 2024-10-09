/******************************************************************************/
/*  exciton.cpp                                                               */
/*        exciton model for precompound reaction                              */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "physicalconstant.h"
#include "structur.h"
#include "nucleus.h"
#include "statmodel.h"
#include "exciton.h"
#include "global.h"
#include "output.h"
#include "terminate.h"

static void   preqParameterSetting (System *, Preeq *, SPdens *);
static void   preqExcitonCalc (System *, Pdata *, double, Exconf *, Preeq *, Transmission **, Transmission **, double **);
static void   preqTransitionRate (System *, Pdata *, Exconf *, Preeq *, Transmission **);
static void   preqAddPopulation (double, Pdata *, Exconf *, Preeq *, Transmission **, Transmission **, double **);
static void   preqPrepareExcitonParameter (System *, Exconf *, Preeq *);
static double preqEmissionRate (Pdata *,Transmission **, Preeq *,Exconf *);
static void   preqLifeTime (const int, const int, Preeq *, Exconf *);
static void   preqLevelPopulation (const int, Transmission *, double *, Nucleus *);
static void   preqAllocateMemory (void);
static void   preqDeleteAllocated (void);
static void   preqClearMemory (void);

static inline double preqTransmissionAverage(int, int, double *);

static double  **wkc;  // continuum particle emission rate
static double  **pop;  // occupation probability
static double  **wk;   // integrated emission rate
static double  **tau;  // lifetime of exciton state
static double  **tap;  // lifetime only for exchange interaction
static double  **tzp;  // lambda z+
static double  **tnp;  // lambda n+
static double  **tz0;  // lambda z0
static double  **tn0;  // lambda n0

#undef LEVEL_POPULATION
#undef DEBUG

//#define FKK_SPIN_DISTRIBUTION
#undef FKK_SPIN_DISTRIBUTION
#ifdef FKK_SPIN_DISTRIBUTION
static double  preqSpinDistribution (double, double);
#endif

/**********************************************************/
/*      Exciton Model for Pre-Equilibrium Reaction        */
/**********************************************************/
void preqExcitonModel(System *sys, Pdata *pdt, Transmission **tc, Transmission **td, Spectra *spc)
{
  Preeq    prq;               // preequilibrium parameters
  Exconf   exc;               // exciton configuration
  SPdens   spd[MAX_CHANNEL];  // single-particle state densities

  /*** Clear spc array */
  spc->memclear("pe");


//---------------------------------------
//      Memory Allocation

  preqAllocateMemory();


//---------------------------------------
//      Parameter Setting

  preqParameterSetting(sys,&prq,spd);


//---------------------------------------
//      Main Calculation Part

  /*** photo induced reaction case */
  if(sys->incident.pid == gammaray){

    double sigcn = crx.reaction;
    double sigQDB = sigcn * photoQDratio();
    double sigGDR = sigcn * (1.0 - photoQDratio());

    /*** Total CN formation cross section for GDR, scaled by phase space factor
         R = omega_B(1p,1h)/omega(1p,1h) = B/(E-delta) */
    double bn = ncl[0].cdt[neutron].binding_energy;
    double e0 = prq.ex_total - spd[neutron].pairing;
    double r  = (e0 > 0.0) ? bn/e0 : 0.0;
    if((0.0 < r) && (r <= 1.0)) sigGDR *= r;

    /*** Initial number of particle/holes
         GDR reaction: start with (1p,1h) in n-shell
         quasi-deuteron reaction: (1p,1h) in both shells */
    Exconf exGDR(0,0,1,1);
    Exconf exQDB(1,1,1,1);

    /*** GDR */
    preqExcitonCalc(sys,pdt,sigGDR,&exGDR,&prq,tc,td,spc->pe);

    /*** QD */
    preqExcitonCalc(sys,pdt,sigQDB,&exQDB,&prq,tc,td,spc->pe);


//    double sigcn = crx.reaction;
//    Exconf ex0(0,0,2,1);
//    preqExcitonCalc(sys,pdt,sigcn,&ex0,&prq,tc,td,spc->pe);
  }
  else{
    double sigcn = crx.reaction;

    /*** Initial number of particle/holes ((z+n)p, 0h) */
    Exconf ex0(sys->incident.za.getZ(),0,sys->incident.za.getN(),0);

    preqExcitonCalc(sys,pdt,sigcn,&ex0,&prq,tc,td,spc->pe);
  }


//---------------------------------------
//      Free Allocated Memory

  preqDeleteAllocated();
}


/**********************************************************/
/*      Set Model Parameters in Each Channel              */
/**********************************************************/
void preqParameterSetting(System *sys, Preeq *prq, SPdens *spd)
{
  /*** Total excitation energy */
  prq->ex_total = sys->ex_total;

  /*** Parameters for state density */
  for(int id=0 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;
    int i = ncl[0].cdt[id].next;

    /*** Default single particle state density */
    spd[id].gz         = preqSingleStateDensity(id, ncl[i].za.getZ());
    spd[id].gn         = preqSingleStateDensity(id, ncl[i].za.getN());

    /*** Pairing energy */
    spd[id].pairing    = preqPairingDelta(ncl[i].za.getZ(),ncl[i].za.getA());

    /*** Potential well depth */
    spd[id].well_depth = preqPotentialDepth(sys->lab_energy,sys->incident.za.getZ(),ncl[i].za.getA());
  }
  prq->spd = spd;
}


/**********************************************************/
/*      Exciton Model Main Calculation                    */
/**********************************************************/
void preqExcitonCalc(System *sys, Pdata *pdt, double sigcn, Exconf *ex0, Preeq *prq, Transmission **tc, Transmission **td, double **spc)
{
  /*** Transition Rates for Each Configuration */
  preqClearMemory();
  preqTransitionRate(sys,pdt,ex0,prq,tc);

  /*** Ocupation Probability for Each Configuration */
  preqPOccupation(pop,tau,tap,tzp,tnp,tz0,tn0);

  /*** Add to Population */
  preqAddPopulation(sigcn,pdt,ex0,prq,tc,td,spc);

  /*** Exciton Model Parameter Output */
  if(prn.preeqparm) preqPrepareExcitonParameter(sys,ex0,prq);
}


/**********************************************************/
/*      Decay and Emission Rates at Each Configuration    */
/**********************************************************/
void preqTransitionRate(System *sys, Pdata *pdt, Exconf *ex0, Preeq *prq, Transmission **tc)
{
  Exconf exc; // exciton configuration

  /*** Index i,j set to the number of holes */
  for(int i=0 ; i<MAX_PREEQ ; i++){
    exc.zp = ex0->zp + i;
    exc.zh = ex0->zh + i;

    for(int j=0 ; j<MAX_PREEQ ; j++){
      exc.np = ex0->np + j;
      exc.nh = ex0->nh + j;

      /*** Total exciton state density omega(p,h,Ex) */
      prq->omega_total = preqStateDensity(prq->ex_total,prq->spd,&exc);
      if(prq->omega_total == 0.0) continue;

      /*** Particle emission rate */
      wk[i][j] = preqEmissionRate(pdt,tc,prq,&exc);

      /*** Effective squared matrix elements */
      preqMSquare(sys->target.getA()+sys->incident.za.getA(),exc.sum(),prq);

      /*** Lambdas for all configurations */
      preqLifeTime(i,j,prq,&exc);
    }
  }
}


/**********************************************************/
/*      Decay and Emission Rates at Each Configuration    */
/**********************************************************/
void preqAddPopulation(double sigcn, Pdata *pdt, Exconf *ex0, Preeq *prq, Transmission **tc, Transmission **td, double **spc)
{
  Exconf exc;

  for(int i=0 ; i<MAX_PREEQ ; i++){
    exc.zp = ex0->zp + i;
    exc.zh = ex0->zh + i;

    for(int j=0 ; j<MAX_PREEQ ; j++){
      exc.np = ex0->np + j;
      exc.nh = ex0->nh + j;

      prq->omega_total = preqStateDensity(prq->ex_total,prq->spd,&exc);
      if(prq->omega_total == 0.0) continue;

      preqEmissionRate(pdt,tc,prq,&exc);

      double taup = tau[i][j] * pop[i][j];

      for(int id=1 ; id<MAX_CHANNEL ; id++){
        if(!ncl[0].cdt[id].status) continue;
        int ip = ncl[0].cdt[id].next;
        for(int k=1 ; k<ncl[ip].ntotal ; k++) spc[id][k] += sigcn * wkc[id][k] * taup;
      }
    }
  }

  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;
    preqExcitonPopulation(id,tc[id],td[id],spc[id]);
  }
}


/**********************************************************/
/*      Copy Parameter in Temporal Strage and Call Output */
/**********************************************************/
void preqPrepareExcitonParameter(System *sys, Exconf *ex0, Preeq *prq)
{
  /*** wkc, tz0, and tn0 are used for temporal space for printing */
  for(int id=0 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;
    wkc[id][0] = prq->spd[id].gn;
    wkc[id][1] = prq->spd[id].gz;
    wkc[id][2] = prq->spd[id].pairing;
    wkc[id][3] = prq->spd[id].well_depth;
  }

  for(int i=0 ; i<MAX_PREEQ ; i++){
    for(int j=0 ; j<MAX_PREEQ ; j++) tz0[i][j] = tn0[i][j] = 0.0;
  }

  Exconf exc;
  for(int i=0 ; i<MAX_PREEQ ; i++){
    exc.zp = ex0->zp + i;
    exc.zh = ex0->zh + i;

    for(int j=0 ; j<MAX_PREEQ ; j++){
      exc.np = ex0->np + j;
      exc.nh = ex0->nh + j;
        
      prq->omega_total = preqStateDensity(prq->ex_total,prq->spd,&exc);
      preqMSquare(sys->target.getA()+sys->incident.za.getA(),exc.sum(),prq);
      tz0[i][j] = prq->m2zz;
      tn0[i][j] = prq->m2nz;
    }
  }

  outExcitonParameter(wkc,tz0,tn0);
}


/**********************************************************/
/*      Particle Emission Rate                            */
/*      calculated in the unit of h-bar [MeV sec]         */
/**********************************************************/
double preqEmissionRate(Pdata *pdt, Transmission **tc, Preeq *prq, Exconf *exc)
{
  double c  = AMUNIT/(HBARSQ*VLIGHTSQ*NORM_FACT*PI*PI);
  double wt = 0.0;

  for(int id=1 ; id<MAX_CHANNEL ; id++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) wkc[id][k] = 0.0;
    if(!ncl[0].cdt[id].status) continue;

    int j = ncl[0].cdt[id].next;

    Exconf res(exc->zp - pdt[id].za.getZ(),exc->zh,
               exc->np - pdt[id].za.getN(),exc->nh);
    if(res.zp < 0 || res.np < 0) continue;

    double rm = ncl[j].mass * pdt[id].mass
              /(ncl[j].mass + pdt[id].mass);

    double x  = c * rm * (pdt[id].spin2+1.0)/prq->omega_total;

    for(int k=1 ; k<ncl[j].ntotal ; k++){
      double omega = preqStateDensity(ncl[j].excitation[k],&prq->spd[id],&res);
      double coll  = preqCollectiveEnhancement(ncl[j].excitation[k],&res);

      wkc[id][k]   = x * tc[id][k].ecms * tc[id][k].sigr * omega * coll;
      wt += wkc[id][k] * ncl[j].de;
#ifdef DEBUG
      if( (res.zp + res.np + res.zh + res.nh) == 2 ){
        std::cout << std::setw(2) << exc->zp << std::setw(2) << exc->zh << std::setw(2) << exc->np << std::setw(2) << exc->nh;
        std::cout << std::setw(4) << id << std::setw(4) << k;
        std::cout << std::setw(11) << tc[id][k].ecms << std::setw(11) << ncl[j].excitation[k];
        std::cout << std::setw(11) << omega << std::setw(11) << coll << std::setw(11) << prq->omega_total;
        std::cout << std::setw(11) << tc[id][k].sigr << std::setw(11) << wkc[id][k]/HBAR << std::endl;
      }
#endif
    }
  }

  return(wt);
}


/**********************************************************/
/*      Lifetimes Tau of Exciton State                    */
/*      calculated in the unit of h-bar [MeV sec]         */
/**********************************************************/
void preqLifeTime(const int i, const int j, Preeq *prq, Exconf *exc)
{
  tzp[i][j] = preqLambdaPlusZ(prq,exc);
  tnp[i][j] = preqLambdaPlusN(prq,exc);
  tz0[i][j] = preqLambdaZeroZ(prq,exc);
  tn0[i][j] = preqLambdaZeroN(prq,exc);

  double d1 = tzp[i][j] + tnp[i][j] + tz0[i][j] + tn0[i][j] + wk[i][j];
  double d2 = tzp[i][j] + tnp[i][j]                         + wk[i][j];
  tau[i][j] = (d1 > 0.0) ? 1.0/d1 : 0.0;
  tap[i][j] = (d2 > 0.0) ? 1.0/d2 : 0.0;
}


/**********************************************************/
/*      Pre-Equilibrium Population Increment              */
/*      Ad-Hoc Angular Momentum Distribution              */
/**********************************************************/
void  preqExcitonPopulation(const int id, Transmission *tc, Transmission *td, double *spc)
{
  int i = ncl[0].cdt[id].next;

  double sdist[MAX_J];

  /*** Distribute pre-equilibrium cross section to each spin state */
  /*** This part emulates compound J-distributions at each excitation */
  for(int k=1 ; k<ncl[i].ncont ; k++){

    for(int j=0 ; j<MAX_J ; j++) sdist[j] = 0.0;

    int jmax0 = 0;
    for(int jcn=0 ; jcn<ncl[0].jmax ; jcn++){
      for(int lp=0 ; lp<=tc[k].lmax ; lp++){

        double t = preqTransmissionAverage(ncl[0].cdt[id].spin2,lp,&tc[k].tran[3*lp]);

        int jmin = std::abs(jcn-lp);
        int jmax =          jcn+lp ; if(jmax > MAX_J-1) jmax = MAX_J-1;
        if(jmax>jmax0) jmax0 = jmax;

        for(int j=jmin ; j<=jmax ; j++){
          sdist[j] += t*(ncl[i].density[k][j].even + ncl[i].density[k][j].odd)
                       *(ncl[0].pop[0][jcn].even   + ncl[0].pop[0][jcn].odd  );
        }
      }
    }

#ifdef FKK_SPIN_DISTRIBUTION
    double ecms = ncl[i].excitation[0];
    double sig2 = preqSpinDistribution(ecms,tc[k].ecms);

    for(int j=0 ; j<=jmax0 ; j++){
      sdist[j] = (j + 0.5)/sig2 * exp(-(j+0.5)*(j+0.5)/(2*sig2));
    }
#endif

    double s = 0.0;
    for(int j=0 ; j<=jmax0 ; j++) s+=sdist[j];
    if(s > 0.0){
      for(int j=0 ; j<=jmax0 ; j++){
        double x = sdist[j]/s * 0.5;
        ncl[i].pop[k][j].even += spc[k] * x * ncl[i].de;
        ncl[i].pop[k][j].odd  += spc[k] * x * ncl[i].de;
      }
    }
  }

  if(ncl[i].ndisc > 1) preqLevelPopulation(id,td,spc,&ncl[i]);
/*
  for(int k=ncl[i].ncont ; k<ncl[i].ntotal ;  k++){
    if(spc[k]==0.0) continue;
    double e1 = ncl[i].excitation[k+1];
    double e2 = ncl[i].excitation[k  ];
    int m=0;
    for(int n=1 ; n<ncl[i].ndisc ; n++){
      if( (e1<=ncl[i].lev[n].energy) && (ncl[i].lev[n].energy<e2) ){
        m++;
      }
    }
    if(m>0){
      for(int n=1 ; n<ncl[i].ndisc ; n++){
        if( (e1<=ncl[i].lev[n].energy) && (ncl[i].lev[n].energy<e2) ){
          ncl[i].lpop[n] += spc[k] * ncl[i].de /(double)m;
        }
      }
    }
  }
*/
}


/**********************************************************/
/*      Pre-Equilibrium Population to the Levels          */
/**********************************************************/
void preqLevelPopulation(const int id, Transmission *td, double *spc, Nucleus *n)
{
  /*** Distribute the sum of spectrum to the open levels,
       weighted by the inverse reaction cross section */

  double spcsum = 0.0;
  for(int k=n->ncont ; k<n->ntotal ;  k++) spcsum += spc[k];
  if(spcsum == 0.0) return;

  int    m = 0;    // number of open levels
  double w = 0.0;  // sum of inverse reaction cross section
  for(int i=1 ; i<n->ndisc ; i++){
    if(td[i].lmax <= 0) continue;
    m ++;
    w += td[i].sigr;
  }

  /*** distribute spcsum to the levels,
       in the case of charged particles, no SigR weight,
       because they are almost zero */
  w = (id == 1) ? w*m : m;

#ifdef LEVEL_POPULATION
  if(m > 0){
    for(int i=1 ; i<n->ndisc ; i++){
      if(td[i].lmax <= 0) continue;
      if(id == 1) n->lpop[i] += spcsum * n->de * td[i].sigr / w;
      else        n->lpop[i] += spcsum * n->de              / w;
    }
  }
#endif
}

/**********************************************************/
/*      Spin-Averaged Transmission Coefficients           */
/**********************************************************/
double preqTransmissionAverage(int spin2, int l, double *tj)
{
  double t=0.0;

  if(spin2==0)       t  = tj[0];
  else if(spin2==1)  t  = ( (l+1.0)*tj[0] + l *tj[1] )/(2.0*l+1.0);
  else if(spin2==2)  t  = ( (2.0*l+3.0)*tj[0]
                           +(2.0*l+1.0)*tj[1]
                           +(2.0*l-1.0)*tj[2] )/(6.0*l+3.0);
  return (t);
}


/**********************************************************/
/*      Total Pre-Equilibrium Emission Cross Section      */
/**********************************************************/
double preqTotalPopulation()
{
  double s = 0.0;

  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;
    int    i  = ncl[0].cdt[id].next;

    for(int k=1 ; k<ncl[i].ncont ; k++){
      for(int j=0 ; j<MAX_J ; j++){
        s += ncl[i].pop[k][j].even + ncl[i].pop[k][j].odd;
      }
    }
    for(int k=1 ;  k<ncl[i].ndisc ; k++) s += ncl[i].lpop[k];
  }

  return(s);
}


/**********************************************************/
/*      Renormalize Pre-Equilibrium Emission              */
/**********************************************************/
void preqRenormalizePopulation(const double f)
{
  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;
    int    i  = ncl[0].cdt[id].next;

    for(int k=1 ; k<ncl[i].ncont ; k++){
      for(int j=0 ; j<ncl[i].jmax ; j++){
        ncl[i].pop[k][j].even *= f;
        ncl[i].pop[k][j].odd  *= f;
      }
    }
    for(int k=1 ;  k<ncl[i].ndisc ; k++) ncl[i].lpop[k] *= f;
  }
}


/**********************************************************/
/*      Pre-Equilibrium Particle Emission Spectra         */
/**********************************************************/
double preqExcitonSpectra(Spectra *spc)
{
  double sig[MAX_CHANNEL];

  crx.preeq = 0.0;
  spc->memclear("pe");

  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;
    int    i  = ncl[0].cdt[id].next;
    double e0 = ncl[0].excitation[0] - ncl[0].cdt[id].binding_energy;

    sig[id]=0.0;
    for(int k=1 ; k<ncl[i].ncont ; k++){
      double x=0.0;
      for(int j=0 ; j<MAX_J ; j++){
        x += ncl[i].pop[k][j].even + ncl[i].pop[k][j].odd;
      }
      spc->pe[id][k] = x/ncl[i].de;
      sig[id] += x;
    }

    for(int k=1 ;  k<ncl[i].ndisc ; k++){
      double e1 = e0 - ncl[i].lev[k].energy;
      int    kp = 0;
      for(kp=0 ; kp<MAX_ENERGY_BIN ; kp++){
        if( (ncl[i].de*kp<= e1) && (e1 < ncl[i].de*(kp+1)) ){
           break;
        }
      }
      if(kp < MAX_ENERGY_BIN){
        spc->pe[id][kp] += ncl[i].lpop[k]/ncl[i].de;
        sig[id]         += ncl[i].lpop[k];
      }
    }
    crx.preeq += sig[id];
  }

  return(crx.preeq);
}


/**********************************************************/
/*      Allocate Working Arrays                           */
/**********************************************************/
void preqAllocateMemory()
{
  try{
    wkc = new double * [MAX_CHANNEL];
    for(int j=0 ; j<MAX_CHANNEL ; j++) wkc[j] = new double [MAX_ENERGY_BIN];

    pop = new double * [MAX_PREEQ];
    wk  = new double * [MAX_PREEQ];
    tzp = new double * [MAX_PREEQ];
    tnp = new double * [MAX_PREEQ];
    tz0 = new double * [MAX_PREEQ];
    tn0 = new double * [MAX_PREEQ];
    tau = new double * [MAX_PREEQ];
    tap = new double * [MAX_PREEQ];
    for(int j=0 ; j<MAX_PREEQ ; j++){
      pop[j] = new double [MAX_PREEQ];
      wk [j] = new double [MAX_PREEQ];
      tzp[j] = new double [MAX_PREEQ];
      tnp[j] = new double [MAX_PREEQ];
      tz0[j] = new double [MAX_PREEQ];
      tn0[j] = new double [MAX_PREEQ];
      tau[j] = new double [MAX_PREEQ];
      tap[j] = new double [MAX_PREEQ];
    }
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("preqAllocateMemory");
  }
}


/**********************************************************/
/*      Clear Working Arrays                              */
/**********************************************************/
void preqClearMemory()
{
  for(int i=0 ; i<MAX_CHANNEL ; i++){
    for(int j=0 ; j<MAX_ENERGY_BIN ; j++) wkc[i][j] = 0.0;
  }

  for(int i=0 ; i<MAX_PREEQ ; i++){
    for(int j=0 ; j<MAX_PREEQ ; j++){
      pop[i][j] = wk[i][j]  = tau[i][j] = tap[i][j] = 0.0;
      tzp[i][j] = tnp[i][j] = tz0[i][j] = tn0[i][j] = 0.0;
    }
  }
}


/**********************************************************/
/*      Free Allocated T and Spec Array                   */
/**********************************************************/
void preqDeleteAllocated()
{
  for(int j=0 ; j<MAX_CHANNEL ; j++) delete [] wkc[j];
  delete [] wkc;

  for(int j=0 ; j<MAX_PREEQ ; j++){
    delete [] pop[j];
    delete [] wk[j];
    delete [] tzp[j];
    delete [] tnp[j];
    delete [] tz0[j];
    delete [] tn0[j];
    delete [] tau[j];
    delete [] tap[j];
  }
  delete [] pop;
  delete [] wk;
  delete [] tzp;
  delete [] tnp;
  delete [] tz0;
  delete [] tn0;
  delete [] tau;
  delete [] tap;
}


#ifdef FKK_SPIN_DISTRIBUTION
/**********************************************************/
/*      Spin Distribution in Populated Continuum          */
/**********************************************************/
double preqSpinDistribution(double ein, double eout)
{
  const int ne =  6;
  const int nm = 20;
  double emesh[nm];
  static double esig[ne] = {10.0, 15.0, 20.0, 25.0, 30.0, 35.0};

  /*** data for Ti48 */
/*
  static double sigma[ne][nm] ={
     {4.94612, 4.52581, 4.17042, 4.05412, 3.97426, 3.91426, 3.84689, 3.79608, 3.72952, 3.70000,
      3.70000, 3.70000, 3.70000, 3.70000, 3.70000, 3.70000, 3.70000, 3.70000, 3.70000, 3.70000},
     {5.6668 , 4.89715, 4.49396, 4.34464, 4.26521, 4.19036, 4.11129, 4.01214, 3.9301 , 3.81121,
      3.70686, 3.61192, 3.54938, 3.46764, 3.50000, 3.50000, 3.50000, 3.50000, 3.50000, 3.50000},
     {7.61267, 5.39837, 4.69063, 4.41014, 4.24143, 4.07311, 3.84503, 3.60981, 3.40308, 3.40000,
      3.40000, 3.40000, 3.40000, 3.40000, 3.40000, 3.40000, 3.40000, 3.40000, 3.40000, 3.40000},
     {9.54057, 8.20322, 6.97566, 5.54943, 4.78723, 4.46325, 4.27775, 4.10757, 3.87651, 3.59899,
      3.32878, 3.12909, 3.10000, 3.10000, 3.10000, 3.10000, 3.10000, 3.10000, 3.10000, 3.10000},
     {10.6726, 9.44205, 8.78424, 8.27432, 7.24602, 5.93951, 4.99992, 4.55558, 4.33381, 4.15854,
      3.94078, 3.67065, 3.32918, 3.08303, 3.10000, 3.10000, 3.10000, 3.10000, 3.10000, 3.10000},
     {10.3983, 9.78949, 10.0803, 9.7471 , 9.00223, 8.44203, 7.69339, 6.53299, 5.31727, 4.69304,
      4.40937, 4.22281, 4.02479, 3.7522 , 3.40316, 3.1018 , 2.89923, 2.90000, 2.90000, 2.90000}};
*/
  /*** data for Ir191 */
/*
  static double sigma[ne][nm] ={
    {6.29812, 5.87556, 5.60913, 5.43508, 5.30576, 5.23844, 5.19664, 5.16545, 5.10267, 5.70000,
     5.70000, 5.70000, 5.70000, 5.70000, 5.70000, 5.70000, 5.70000, 5.70000, 5.70000, 5.70000},
    {9.61283, 8.6076 , 7.69191, 6.88264, 6.2222 , 5.79778, 5.55587, 5.43764, 5.39826, 5.40286,
     5.42323, 5.44247, 5.40281, 5.32062, 5.30000, 5.30000, 5.30000, 5.30000, 5.30000, 5.30000},
    {12.2362, 10.8465, 8.85537, 7.27167, 6.17139, 5.56054, 5.33419, 5.30138, 5.22723, 5.20000,
     5.20000, 5.20000, 5.20000, 5.20000, 5.20000, 5.20000, 5.20000, 5.20000, 5.20000, 5.20000},
    {16.0706, 14.5435, 12.6773, 11.0204, 9.34873, 7.75666, 6.52749, 5.75037, 5.36394, 5.23983,
     5.16567, 4.92705, 4.90000, 4.90000, 4.90000, 4.90000, 4.90000, 4.90000, 4.90000, 4.90000},
    {19.5191, 17.4217, 16.2695, 15.0362, 13.3566, 11.6285, 10.0707, 8.41097, 7.02787, 6.04529,
     5.47471, 5.22899, 5.12851, 4.94994, 4.90000, 4.90000, 4.90000, 4.90000, 4.90000, 4.90000},
    {21.5658, 20.3883, 19.3842, 18.0372, 16.9426, 15.8209, 14.2163, 12.4219, 10.8042, 9.16168,
     7.61702, 6.4267 , 5.6695 , 5.27737, 5.11846, 4.97934, 4.66004, 4.60000, 4.60000, 4.60000}};
*/
  /*** data for Ir193 */
/*
  static double sigma[ne][nm] ={
    {6.45457, 6.01722, 5.7562 , 5.60856, 5.52915, 5.53507, 5.60666, 5.66302, 5.69336, 5.70000,
     5.70000, 5.70000, 5.70000, 5.70000, 5.70000, 5.70000, 5.70000, 5.70000, 5.70000, 5.70000},
    {9.88524, 8.72456, 7.70421, 6.82612, 6.11792, 5.66874, 5.41087, 5.27068, 5.19763, 5.14898,
     5.09104, 4.99966, 4.84852, 4.72548, 4.70000, 4.70000, 4.70000, 4.70000, 4.70000, 4.70000},
    {12.4718, 11.1775, 9.08253, 7.27877, 6.07349, 5.43502, 5.14683, 4.98372, 4.73213, 4.70000,
     4.70000, 4.70000, 4.70000, 4.70000, 4.70000, 4.70000, 4.70000, 4.70000, 4.70000, 4.70000},
    {16.2325, 14.8403, 12.988 , 11.3411, 9.62599, 7.84253, 6.46241, 5.63577, 5.20739, 4.97739,
     4.72278, 4.37075, 4.30000, 4.30000, 4.30000, 4.30000, 4.30000, 4.30000, 4.30000, 4.30000},
    {19.4832, 17.5301, 16.4231, 15.3259, 13.6751, 11.9577, 10.3981, 8.60205, 7.0178 , 5.94773,
     5.34376, 5.01562, 4.76449, 4.40341, 4.40000, 4.40000, 4.40000, 4.40000, 4.40000, 4.40000},
    {21.0782, 20.1689, 19.3755, 18.0689, 17.1246, 16.0806, 14.5732, 12.7645, 11.1492, 9.44293,
     7.69411, 6.35831, 5.55102, 5.10142, 4.82158, 4.49316, 4.07305, 4.00000, 4.00000, 4.00000}};
*/
  /*** data for U238 */

  static double sigma[ne][nm] ={
    {7.34803, 6.18125, 5.88896, 5.7306 , 5.74426, 5.86739, 6.04142, 6.19847, 6.20506, 6.20000,
     6.20000, 6.20000, 6.20000, 6.20000, 6.20000, 6.20000, 6.20000, 6.20000, 6.20000, 6.20000},
    {10.902 , 9.5573 , 8.51022, 7.39234, 6.60081, 6.08198, 5.74622, 5.5702 , 5.52527, 5.57587,
     5.68646, 5.81751, 5.90119, 5.8215 , 5.80000, 5.80000, 5.80000, 5.80000, 5.80000, 5.80000},
    {13.1021, 12.3896, 11.5142, 10.631 , 9.84031, 9.02356, 8.14444, 7.33029, 6.67214, 6.18081,
     5.8419 , 5.62804, 5.51428, 5.47933, 5.50265, 5.56628, 5.50000, 5.50000, 5.50000, 5.50000},
    {16.4085, 14.5985, 12.8082, 11.033 , 9.56834, 8.15703, 6.81883, 5.92264, 5.51098, 5.4259 ,
     5.47728, 5.39691, 5.40000, 5.40000, 5.40000, 5.40000, 5.40000, 5.40000, 5.40000, 5.40000},
    {19.9834, 17.8214, 16.3039, 14.9384, 13.3437, 11.6879, 10.1224, 8.73502, 7.35124, 6.25145,
     5.63924, 5.40432, 5.38634, 5.39488, 5.40000, 5.40000, 5.40000, 5.40000, 5.40000, 5.40000},
    {22.5612, 20.7902, 19.5887, 18.3384, 16.9733, 15.7226, 14.1647, 12.4242, 10.7686, 9.36118,
     7.97722, 6.70136, 5.85954, 5.46075, 5.34639, 5.32269, 5.14016, 5.20000, 5.20000, 5.20000}};


  int k1 = 0, k2 = 0;
  if(     ein < esig[0   ]) k1 = k2 = 0;
  else if(ein > esig[ne-1]) k1 = k2 = ne-1;
  else{
    for(int i=0 ; i<ne-1 ; i++){
      if( (esig[i] <= ein) && (ein < esig[i+1]) ){
        k1 = i;
        k2 = i+1;
        break;
      }
    }
  }

//if(k1 <= 1){
  if(k1 <= 2){
    for(int j=0 ; j<nm ; j++) emesh[j] = j + 0.5;
  }
  else{
    for(int j=0 ; j<nm ; j++) emesh[j] = 2*j + 1.0;
  }

  int n1 = 0, n2 = 0;
  if(     eout < emesh[0   ]) n1 = n2 = 0;
  else if(eout > emesh[nm-1]) n1 = n2 = nm-1;
  else{
    for(int j=0 ; j<nm-1 ; j++){
      if( (emesh[j] <= eout) && (eout < emesh[j+1]) ){
        n1 = j;
        n2 = j+1;
        break;
      }
    }
  }

  double s1 = 0.0, s2 = 0.0;
  if((n1 == 0) && (n2 == 0)){
    s1 = sigma[k1][0];
    s2 = sigma[k2][0];
  }
  else if( (n1 == nm-1) && (n2 == nm-1) ){
    s1 = sigma[k1][nm-1];
    s2 = sigma[k2][nm-1];
  }
  else{
    s1 = (sigma[k1][n2]-sigma[k1][n1])/(emesh[n2]-emesh[n1])*(eout-emesh[n1])+sigma[k1][n2];
    s2 = (sigma[k2][n2]-sigma[k2][n1])/(emesh[n2]-emesh[n1])*(eout-emesh[n1])+sigma[k2][n2];
  }

  double sigma2 = (k2 == k1) ? s1 : (s2-s1)/(esig[k2]-esig[k1])*(ein-esig[k1])+s1;

  return(sigma2);
}
#endif

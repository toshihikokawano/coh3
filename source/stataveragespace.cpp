/******************************************************************************/
/*  stataveragespace.cpp                                                      */
/*        caclculated D0, and adjust level density to reproduce given D0      */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "nucleus.h"
#include "levden.h"
#include "outext.h"
#include "parameter.h"
#include "terminate.h"

static const double dummyincident = 0.001; // 1 keV added to Sn
static const double d0range       = 0.5;   // a-parameter search range, +/- 50%
static const int    max_iter      = 100;

/**********************************************************/
/*      Calculate D0 from Level Density                   */
/**********************************************************/
double statCalculateD0(const double sp, Nucleus *n0, Nucleus *n1)
{
  double d0 = 0.0;

  double ex = sp + dummyincident;
  double r  = ldLevelDensity(ex,(double)n0->za.getA(),&n0->ldp);
  int    i0 = (int)(2.0*n1->lev[0].spin);

  double fs = parmGetFactor(parmSPIN,n0->za); if(fs < 0.0) fs = 1.0;
  double fd = parmGetFactor(parmPDST,n0->za); if((fd < 0.0) || (fd >2.0)) fd = 1.0;

  double pd = ldParityDistributionEdepend(1,ex,fd);

  /*** when target spin is zero */
  if(i0 == 0){
    double sd0 = ldSpinDistribution(       0.5,n0->ldp.spin_cutoff,n0->ldp.sigma0,fs);
    d0 = 1.0 / (r * sd0 * pd);
  }
  /*** otherwise, sum level densities */
  else{
    double sd1 = ldSpinDistribution((i0+1)*0.5,n0->ldp.spin_cutoff,n0->ldp.sigma0,fs);
    double sd2 = ldSpinDistribution((i0-1)*0.5,n0->ldp.spin_cutoff,n0->ldp.sigma0,fs);
    d0 = 1.0 / (r * (sd1+sd2) * pd);
  }

  return(d0);
}    


/**********************************************************/
/*      Adjust Level Density Parameter for Given D0       */
/**********************************************************/
void statAdjustD0(const double sp, Nucleus *n0, Nucleus *n1)
{
  int parID = adj.getindex(parmD0SR);
  if(parID < 0) return;

  double d0exp = parmGetFactor(parmD0SR);

  /*** save level density parameter for CN */
  double a0 = n0->ldp.a;
  double f0 = 1.0 - d0range;
  double f1 = 1.0 + d0range;
  double f2 = (f0 + f1) / 2.0;

  /*** binary search for the a-parameter */
  int k = 0;
  while(k < max_iter){
    n0->ldp.a = a0 * f2;
    ldTGconnect((double)n0->za.getA(), n0->lev[n0->ndisc-1].energy, n0->ndisc, &n0->ldp);

    double d0  = statCalculateD0(sp,n0,n1);
    if(fabs(d0/d0exp - 1.0) < 1e-5) break;

    if(d0 > d0exp) f0 = f2;
    else           f1 = f2;
    f2 = (f0 + f1) / 2.0;
    k ++;
  }
  if(k == max_iter){
    message << "D0 search faild, parameter out of +/- " << d0range;
    cohTerminateCode("statAdjustD0");
  }

  /*** save factor, change parameter name */  
  adj.parm[parID].type = parmLD;
  adj.parm[parID].za   = n0->za;
  adj.parm[parID].setfactor(f2);

  /*** recalculate level density with adjusted parameter */
  statSetupLevelDensity(n0,&n0->ldp);
}


/**********************************************************/
/*      Compound Decay Width                              */
/**********************************************************/
void    statDecayWidth(System *sys, Pdata *pdt, Transmission **tc, Transmission **td, double **tg, Spectra *spc)
{
  /*** calculate only for neutron below 500 keV */
  if(sys->lab_energy > 0.5) return;
  if(sys->incident.pid != neutron) return;

  const int c0 = 0;
  const int k0 = 0;

  Nucleus *n0 = &ncl[c0];
  double w[MAX_CHANNEL+1];

  /*** target spin and parity */
  int j0  = (int)(2.0*ncl[sys->target_id].lev[sys->target_level].spin);
  int pcn = ncl[sys->target_id].lev[sys->target_level].parity; 

  /*** calculate transmission coefficients */
  statStoreContinuumTransmission(c0,n0->excitation[k0],pdt,tc);
  statStoreContinuumGammaTransmission(k0,tg,n0);
  statStoreDiscreteTransmission(c0,n0->excitation[k0],pdt,td);

  /*** determine level density at neutron separation energy */
  double sn = n0->cdt[sys->inc_id].binding_energy;
  int kx = 0;
  for(int k=0 ; k<n0->ntotal-1 ; k++){
    if((ncl->excitation[k+1] <= sn) && (sn < ncl->excitation[k])){ kx = k; break; }
  }
  double e0 = ncl->excitation[kx+1];
  double e1 = ncl->excitation[kx  ];

  /*** for all compound nucleus spins */
  int jmin = std::abs(j0 - 1);
  int jmax =          j0 + 1;

  for(int jcn=jmin ; jcn<=jmax ; jcn+=2){
    for(int id=0 ; id<MAX_CHANNEL+1 ; id++) w[id] = 0.0;

    w[0] = specTransitionGamma(sumall,k0,pcn,jcn,tg,0.0,n0,spc->cn[0],spc->dp[0]);

    for(int id=1 ; id<MAX_CHANNEL ; id++){
      if(!n0->cdt[id].status) continue;
      Nucleus *n1 = &ncl[n0->cdt[id].next];
      w[id] = specTransitionParticle(sumall,k0,id,pcn,jcn,tc[id],td[id],0.0,n0,n1,spc->cn[id],spc->dp[id]);
    }

    if(n0->fissile) w[MAX_CHANNEL] = specTransitionFission(sumall,k0,pcn,jcn,0.0,n0);

    /*** calculate rho = 1/D0 for each Jcn */
    double r0 = (pcn > 0) ? ncl->density[kx+1][jcn/2].even : ncl->density[kx+1][jcn/2].odd;
    double r1 = (pcn > 0) ? ncl->density[kx  ][jcn/2].even : ncl->density[kx  ][jcn/2].odd;
    double r  = (r1 - r0)/(e1 - e0) * (sn -e0) + r0;

    /*** <Gamma> = D0/2pi int T dE */
    for(int id=0 ; id<MAX_CHANNEL+1 ; id++) w[id] /= r * PI2;

    /*** print out decay width */
    extDecayWidth(pcn*jcn,w);
  }
}



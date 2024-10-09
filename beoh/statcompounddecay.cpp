/******************************************************************************/
/*  statcompounddecay.cpp                                                     */
/*        decay of compound nucleus, sum all transmission coefficients        */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "structur.h"
#include "statmodel.h"
#include "omcalc.h"
#include "optical.h"
#include "nucleus.h"
#include "output.h"
#include "global.h"

/*** skip loop if current population is less than this value */
static const double gPopulationCut = 1.0e-20;

static void specPopulationBypass (const int, double, Nucleus *, double *, double *);


/**********************************************************/
/*      Particle Emission from CNparent to CNdaughter     */
/*      -----                                             */
/*             Usual Hauser-Feshhach calculation          */
/**********************************************************/
void    specCompoundDecay(const int c0, const int k0, Transmission **tc, Transmission **td, double **tg, Spectra *spc)
{
  bool tstat[MAX_CHANNEL];
  int  i0 = (int)(2.0*halfint(ncl[c0].lev[0].spin));

  /*** Check if channel is closed */
  tstat[gammaray] = true;
  for(int id=1 ; id<MAX_CHANNEL ; id++){
    tstat[id] = ncl[c0].cdt[id].status;
    if(tstat[id] && (ncl[c0].excitation[k0] < ncl[c0].cdt[id].binding_energy)){
      tstat[id] = false;
    }
  }

  double sigreac = 0.0;

  /** For photo-induced reaction, spc.dp still contains PE at the top bin, add it to sigreac */
  if((c0 == 0) && (k0 == 0)){
    for(int id=0 ; id<MAX_CHANNEL ; id++){
      if(!tstat[id]) continue;
      Nucleus *n1 = &ncl[ncl[c0].cdt[id].next];
      for(int k=0 ; k<n1->ntotal ; k++) sigreac += spc->dp[id][k];
    }
    sigreac *= ncl[c0].de;
  }

  /*** Loop over CN J and Pariry */
  for(int j0=i0 ; j0<=ncl[c0].jmax*2+i0 ; j0+=2){
    int jdx = (j0-i0)/2;
    for(int p0=-1 ; p0<=1 ; p0+=2){
      double pop = (p0 == 1) ? ncl[c0].pop[k0][jdx].even 
                             : ncl[c0].pop[k0][jdx].odd;
      if( pop <= gPopulationCut ) continue;

      sigreac += pop;
      double tsum = 0.0;

      /*** Sum over gamma and particle decay */
      tsum = specTransmissionSum(sumall,tstat,k0,p0,j0,tg,tc,td,0.0,&ncl[c0],spc);

      /*** If no transition, add population to discrete levels */
      if((tsum == 0.0) && (pop > 0.0)){
        specPopulationBypass(k0,pop,&ncl[c0],spc->cn[gammaray],spc->dp[gammaray]);
        continue;
      }

      /*** Population add to each resudual CN */
      double dp = (tsum == 0.0) ? 0.0 : pop/tsum;

      specTransmissionSum(hauser,tstat,k0,p0,j0,tg,tc,td,dp,&ncl[c0],spc);
    }
  }

  /*** Fission probability at each bin is stored at the top of delta-pop */
  if(ctl.exclusive && ctl.fission){
    specLostPopFraction(c0,k0,sigreac,tstat,spc->dp);

    /*** avoid unrealistic fission population in the no-particle emission case */
    if(k0 == ncl[c0].ncont - 1) ncl[c0].popfis[k0] = 0.0;
  }
}


/**********************************************************/
/*      Sum All Transmission Coefficients from CNparent   */
/**********************************************************/
double  specTransmissionSum(const Statcalcmode mode, bool *tstat, const int k0, const int pcn, const int jcn, double **tg, Transmission **tc, Transmission **td, const double dpop, Nucleus *n0, Spectra *spc)
{
  double tsum = 0.0, tsumf = 0.0, tsumg = 0.0, tsumn = 0.0;

  /*** Sum up gamma-ray transmissions */
  if(tstat[gammaray]){
    tsumg = specTransitionGamma(mode,k0,pcn,jcn,tg,dpop,n0,spc->cn[gammaray],spc->dp[gammaray]);
  }

  /*** Sum up fission transmissions if fissile */
  if(n0->fissile){
    tsumf = specTransitionFission(mode,k0,pcn,jcn,dpop,n0);
  }

  /*** Sum up particle transmissions */
  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!tstat[id]) continue;
    Nucleus *n1 = &ncl[n0->cdt[id].next];
    tsumn += specTransitionParticle(mode,k0,id,pcn,jcn,tc[id],td[id],dpop,n0,n1,spc->cn[id],spc->dp[id]);
  }

  if(mode == sumall){
    /*** Avoid strange fission-only case */
    if((tsumn == 0.0) && (tsumg == 0.0) && (tsumf > 0.0)) tsum = 0.0;
    else tsum = tsumn + tsumg + tsumf;
  }
  
  // Call the analysis routine for potential print out of neutron-gamma competition
//ngCompetitionAnalysis(n0,jcn,pcn,k0,tsumg,tsumn,tsumf,tsum);

  return(tsum);
}


/**********************************************************/
/*      Avoid Population Trap for No Transition Case      */
/*      -----                                             */
/*             In case J in the continuum is so high,     */
/*             there is no transition to the discrete     */
/*             levels due to low multipolarity.           */
/*             The trapped flux is distributed  to        */
/*             the discrete states                        */
/**********************************************************/
void    specPopulationBypass(const int k0, double dp, Nucleus *n, double *spc, double *dlp)
{
  bool addtospectrum = false;
  bool weighteddist  = false;

  double sum = 0.0;
  for(int i=0 ; i<n->ndisc ; i++){
    /*** distribute dPop weighted by 2J+1 */
    if(weighteddist) sum += 2.0*n->lev[i].spin + 1.0;
    /*** distribute dPop evenly */
    else             sum += 1.0;
  }
  dp = dp / sum;

  for(int i=0 ; i<n->ndisc ; i++){
    if(weighteddist) n->lpop[i] += dp * (2.0*n->lev[i].spin + 1.0);
    else             n->lpop[i] += dp;

    /*** Since this gamma-ray is not real,
         we may not add this to the gamma-ray spectrum */
    if(addtospectrum){
      double eg = n->excitation[k0] - n->lev[i].energy;
      int    kg = specFindEnergyBin(eg,n->de);
      if(kg >= 0){
        double y = dp/n->de;
        spc[kg   ] += y;
        dlp[kg+k0] += y;
      }
    }
  }
}


/**********************************************************/
/*      Sum All Population Increments Other Than Fission  */
/**********************************************************/
void    specLostPopFraction(const int c0, const int k0, double sigreac, bool *tstat, double **dp)
{
  const double eps = 1e-8; // round-off error tolerance

  /*** Total decay probabilities at each bin
       sigreac :
          if at the top, total reaction cross section = sgimaR + sum sigmaDI
          else population at k0 in the parent CN
       sigresi : popuation transferred to residual nuclei
       sigcn - sigres = partial fission cross section

       Note that this procedure is not so exact at low energies, 
       since discrete transitions are binned. */

  dp[0][0] = 0.0;
  ncl[c0].popfis[k0] = 0.0;

  if(sigreac <= 0.0) return;

  double sigresi = 0.0;

  for(int id=0 ; id<MAX_CHANNEL ; id++){
    if(!tstat[id]) continue;
    Nucleus *n1 = &ncl[ncl[c0].cdt[id].next];
    for(int k=0 ; k<n1->ntotal ; k++) sigresi += dp[id][k];
  }
  sigresi *= ncl[c0].de;

  double xp = fabs(1.0 - sigresi/sigreac);

  double sigfiss = (xp < eps) ? 0.0 : sigreac - sigresi;
  if(sigfiss < 0.0) sigfiss = 0.0;

  /*** Fission probability at each bin is stored at the top of delta-pop */
  dp[0][0] = sigfiss/sigreac;

  /*** population lost by fission stored in popfis array */
  ncl[c0].popfis[k0] = sigfiss;
/*
  std::cout << std::setw( 5) << ncl[c0].za.getA();
  std::cout << std::setw(12) << ncl[c0].excitation[k0];
  std::cout << std::setw(12) << ncl[c0].popfis[k0];
  std::cout << std::setw(12) << dp[0][0] << std::endl;
*/
}


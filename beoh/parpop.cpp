/******************************************************************************/
/*  parpop.cpp                                                                */
/*        particle population                                                 */
/******************************************************************************/

#include <iostream>

#include "structur.h"
#include "statmodel.h"

static double specTransitionParticleCont (const Statcalcmode, const int, const int, const int, const int, Transmission *, const double, Nucleus *, Nucleus *);
static double specTransitionParticleDisc (const Statcalcmode, const int, const int, const int, const int, Transmission *, const double, Nucleus *, Nucleus *);


/**********************************************************/
/*      Particle Transition                               */
/**********************************************************/
double specTransitionParticle
( const Statcalcmode mode,         // calc. control for sum T/Moldauer/HF
  const int     k0,                // index of parent excited state
  const int     id,                // index of emitting particle
  const int     pcn,               // parent state parity (-1/+1)
  const int     jcn,               // parent spate spin (doubled)
  Transmission  *tc,               // transmission coeff. for continuum
  Transmission  *td,               // transmission coeff. for discrete
  const double  q,                 // pop(parent)/Sigma T, if mode=HF
  Nucleus       *n0,               // nucleus information for the parent
  Nucleus       *n1,               // nucleus information for the daughter
  double        *spc,              // particle emission spectra for this channel
  double        *dlp)              // population increment
{
  double x, y, sum = 0.0;
  int    kp;
  Transmission *tp;

  /*** Continuum to continuum : loop over daughter excitation */
  for(int k1=k0 ; k1<n1->ncont ; k1++){
    kp = k1-k0;
    tp = &tc[kp];
    if(tp->lmax <= 0) continue;
    if(n1->excitation[k1] == 0.0) continue;

    x = specTransitionParticleCont(mode,k1,id,pcn,jcn,tp,q,n0,n1);

    if( (mode == hauser) || (mode == fluctuation) ){
      y = x/n1->de;
      spc[kp] += y;
      dlp[k1] += y;
    }
    sum += x;
  }

  /*** Continuum to levels: Loop over daughter excitation */
  for(int k1=0 ; k1<n1->ndisc ; k1++){

    if(n1->lev[k1].flag == 1) continue;

    tp = &td[k1];
    if(tp->lmax <= 0) continue;

    x = specTransitionParticleDisc(mode,k1,id,pcn,jcn,tp,q,n0,n1);

    if( (mode == hauser) || (mode == fluctuation) ){
      int kp = specFindEnergyBin(tp->ecms,n1->de);
      if(kp >= 0){
        y = x/n1->de;
        spc[kp]    += y;
        dlp[kp+k0] += y;
      }
    }
    sum += x;
  }

  return(sum);
}


/**********************************************************/
/*      Particle Transition : From Level to Levels        */
/**********************************************************/
double specLevelTransitionParticle(const Statcalcmode mode, const int k0, const int id, const int p0, const int j0, Transmission *tc, Transmission *td, const double q, Nucleus *n0, Nucleus *n1, double *spc, double *dlp)
{
  double x, y, sum = 0.0;
  Transmission *tp;

  /*** Level to continuum : loop over daughter excitation */
  for(int k1=0 ; k1<n1->ncont ; k1++){
    tp = &tc[k1];
    if(tp->lmax <= 0) continue;

    x = specTransitionParticleCont(mode,k1,id,p0,j0,tp,q,n0,n1);

    if(mode == hauser){
      int kp = specFindEnergyBin(tp->ecms,n1->de);
      if(kp >=0 ){
        y = x/n1->de;
        spc[kp] += y;
        dlp[k1] += y;
      }
    }
    sum += x;
  }

  /*** Level to level : loop over discrete levels in the daughter */
  for(int k1=0 ; k1<n1->ndisc ; k1++){

    if(n1->lev[k1].flag == 1) continue;

    tp = &td[k1];
    if(tp->lmax <= 0) continue;

    x = specTransitionParticleDisc(mode,k1,id,p0,j0,tp,q,n0,n1);

    if(mode == hauser){
      int kp = specFindEnergyBin(tp->ecms,n1->de);
      if(kp >=0 ){
        y = x/n1->de;
        spc[kp] += y;
        dlp[kp+k0] += y;
      }
    }
    sum += x;
  }

   return(sum);
}


/**********************************************************/
/*      Particle Transition : Continuum to Continuum      */
/**********************************************************/
double specTransitionParticleCont(const Statcalcmode mode, const int k1, const int id, const int pcn, const int jcn, Transmission *tc, const double q, Nucleus *n0, Nucleus *n1)
{
  int    jtop = 0;
  double x0 = 0.0, x1 = 0.0, dx = 0.0, sum = 0.0;

  /*** Loop over particle angular momentum */
  for(int lp=0 ; lp<=2*tc->lmax ; lp+=2){
    int pl = ((lp/2)%2 == 1) ? -1 : 1;
    for(int sp=n0->cdt[id].spin2 ; sp>=-n0->cdt[id].spin2 ; sp-=2){
      int jp = lp + sp;  if(jp < 0) continue;

      double tj  = tc->tran[tj_index(lp,sp,n0->cdt[id].spin2)];
      double wfc = (mode == fluctuation) ? statMoldauer(-1,0,0,0,tj,0.0) : 1.0;

      int j1min = std::abs(jcn-jp);
      int j1max =          jcn+jp ; if(j1max > 2*(MAX_J-1)) j1max=2*(MAX_J-1);

      if(j1max > jtop) jtop = j1max;

      for(int j1=j1min ; j1<=j1max ; j1+=2){

        int idx = (j1-(int)(2.0*halfint(n1->lev[0].spin)))/2;

        if(pl*pcn == 1){
          dx = n1->density[k1][idx].even*n1->de;
          x0 = tj*dx;
          if(mode == wfactor){
            statWidthFluctuation(tj,dx,0.0);
          }else if( (mode == hauser) || (mode == fluctuation) ){
            n1->pop[k1][idx].even += (x1=q*x0*wfc);
          }
        }else{
          dx = n1->density[k1][idx].odd*n1->de;
          x0 = tj*dx;
          if(mode == wfactor){
            statWidthFluctuation(tj,dx,0.0);
          }else if( (mode == hauser) || (mode == fluctuation) ){
            n1->pop[k1][idx].odd  += (x1=q*x0*wfc);
          }
        }
        if(mode == sumall) sum += x0;
        else if( (mode == hauser) || (mode == fluctuation) ) sum += x1;
      }
    }
  }
  if(n1->jmax < jtop/2) n1->jmax = jtop/2;
  return(sum);
}


/**********************************************************/
/*      Particle Transition : Continuum to Levels         */
/**********************************************************/
double specTransitionParticleDisc(const Statcalcmode mode, const int k1, const int id, const int pcn, const int jcn, Transmission *td, const double q, Nucleus *n0, Nucleus *n1)
{
  double x1 = 0.0, sum = 0.0;

  int j1 = (int)(2.0*n1->lev[k1].spin);

  /*** Loop over particle angular momentum */
  for(int lp=0 ; lp<=2*td->lmax ; lp+=2){
    int pl = ((lp/2)%2 == 1) ? -1 : 1;
    if(pl != pcn*n1->lev[k1].parity) continue;

    for(int sp=n0->cdt[id].spin2 ; sp>=-n0->cdt[id].spin2 ; sp-=2){
      int jp = lp + sp;  if(jp < 0) continue;

      int j1min = std::abs(jcn-jp);
      int j1max =          jcn+jp ;

      if(j1 >= j1min && j1 <= j1max){

        double tj  = td->tran[tj_index(lp,sp,n0->cdt[id].spin2)];
        double wfc = (mode == fluctuation) ? statMoldauer(lp,jp,id,k1,tj,0.0) : 1.0;

        if(mode == wfactor){
          statWidthFluctuation(tj,1.0,0.0);
        }else if( (mode == hauser) || (mode == fluctuation) ){
          n1->lpop[k1] += (x1=tj*q*wfc);
        }
        if(mode == sumall) sum += tj;
        else if( (mode == hauser) || (mode == fluctuation) ) sum += x1;
      }
    }
  }

  return(sum);
}



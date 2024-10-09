/******************************************************************************/
/*  gampop.cpp                                                                */
/*        gamma ray population                                                */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "statmodel.h"
#include "levden.h"

static double specTransitionGammaCont (const Statcalcmode, const int, const int, const int, const int, const int, const double, const double, Nucleus *);
static double specTransitionGammaDisc (const Statcalcmode, const int, const int, const int, const int, const int, const double, const double, Nucleus *);

/*  Parity and L of Gamma-ray E1(-), M1(+), E2(+), M2(-)... */
inline int gammaL(const int m){ return (m/2 + 1); }
inline int gammaParity(const int m){
  int pm = (m%2 == 0) ? 1 : -1;  // E:even, M:odd parity
  int l  = gammaL(m);            // L
  int pl = (l%2 == 0) ? 1 : -1;  // parity of L
  return(pm * pl);
}


/**********************************************************/
/*      Gamma-ray Transition                              */
/**********************************************************/
double specTransitionGamma
( const Statcalcmode mode,         // calc. control for sum T/Moldauer/HF
  const int     k0,                // index of parent excited state
  const int     p0,                // parent state parity (-1/+1)
  const int     j0,                // parent spate spin (doubled)
  double        **tg,              // transmission coeff. for gamma-ray
  const double  q,                 // pop(parent)/Sigma T, if mode=HF
  Nucleus       *n,                // nucleus information
  double        *spc,              // gamma-ray spectrum
  double        *dlp)              // poplation increment
{
  /*** continuum to continuum */
  double x,y,sum = 0.0,td[MAX_MULTIPOL];
  int    kg;

  for(int k1=k0+1 ; k1<n->ncont ; k1++){

    x = 0.0;
    for(int m=0 ; m<MAX_MULTIPOL ; m++){
      int pg = gammaParity(m);
      int l  = gammaL(m);
      x += specTransitionGammaCont(mode,pg,k1,p0,j0,l,tg[k1][m],q,n);
    }

    if( (mode == hauser) || (mode == fluctuation) ){
      kg = k1-k0;
      y  = x/n->de;
      spc[kg] += y;
      dlp[k1] += y;
    }
    sum += x;
  }

  /*** continuum to level */
  for(int i=0 ; i<n->ndisc ; i++){

    if(n->lev[i].flag == 1) continue;

    double eg = n->excitation[k0] - n->lev[i].energy;
    if(eg < 0.0) continue;

    statStoreDiscreteGammaTransmission(eg,td,n);

    x = 0.0;
    for(int m=0 ; m<MAX_MULTIPOL ; m++){
      int pg = gammaParity(m);
      int l  = gammaL(m);
      x += specTransitionGammaDisc(mode,pg,i,p0,j0,l,td[m],q,n);
    }

    if( (mode == hauser) || (mode == fluctuation) ){
      kg = specFindEnergyBin(eg,n->de);
      if(kg >= 0){
        y  = x/n->de;
        spc[kg   ] += y;
        dlp[kg+k0] += y;
      }
    }
    sum += x;
  }

  return(sum);
}


/**********************************************************/
/*      Gamma-ray Transition : From Level to Levels       */
/**********************************************************/
double  specLevelTransitionGamma(const Statcalcmode mode, const int k0, const int p0, const int j0, const double q, Nucleus *n)
{
  double x = 0.0,td[MAX_MULTIPOL];

  for(int k1=0 ; k1<k0-1 ; k1++){

    double eg = n->lev[k0].energy - n->lev[k1].energy;
    if(eg < 0.0) continue;

    statStoreDiscreteGammaTransmission(eg,td,n);

    x = 0.0;
    for(int m=0 ; m<MAX_MULTIPOL ; m++){
      int pg = gammaParity(m);
      int l  = gammaL(m);
      x += specTransitionGammaDisc(mode,pg,k1,p0,j0,l,td[m],q,n);
    }
  }

  return(x);
}

 
/**********************************************************/
/*      Gamma-ray Transition : Continuum to Continuum     */
/**********************************************************/
double specTransitionGammaCont(Statcalcmode mode, int pg, int k, int p0, int j0, int l,
                double tg, double q, Nucleus *n)
{
  double  x0=0.0,x1=0.0,dx=0.0,sum=0.0,wfc=1.0;
  int j1min = std::abs(j0-2*l);
  int j1max =          j0+2*l ;

  if(tg == 0.0) return(sum);

  if(j1max > 2*MAX_J-1) j1max=2*MAX_J-1;
  if(n->jmax < j1max/2) n->jmax=j1max/2;

  wfc = (mode == fluctuation) ? statMoldauer(-1,0,0,0,tg,0.0) : 1.0;

  for(int j1=j1min ; j1<=j1max ; j1+=2){
    if( (j0 == 0) && (j1 == 0) ) continue; // gamma decay prohibited for 0->0
    int idx = (j1-(int)(2.0*halfint(n->lev[0].spin)))/2;

    if(pg*p0 == 1){
      dx = n->density[k][idx].even*n->de;
      x0 = tg*dx;
      if(mode == wfactor){
        statWidthFluctuation(tg,dx,0.0);
      }else if( (mode == hauser) || (mode == fluctuation) ){
        n->pop[k][idx].even += (x1=x0*q*wfc);
      }
    }else{
      dx = n->density[k][idx].odd *n->de;
      x0 = tg*dx;
      if(mode == wfactor){
        statWidthFluctuation(tg,dx,0.0);
      }else if( (mode == hauser) || (mode == fluctuation) ){
        n->pop[k][idx].odd  += (x1=x0*q*wfc);
      }
    }
    if(mode == sumall) sum += x0;
    else if( (mode == hauser) || (mode == fluctuation) ) sum += x1;
  }

  return(sum);
}


/**********************************************************/
/*      Gamma-ray Transition : Continuum to Levels        */
/**********************************************************/
double specTransitionGammaDisc(const Statcalcmode mode, const int pg, const int i, const int p0, const int j0, const int l, const double tg, const double q, Nucleus *n)
{
  double x1=0.0,sum=0.0,wfc=1.0;
  double s1min = fabs(j0/2.0-l);
  double s1max =      j0/2.0+l ;

  if(tg == 0.0) return(sum);
  if( (j0 == 0) && ((int)n->lev[i].spin == 0) ) return (sum);

  wfc = (mode == fluctuation) ? statMoldauer(-1,0,0,0,tg,0.0) : 1.0;

  if(n->lev[i].spin >=s1min && n->lev[i].spin<=s1max){
    if(pg == p0*n->lev[i].parity){
      if(mode == wfactor){
        statWidthFluctuation(tg,1.0,0.0);
      }else if( (mode == hauser) || (mode == fluctuation) ){
        n->lpop[i] += (x1=tg*q*wfc);
      }
      if(mode == sumall) sum += tg;
      else if( (mode == hauser) || (mode == fluctuation) ) sum += x1;
    }
  }

  return(sum);
}


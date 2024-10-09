/******************************************************************************/
/*  beohprofile.cpp                                                           */
/*        calculate strength distribution after beta-decay                    */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "beoh.h"
#include "statmodel.h"
#include "beohstat.h"
#include "terminate.h"

static bool   beohStrengthENSDF (Beta *);
static void   beohGaussianBroadening (Beta *, BetaProfile *);
static double beohDiscreteStrength (Beta *, BetaProfile *);
static double beohDiscreteRemap (BetaProfile *);
static double beohDiscreteScale (const double, BetaProfile *);
static double beohRenormalizeContinuum (double, BetaProfile *);
static double beohRenormalizeDiscrete (double, BetaProfile *);


/*** if ENSDF sum is more than this number, we use ENSDF only */
static const double THRESHOLD_COMPLETENESS = 0.90 ;

/*** continuum beta strength zero if below this value */
static const double THRESHOLD_STRENGTH = 1.0e-10;

/*** default beta strength function Gaussian broadening resolution */
static const double PROFILE_BROADENING_WIDTH = 0.1;

/*** highest discrete level energy */
static double elevhigh = 0.0;

/*** scaling factor for low lying GT discrete transitions */
static double pfactor = 1.0;



/***********************************************************/
/*      Read in GT-Strength and ENSDF, and Merge Them      */
/***********************************************************/
void beohStrengthProfile(Beta *gts, Beta *ens, BetaProfile *bpf)
{
  elevhigh = bpf->lev[bpf->ndisc-1].energy;
  double td = 0.0, tc = 0.0;

  /*** if ENSDF is thought to be complete, use all the data */
  if(beohStrengthENSDF(ens)){
    beohGaussianBroadening(ens,bpf);
    td = beohDiscreteStrength(ens,bpf);
    tc = beohRenormalizeContinuum(td,bpf);
    bpf->source = "   ENSDF   ";
  }
  else{
    beohGaussianBroadening(gts,bpf);

    /*** if ENSDF not found, use Gamow-Teller only */
    if(ens->nstate == 0){
      td = beohDiscreteRemap(bpf);
      bpf->source = "GamowTeller";
    }
    /*** mix ENSDF and Gamow-Teller */
    else{
      td = beohDiscreteStrength(ens,bpf);
      bpf->source = "   Mixed   ";
    }
    if(pfactor > 1.0) td = beohDiscreteScale(pfactor,bpf);
    tc = beohRenormalizeContinuum(td,bpf);
  }

  /*** renormalize again, if no continuum */
  if( (tc == 0.0) && ((tc+td)<1.0) ){
    td = beohRenormalizeDiscrete(td,bpf);
  }

  bpf->tcont = tc;
  bpf->tdisc = td;
}


/***********************************************************/
/*      Check If ENSDF Is Complete                         */
/***********************************************************/
bool beohStrengthENSDF(Beta *ensdf)
{
  bool c = false;
  /*** get a fraction of decay to discrete levels from ENSDF */
  double x = 0.0;
  for(int i=0 ; i<ensdf->nstate ; i++) x += ensdf->br[i].getR();

  if(x > THRESHOLD_COMPLETENESS){
    x = 1.0/x;
    for(int i=0 ; i<ensdf->nstate ; i++) ensdf->br[i].scaleR(x);
    c = true;
  }

  return(c);
}


/***********************************************************/
/*      Gaussian Broadening Beta Strength Data             */
/***********************************************************/
void beohGaussianBroadening(Beta *b, BetaProfile *p)
{
  /*** Gaussian broadening for all beta transitions,
       generate normalized continuum final distribution */

  double g = PROFILE_BROADENING_WIDTH;
  double c = 1.0 / (sqrt(PI2)*g);

  if(b->width > 0.0) g = b->width;

  double s = 0.0;
  /*** start k=1 to skip zero energy transition */
  for(int k=1 ; k<p->ntotal ; k++){
    p->rcont[k] = 0.0;
    for(int j=0 ; j<b->nstate ; j++){
      double d = b->br[j].getE() - p->excitation[k];
      p->rcont[k] += c * b->br[j].getR() * exp( -d*d/(g*g) );
    }

    if(p->rcont[k] < THRESHOLD_STRENGTH) p->rcont[k] = 0.0;

    s += p->rcont[k];
  }
  if(s == 0.0){
    message << "beta-decay strength zero";
    cohTerminateCode("beohGaussianBroadening");
  }

  s = 1.0/s;
  for(int k=0 ; k<p->ntotal ; k++) p->rcont[k] *= s;
}


/***********************************************************/
/*      Discrete Beta Strength for Same Level Scheme       */
/***********************************************************/
double beohDiscreteStrength(Beta *b, BetaProfile *p)
{
  /*** first, calculate sum of discrete transition strengths
       that are above Elevhigh (cont/disc boundary) */
  for(int k=0 ; k<b->nstate ;  k++){
    if(b->br[k].getE() > elevhigh) break;
  }

  /*** find the same energy level */
  for(int k=0 ; k<b->nstate ;  k++){
    if(b->br[k].getE() > elevhigh) break;

    /*** RIPL and ENSDF level energies should be consistent */
    bool found = false;
    for(int i=0 ; i<p->ndisc ; i++){
      if( fabs(b->br[k].getE() - p->lev[i].energy) < 0.001 ){
        p->rdisc[i] = b->br[k].getR();
        found = true;
        break;
      }
    }
    /*** if not, force feed to the nearest level */
    if(!found){
      int i0 = 0;
      double d0 = fabs(b->br[k].getE() - p->lev[i0].energy);
      for(int i1=1 ; i1<p->ndisc ; i1++){
        double d1 = fabs(b->br[k].getE() - p->lev[i1].energy);
        if( d1<d0 ){
          d0 = d1;
          i0 = i1;
        }
      }
      p->rdisc[i0] = b->br[k].getR();
    }
  }

  /*** sum strength in the discrete levels */
  double t = 0.0;
  for(int i=0 ; i<p->ndisc ; i++) t += p->rdisc[i];

  return(t);
}


/***********************************************************/
/*      Renormalize Discrete Strength                      */
/***********************************************************/
double beohRenormalizeDiscrete(double t, BetaProfile *p)
{
  double s = 0.0;
  if(t > 0.0){
    t = 1.0/t;
    for(int i=0 ; i<p->ndisc ; i++){
      p->rdisc[i] *= t;
      s += p->rdisc[i];
    }
  }

  return(s);
}


/***********************************************************/
/*      Renormalize Continuum Strength                     */
/***********************************************************/
double beohRenormalizeContinuum(double t, BetaProfile *p)
{
  /*** renormalize continuum part to be 1 - sum(discrete) */
  double s = 0.0;
  for(int k=0 ; k<p->ncont ; k++) s += p->rcont[k];

  if(s > 0.0){
    s = (1.0-t)/s;
    for(int k=0 ; k<p->ncont ; k++) p->rcont[k] *= s;
  }

  /*** clear overlapped region */
  for(int k=p->ncont ; k<p->ntotal ; k++) p->rcont[k] = 0.0;

  s = 0.0;
  for(int k=0 ; k<p->ntotal ; k++) s += p->rcont[k];

  return(s);
}


/***********************************************************/
/*      Remap Beta Strength Data onto Discrete Levels      */
/***********************************************************/
double beohDiscreteRemap(BetaProfile *p)
{
  /*** distribute the strength in the discrete region
       onto the nearest discerete levels,
       regardless the spin and parity.
       when no final state exists in the energy bin,
       this flux is distributed to all other level population. */

  double nolevpop = 0.0, levpop = 0.0;
  /*** scan discrete region */
  for(int k=p->ncont ; k<p->ntotal ;  k++){
    if(p->rcont[k]==0.0) continue;

    double e1 = p->excitation[k+1];
    double e2 = p->excitation[k  ];
    int    m  = 0; // number of levels in [E1,E2]
    for(int i=1 ; i<p->ndisc ; i++){
      if( (e1 <= p->lev[i].energy) && (p->lev[i].energy < e2) ){
        m++;
        }
    }

    /*** when levels exist between E1 and E2 */
    if(m>0){
      for(int i=1 ; i<p->ndisc ; i++){
        if( (e1 <= p->lev[i].energy) && (p->lev[i].energy < e2) ){
          p->rdisc[i] += p->rcont[k] / (double)m;
        }
      }
    }
    /*** no level case */
    else if(m==0) nolevpop += p->rcont[k];

    levpop += p->rcont[k];
  }

  /*** when no-level bins are found, scale-up others */
  if(levpop>nolevpop){
    double scaling = levpop/(levpop-nolevpop);
    for(int i=0 ; i<p->ndisc ; i++) p->rdisc[i] *= scaling;
  }

  double t = 0.0;
  for(int i=0 ; i<p->ndisc ; i++) t += p->rdisc[i];
  
  return(t);
}


/***********************************************************/
/*      Renormalize Beta Strength in Discrete Levels       */
/***********************************************************/
double beohDiscreteScale(const double f, BetaProfile *p)
{
  double s = 0.0;

  for(int i=0 ; i<p->ndisc ; i++) s += f*p->rdisc[i];

  /*** if scaled p->rcont[] exceed 1.0, set zero in continuum,
       and remap the 1.0 strenght on the discrete states */
  if(s>1.0){
    for(int k=0 ; k<p->ncont ; k++) p->rcont[k] = 0.0;
    for(int i=0 ; i<p->ndisc ; i++) p->rdisc[i] = f * p->rcont[i] / s;
  }else{
    for(int k=0 ; k<p->ncont ; k++) p->rcont[k] *= 1-s;
    for(int i=0 ; i<p->ndisc ; i++) p->rdisc[i] *= f;
  }

  s = 0.0;
  for(int i=0 ; i<p->ndisc ; i++) s += p->rdisc[i];

  return(s);
}


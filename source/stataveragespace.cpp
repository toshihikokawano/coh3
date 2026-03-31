/******************************************************************************/
/*  stataveragespace.cpp                                                      */
/*        caclculated D0, and adjust level density to reproduce given D0      */
/******************************************************************************/

#include <iostream>
#include <cmath>
#include <complex>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "nucleus.h"
#include "levden.h"
#include "outext.h"
#include "global.h"
#include "parameter.h"
#include "terminate.h"

static double statPenetrationFactor(const int, const double);

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
/*      Compound Decay Width in URR format                */
/**********************************************************/
void    statDecayWidth(System *sys, Transmission *tin, Transmission **tc, Transmission **td, double **tg, Spectra *spc)
{
  if(sys->incident.pid != neutron){
    message << "decay width calculation only for neutron incident reactions";
    cohTerminateCode("statDecayWidth");
    return;
  }

  if(ncl[sys->target_id].ncont > 0){
    message << "no decay width calculation when continuum inelastic exists";
    cohTerminateCode("statDecayWidth");
    return;
  }

  const int c0 = 0;
  const int k0 = 0;
  const double eps = 1e-6; // cut-off for the inelastic scattering width

  Nucleus *n0 = &ncl[c0];

  int lmax = td[neutron][sys->target_level].lmax;
  int jmax = ncl[c0].jmax;

  /*** memory allocation */
  double ***gamin = new double ** [jmax+1]; // Gamma for incident channel  [J][Pi][L]
  double ***gamma = new double ** [jmax+1]; // Gamma for decay channels    [J][Pi][channel]
  double  **rho   = new double  * [jmax+1]; // level density

  for(int j=0 ; j<=jmax ; j++){
    gamma[j] = new double * [2];
    gamin[j] = new double * [2];
    rho[j]   = new double [2];
    for(int p=0 ; p<2 ; p++){
      gamma[j][p] = new double [MAX_CHANNEL+2];  // [MC] = fission, [MC+1] = inel
      gamin[j][p] = new double [lmax+1];
      for(int c=0 ; c<MAX_CHANNEL+2 ; c++) gamma[j][p][c] = 0.0;
      for(int l=0 ; l<=lmax ; l++) gamin[j][p][l] = 0.0;
    }
  }

  /*** copy channel status to tstat array */
  bool tstat[MAX_CHANNEL+2];
  for(int id=0 ; id<MAX_CHANNEL ; id++) tstat[id] = n0->cdt[id].status;
  tstat[gammaray] = true;
  tstat[MAX_CHANNEL] = n0->fissile;
  tstat[MAX_CHANNEL+1] = false;

  /*** calculate transmission coefficients */
  statStoreContinuumGammaTransmission(k0,tg,n0);

  /*** projectile, target GS, and target spin and parity */
  int spinproj = ncl[c0].cdt[sys->incident.pid].spin2;                       // projectile spin
  int spintarg = (int)(2.0*ncl[sys->target_id].lev[sys->target_level].spin); // target spin
  int paritarg = ncl[sys->target_id].lev[sys->target_level].parity;          // parget parity
  int spings   = (int)(2.0*halfint(ncl[c0].lev[0].spin));                    // CN ground state spin


  /*** for entrance channel partial waves */
  for(int j=0 ; j<=jmax ; j++){
    int j0 = 2*j + spings;
    for(int p=0 ; p<2 ; p++){
      int p0 = (p == 0) ? 1 : -1; // first even then odd

      for(int l=0 ; l<=lmax ; l++){
        int lp0 = 2*l;
        if(parity(lp0) != p0*paritarg) continue;

        for(int sp0=spinproj ; sp0>=-spinproj ; sp0-=2){
          int jp0 = lp0 + sp0;
          if(jp0 < 0) continue;

          int j0min = std::abs(jp0 - spintarg);
          int j0max =          jp0 + spintarg;

          /*** sum incident channel transmission to make it l-dependent only */
          if( (j0min <= j0) && (j0 <= j0max) ){
            gamin[j][p][l] += tin->tran[tj_index(lp0,sp0,spinproj)];
          }
        }
      }
    }
  }


  bool firstcall = true; // to print header

  /*** for exit channels, calculate decay widths from CN, which depend on JPi only */
  for(int j=0 ; j<=jmax ; j++){
    int j0 = 2*j + spings;
    for(int p=0 ; p<2 ; p++){
      int p0 = (p == 0) ? 1 : -1; // even then odd

      /*** sum gamma-ray transmission coefficient */
      gamma[j][p][gammaray] = specTransitionGamma(sumall,k0,p0,j0,tg,0.0,n0,spc->cn[gammaray],spc->dp[gammaray]);

      /*** sum fission transmissions if fissile */
      if(n0->fissile) gamma[j][p][MAX_CHANNEL] = specTransitionFission(sumall,k0,p0,j0,0.0,n0);

      /*** particle widths */
      for(int id=1 ; id<MAX_CHANNEL ; id++){
        if(!tstat[id]) continue;
        Nucleus *n1 = &ncl[n0->cdt[id].next];
        gamma[j][p][id] = specTransitionParticle(sumall,k0,id,p0,j0,tc[id],td[id],0.0,n0,n1,spc->cn[id],spc->dp[id]);
      }

      /*** if inelastic channels open, subtract elastic scattering */
      if(td[neutron][1].lmax > 0){
        double g0 = 0.0;
        for(int l=0 ; l<=lmax ; l++) g0 += gamin[j][p][l];
        gamma[j][p][MAX_CHANNEL+1] = gamma[j][p][neutron] - g0;
        if(std::abs(gamma[j][p][MAX_CHANNEL+1]) < eps) gamma[j][p][MAX_CHANNEL+1] = 0.0;
      }

      /*** rho = 1/D for each Jcn */
      rho[j][p] = (p0 > 0) ? n0->density[k0][j].even : n0->density[k0][j].odd;

      /*** <Gamma> = D0/2pi int T dE */
      for(int id=0 ; id<MAX_CHANNEL+2 ; id++){
        gamma[j][p][id] = (rho[j][p] > 0.0) ? gamma[j][p][id] / rho[j][p] / PI2 : 0.0;
      }
    }
  }

  /*** for neutron, include penetration factor */
  double alpha = sys->wave_number * 1.2 * pow((double)sys->target.getA(),1/3.0);

  /*** for each of incident channel L */
  for(int l=0 ; l<=lmax ; l++){
    int lp0 = 2*l;

    double pfact = alpha / statPenetrationFactor(l,alpha) / sqrt(sys->lab_energy * 1e+6);

    /*** determine J range when intrinsic spin sum is not taken */
    int j0min = 2*jmax;
    int j0max = 0;

    for(int sp0=spinproj ; sp0>=-spinproj ; sp0-=2){
      int jp0 = lp0 + sp0;
      if(jp0 < 0) continue;

      int j1min = std::abs(jp0 - spintarg);
      int j1max =          jp0 + spintarg;

      for(int j=j1min ; j<=j1max ; j+=2){
        for(int p0=-1 ; p0<=1 ; p0+=2){
          if(parity(lp0) != p0*paritarg) continue;
          if(j > j0max) j0max = j;
          if(j < j0min) j0min = j;
        }
      }
    }

    for(int j2=j0min ; j2<=j0max ; j2+=2){
      int j  = j2/2;
      for(int p=0 ; p<2 ; p++){
        int p0 = (p == 0) ? 1 : -1; // even then odd
        if(parity(lp0) != p0*paritarg) continue;

        /*** replace neutron decay channel by incident channel */
        gamma[j][p][neutron] = (rho[j][p] > 0.0) ? gamin[j][p][l] * pfact / rho[j][p] / PI2 : 0.0;

        /*** print out decay width */
        extDecayWidth(firstcall,lp0,j2,p0,rho[j][p],gamma[j][p]);
        firstcall = false;
      }
    }
  }

  /*** free allocated memory */
  for(int j=0 ; j<=jmax ; j++){
    for(int p=0 ; p<2 ; p++){
      delete [] gamma[j][p];
      delete [] gamin[j][p];
    }
    delete [] gamma[j];
    delete [] gamin[j];
    delete [] rho[j];
  }
  delete [] gamma;
  delete [] gamin;
  delete [] rho;
}


/**********************************************************/
/*      Penetrability at neutron-incident energy          */
/**********************************************************/
static double statPenetrationFactor(const int l, const double a)
{
  std::complex<double> h,d;

  double g0 = cos(a);
  double f0 = sin(a);
  double g1 =  f0 + g0/a;
  double f1 = -g0 + f0/a;

  if(l == 0){
    h = std::complex<double>( g0,f0);
    d = std::complex<double>(-f0,g0);
  }
  else if(l == 1){
    h = std::complex<double>(g1,f1);
    d = std::complex<double>(-f1-g0/(a*a),g1-f0/(a*a));
  }
  else{
    double g2 = 0.0, f2 = 0.0;
    for(int k=2 ; k<=l ; k++){
      g2 = (2*k-1)/a * g1 - g0;  g0 = g1;  g1 = g2;
      f2 = (2*k-1)/a * f1 - f0;  f0 = f1;  f1 = f2;
    }
    h = std::complex<double>(g2,f2);
    d = std::complex<double>(g0-l/a*g1, f0-l/a*f1);
  }
 
  double pen = a * (d / h).imag();

  if(pen < 0.0){
    double p = 0.0;
    switch(l){
    case 0: p = a; break;
    case 1: p = pow(a,3) / (  1.0 +      a*a); break;
    case 2: p = pow(a,5) / (  9.0 +  3.0*a*a +     pow(a,4)); break;
    case 3: p = pow(a,7) / (225.0 + 45.0*a*a + 6.0*pow(a,4) + pow(a,6)); break;
    default: break;
    }                             
    pen = p;
  }

  return(pen);
}

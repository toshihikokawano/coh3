/******************************************************************************/
/*  statdecay.cpp                                                             */
/*        decay of initial compound nucleus, including width fluctuation      */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "structur.h"
#include "statmodel.h"
#include "omcalc.h"
#include "optical.h"
#include "nucleus.h"
#include "output.h"
#include "outformat.h"
#include "global.h"
#include "ewtrans.h"
#include "etc.h"


static int  specCNDecaySpherical (const int, const int, const int, const int, const double, Transmission *, Transmission **, Transmission **,  double **, Spectra *);
static int  specCNDecayDeformed (const int, const int, const int, const int, const double, Transmission **, Transmission **, double **, Spectra *, int *);
static void specCNDecayPreProcess (System *, const int);
static void specCNDecayPostProcess (Transmission **, Spectra *);
static void specContinuumAngularDistribution (const int, const int, Nucleus *);
static void specContinuumLegendreAllocate (void);
static void specContinuumLegendreFree (void);

static bool tstat[MAX_CHANNEL];
static int  targid = 0, incid = 0, targlev = 0;
static int  spinproj = 0, spings = 0, spintarg = 0, paritarg = 0;

static double ***cleg; // Legendre coefficients in continuum, optional calc

#undef SelectJPi
//#define SelectJPi 0


/**********************************************************/
/*      Compound Elastic Scattering, Witdth Fluctuation   */
/*      -----                                             */
/*             For decay of the first compound            */
/**********************************************************/
void specCompoundElastic
( System       *sys,               // system parameters
  double       renorm,             // re-normalizatoin of total CN cross section
  Transmission *tin,               // transmission coefficient for entrance channel
  Transmission **tc,               // transmission in the continuum region
  Transmission **td,               // transmission for the discrete levels
  double       **tg,               // gamma-ray transmission
  Direct       *dir,               // discrete levels for which direct c.s. given
  Spectra      *spc)               // particle emission spectra
{
  const int c0 = 0;
  const int k0 = 0;

  /*** Store some global variables */
  specCNDecayPreProcess(sys,c0);

  double c1 = NORM_FACT*PI/(sys->wave_number*sys->wave_number) * renorm;
  double c2 = c1 /((spintarg+1.0)*(spinproj+1.0));

  /*** Copy channel status to tstat array */
  for(int id=0 ; id<MAX_CHANNEL ; id++) tstat[id] = ncl[c0].cdt[id].status;
  tstat[gammaray] = true;

  /*** Engelbrecht-Weidenmueller transformation */
  if(ctl.ewtransform){
    /*** Preset coupled-channels parameters */
    ccPreset(sys->cms_energy,sys->excitation,&sys->incident,&sys->target,sys->reduced_mass,dir);
    
    /*** skip coupled-levels in the summation by setting lmax = 0 */
    for(int k=0 ; k<ncl[ncl[c0].cdt[incid].next].ndisc ; k++)
      if(sys->bandidx[k] >= 0) td[incid][k].lmax = 0;

    /*** Loop over CN J and Parity*/
    int jmax = ncl[c0].jmax*2 + spings;
    for(int j0=spings ; j0<=jmax ; j0+=2){
      double c3 = c2 * (j0+1.0);
      int f = 0;
      for(int p0=-1 ; p0<=1 ; p0+=2){
        f = specCNDecayDeformed(c0,k0,j0,p0,c3,tc,td,tg,spc,sys->bandidx);
      }
      if(f < 0) break;
    }
    ccCleanUp();
  }

  /*** Hauser-Feshbach-Moldauer calculation */
  else{

    /*** Loop over CN J and Parity*/
    int jmax = ncl[c0].jmax*2 + spings;
    for(int j0=spings ; j0<=jmax ; j0+=2){
      double c3 = c2 * (j0+1.0);
      for(int p0=-1 ; p0<=1 ; p0+=2){
        specCNDecaySpherical(c0,k0,j0,p0,c3,tin,tc,td,tg,spc);
      }
    }
  }

  if(ctl.fluctuation) statWidthFluctuationReset(0.0);

  /*** Post processing */
  specCNDecayPostProcess(td,spc);
}


/**********************************************************/
/*      CN Decay for Given JPi, Spherical Case            */
/**********************************************************/
int specCNDecaySpherical
(const int c0,      // CN index
 const int k0,      // CN initial excitation energy index
 const int j0,      // CN spin
 const int p0,      // CN parity
 const double g,    // spin statistical factor
 Transmission *tin, Transmission **tc, Transmission **td, double **tg, Spectra *spc)
{
  /*** denominator summation over gamma and particle decay */
  double tsum = specTransmissionSum(sumall,tstat,k0,p0,j0,tg,tc,td,0.0,&ncl[c0],spc);
  if(tsum == 0.0) return -1;

  /*** Width fluctuation correction for the first CN */
  if(ctl.fluctuation){
    statWidthFluctuationReset(tsum);
    specTransmissionSum(wfactor,tstat,k0,p0,j0,tg,tc,td,0.0,&ncl[c0],spc);
  }

  /*** For entrance channel Tlj coupled to JP */
  for(int lp0=0 ; lp0<=tin->lmax*2 ; lp0+=2){
    if(parity(lp0) != p0*paritarg) continue;

#ifdef SelectJPi
    if(lp0 != SelectJPi*2) continue;
#endif

    for(int sp0=spinproj ; sp0>=-spinproj ; sp0-=2){
      int jp0 = lp0 + sp0;
      if(jp0 < 0) continue;

      if( std::abs(jp0-spintarg) > j0 || j0 > (jp0+spintarg) ) continue;

      double tlj  = tin->tran[tj_index(lp0,sp0,spinproj)];
      if(tlj == 0.0) continue;

      if(ctl.fluctuation) statWidthFluctuationSet(lp0,jp0,incid,targlev,tlj);

      /*** Population add to each residual CN */
      double pini = g * tlj;
      double dp = pini/tsum;

      Statcalcmode scm = (ctl.fluctuation) ? fluctuation : hauser;
      specTransmissionSum(scm,tstat,k0,p0,j0,tg,tc,td,dp,&ncl[c0],spc);

      if(prn.angdist){
        specLegendreCoefficient(targid,incid,targlev,dp,p0,j0,lp0,jp0,td);
        if(opt.continuumangdist) specLegendreCoefficientContinuum(dp,p0,j0,lp0,jp0,tc,cleg);
      }

      int jdx = (j0-spings)/2;
      if(p0 > 0) ncl[c0].pop[k0][jdx].even += pini;
      else       ncl[c0].pop[k0][jdx].odd  += pini;
    }
  }

  return 0;
}


/**********************************************************/
/*      CN Decay for Deformed Case                        */
/*      --------                                          */
/*      width fluctuation always included                 */
/**********************************************************/
int specCNDecayDeformed
(const int c0,      // CN index
 const int k0,      // CN initial excitation energy index
 const int j0,      // CN spin
 const int p0,      // CN parity
 const double g,    // spin statistical factor
 Transmission **tc, Transmission **td, double **tg, Spectra *spc, int *bidx)
{
  int c1 = ncl[c0].cdt[incid].next;

  /*** calculate penetration matrix P, and diagonalize */
  int m = ccHermiteMatrix(j0,p0); if(m < 0) return -1;

  /*** Sum of the eigenvalues */
  double tsum = 0.0;
  for(int i=0 ; i<m ; i++) tsum += gHermiteEigenvalue[i];

  /*** sum over gamma and particle decay, except for the coupled levels
       this is distinguished by lmax = 0 */
  tsum += specTransmissionSum(sumall,tstat,k0,p0,j0,tg,tc,td,0.0,&ncl[c0],spc);
  if(tsum == 0.0) return 0;

  /*** width fluctuation factor in the channel space */
  statWidthFluctuationReset(tsum);
  for(int j=0 ; j<m ; j++) statWidthFluctuation(gHermiteEigenvalue[j],1.0,0.0);
  specTransmissionSum(wfactor,tstat,k0,p0,j0,tg,tc,td,0.0,&ncl[c0],spc);

  /*** average cross section in the channel space, coupled levels only */
  for(int i=0 ; i<m ; i++){
    double dp = g*gHermiteEigenvalue[i]/tsum;

    /*** channel-space cross sections for coupled channels */
    for(int j=0 ; j<=i ; j++){
      double ef = statElasticEnhancement(gHermiteEigenvalue[j]);
      double w1 = statMoldauerCC(gHermiteEigenvalue[i],gHermiteEigenvalue[j]);
      gChannelSigma[i*(i+1)/2 + j] = gHermiteEigenvalue[j] * w1 * dp * ((i == j) ? ef : 1.0);
    }
  }

  /*** cross sections for coupled levels,  looking at in-coming channel only */
  std::complex<double> x1(0.0,0.0), x2(0.0,0.0), xp(0.0,0.0);

  for(int a=0 ; a<m ; a++){ if(gCdt[a].level != gCol.target_index) continue;

    /*** outgoing channels in the coupled levels */
    for(int b=0 ; b<m ; b++){

      /*** back-transformation of Engelbrecht-Weidenmueller, by Moldauer */
      double pop = 0.0;

      std::complex<double> u1(0.0,0.0), u2(0.0,0.0);

      for(int c=0 ; c<m ; c++){
        int ca = a*m + c;
        int cb = b*m + c;

        pop += norm(gHermiteEigenvector[ca]) * norm(gHermiteEigenvector[cb])
             * gChannelSigma[c*(c+1)/2+c];

        for(int d=0 ; d<m ; d++){ if(c == d) continue;

          int da = a*m + d;
          int db = b*m + d;

          int cd = (c >= d) ? c*(c+1)/2+d : d*(d+1)/2+c;

          /*** <Saa Sbb^*> by Moldauer */
          double sm = statMoldauerCorrelation(gHermiteEigenvalue[c],gHermiteEigenvalue[d]);
          /*** phase of Saa Sbb^* */
          double th = gSPhase[c] - gSPhase[d];
          xp = std::complex<double>(cos(th),sin(th));
          
          x1 = conj(gHermiteEigenvector[ca]) * conj(gHermiteEigenvector[db])
               * (  gHermiteEigenvector[ca]  *      gHermiteEigenvector[db]
                  + gHermiteEigenvector[da]  *      gHermiteEigenvector[cb] );

          x2 = conj(gHermiteEigenvector[ca]) * conj(gHermiteEigenvector[cb])
                  * gHermiteEigenvector[da]  *      gHermiteEigenvector[db];
          
          u1 += x1 * gChannelSigma[cd];
          u2 += (x2 * sm * xp) * gChannelSigma[cd];
        }
      }
      pop += u1.real() + u2.real();

      /*** add to discrete level population */
      for(int k=0 ; k<ncl[c1].ndisc ; k++){
        if(bidx[k] == gCdt[b].level){
          ncl[c1].lpop[k] += pop; break;
        }
      }
    }
  }

  /*** cross sections for uncoupled levels */
  double pini = 0.0;
  for(int a=0 ; a<m ; a++){ if(gCdt[a].level != gCol.target_index) continue;

    for(int c=0 ; c<m ; c++){
      /*** |Uca|^2 x p[a] / sum_k p[k] */
      double dp = g*gHermiteEigenvalue[c]/tsum * norm(gHermiteEigenvector[a*m + c]);
      if(dp != 0.0){
        statWidthFluctuationSetT0(gHermiteEigenvalue[c]);
        specTransmissionSum(fluctuation,tstat,k0,p0,j0,tg,tc,td,dp,&ncl[c0],spc);
        pini += dp;
      }
    }
  }

  /*** initial population, cross section sum over exit channel */
  int jdx = (j0-spings)/2;
  if(p0 > 0) ncl[c0].pop[k0][jdx].even += pini;
  else       ncl[c0].pop[k0][jdx].odd  += pini;

  return 0;
}


/**********************************************************/
/*      Compound Elastic Scattering, Pre-Process          */
/**********************************************************/
void    specCNDecayPreProcess(System *sys, const int c0)
{
  targid   = sys->target_id;
  incid    = sys->inc_id;
  targlev  = sys->target_level;

  spinproj = ncl[c0].cdt[incid].spin2;                 // projectile spin
  spings   = (int)(2.0*halfint(ncl[c0].lev[0].spin));  // CN ground state spin
  spintarg = (int)(2.0*ncl[targid].lev[targlev].spin); // target state spin
  paritarg = ncl[targid].lev[targlev].parity;          // target state parity

  if(opt.continuumangdist) specContinuumLegendreAllocate();
}


/**********************************************************/
/*      Compound Elastic Scattering, Post-Process         */
/**********************************************************/
void    specCNDecayPostProcess(Transmission **td, Spectra *spc)
{
  const int c0 = 0;

  /*** fission probability for binary reaction */
  if(ctl.exclusive && ctl.fission){
    specLostPopFraction(c0,0,crx.reaction + crx.totaldir,tstat,spc->dp);
  }

  /*** Compound elastic/inelastic, discrete population cross sections before cascading */
  crx.compound = ncl[targid].lpop[targlev];
  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!tstat[id]) continue;
    int c1 = ncl[c0].cdt[id].next;
    for(int k=0 ; k<ncl[c1].ndisc ; k++) crx.levexcite[id][k] = ncl[c1].lpop[k];
  }

  /*** Reconstruct differential cross sections from Legendre coefficients */
  if(prn.angdist){
    /*** First, output angular distributions for the inelastic scattering */
    int n = statInelAngularDistribution(incid,ncl[c0].jmax,td[incid]);
    outAngularDistribution(3,0,n,ANGLE_STEP,&ncl[targid].za);
    outLegendreCoefficient(3,0,n,incid,&ncl[targid].za,crx.legcoef);

    for(int id=1 ; id<MAX_CHANNEL ; id++){
      int c1 = ncl[c0].cdt[id].next;
      if(!tstat[id]) continue;
      if(id == incid) continue;

      if(!ctl.fileout){
        /*** since ctl.angdist will be overwritten, 
             this is unable when fileout option is ON */
        n = statAngularDistribution(id,ncl[c0].jmax,td[id]);
        outAngularDistribution(3,0,n,ANGLE_STEP,&ncl[c1].za);
      }

      for(n=0 ; n<MAX_ANGDISTLEVELS ; n++) if(td[id][n].lmax <= 0) break;
      outLegendreCoefficient(3,0,n,id,&ncl[c1].za,crx.legcoef);
    }

    /*** angular distribution in continuum, optional calculation */
    if(opt.continuumangdist){
      for(int id=1 ; id<MAX_CHANNEL ; id++){
        int c1 = ncl[c0].cdt[id].next;
        if(!tstat[id]) continue;

        outLegendreCoefficient(4,0,ncl[c1].ncont,id,&ncl[c1].za,cleg);
        specContinuumAngularDistribution(id,ncl[c0].jmax,&ncl[c1]);
      }
      specContinuumLegendreFree();
    }
  }
}


/**********************************************************/
/*      Pring Continuum Angular Distribution              */
/**********************************************************/
void specContinuumAngularDistribution(const int id, const int jmax, Nucleus *n)
{
  for(int k=0 ; k<n->ncont ; k++){
    for(int i=0 ; i<MAX_ANGDIST ; i++){
      double x = 0.0;
      for(int l=0 ; l<jmax ; l+=2) x += cleg[id][k][l] * legendre(l,crx.costh[i]);
      std::cout << std::setw(11) << cleg[id][k][MAX_J];
      std::cout << std::setw(11) << crx.theta[i];
      std::cout << std::setw(11) << x << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}


/**********************************************************/
/*      Prepare Legendre Coeff. in Continuum              */
/**********************************************************/
void specContinuumLegendreAllocate()
{
  cleg = new double ** [MAX_CHANNEL];
  for(int id=0 ; id<MAX_CHANNEL ; id++){
    cleg[id] = new double * [MAX_ENERGY_BIN];
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
      cleg[id][k] = new double [MAX_J+1]; // store excitation energy at MAX_J
      for(int j=0 ; j<MAX_J+1 ; j++) cleg[id][k][j] = 0.0;
    }
  }

  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;
    int c1 = ncl[0].cdt[id].next;
    for(int k=0 ; k<ncl[c1].ncont ; k++){
      cleg[id][k][MAX_J] = ncl[c1].excitation[k];
    }
  }
}

void specContinuumLegendreFree()
{
  for(int id=0 ; id<MAX_CHANNEL ; id++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++)  delete [] cleg[id][k];
    delete cleg[id];
  }
  delete [] cleg;
}




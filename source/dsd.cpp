/******************************************************************************/
/*  dsd.cpp                                                                   */
/*        direct/semidirect nucleon capture calculation                       */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "structur.h"
#include "optical.h"
#include "dsd.h"
#include "nucleus.h"
#include "coupling.h"
#include "etc.h"
#include "global.h"
#include "output.h"
#include "statmodel.h"
#include "parameter.h"
#include "terminate.h"


// gamma-ray spectrum broadening witdth, to be 100 keV
static const double gBroadeningWidth = 0.0;

static bool deformedDSD = false;

static void dsdCalc (const int, CCdata *, Bstate *, std::complex<double> *, double *);
static void dsdCalcHF (CCdata *, Bstate *, std::complex<double> *, double *, HFInterface *);
static bool dsdAmplitudeSpherical (int, int, double *, CCdata *, std::complex<double> *, Bstate *);
static bool dsdAmplitudeDeformed (const int, const int, const int, const int, double *, CCdata *, std::complex<double> *, Bstate *);
static void dsdSetGDRParameter (const double, const double, GDR *);
static void dsdSetConstant (const Particle, const int, const int, Optical *, CCdata *, const double);
static void dsdSpecBroadening (double *, double *, const int, const double, const double);
static inline void dsdHOWaveFunction (const double, const double, Bstate *);
static inline std::complex<double> phasefactor(const int);


static double  const0 = 0.0, const1 = 0.0, const2 = 0.0;
static double  const30 = 0.0, const31 = 0.0;
static double  integwidth = 0.0;
static int     ninteg = 0;

static std::complex<double> *f1=nullptr, *f2=nullptr;
static double  sig[3];
static GDR     gdr0,gdr1;

extern double *fact;


/**********************************************************/
/*      DSD Calculation for Nucleon Capture               */
/**********************************************************/
int     dsdDirectCaptureModel
(const int     incid,         // index for incident channel
 const double  energy,        // LAB incident energy
 Pdata        *proj,          // particle data for projectile
 ZAnumber     *targ,          // target ZA
 const double  mu,            // reduced mass for Targ+Proj system
 Dcapt        *dsd,           // DSD parameters
 const double  beta2,         // deformation parameter for double-humped GDR
 double       *gspec)         // gamma-ray energy spectrum
{
  Optical      omp;
  Potential    pot;
  Wavefunc     wfn;
  Bstate      *bst = nullptr;
  CCdata       cdt;
  LevelData    lev;
  std::complex<double> *dwave = nullptr;
  double      *dsdspec = nullptr;

  if((proj->pid != neutron) && (proj->pid != proton)){
    ctl.dsd = false;
    return(0);
  }

//---------------------------------------
//     Memory Allocation

  try{
    pot.memalloc(MAX_POINTS);
    wfn.memalloc(MAX_L,MAX_POINTS);

    bst            = new Bstate [MAX_SPLEVEL];
    for(int i=0 ; i<MAX_SPLEVEL ; i++){
      bst[i].wave = new double [MAX_POINTS_BW];
    }
    f1             = new std::complex<double> [MAX_POINTS];
    f2             = new std::complex<double> [MAX_POINTS];
    dwave          = new std::complex<double> [MAX_POINTS*MAX_L*3];
    dsdspec        = new double          [MAX_ENERGY_BIN];
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("dsdDirectCaptureModel");
  }

  int zzprod = targ->getZ()*proj->za.getZ();
  ninteg = MAX_POINTS_BW-1;

  if(dsd->nonloc == 0.0) dsd->nonloc = 0.85;


//---------------------------------------
//      Calc. All Bound State and Bound Potential

  double binding = ncl[0].cdt[incid].binding_energy;

  /*** set bound state potential */
  /*** energy set to 1eV to avoid zero-div for proton case */
  cdt.lev = &lev;
  omSetEnergy(1e-06, zzprod, mu, &lev);

  unsigned int bindex = find_omp("bound");
  omSetOmp(bindex,0.0,targ->getZ(),targ->getA(),proj->za.getZ(),proj->za.getA(),&omp);

  pot.rad_match = 30.0;
  pot.width     = integwidth = INTEG_WIDTH;
  pot.n_match   = (int)(pot.rad_match / integwidth);
  if(pot.n_match > MAX_POINTS-2) pot.n_match = MAX_POINTS-2;
  pot.setRho(lev.wave_number);

  omPotentialFixedLength(zzprod,&omp,mu,&pot);
  
  /*** calc. single particle wavefunctions */
  ninteg = bstateWavefunc(proj,targ,binding,mu,dsd->nonloc,dsd->bmax,dsd->fshift,&omp,&pot,bst);


//---------------------------------------
//      Entrance Channel Optical Potential
  omSetEnergy(energy, zzprod, mu, &lev);

  omSetOmp(proj->omp,lev.energy,targ->getZ(),targ->getA(),proj->za.getZ(),proj->za.getA(),&omp);

  pot.width    = integwidth = INTEG_WIDTH;
  pot.n_match  = ninteg-4;
  pot.rho_match = pot.rad_match*sqrt(cdt.lev->wavesq);

  omPotentialFixedLength(zzprod,&omp,mu,&pot);

  /*** Store entrance channel distorted waves */
  dwStoreDistortedWave(proj->spin2,dsd->nonloc,&cdt,&pot,&wfn,dwave);


//---------------------------------------
//      GDR parameters

  dsdSetGDRParameter(beta2,(double)targ->getA(),dsd->gdr);


//---------------------------------------
//      V0 systematics

  if(dsd->v1 <= 0.0){
    dsd->v1 = 1.0/(4e-5*targ->getA()+0.0014);
    dsd->w1 = 0.0;
  }

  double v1s = parmGetFactor(parmDSDV);

  dsdFormfactor(VIB_FORMFACTOR,&omp,&pot,dsd->v1*v1s,dsd->w1,f1,f2);
  cdt.lmax = ncl[0].jmax;


//---------------------------------------
//      DSD calculation

  dsdSetConstant(proj->pid,targ->getZ(),targ->getA(),&omp,&cdt,dsd->v1*v1s);

  /*** reset cross section array */
  sig[0] = sig[1] = sig[2] = 0.0;

  int nlev = 0;
  for(nlev=0 ; nlev<MAX_SPLEVEL ; nlev++) if(bst[nlev].bind > 0.0) break;

  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) dsdspec[k] = 0.0;

  if(dsd->hf) dsdCalcHF(&cdt,bst,dwave,dsdspec,dsd->hfsp);
  else        dsdCalc(nlev,&cdt,bst,dwave,dsdspec);

  dsdSpecBroadening(dsdspec,gspec,ncl[0].ncont,ncl[0].de,gBroadeningWidth);


//---------------------------------------
//      Output resulsts, set cross section

  if(prn.xsection) outDSD(sig[0],sig[1],sig[2],dsd);
  crx.dsd = sig[0];

  /*** DSD final state is assumed to be the ground state */
  ncl[0].lpop[0] += crx.dsd;


//---------------------------------------
//      Free Allocated Memory

  for(int i=0 ; i<MAX_SPLEVEL ; i++){
    delete []  bst[i].wave;
  }
  delete [] bst;
  delete [] f1;
  delete [] f2;
  delete [] dwave;
  delete [] dsdspec;

  return(0);
}


/**********************************************************/
/*      Direct /Semi-Direct Capture Calculation           */
/**********************************************************/
void    dsdCalc(const int nlev, CCdata *cdt, Bstate *bst, std::complex<double> *dwave, double *gspec)
{
  double  sigc[3], siglev[3];
  std::complex<double> *wf;

  const int spin2 = 1;
  int smax = spin2 + 1;

  /*** number of bound orbit */
  for(int ip=0 ; ip<nlev ; ip++){
    if(bst[ip].w == 0.0) continue;

    double eg  = cdt->lev->energy - bst[ip].bind;
    double kg  = eg/HBAR/VLIGHT;
    double c   = kg*kg*kg*NORM_FACT;

    /*** sum cross sections for different (l,j) */
    for(int i=0 ; i<3 ; i++) siglev[i] = 0.0;
    for(int l=0 ; l<=cdt->lmax ; l++){

      for(int s=0 ; s<smax ; s++){
        int ss = spin2 - 2*s;
        int j2 = 2*l + ss;
        if(j2 < 0) continue;

        wf = &dwave[(3*l+s)*MAX_POINTS];

        bool p = false;
        if(deformedDSD) p = dsdAmplitudeDeformed(nlev,ip,l*2,j2,sigc,cdt,wf,bst);
        else            p = dsdAmplitudeSpherical(l*2,j2,sigc,cdt,wf,&bst[ip]);

        if(p){
          for(int i=0 ; i<3 ; i++) siglev[i] += c*sigc[i];
        }
      }
    }

    for(int i=0 ; i<3 ; i++) sig[i] += siglev[i];

    int k = specFindEnergyBin(eg,ncl[0].de);
    if(k >= 0) gspec[k] += siglev[0]/ncl[0].de;
  }
}


/**********************************************************/
/*      Direct /Semi-Direct for Hartree-Fock Results      */
/**********************************************************/
void    dsdCalcHF(CCdata *cdt, Bstate *bst, std::complex<double> *dwave, double *gspec, HFInterface *hf)
{
  const double ucut = 1.0e-8;

  Bstate  bho;
  bho.wave = new double [MAX_POINTS_BW];

  /*** for each Hartree-Fock state */
  for(int ih=0 ; ih<hf->n ; ih++){
    if(hf->state[ih].energy > 0.0) break;

    /*** u = sqrt(1-v^2) */
    double u = sqrt(1.0 - hf->state[ih].v2);
    if(u < ucut) continue;

    /*** L and J in the spherical nucleus */
    int nlev = 0;
    for(int l=0 ; l<=hf->state[ih].lmax ; l++){
      for(int s=-1 ; s<=1 ; s+=2){
        int j2 = 2*l+s; if(j2 < 0) continue;

        if(nlev < MAX_SPLEVEL){
          for(int i=0 ; i<=ninteg ; i++) bst[nlev].wave[i] = 0.0;

          /*** look for the expanded states with the same L,J */
          bool flag = false;
          for(int is=0 ; is<hf->state[ih].m ; is++){
            if( (hf->state[ih].q[is].l == l) && (hf->state[ih].q[is].j == j2) ){

              bho.n  = hf->state[ih].q[is].n;
              bho.l  = l;
              bho.j2 = j2;
              flag   = true;

              dsdHOWaveFunction(hf->b,integwidth,&bho);

              for(int i=0 ; i<=ninteg ; i++){
                bst[nlev].wave[i] += bho.wave[i] * hf->state[ih].q[is].c * u;
              }
            }
          }
          if(flag){
            bst[nlev].bind = hf->state[ih].energy;
            bst[nlev].k2   = hf->state[ih].k2;
            bst[nlev].l    = l;
            bst[nlev].j2   = j2;
            bst[nlev].w    = sqrt(2.0/(j2+1.0));
            nlev++;
          }
        }
      }
    }
    dsdCalc(nlev,cdt,bst,dwave,gspec);
  }

  delete [] bho.wave;
}


/**********************************************************/
/*      DSD Amplitude for Spherical Nucleus               */
/**********************************************************/
bool    dsdAmplitudeSpherical(int l2, int j2, double *x, CCdata *cdt, std::complex<double> *wf, Bstate *bw)
{
  double  z;
  std::complex<double> tdir,tcol,tall;

  x[0] = x[1] = x[2] = 0.0;
  if( (z = z_coefficient(l2,j2,2*bw->l,bw->j2,1,2)) == 0.0 ) return(false);

  std::complex<double> rmtd = dsdRadialIntegral(ninteg,integwidth,f1,bw->wave,wf);
  std::complex<double> rmtc = dsdRadialIntegral(ninteg,integwidth,f2,bw->wave,wf);
  std::complex<double> cnu  = std::complex<double>(const30,0.0)/std::complex<double>(cdt->lev->energy-gdr0.getEnergy()-bw->bind,gdr0.getWidth()/2.0);

  tdir = rmtd       * z * const1 * bw->w;
  tcol = rmtc * cnu * z * const2 * bw->w;
  tall = tdir + tcol;

  x[0] = 0.5*const0*norm(tall);  /* total       */
  x[1] = 0.5*const0*norm(tdir);  /* direct      */
  x[2] = 0.5*const0*norm(tcol);  /* semi-direct */

  return(true);
}


/**********************************************************/
/*      DSD Amplitude for Deformed Nucleus                */
/**********************************************************/
bool    dsdAmplitudeDeformed(const int nlev, const int ip, const int l2, const int j2, double *x, CCdata *cdt, std::complex<double> *wf, Bstate *bw)
{
  std::complex<double> tdir,tcol,tall,cnu,cnu0,cnu1,cnu2;

  x[0] = x[1] = x[2] = 0.0;

  double z1 = z_coefficient(l2,j2,2*bw[ip].l,bw[ip].j2,1,2)/sqrt(bw[ip].j2+1.0);

  std::complex<double> rmtd(0.0,0.0);
  if(z1 != 0.0) rmtd = dsdRadialIntegral(ninteg,integwidth,f1,bw[ip].wave,wf);

  tdir = rmtd * phasefactor(bw[ip].l) * z1 * const1;

  double cg00 = clebsh_gordan(2,j2, 0,bw[ip].k2  ,bw[ip].j2);
  double cg01 = clebsh_gordan(2,j2, 2,bw[ip].k2-2,bw[ip].j2);
  double cg02 = clebsh_gordan(2,j2,-2,bw[ip].k2+2,bw[ip].j2);

  tcol = std::complex<double>(0.0,0.0);
  for(int iq=0 ; iq<nlev ; iq++){

    double z2 = z_coefficient(l2,j2,2*bw[iq].l,bw[iq].j2,1,2)/sqrt(bw[iq].j2+1.0);
    if(z2 == 0.0) continue;

    std::complex<double> rmtc(0.0,0.0);
    rmtc = dsdRadialIntegral(ninteg,integwidth,f2,bw[iq].wave,wf);
    rmtc *= phasefactor(bw[iq].l);

    cnu1 =
    cnu2 = std::complex<double>(const31,0.0)
          /std::complex<double>(cdt->lev->energy-gdr1.getEnergy()-bw[iq].bind,gdr1.getWidth()/2.0);
    cnu0 = std::complex<double>(const30,0.0)
          /std::complex<double>(cdt->lev->energy-gdr0.getEnergy()-bw[iq].bind,gdr0.getWidth()/2.0);

    double cg10 = clebsh_gordan(2,j2, 0,bw[iq].k2  ,bw[iq].j2)*cg00;
    double cg11 = clebsh_gordan(2,j2, 2,bw[iq].k2-2,bw[iq].j2)*cg01;
    double cg12 = clebsh_gordan(2,j2,-2,bw[iq].k2+2,bw[iq].j2)*cg02;

    cnu = cg10*cnu0 + cg11*cnu1 + cg12*cnu2;
    tcol += rmtc*cnu * z2 * const2;
  }

  tall = tdir + tcol;

  x[0] = const0*norm(tall); /* total       */
  x[1] = const0*norm(tdir); /* direct      */
  x[2] = const0*norm(tcol); /* semi-direct */

  return(true);
}


/**********************************************************/
/*      Local GDR Parameters                              */
/**********************************************************/
void dsdSetGDRParameter(const double beta2, const double targA, GDR *gdr)
{
  /*** for spherical calculation, enforce single humped GDR */
  if(!deformedDSD){
    gdrE1(targA,beta2,&gdr0);
    gdr1.clear();
  }
  else{
    /*** copy GDR parameter to local variables */
    std::string XL = "E1";
    gdr0.setGDR(XL,gdr[0].getEnergy(),gdr[0].getWidth(),gdr[0].getSigma(),gdr[0].getProfile());

    if(beta2 > 0.0 && gdr[1].getXL() == XL){
      gdr1.setGDR(XL,gdr[1].getEnergy(),gdr[1].getWidth(),gdr[1].getSigma(),gdr[1].getProfile());
    }
    else gdr1.clear();
  }
}


/**********************************************************/
/*      Local Pre-Calculated Constants                    */
/**********************************************************/
void dsdSetConstant(const Particle p, const int targZ, const int targA, Optical *omp, CCdata *cdt, const double v1)
{
  double const3;

  double numZ = (double)targZ;
  double numA = (double)targA;
  double numN = (double)(targA - targZ);

  double e2  = PERMITTIV * COULOMBSQ;
  double msr = 0.2*(3.0*omp->R0*omp->R0 + 7.0*PI*PI*omp->a0*omp->a0);

  double k   = sqrt(cdt->lev->wavesq);

  const0  = 16*PI/9.0 * cdt->lev->reduced_mass * AMUNIT/VLIGHTSQ/HBARSQ / (k*k*k);
  const1  = sqrt(e2) * ((p==neutron) ? -numZ/numA : numN/numA);
  const2  = ((p==neutron) ? -1 :1) * sqrt(e2)
           * (numN*numN*numZ*numZ) / (numA*numA*numA) / 2.0;
  if(VIB_FORMFACTOR != 0) const2 *= 3.0/msr;
  const3  = (numA*numA)/(numN*numN*numZ*numZ) * HBAR*VLIGHT / e2 /(PI4*PI);

  if(v1 == 0.0){ const30 = const31 = 0.0; }
  else{
    if(deformedDSD){
      const30 = const3*gdr0.getSigma()/NORM_FACT *PI_2*gdr0.getWidth()/gdr0.getEnergy();
      const31 = const3*gdr1.getSigma()/NORM_FACT *PI_2*gdr1.getWidth()/gdr1.getEnergy();
    }
    else{
      const30 = const3*gdr0.getSigma()/NORM_FACT *PI_2*gdr0.getWidth()/gdr0.getEnergy();
      const31 = 0.0;
    }
  }
}


/**********************************************************/
/*     Add to Spectrum                                    */
/**********************************************************/
void   dsdSpecBroadening(double *ds, double *gs, const int ne, const double de, const double w)
{
  if(w == 0.0){
    for(int k=0 ; k<ne ; k++) gs[k] += ds[k];
  }
  else{
    /*** bin number for 3-sigma width */
    int dk = (int)(2.0*w/de);
  
    for(int k1=0 ; k1<ne ; k1++){
      if(ds[k1] == 0.0) continue;

      double e1 = k1*de;
      double sumg = 0.0;
      for(int k2=k1-dk ; k2<k1+dk ; k2++){
        if( (k2 < 0) || (k2 >= ne) ) continue;
        double e2 = k2*de;
        double wg = 1.0/((e1-e2)*(e1-e2) + w*w);
        sumg += wg;
      }

      if(sumg > 0.0){
        for(int k2=k1-dk ; k2<k1+dk ; k2++){
          if( (k2 < 0) || (k2 >= ne) ) continue;
          double e2 = k2*de;
          double wg = 1.0/((e1-e2)*(e1-e2) + w*w) / sumg;
          gs[k2] += wg * ds[k1];
        }
      }
    }
  }
}


/**********************************************************/
/*     Bound Wave Calculation (harmonic oscillator)       */
/**********************************************************/
inline void dsdHOWaveFunction(const double b, const double d, Bstate *bst)
{
  double r,c,p,q;
  int n = bst->n;
  int l = bst->l;

  for(int i=0 ; i<=ninteg ; i++){
    r = i*d;
    p = r/b;
    c = sqrt(2.0*exp(fact[n])/gam((l+n)+1.5)/b)/b;
    q = (p==0.0 && l==0) ? 1.0 : pow(p,(double) l);
    bst->wave[i] = c * q * laguerre(n,(double)l+0.5,p*p) * exp(-p*p/2.0)*r;
  }
}


/**********************************************************/
/*     Complex Phase Factor                               */
/**********************************************************/
inline std::complex<double> phasefactor(const int l)
{
  double preal, pimag;

  switch((l+1)%4){
  case  0: preal =  1.0; pimag =  0.0; break;
  case  1: preal =  0.0; pimag = -1.0; break;
  case  2: preal = -1.0; pimag =  0.0; break;
  case  3: preal =  0.0; pimag =  1.0; break;
  default: preal =  1.0; pimag =  0.0; break;
  }

  return( std::complex<double>(preal,pimag) );
}


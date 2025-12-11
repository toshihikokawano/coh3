/******************************************************************************/
/*  eclddx.cpp                                                                */
/*        calculate double differential cross sections for ENDF-6 format      */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "physicalconstant.h"
#include "structur.h"
#include "nucleus.h"
#include "eclipse.h"
#include "ddxkalbach.h"
#include "terminate.h"

static void   ddxMainExclusive (const int, const int, const int, const double, double **, double **, EXSpectra *);
static int    ddxEnergyArray (const int, const int, const double, double *, EXSpectra *);
static double ddxDiscreteFraction (const int, EXSpectra *);
static int    ddxDiscreteGamma (int, double *, double *, const double, EXSpectra *);
static void   ddxNormalizeContGamma (const int, const int, const double, double *, double **);
static void   ddxContGammaMultiplicity (const int, double *, int, double *, double *, EXSpectra *, Nucleus *);


static double   ***cleg;      // Legendre expansion coefficients
static double    **epar;      // secondary particle energies
static int        *kmax;      // number of secondary particle energy points

static bool checkbalance = false;

/**********************************************************/
/*      Double Differential Cross Sections                */
/**********************************************************/
void eclDDX(const int km, System *sys, double **np, double **fmsd, EXSpectra *dat)
{
  try{
    /*** cleg[Particle][LegOrder][Energy] */
    cleg = new double ** [sys->max_channel];
    for(int c=0 ; c<sys->max_channel ; c++){
      cleg[c] = new double * [NLEG];
      for(int l=0 ; l<NLEG ; l++)  cleg[c][l] = new double [km];
    }

    /*** kmax[Energy], epar[Particle][Energy] */
    kmax  = new int      [sys->max_channel];
    epar  = new double * [sys->max_channel];
    for(int c=0 ; c<sys->max_channel ; c++){
      epar[c] = new double [NCNT];
    }
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("eclDDX");
  }

  /*** Generate DDX from Kalbach systematics */
  ddxKalbachSetParm(sys->cms_energy,sys->compound,sys->incident.za);
  ddxMainExclusive(sys->max_compound,km,sys->max_channel,sys->energy_bin,np,fmsd,dat);


  for(int c=0 ; c<sys->max_channel ; c++){
    delete [] epar[c];

    for(int l=0 ; l<NLEG ; l++)  delete [] cleg[c][l];
    delete [] cleg[c];
  }
  delete [] cleg;
  delete [] epar;
  delete [] kmax;
}


/**********************************************************/
/*      Main Part for DDX and Gamma-Ray Spectra           */
/**********************************************************/
void ddxMainExclusive(const int nm, const int km, const int cm, const double de, double **np, double **fmsd, EXSpectra *dat)
{
  double *gx, *gy;
  double *ctmp = new double [NLEG];
  int    *nleg, ngam;

  nleg = new int [cm];
  gx = new double [NGLN];
  gy = new double [NGLN];

  /*** for each of compound nuclei */
  for(int n=0 ; n<nm ; n++){

    for(int c=0 ; c<cm ; c++){
      kmax[c] = 0;
      nleg[c] = NLEG;
    }
    ngam = 0;

    for(int i=0 ; i<NGLN ; i++){ gx[i] = gy[i] = 0.0; }

    /*** find non-zero element from the high-side */
    int  km2  = eclFindKmax(km,cm,dat[n].spec);

    /*** for each of particle emission channels */
    for(int c=0 ; c<cm ; c++){

      if(km2 <= 0) continue;
      if( (c > 0) && (np[n][c] == 0.0) ) continue;

      /*** outgoing energy array */
      kmax[c] = ddxEnergyArray(c,km2,de,epar[c],&dat[n]);

      for(int l=0 ; l<NLEG ; l++){
        for(int k=0 ; k<kmax[c] ; k++) cleg[c][l][k] = 0.0;
      }

      for(int k=0 ; k<kmax[c]-2 ; k++){
        if( (c == 0) || (k > dat[n].getNc()) ) cleg[c][0][k+1] = dat[n].spec[c][k];
        else{
          ddxKalbach(c,epar[c][k+1],fmsd[c][k],NLEG,ctmp);

          /*** convert into ENDF-6 data, F(E') = 2 f/(2L+1) x (ds/dE) */
          for(int l=0 ; l<NLEG ; l++) cleg[c][l][k+1] = ctmp[l] / (l+0.5) * dat[n].spec[c][k];
        }
      }

      double frac = 0.0;
      if(c == 0){
        nleg[c] = 1;
        frac = ddxDiscreteFraction(km2,&dat[n]);
        ngam = ddxDiscreteGamma(n,gx,gy,frac,dat);
        if(ngam == 0) frac = 0.0;
        //if(n == 0 && npgam > 0) kmax[c] = ddxPrimaryGamma(kmax[c],epar[c]);
      }

      ddxNormalizeContGamma(kmax[c],nleg[c],frac,epar[c],cleg[c]);
    }

    ddxContGammaMultiplicity(cm,np[n],ngam,gx,gy,&dat[n],&ncl[n]);

    eclOutNucleusHead(n,cm,np[n]);

    if(checkbalance){
      eclCheckEnergyBalance(cm,np[n],kmax,epar,ngam,gx,gy,cleg,&dat[n],&ncl[n]);
      eclTotalEnergySpectrum(cm,km,crx.prod[n].xsec,np[n],&dat[n],&ncl[n]);
    }
    else{
      /*** check if no particle spectrum */
      if(n != 0){
        bool zero = true;
        for(int c=1 ; c<cm ; c++){
          for(int k=0 ; k<kmax[c] ; k++){
            for(int l=0 ; l<nleg[c] ; l++) if(cleg[c][l][k] != 0.0) zero = false;
            if(!zero) break;
          }
        }
        if(zero){
          for(int c=0 ; c<cm ; c++) np[n][c] = 0.0;
        }
      }

      for(int c=0 ; c<cm ; c++){
        eclOutChannelHead(c,np[n][c]);
        if(np[n][c] == 0.0) continue;
        if(kmax[c] > 0){
          if(c == 0) eclOutGammaLine(ngam,nleg[c],gx,gy);
          eclOutLegCoeff(kmax[c],nleg[c],epar[c],cleg[c]);
        }        
        else{
          eclOutLegCoeff(-1,0,epar[c],cleg[c]);
        }
      }
    }

  }

  delete [] ctmp;
  delete [] nleg;
  delete [] gx;
  delete [] gy;
}


/**********************************************************/
/*      Set Up Outgoing Secondary Particle Energies       */
/**********************************************************/
int ddxEnergyArray(const int c, const int km2, const double de, double *ep, EXSpectra *n)
{
  ep[0] = 0.0;
  ep[1] = de/4.0;

  double epmax = n->getEm();

  /*** check if this is a binary reaction (n,n'), (n,p), (n,alpha),
       and remove discrete level population part */

  if( n->isBinary() && (c != 0) ) epmax -= n->getEl();
  if(epmax <= ep[1]) return (0);

  /*** generate energy bins for spectrum */
  bool emx = true;
  int  km3 = km2+1;
  for(int k=1 ; k<km2 ; k++){
    ep[k+1] = de*k;
    if(ep[k+1] >= epmax){
      ep[k+1] = epmax;
      km3 = k+2;
      emx = false;
      break;
    }
  }
  /*** if km2 range does not cover the Epmax, add the max point */
  if(emx) ep[km3++] = epmax;

  return(km3);
}


/**********************************************************/
/*      Renormalize Legendre Coefficients                 */
/**********************************************************/
void ddxNormalizeContGamma(const int km, const int nleg, const double f, double *ep, double **cl)
{
  if(km <= 1) return;

 /*** replace the Km-2 point to adjust the highest bin integral */
  double d1 = ep[km-2]-ep[km-3];
  double d2 = ep[km-1]-ep[km-2];
  double p  = 1.5 * d1/(d1+d2);

  for(int l=0 ; l<nleg ; l++) cl[l][km-2] *= p;

  /*** add all trapezoid pieces */
  double s = 0.0;
  for(int k=0 ; k<km-1 ; k++) s += (cl[0][k+1] + cl[0][k])*(ep[k+1] - ep[k])* 0.5;

  if(s > 0.0) s = (1.0-f)/s;

  /*** renormalize */
  for(int k=0 ; k<km ; k++){
    for(int l=0 ; l<nleg ; l++) cl[l][k] *= s;
  }
}


/**********************************************************/
/*      Fix Gamma Multiplicity                            */
/**********************************************************/
void ddxContGammaMultiplicity(const int cm, double *np, int ngam, double *gx, double *gy, EXSpectra *dat, Nucleus *n)
{
  double *sum = new double [cm];
  double *tot = new double [cm];

  /*** integrate continuum spectra, sum(cleg[0] + gy) are normalized to unity */
  for(int c=0 ; c<cm ; c++){
    sum[c] = tot[c] = 0.0;
    for(int k=0 ; k<kmax[c]-1 ; k++){
      double e0 = epar[c][k  ];
      double e1 = epar[c][k+1];
      double s0 = cleg[c][0][k  ];
      double s1 = cleg[c][0][k+1];

      sum[c] += (s0             + s1            ) * (e1-e0) / 2.0;
      tot[c] += (s0*(e1+2.0*e0) + s1*(e0+2.0*e1)) * (e1-e0) / 6.0;
    }
  }

  /*** add discrete gamma energies */
  for(int k=0 ; k<ngam ; k++){
    sum[0] += gy[k];
    tot[0] += gx[k] * gy[k];
  }
  double gtot1 = tot[0]; // total gamma-ray energy from spectrum

  /*** available energy for this reaction */
  double emax = dat->getEm();

  /*** energy taken by emitted particles */
  double ptot = 0.0;
  for(int c=1 ; c<cm ; c++){
    ptot += tot[c] * np[c];
    if(dat->isBinary()){
      /*** particle energy = Emax - Elev, then Elev is released by gammas,
           so the total energy taken by the particle + gamma = Emax */
      double plev = 0.0;
      for(int k=0 ; k<n->ndisc ; k++) plev += dat->lpop[k].particle;
      ptot += emax * plev * np[c];
    }
  }
  double gtot2 = emax - ptot; // energy available to continuum gammas

  np[0] = (gtot1 > 0.0) ? gtot2/gtot1 : 0.0;
  /*** if no spec given, set multiplicity zero */
  if(np[0] < 0.0) np[0] = 0.0;

/*
  std::cout << " Em" << std::setw(14) << emax;
  std::cout << " Pt" << std::setw(14) << ptot;
  std::cout << " G1" << std::setw(14) << gtot1;
  std::cout << " G2" << std::setw(14) << gtot2 << std::endl;
  std::cout << "#    Sum.Spec      Etotal        Multipl       Mul x Eave" << std::endl;
  for(int c=0 ; c<cm ; c++){
    std::cout << std::setw( 3) << c;
    std::cout << std::setw(14) << sum[c];
    std::cout << std::setw(14) << tot[c];
    std::cout << std::setw(14) << np[c];
    std::cout << std::setw(14) << np[c] * tot[c] << std::endl;
  }
*/
  delete [] sum;
  delete [] tot;
}


/**********************************************************/
/*      Discrete Gamma Fraction                           */
/**********************************************************/
double ddxDiscreteFraction(const int km, EXSpectra *d)
{
  double sc = 0.0, sd = 0.0, st = 0.0;

  for(int k=0 ; k<km ; k++){
    sc += d->spec[0][k];   // continuum gamma-ray spectrum
    sd += d->glin[k];      // spectrum due to dicrete transition only
  }
  st = sc + sd;

  if(st == 0.0) return(0.0);
  else          return(sd / st);
}


/**********************************************************/
/*      Discrete Gamma-Rays                               */
/**********************************************************/
int ddxDiscreteGamma(int n, double *x, double *y, const double f, EXSpectra *dat)
{
  double eg, dp, pop[NLEV];
  int    i0, i1, k = 0, nlev = dat[n].getNl();

  if(f == 0.0) return(0);

  /*** lpop.particle: probability of level population directly by a particle
       lpop.gamma: all other level population, such as continuum transition */
  for(i0=0 ; i0<nlev ; i0++){
    pop[i0] = (dat[n].isBinary()) ? dat[n].lpop[i0].gamma
                                  : dat[n].lpop[i0].particle + dat[n].lpop[i0].gamma;
  }

  /*** calculate discrete gamma lines by gamma-cascade */
  for(i0=nlev-1 ; i0>0 ; i0--){
    for(int j=0 ; j<ncl[n].lev[i0].ngamma ; j++){
      i1 = ncl[n].lev[i0].fstate[j];
      eg = ncl[n].lev[i0].energy - ncl[n].lev[i1].energy;
      dp = pop[i0] * ncl[n].lev[i0].branch[j];

      pop[i1] += dp;

      if(dp == 0.0) continue;

      if(k < NGLN){
        x[k] = eg;
        y[k] = dp;
      }
      k++;
    }
  }

  int ng = k;
  if(ng >= NGLN) ng = NGLN-1;

  /*** sort the gamma lines by energies */
  for(int j=0 ; j<ng ; j++){
    int k = j;
    for(int i=j ; i<ng ; i++){
      if(x[i] > x[k]) k = i;
    }
    eg = x[j];  x[j] = x[k];  x[k] = eg;
    dp = y[j];  y[j] = y[k];  y[k] = dp;
  }

  /*** normalize discrete lines relative to continuum */
  double s = 0.0;
  for(int j=0 ; j<ng ; j++) s += y[j];
  if(s == 0.0) return(0);

  s =  f/s;
  for(int j=0 ; j<ng ; j++) y[j] *= s;

  return(ng);
}



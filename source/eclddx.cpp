/******************************************************************************/
/*  eclddx.cpp                                                                */
/*        calculate double differential cross sections for ENDF-6 format      */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "nucleus.h"
#include "eclipse.h"
#include "terminate.h"
#include "polysq.h"

static void   ddxMainProc (const int, const int, const int, const double, double **, double **, EXSpectra *);
static int    ddxEnergyArray (const int, const int, const double, double *, EXSpectra *);
static double ddxDiscreteFraction (const int, EXSpectra *);
static int    ddxDiscreteGamma (int, double *, double *, const double, EXSpectra *);
static void   ddxCheckEnergyBalance (const int, double *, int, double *, double *, EXSpectra *);
static void   ddxLegNormalize (const int, const int, const double, double *, double **);
static void   ddxKalbach (const int, const double, const double, const double, double *, double *);
static void   ddxKalbachSetParm (void);
static double ddxKalbachAfactor (const int, const int, const double, const double, const double);


static double ecms = 0.0;
static ZAnumber comp(0,0), proj(0,0);

static double   ***cleg;      // Legendre expansion coefficients
static double    **epar;      // secondary particle energies
static int        *kmax;      // number of secondary particle energy points

static double  sa=0.0, sb[7]; // separation energies
static int     incid = 0;     // incident particle ID

/**********************************************************/
/*      Double Differential Cross Sections                */
/**********************************************************/
void eclDDX(const int nm, const int km, const int cm, System *sys, double **np, double **fmsd, EXSpectra *dat)
{
  try{
    /*** cleg[Particle][LegOrder][Energy] */
    cleg = new double ** [cm];
    for(int c=0 ; c<cm ; c++){
      cleg[c] = new double * [NLEG];
      for(int l=0 ; l<NLEG ; l++)  cleg[c][l] = new double [km];
    }

    /*** kmax[Energy], epar[Particle][Energy] */
    kmax  = new int      [cm];
    epar  = new double * [cm];
    for(int c=0 ; c<cm ; c++){
      epar[c] = new double [NCNT];
    }
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("eclDDX");
  }

  /*** Generate DDX from Kalbach systematics */
  comp = sys->compound;
  proj = sys->incident.za;
  ecms = sys->cms_energy;

  ddxKalbachSetParm();
  ddxMainProc(nm,km,cm,sys->energy_bin,np,fmsd,dat);


  for(int c=0 ; c<cm ; c++){
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
void ddxMainProc(const int nm, const int km, const int cm, const double de, double **np, double **fmsd, EXSpectra *dat)
{
  double ddx[NDDX];        // double differential cross section, particle energies
  double *eg, *gx;         // discrete gamma lines
  double ctmp[NLEG];
  int    nleg[cm], ngam;

  eg = new double [NGLN];
  gx = new double [NGLN];


  for(int n=0 ; n<nm ; n++){

    for(int c=0 ; c<cm ; c++){
      kmax[c] = 0;
      nleg[c] = NLEG;
    }
    ngam = 0;

    for(int i=0 ; i<NGLN ; i++){ eg[i] = gx[i] = 0.0; }

    /*** find non-zero element from the high-side */
    int  km2  = eclFindKmax(km,cm,dat[n].spec);
      
    for(int c=0 ; c<cm ; c++){

      if(km2 <= 0) continue;
      if( (c > 0) && (np[n][c] == 0.0) ) continue;

      /*** outgoing energy array */
      kmax[c] = ddxEnergyArray(c,km2,de,epar[c],&dat[n]);

      for(int l=0 ; l<NLEG ; l++){
        for(int k=0 ; k<kmax[c] ; k++) cleg[c][l][k] = 0.0;
      }

      for(int k=0 ; k<kmax[c]-2 ; k++){
        if( (c == 0) || (k > dat[c].getNc()) ) cleg[c][0][k+1] = dat[n].spec[c][k];
        else{
          ddxKalbach(c,epar[c][k+1],fmsd[c][k],dat[n].spec[c][k],ddx,ctmp);
          for(int l=0 ; l<NLEG ; l++) cleg[c][l][k+1] = ctmp[l];
        }
      }

      double frac = 0.0;
      if(c == 0){
        nleg[c] = 1;
        frac = ddxDiscreteFraction(km2,&dat[n]);
        ngam = ddxDiscreteGamma(n,eg,gx,frac,dat);
        if(ngam == 0) frac = 0.0;
        //if(n==0 && npgam>0) kmax[c] = ddxPrimaryGamma(kmax[c],epar[c]);
      }

      ddxLegNormalize(kmax[c],nleg[c],frac,epar[c],cleg[c]);
    }

    ddxCheckEnergyBalance(cm,np[n],ngam,eg,gx,&dat[n]);

    eclOutNucleusHead(n,cm,np);

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
        if(c == 0) eclOutGammaLine(ngam,nleg[c],eg,gx);
        eclOutLegCoeff(kmax[c],nleg[c],epar[c],cleg[c]);
      }        
      else{
        eclOutLegCoeff(-1,0,epar[c],cleg[c]);
      }
    }
  }

  delete [] eg;
  delete [] gx;
}


/**********************************************************/
/*      Check Energy Balance and Fix Gamma Multiplicity   */
/**********************************************************/
void ddxCheckEnergyBalance(const int cm, double *np, int ngam, double *eg, double *gx, EXSpectra *n)
{
  double ptot = 0.0, gtot1 = 0.0, gtot2 = 0.0;
  double *sum = new double [cm];
  double *tot = new double [cm];

  double emax = n->getEm();

  /*** integrate continuum spectra */
  for(int c=0 ; c<cm ; c++){
    sum[c] = tot[c] = 0.0;
    for(int k=0 ; k<kmax[c]-1 ; k++){
      double e0 = epar[c][k  ];
      double e1 = epar[c][k+1];
      double s0 = cleg[c][0][k  ];
      double s1 = cleg[c][0][k+1];
      sum[c] += (s0             + s1            )*(e1-e0)/2.0;
      tot[c] += (s0*(e1+2.0*e0) + s1*(e0+2.0*e1))*(e1-e0)/6.0;
    }
  }

  /*** add discrete gamma energies */
  for(int k=0 ; k<ngam ; k++){
    sum[0] += gx[k];
    tot[0] += gx[k]*eg[k];
  }
  gtot1 = tot[0]; // total gamma-ray energy from spectrum

  for(int c=1 ; c<cm ; c++) ptot += tot[c] * np[c];
  gtot2 = emax - ptot; // all available gamma-ray energy

  np[0] = (gtot1 > 0.0) ? gtot2/gtot1 : 0.0;
  /*** if no spec given, set multiplicity zero */
  if(np[0] < 0.0) np[0] = 0.0;

/*
  cout <<" Em " <<setw(14)<<emax<< " " <<setw(14)<<ptot<<setw(14)<<gtot1<< " " <<setw(14)<<gtot2<< endl;
  cout <<" Mg " <<setw(14)<<np[0]<< endl;
  for(int c=0 ; c<cm ; c++)
    cout << " Intg" << setw(3) << c << setw(14) << sum[c] << setw(14) << tot[c] << endl;
*/

  delete [] sum;
  delete [] tot;
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

  if( n->getGcon() && (c != 0) ) epmax -= n->getEl();
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
/*      Discrete Gamma Fraction                           */
/**********************************************************/
double ddxDiscreteFraction(const int km, EXSpectra *d)
{
  double sc = 0.0, sd = 0.0, st = 0.0;

  for(int k=0 ; k<km ; k++){
    sc += d->spec[0][k];
    sd += d->glin[k];
  }
  st = sc+sd;

  if(st == 0.0) return(0.0);
  else          return(sd/st);
}


/**********************************************************/
/*      Discrete Gamma-Rays                               */
/**********************************************************/
int ddxDiscreteGamma(int n, double *e, double *x, const double f, EXSpectra *dat)
{
  double eg, dp, pop[NLEV];
  int    i0, i1, k = 0, nlev = dat[n].getNl();

  if(f == 0.0) return(0);

  /*** if binary reaction is separated, discrete gammas exclude
       fraction populated by particle transition */
  for(i0=0 ; i0<nlev ; i0++){
    pop[i0] = (dat[n].getGcon()) ? dat[n].lpop[i0].gamma :
                                   dat[n].lpop[i0].gamma + dat[n].lpop[i0].particle;
  }

  /*** calculate gamma-cascade, populated by continuum gamma only */
  for(i0=nlev-1 ; i0>0 ; i0--){
    for(int j=0 ; j<ncl[n].lev[i0].ngamma ; j++){
      i1 = ncl[n].lev[i0].fstate[j];
      eg = ncl[n].lev[i0].energy - ncl[n].lev[i1].energy;
      dp = pop[i0] * ncl[n].lev[i0].branch[j];
      pop[i1] += dp;

      if(dp == 0.0) continue;

      if(k < NGLN){
        e[k] = eg;
        x[k] = dp;
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
      if(e[i] > e[k]) k = i;
    }
    eg = e[j];  e[j] = e[k];  e[k] = eg;
    dp = x[j];  x[j] = x[k];  x[k] = dp;
  }

  /*** normalize discrete lines */
  double s = 0.0;
  for(int j=0 ; j<ng ; j++) s += x[j];
  if(s == 0.0) return(0);

  s = f/s;
  for(int j=0 ; j<ng ; j++) x[j] *= s;

  return(ng);
}


/**********************************************************/
/*      Renormalize Legendre Coefficients                 */
/**********************************************************/
void ddxLegNormalize(const int km, const int nleg, const double f, double *ep, double **cl)
{
  if(km<=1) return;

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
/*      Kalbach Systematics for DDX                       */
/**********************************************************/
void ddxKalbach(const int c, const double ep, const double f, const double p, double *ddx, double *cl)
{
  double ang[NDDX];

  double dt = 180.0/(NDDX-1);

  /*** Kalbach systematics a-parameter */
  double a = ddxKalbachAfactor(incid,c,ep,sa,sb[c]);
  double x = a/(2.0*sinh(a));

  /*** calculate angular distribution */
  for(int i=0 ; i<NDDX ; i++){
    ang[i] = (double)i*dt;
    double q = ang[i]/180.0 * PI;
    ddx[i] = x * (cosh(a*cos(q)) + f*sinh(a*cos(q)));
  }

  /*** calculate Legendre expansion coefficients by LSQ method */
  LSQLegendre(false,NDDX,NLEG,ang,ddx,cl);

  /*** convert into ENDF-6 data, F(E') = 2 f/(2L+1) x (ds/dE) */
  for(int l=0 ; l<NLEG ; l++) cl[l] = cl[l] / (1.0*l+0.5) * p;
}


/**********************************************************/
/*      Separation Energies Defined in Kalbach Syst.      */
/**********************************************************/
void ddxKalbachSetParm(void)
{
  double   ib = 0.0;
  ZAnumber ejec(0,0),resd(0,0);

  sb[0] = 0.0;

  for(int c=1 ; c<7 ; c++){

    switch(c){
    case  1: ejec.setZA(0,1); ib =  0.0  ;  break;
    case  2: ejec.setZA(1,1); ib =  0.0  ;  break;
    case  3: ejec.setZA(2,4); ib = 28.296;  break;
    case  4: ejec.setZA(1,2); ib =  2.225;  break;
    case  5: ejec.setZA(1,3); ib =  7.718;  break;
    case  6: ejec.setZA(2,3); ib =  8.482;  break;
    default: ejec.setZA(0,0); ib =  0.0  ;  break;
    }

    resd = comp - ejec;
    if( (proj.getZ() == ejec.getZ()) && (proj.getA() == ejec.getA()) ) incid = c;

    double xc = (double)(comp.getN()-comp.getZ())*(double)(comp.getN()-comp.getZ());
    double xb = (double)(resd.getN()-resd.getZ())*(double)(resd.getN()-resd.getZ());
    double yc = pow((double)comp.getA(),1.0/3.0);
    double yb = pow((double)resd.getA(),1.0/3.0);
    double zc = (double)comp.getZ()*(double)comp.getZ();
    double zb = (double)resd.getZ()*(double)resd.getZ();

    sb[c] = 15.68 * (comp.getA()    - resd.getA()   )
          - 28.07 * (xc/comp.getA() - xb/resd.getA())
          - 18.56 * (yc*yc - yb*yb)
          + 33.22 * (xc/pow(yc,4.0) - xb/pow(yb,4.0))
          - 0.717 * (zc/yc - zb/yb)
          + 1.211 * (zc/comp.getA() - zb/resd.getA())
          - ib;
  }

  sa = sb[incid];
}


/**********************************************************/
/*      a(E) Factor in Kalbach Systematics                */
/**********************************************************/
double ddxKalbachAfactor(const int a, const int b, const double eout, const double sa, const double sb)
{
  const double et1 = 130.0;
  const double et3 =  41.0;
  const double c1 = 0.04, c2 = 1.8e-06, c3 = 6.7e-07;

  double ea = ecms + sa;
  double eb = eout + sb;

  double e1 = fmin(ea,et1);
  double e3 = fmin(ea,et3);

  double x1 = e1*eb/ea;
  double x3 = e3*eb/ea;

  double ma = (a == 3) ? 0.0 : 1.0;
  double mb = (b == 3) ? 2.0 : ((b == 1) ? 0.5 : 1.0);

  return( c1*x1 + c2*pow(x1,3.0) + c3*pow(x3,4.0)*ma*mb );
}

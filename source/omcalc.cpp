/******************************************************************************/
/*  omcalc.cpp                                                                */
/*        spherical optical model calculation                                 */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "structur.h"
#include "optical.h"
#include "omcalc.h"
#include "global.h"
#include "output.h"
#include "omoutput.h"
#include "etc.h"
#include "coupling.h"
#include "terminate.h"

#undef PartialCrossSection
#undef NormalizedWavefunction
#undef omLegendre


#ifdef PartialCrossSection
static const int LPartial = 4;
static double parSigE[LPartial], parSigR[LPartial], parSigT[LPartial];
#endif


#ifdef omLegendre
static void omLegendreCoefficient (const int, double, CrossSection *, std::complex<double> *);
#endif

#ifdef NormalizedWavefunction
void omWaveFunction(const int, const double, const double, const double, std::complex<double>, LevelData *, Potential *, Wavefunc *);
#endif

/**********************************************************/
/*     Entrance channel calculation                       */
/*         usual optical model calculation for incident   */
/*         particle and determine max-J                   */
/**********************************************************/
int omCalc
(const double  energy,        // CMS incident energy
 Pdata        *proj,          // particle data for projectile
 ZAnumber     *targ,          // target ZA
 const double  mu,            // reduced mass for Targ+Proj system
 double       *tran,          // transmission coefficients (output)
 CrossSection *crx)           // calculated cross sections (output)
{
  Optical      omp;
  Potential    pot;
  Wavefunc     wfn;
  CCdata       cdt;
  LevelData    lev;
  std::complex<double> *smat = nullptr;

  int          l,zzprod=0;
  double       cx[3], ce;

  if(energy <= 0.0) return(-1);
  if((proj->za.getZ() > 0) && (energy < ECUT_CHARGED)) return(-1);

//---------------------------------------
//      Memory Allocation

  try{
    pot.memalloc(MAX_POINTS);              // optical potential
    wfn.memalloc(MAX_L,MAX_POINTS);        // wave functions (internal / external)
    smat = new std::complex<double> [MAX_L*3];  // S-matrix elements
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("omCalc");
  }

  for(int j=0 ; j<3*MAX_L ; j++) tran[j] = 0.0;

  zzprod = targ->getZ() * proj->za.getZ();

 
//---------------------------------------
//      Energy and Coulomb Parameter

  cdt.lev = &lev;
  omSetEnergy(energy, zzprod, mu, &lev);

//---------------------------------------
//      Setup Optical Potential Geometry

  omSetOmp(proj->omp,lev.energy,targ->getZ(),targ->getA(),proj->za.getZ(),proj->za.getA(),&omp);
  if(prn.system && ctl.entrance) outOMP(0,&omp);

  omPotentialForm(zzprod, &omp, &lev, &pot);

//---------------------------------------
//      Main Calculation

  /***  free space wave function */
  cdt.lmax = omExternalFunction(0,pot.rho_match,lev.coulomb.eta,&wfn);

#ifdef  ParticalCrossSection
 for(int l=0 ; l<LPartial ; l++) parSigE[l] = parSigR[l] = parSigT[l] = 0.0;
#endif


  /***  internal wave function */
  std::complex<double> a, b;
  double  cr, ci, dr, di, pspin = proj->spin2/2.0;
  int     jmax   = proj->spin2+1;

  crx->elastic = crx->reaction = crx->total = 0.0;
  for(l=0 ; l<=cdt.lmax ; l++){
    a = wfn.extderiv[l] / wfn.external[l] * pot.rho_match;
    b = conj(wfn.external[l]) / wfn.external[l];

    for(int j=0 ; j<jmax ; j++){
      double xj = l + pspin - (double)j;
      if(xj < fabs(l-pspin)) continue;
      omInternalFunction(pot.n_match,pot.width,lev.wavesq,(double)l,pspin,xj,&pot,&wfn);

      int index = l*3+j;
      smat[index] = omSmatrix(pot.n_match,pot.width,pot.rad_match,a,b,&wfn);
      tran[index] = 1.0-norm(smat[index]); if(tran[index] < 0.0) tran[index] = 0.0;

#ifdef NormalizedWavefunction
      if(ctl.entrance) omWaveFunction(l,pspin,xj,mu,smat[index],&lev,&pot,&wfn);
#endif
    }

    /***  reaction cross section */
    cx[0] = crx->reaction;
    cx[1] = 0.0;
    switch(proj->pid){
    case neutron :
      cr = 1.0 - smat[l*3  ].real();
      ci =       smat[l*3  ].imag();
      dr = 1.0 - smat[l*3+1].real();
      di =       smat[l*3+1].imag();
      ce = (l+1.0)*(cr*cr + ci*ci) + l*(dr*dr + di*di);
      crx->elastic += ce;
#ifdef PartialCrossSection
      if(l < LPartial) parSigE[l] = ce;
#endif
      // fall through
    case proton :
    case triton :
    case helion :
      crx->reaction += (cx[1]=    (l+1)*tran[l*3  ]
                                 + l   *tran[l*3+1]);
      break;
    case alpha :
      crx->reaction += (cx[1]=  (2*l+1)*tran[l*3  ]);
      break;
    case deuteron :
      crx->reaction += (cx[1]=( (2*l+3)*tran[l*3  ]
                               +(2*l+1)*tran[l*3+1]
                               +(2*l-1)*tran[l*3+2])/3.0);
      break;
    default:
      break;
    }

    if(proj->pid == userdef){
      crx->reaction += (cx[1]=  (2*l+1)*tran[l*3  ]);
    }


#ifdef PartialCrossSection
      if(l < LPartial) parSigR[l] = cx[1];
#endif

    if((cx[0] != 0.0) && ((cx[1]/cx[0]) < CRIT_LCUT)) break;
  }
  if(l == MAX_L-1){
//    std::cerr << "reaction cross section not converge" << std::endl;
    cdt.lmax=l;
  }
  else if(l<cdt.lmax) cdt.lmax=l;

//---------------------------------------
//      Scattering Amplitude

  /*** only if print opt is given, and for entrance channel */
  if(prn.angdist && ctl.entrance){
    double w = sqrt(lev.wavesq);

    for(int i=0 ; i<MAX_ANGDIST ; i++){
      double *p = crx->angdist[0]; // point elastic data entry

      switch(proj->pid){
      case neutron  :
        p[i] = omScatteringAmplitudeN(cdt.lmax,cx,crx->costh[i],w,smat);
        break;
      case proton   :
      case triton   :
      case helion   :
        p[i] = omScatteringAmplitudeP(cdt.lmax,cx,crx->costh[i],w,lev.coulomb,smat);
        break;
      case alpha    :
        p[i] = omScatteringAmplitudeA(cdt.lmax,cx,crx->costh[i],w,lev.coulomb,smat);
        break;
      case deuteron :
        p[i] = omScatteringAmplitudeD(cdt.lmax,cx,crx->costh[i],w,lev.coulomb,smat);
        break;
      default       :
        break;
      }
      p[i] *= NORM_FACT;
    }
    /*** print out every ANGLE_STEP angle points */
    outAngularDistribution(0,0,0,ANGLE_STEP,targ);
  }

//---------------------------------------
//      Strength Function if E < 100keV

  if(prn.system && energy < 0.1 && ctl.entrance && proj->pid == neutron){
    omStrengthFunction(cdt.lmax,targ->getA(),lev.wavesq,energy,&wfn,smat);
  }

//---------------------------------------
//      Output Cross Sections

  double kfact = NORM_FACT*PI/lev.wavesq;
  crx->elastic  *= kfact;
  crx->reaction *= kfact;
  crx->total  = crx->elastic + crx->reaction;

  if(ctl.entrance){
    if(prn.xsection){
      if(proj->pid == neutron) outCrossSection(0,crx->total,crx->elastic,crx->reaction);
      else                     outCrossSection(1,crx->total,crx->elastic,crx->reaction);
    }
    if(prn.transmission) outTransmission(cdt.lmax,proj->spin2,NORM_FACT*PI/lev.wavesq,tran);
    if(prn.smatrix) outSmatrix(cdt.lmax,proj->spin2,smat);


#ifdef PartialCrossSection
    std::cout << std::setw(12) << energy;
    for(int l=0 ; l<LPartial ; l++){
      parSigE[l] *= kfact;
      parSigR[l] *= kfact;
      parSigT[l]  = parSigE[l] + parSigR[l];
      std::cout << std::setw(12) << parSigT[l];
      std::cout << std::setw(12) << parSigE[l];
      std::cout << std::setw(12) << parSigR[l];
    }
    std::cout << std::endl;
#endif

#ifdef omLegendre
    if(ctl.entrance && prn.angdist && proj->pid == neutron){
      omLegendreCoefficient(cdt.lmax, kfact, crx, smat);
    }
#endif
  }


//---------------------------------------
//      Free Allocated Memory

  delete [] smat;

  return(cdt.lmax);
}


/**********************************************************/
/*      Legendre Coefficient from S Matrix Elements       */
/**********************************************************/
#ifdef omLegendre
void omLegendreCoefficient(const int lmax, double x, CrossSection *crx, std::complex<double> *smat)
{
  if(lmax == 0) return;

  double *p;
  p = new double [MAX_L*2];

  for(int l=0 ; l<2*MAX_L ; l++) p[l] = 0.0;

  x = x/(8*PI);
  int smax = 2;
  for(int l=0 ; l<2*lmax ; l++){

    for(int l0=0 ; l0<lmax ; l0++){
      for(int s0=0 ; s0<=smax ; s0+=2){
        int j0 = l0*2 + 1 - s0; if(j0 < 0) continue;
        int idx0 = 3*l0+s0/2;

        for(int l1=0 ; l1<lmax ; l1++){
          for(int s1=0 ; s1<=smax ; s1+=2){
            int j1 = l1*2 + 1 - s1; if(j1 < 0) continue;
            int idx1 = 3*l1+s1/2;

            /*** Z-coefficient */
            double z = z_coefficient(2*l0,j0,2*l1,j1,1,2*l);
            if(z == 0.0) continue;

            /*** Re{(1-S0)(1-S1)^*} */
            double s =  (1.0-smat[idx0].real())*(1.0-smat[idx1].real())
                       +     smat[idx0].imag() *     smat[idx1].imag();
            p[l] += x*z*z*s;
          }
        }
      }
    }
  }

  for(int l=0 ; l<2*lmax ; l++)  std::cout << l<<" " << p[l] * 4*PI / (2*l+1) << std::endl;

  for(int i=0 ; i<MAX_ANGDIST ; i++){

    double y = 0.0;
    for(int l=0 ; l<2*lmax ; l++){
      y += p[l] * legendre(l,crx->costh[i]);
    }
    std::cout << std::setw(12) << crx->theta[i] << std::setw(12) << y << std::endl;
  }

  delete [] p;
}
#endif


#ifdef NormalizedWavefunction
void omWaveFunction(const int l, const double xs, const double xj, const double mu, std::complex<double>smat, LevelData *lev, Potential *pot, Wavefunc *wfn)
{
  if(l != 3) return;

  std::complex<double>c = omNormalizationFactor(l,pot->n_match-3,lev->coulomb,lev->coulomb_scat0,smat,wfn);
  omNormalization(pot->n_match,0.0,c,pot,wfn);

  double k = sqrt(lev->wavesq);
  double x = 2.0*mu*AMUNIT/VLIGHTSQ/HBARSQ;
  double xl_term = l*(l+1.0);
  double so_term = xj*(xj+1.0)-xl_term-xs*(xs+1.0);

  std::cout << "#  L = " << l << "  J = " << std::setprecision(1) << std::setw(5) << xj;
  std::cout << "           S = " << std::setprecision(4) << std::setw(12) << smat.real() << std::setw(12) << smat.imag();
  std::cout << "        T = " << std::setprecision(4) << std::setw(12) << 1.0 - norm(smat) << std::endl;
  std::cout << "# R[fm]       Real        Imag        Norm      ";
  std::cout << " V            W           V+L+SO      W+SO" << std::endl;

  for(int ir=0 ; ir<pot->n_match ; ir++){
    std::complex<double> q = xl_term*pot->r2inv[ir] - pot->mean_field[ir] - so_term*pot->spin_orbit[ir];
    double r = pot->width*ir;

    std::cout << std::setprecision(4);
    std::cout << std::setw(12) << r;
    std::cout << std::setw(12) << wfn->internal[ir].real() << std::setw(12) <<wfn->internal[ir].imag();
    std::cout << std::setw(12) << norm(wfn->internal[ir]);

    std::cout << std::setw(12) << -pot->mean_field[ir].real() / x;
    std::cout << std::setw(12) << -pot->mean_field[ir].imag() / x;
    std::cout << std::setw(12) << q.real() / x;
    std::cout << std::setw(12) << q.imag() / x;

    double d = arg(smat);
    std::cout << std::setw(12) << 0.5*(sin(k*r) + std::abs(smat)*sin(k*r + d));
    std::cout << std::setw(12) << 0.5*(cos(k*r) - std::abs(smat)*cos(k*r + d));

    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;
/*
  double s = 0.0;
  for(int ir=0 ; ir<pot->n_match ; ir++){
    s += std::abs(wfn->internal[ir]) * pot->mean_field[ir].imag() * pot->width;
  }
  std::cout << std::setprecision(4);
  std::cout << std::setw(12) << lev->energy << std::setw(12) << 1.0 - norm(smat) << std::setw(12) << s << std::endl;
*/
}
#endif

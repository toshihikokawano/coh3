/******************************************************************************/
/*  omcalc.cpp                                                                */
/*        spherical optical model calculation                                 */
/******************************************************************************/

#include <iostream>
#include <cmath>
#include <cstdio>

#include "structur.h"
#include "optical.h"
#include "omcalc.h"
#include "terminate.h"


/**********************************************************/
/*     Entrance channel calculation                       */
/*         usual optical model calculation for incident   */
/*         particle and determine max-J                   */
/**********************************************************/
int omCalc
(double        energy,        // CMS incident energy
 Pdata        *proj,          // particle data for projectile
 ZAnumber     *targ,          // target ZA
 double        mu,            // reduced mass for Targ+Proj system
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
  double       cx[3];

  if(energy <= 0.0) return(-1);
  if((proj->za.getZ() > 0) && (energy < ECUT_CHARGED)) return(-1);

//---------------------------------------
//      Memory Allocation

  try{
    pot.memalloc(MAX_POINTS);
    wfn.memalloc(MAX_L,MAX_POINTS);
    smat = new std::complex<double> [MAX_L*3];
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("omCalc");
  }

  for(int j=0 ; j<3*MAX_L ; j++) tran[j] = 0.0;

  zzprod = targ->getZ()*proj->za.getZ();


//---------------------------------------
//      Energy and Coulomb Parameter

  cdt.lev = &lev;
  omSetEnergy(energy, zzprod, mu, &lev);

//---------------------------------------
//      Setup Optical Potential Geometry

  omSetOmp(proj->omp,lev.energy,targ->getZ(),targ->getA(),proj->za.getZ(),proj->za.getA(),&omp);
  pot.width = INTEG_WIDTH;
  omPotentialForm(zzprod, &omp, &lev, &pot);

//---------------------------------------
//      Main Calculation

  /***  free space wave function */
  cdt.lmax = omExternalFunction(0,pot.rho_match,lev.coulomb.eta,&wfn);

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
      if(xj<fabs(l-pspin)) continue;
      omInternalFunction(pot.n_match,pot.width,lev.wavesq,(double)l,pspin,xj,&pot,&wfn);

      int index = l*3+j;
      smat[index] = omSmatrix(pot.n_match,pot.width,pot.rad_match,a,b,&wfn);
      tran[index] = 1.0-norm(smat[index]); if(tran[index] < 0.0) tran[index] = 0.0;
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
      crx->elastic += (l+1.0)*(cr*cr + ci*ci) + l*(dr*dr + di*di);
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

    if((cx[0] != 0.0) && ((cx[1]/cx[0]) < CRIT_LCUT)) break;
  }
  if(l==MAX_L-1){
//    cerr << "reaction cross section not converge" << endl;
    cdt.lmax=l;
  }
  else if(l<cdt.lmax) cdt.lmax=l;

  double kfact = NORM_FACT*PI/cdt.lev->wavesq;
  crx->elastic  *= kfact;
  crx->reaction *= kfact;
  crx->total  = crx->elastic + crx->reaction;


//---------------------------------------
//      Free Allocated Memory

  delete [] smat;

  return(cdt.lmax);
}


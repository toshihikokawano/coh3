/******************************************************************************/
/*  FRDMmain.cpp                                                              */
/*        Finite Range Droplet Model calculation                              */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "constant.h"
#include "FRDM.h"
#include "global.h"
#include "terminate.h"

static const int beta_lmax = 10;

/**********************************************************/
/*      FRDM Calculation                                  */
/**********************************************************/
void FRDMmain(MFTSystem *sys, Basis *basis, double *beta, HFInterface *spspec)
{
//---------------------------------------
//      Initial Setting

  double **buf, **alm;

  /*** shperical harmonics expansion coefficients */
  buf = new double * [beta_lmax+1];
  alm = new double * [beta_lmax+1];
  for(int l=0 ; l<=beta_lmax ; l++){
    buf[l] = new double [2*l+1];
    alm[l] = &buf[l][l];
  }


//---------------------------------------
//      Create Nucleus Shape

  /*** 3D shape normalization factor */
  PolarGLPoint gp;
  FRDMShapeNormalization(&gp,sys);

  /*** Gauss-Legendre points */
  FRDMShapeCreateSurface(&gp,sys);
  FRDMShapeStoreCurvature(&gp,sys);

  /*** betas by spherical harmonics expansion */
  if(pmf.FRDMbeta){
    BetaExpansion(beta_lmax,alm,&gp,sys);
    FRDMPrintBeta(beta_lmax,alm,sys->R0);
  }

  /*** copy deformation parameter into beta array */
  for(int l=0 ; l<MAX_LAMBDA ; l++) beta[l] = alm[(l+1)*2][0]/sys->R0;

  /*** create surface at the equi-distant 3D mesh */
  if(pmf.FRDMshape){
    PolarCoordinate sp;
    for(int i = 0 ; i < sp.nt ; i ++){
      for(int j = 0 ; j < sp.np ; j ++){
        sp.r[i][j] = gp.getNormalization() * FRDMShape(sp.t[i], sp.p[j], sys);
      }
    }
    FRDMWrite3DMesh(fileFRDMShape, &sp);
  }


//---------------------------------------
//      Main FRDM Calculation

  /*** energy compoents */
  FRDMEnergy e;

  /*** macro energies */
  e.macro = FRDMMacroEnergy(&gp,sys);
  e.spherical = FRDMMacroEnergySpherical(sys);

  /*** micro energy and total */
  e.micro = FRDMMicroEnergy(sys,basis,spspec);


//---------------------------------------
//      Zero-Point Energy

  if(pmf.FRDMzero) e.zero = FRDMZeropointEnergy(&gp,sys,basis,spspec,&e);

  e.setTotal();


//---------------------------------------
//      Print Macro/Micro Energies
  if(pmf.energy) FRDMPrintTotalEnergy(&e);
  if(pmf.FRDMfileout) FRDMWriteEnergy(fileFRDMEnergy,sys,&e);


  for(int l=0; l<=beta_lmax ; l++)  delete [] buf[l];
  delete [] buf;
  delete [] alm;
}

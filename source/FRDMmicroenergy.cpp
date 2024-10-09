/******************************************************************************/
/*  FRDMmicroenergy.cpp                                                       */
/*        Micro Energy Part in FRDM                                           */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "global.h"
#include "FRDM.h"
#include "eispack.h"
#include "terminate.h"

const int HFGridSize_Z = 50;
const int HFGridSize_R = 25;

static HermitePolynomial  fh;
static LaguerrePolynomial fl;
static Operator  *matH;


static void FRDMEigenvalue (Basis *, OneBodyPotential *, SPEnergy *);
static void FRDMAllocateMemory (Basis *, SPEnergy *);
static void FRDMDeleteAllocated (void);

/**********************************************************/
/*   Microscopic Energy Emac(Z,N,shape)                   */
/**********************************************************/
double FRDMMicroEnergy(MFTSystem *sys, Basis *basis, HFInterface *spspec)
{
  OneBodyPotential pot[2];
  SPEnergy         spe[2];

  /*** Lipkin-Nogami parameters */
  LNParameter bcs[2];

  /*** Strutinsky shell correction parameters */
  SCParameter scp[2];

  /*** re-adjust basis pamaeters, whatever given */
  double homega0 = gCurvature / sys->a13;

  basis->setN0(12);
  basis->q = 1.0;
  basis->b = HBAR*VLIGHT / sqrt(MNEUTRON*AMUNIT * homega0);
  basis->b = 1.0/basis->b; // sqrt{m omega0 / hbar}

  HOCylindricalBasis(basis);

  /*** Hermite and Laguerre Polynomials on the basis */
  HOPolynomials(HFGridSize_R, HFGridSize_Z, basis, &fl, &fh);

  /*** allocate potential arrays, and prepare FRDM potentials */
  if(pot[0].memalloc(fl.getNsize(),fh.getNsize())){
    message << "FRDM allocation of neutron potential";
    cohTerminateCode("FRDMMicroEnergy");
  }
  if(pot[1].memalloc(fl.getNsize(),fh.getNsize())){
    message << "FRDM allocation of proton potential";
    cohTerminateCode("FRDMMicroEnergy");
  }

  FRDMPotential(sys, basis, pot, &fl, &fh);

  /*** allocate arrays */
  FRDMAllocateMemory(basis,spe);

  /*** for neutron and proton shells */
  double emic = 0.0;
  for(int p=0 ; p<=1 ; p++){
    int np = (p == 0) ? sys->getN() : sys->getZ();

    /*** calculate single particle levels */
    FRDMEigenvalue(basis,&pot[p],&spe[p]);

    /*** shell correction by Strutinsky */
    scp[p].energy_average = FRDMShellEnergy(sys,np,&spe[p],&scp[p]);
    scp[p].energy_micro   = FRDMSumSingleParticleEnergy(np,&spe[p]);
    emic += scp[p].energy_micro - scp[p].energy_average;

    /*** solve BCS gap equation */
    bcs[p].np = np;
    FRDMSoleveBCS(sys,scp[p].rho,&bcs[p],&spe[p]);
    bcs[p].energy_micro   = FRDMPairingEnergy(&bcs[p],&spe[p]);
    bcs[p].energy_average = FRDMAveragePairingEnergy(scp[p].rho,&bcs[p]);

    /*** total microscopic shell+pairing correction */
    emic += bcs[p].energy_micro - bcs[p].energy_average;

    /*** expand single particle wavefunctions by spherical harmonics
         this also determines parities of the states */
    if(ctl.expansion || pmf.spstate) HOSphericalExpansion(basis,&spe[p],matH,&spspec[p]);
  }


//---------------------------------------
//      Output Results

  /*** BCS and Shell Correction parameters */
  if(pmf.bcs){
    FRDMPrintBCS(bcs);
    FRDMPrintSC(sys,scp);
  }

  /*** microscopic energy components */
  if(pmf.energy) FRDMPrintMicroEnergy(bcs,scp);

  /*** print 3D potential data on the Gauss mesh */
  if(pmf.potential) FRDMPrintPotential(basis,pot,&fl,&fh);

  /*** single particle spectra */
  if(pmf.spstate) FRDMPrintSPState(false,basis,spe,bcs);
  if(pmf.shiftedstate) FRDMPrintSPState(true,basis,spe,bcs);

  if(pmf.FRDMfileout) FRDMWriteParameters(sys->getZ(),sys->getA(),bcs,scp);


  /*** free allocated memory */
  FRDMDeleteAllocated();

  return(emic);
}


/**********************************************************/
/*      Neutron/Proton Separable Part                     */
/**********************************************************/
void FRDMEigenvalue(Basis *basis, OneBodyPotential *pot, SPEnergy *spe)
{
  /*** loop over each omega block */
  for(int k=0 ; k<basis->nb ; k++){

    /*** Harmonic Oscillator Hamiltonian */
    HOHamiltonian(k,basis,pot,&matH[k],&fl,&fh);

    /*** diagonalize Hamiltonian */
    EISPACKUnitaryDiag(basis->ndata[k],matH[k].matrix,matH[k].value,matH[k].vector);

    /*** copy eigenvector and eigenvalues */
    for(int i=0 ; i<basis->ndata[k] ; i++){
      spe->energy[basis->index[k]+i] = matH[k].value[i];
    }
  }

  /*** sort eigenvalues in the assending order, and determine pointers */
  HOSortEigenvalue(basis,spe->energy,spe->index,spe->block,spe->subindex);
}


/**********************************************************/
/*      Allocate Large Arrays                             */
/**********************************************************/
void FRDMAllocateMemory(Basis *basis, SPEnergy *spe)
{
  try{
    /*** allocate operator H, Nc blocks, and each contains NxN matrix */
    matH = new Operator [basis->nb];
  }
  catch(std::bad_alloc &e){
    message << "MFT memory allocation error for Vector data " << e.what();
    cohTerminateCode("FRDMAllocateMemory");
  }

  int icon = 0;
  for(int k=0 ; k<basis->nb ; k++){
    icon += matH[k].memalloc(basis->ndata[k]);
  }
  icon += spe[0].memalloc(basis->nv);
  icon += spe[1].memalloc(basis->nv);
  if(icon > 0){
    message << "memory allocation error for Operator";
    cohTerminateCode("FRDMAllocateMemory");
  }
}


/**********************************************************/
/*      Delete Allocated Arrays                           */
/**********************************************************/
void FRDMDeleteAllocated(void)
{
  fl.memfree();
  fh.memfree();
  delete [] matH;
}


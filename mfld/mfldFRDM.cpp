/******************************************************************************/
/*  mfldFRDM.cpp                                                              */
/*        Micro Energy Part in FRDM                                           */
/******************************************************************************/

#include <iostream>
using namespace std;

#include "global.h"
#include "FRDM.h"
#include "mfld.h"
#include "coupling.h"
#include "eispack.h"
#include "terminate.h"


const int HFGridSize_Z  = 50;              // cylindrical coordinate grid Z-size
const int HFGridSize_R  = 25;              // cylindrical coordinate grid R-size

static HermitePolynomial  fh;              // z-direction solution of schroedinger equation
static LaguerrePolynomial fl;              // r-direction
static SPEnergy           spe[2];          // energies of single particle states
static Operator          *matH[2];         // cylindrical harmonic oscillator hamiltonian
static double            *v2[2], *matV[2]; // v^2 of BCS solution, and density matrix
static double            *membuf;          // memory buffer

static void FRDMSurfaceEnergyParameter(MFTSystem *);
static void FRDMSingleParticle (MFTSystem *, Basis *, SPEnergy *, HFInterface *);
static void FRDMEigenvalue (Basis *, OneBodyPotential *, SPEnergy *, Operator *);
static void FRDMAllocateMemory (Basis *);
static void FRDMDeleteAllocated (Basis *);


/**********************************************************/
/*     Mean-Field Calculation                             */
/**********************************************************/
void FRDMCalc(MFTSystem *sys, HFInterface *spspec)
{
  /*** FRDM re-adjust basis parameters, whatever given */
  double homega0 = gCurvature / sys->a13;

  /*** create harmonic oscillator basis */
  Basis basis;

  basis.setN0(12);
  basis.q = 1.0;
  basis.b = HBAR*VLIGHT / sqrt(MNEUTRON*AMUNIT * homega0);
  basis.b = 1.0/basis.b; // sqrt{m omega0 / hbar}

  /*** Harmonic Oscillator basis (should be before memory allocation) */
  HOCylindricalBasis(&basis);

  /*** allocate memory */
  FRDMAllocateMemory(&basis);

  /*** Hermite and Laguerre Polynomials on the basis */
  HOPolynomials(HFGridSize_R, HFGridSize_Z, &basis, &fl, &fh);

  /*** calculate single particle spectra */
  FRDMSurfaceEnergyParameter(sys);
  FRDMSingleParticle(sys,&basis,spe,spspec);

  /*** free allocated memory */
  FRDMDeleteAllocated(&basis);

  message << "FRDM produced " << spspec[0].n << " states in N and " << spspec[1].n << " states in Z";
  cohNotice("FRDMCalc");
}


/**********************************************************/
/*   Approximated FRDM Bs Parameter from Liquid Drop      */
/*   --------                                             */
/*   assuming deformation is only eps2, calculated Bs are */
/*   fitted by the polynomial                             */
/**********************************************************/
void FRDMSurfaceEnergyParameter(MFTSystem *sys)
{
  const double a = 0.184681;  // +/- 0.0004879    (0.2642%)
  const double b = 0.0341609; // +/- 0.001104     (3.233%)

  double eps = sys->getEps(2);

  sys->bs = 1.0 + a * eps*eps + b * eps*eps*eps;
}


/**********************************************************/
/*   FRDM Single Particle State Calculation               */
/**********************************************************/
void FRDMSingleParticle(MFTSystem *sys, Basis *basis, SPEnergy *spe, HFInterface *spspec)
{
  OneBodyPotential pot[2];  // FRDM potential
  LNParameter      bcs[2];  // BCS parameters
  SCParameter      scp[2];  // shell correction parameters

  for(int p=0 ; p<=1 ; p++){
    pot[p].memalloc(fl.getNsize(),fh.getNsize());
  }

  FRDMPotential(sys,basis,pot,&fl,&fh);

  for(int p=0 ; p<=1 ; p++){
    int np = (p == 0) ? sys->getN() : sys->getZ();

    /*** calculate single particle levels */
    FRDMEigenvalue(basis,&pot[p],&spe[p],matH[p]);

    /*** shell correction by Strutinsky */
    scp[p].energy_average = FRDMShellEnergy(sys,np,&spe[p],&scp[p]);
    scp[p].energy_micro   = FRDMSumSingleParticleEnergy(np,&spe[p]);

    /*** solve BCS gap equation */
    bcs[p].np = np;
    FRDMSoleveBCS(sys,scp[p].rho,&bcs[p],&spe[p]);
    bcs[p].energy_micro   = FRDMPairingEnergy(&bcs[p],&spe[p]);
    bcs[p].energy_average = FRDMAveragePairingEnergy(scp[p].rho,&bcs[p]);

    for(int k=0 ; k<basis->nv ; k++) v2[p][spe->index[k]] = spe->v2[k];
    HFDensityMatrix(basis,v2[p],matV[p],matH[p]);

    /*** estimate parity of the state and copy data */
    HOSphericalExpansion(basis,&spe[p],matH[p],&spspec[p]);

//      for(int k=0 ; k<basis->nv ; k++){
//        if( (bcs[p].n1 <= k) && (k<=bcs[p].n2) ){
//          double dx = (4.0*bcs[p].lambda2 - bcs[p].strength)*spspec[p].state[k].v2;
//          spspec[p].state[k].energy += dx;
//        }
//      }

    spspec[p].lambda = bcs[p].chemical_potential;
    spspec[p].delta  = bcs[p].pairing_gap;
  }

  /*** deformation parameters */
  if(pmf.deformation) FRDMPrintShape(sys);

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


  for(int p=0 ; p<=1 ; p++){
    pot[p].memfree();
  }
}


/**********************************************************/
/*      Neutron/Proton Separable Part                     */
/**********************************************************/
void FRDMEigenvalue(Basis *basis, OneBodyPotential *pot, SPEnergy *spe, Operator *matH)
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
/*      Allocate Memory                                   */
/**********************************************************/
void FRDMAllocateMemory(Basis *basis)
{
  const int ndim = basis->nv;
  try{
    /*** single particle spectra */
    spe[neutron].memalloc(ndim);
    spe[proton ].memalloc(ndim);

    /*** Harmonic Oscillator basis operator */
    matH[neutron] = new Operator [basis->nb];
    matH[proton ] = new Operator [basis->nb];
    for(int k=0 ; k<basis->nb ; k++){
      matH[neutron][k].memalloc(basis->ndata[k]);
      matH[proton ][k].memalloc(basis->ndata[k]);
    }

    /*** for other memory strage */
    membuf = new double [2*(ndim + ndim*(ndim+1)/2)];

    int pbuf = 0;
    v2[neutron]   = &membuf[pbuf]; pbuf += ndim;
    v2[proton ]   = &membuf[pbuf]; pbuf += ndim;
    matV[neutron] = &membuf[pbuf]; pbuf += ndim*(ndim+1)/2;
    matV[proton ] = &membuf[pbuf]; pbuf += ndim*(ndim+1)/2;

    /*** factorial */
    factorial_allocate();
  }
  catch(std::bad_alloc &e){
    std::cerr << "ERROR     :memory allocation error" << e.what() << std::endl;
    cohTerminateCode("FRDMAllocateMemory");
  }
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void FRDMDeleteAllocated(Basis *basis)
{
  spe[neutron].memfree();
  spe[proton ].memfree();

  fl.memfree();
  fh.memfree();

  basis->memfree();

  delete [] matH[neutron];
  delete [] matH[proton];

  delete [] membuf;

  factorial_delete();
}



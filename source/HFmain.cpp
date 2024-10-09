/******************************************************************************/
/*  HFmain.cpp                                                                */
/*        Hartree-Fock calculation                                            */
/*        the original Hartree-Fock BCS code in FORTRAN                       */
/*        written by L. Bonneau, rewritten in C++ by T. Kawano in 2014        */
/******************************************************************************/

#include <iostream>

#include "HF.h"
#include "global.h"
#include "eispack.h"
#include "terminate.h"


const int HFGridSize_Z = 50;
const int HFGridSize_R = 25;


static void HFSolve(MFTSystem *, BCSParameter *, HFEnergy *, Basis *, Skyrme *, Moment *);
static void HFEigenvalue(int, BCSParameter *, Basis *, OneBodyPotential *, SPEnergy *);
static inline void HFConstraintControl(int, Moment *);
static void HFAllocateMemory (const int, const int, Basis *);
static void HFDeleteAllocated ();

#define HFBCS_TOPLEVEL

static double *membuf, *v2, *matV;
static HermitePolynomial  fh;
static LaguerrePolynomial fl;
static GridData   vcoul, rho[2], tau[2], d2rho[2], divJ[2], *spd;
static Operator  *matH;
static OneBodyPotential  pot[2];
static SPEnergy   spe[2];


/**********************************************************/
/*      Hartree-Fock Calculation                          */
/**********************************************************/
void HartreeFock(MFTSystem *sys, Basis *basis, Skyrme *skyrme,
                 BCSParameter *bcs, Moment *moment, HFInterface *spspec)
{

//---------------------------------------
//      Initial Setting

  /*** cylindrical harmonic-oscillator basis */
  if(HOCylindricalBasis(basis)){
    message << "MFT basis solution error";
    cohTerminateCode("HartreeFock");
  }

  /*** Hermite and Laguerre Polynomials on the basis */
  HOPolynomials(HFGridSize_R, HFGridSize_Z, basis, &fl, &fh);

  /*** initial potential */
  if(pot[0].memalloc(fl.getNsize(),fh.getNsize())){
    message << "HF allocation of neutron potential";
    cohTerminateCode("HartreeFock");
  }
  if(pot[1].memalloc(fl.getNsize(),fh.getNsize())){
    message << "HF allocation of proton potential";
    cohTerminateCode("HartreeFock");
  }

  HFInitialPotential(ctl.frdm,sys,basis,pot,&fl,&fh);

  /*** allocate arrays */
  HFAllocateMemory(fl.getNsize(),fh.getNsize(),basis);


//---------------------------------------
//      Main Hartree-Fock Loop

  /*** Hartree-Fock iteration */
  HFEnergy e0, e1;
  bool   converge = false;

  e0.setNucleonMass(ENEUTRON * sys->getN() + EPROTON * sys->getZ());

  for(int n=0 ; n<sys->max_iteration ; n++){

    HFConstraintControl(n,moment);

    HFSolve(sys,bcs,&e0,basis,skyrme,moment);

    double de = 100.0;
    if(n > 0){
      de  = fabs(e0.kinetic   /e1.kinetic   - 1.0);
      de += fabs(e0.volume    /e1.volume    - 1.0);
      de += fabs(e0.surface   /e1.surface   - 1.0);
      de += fabs(e0.spinorbit /e1.spinorbit - 1.0);
      de += fabs(e0.coulomb   /e1.coulomb   - 1.0);
    }

    if(pmf.HFmonitor) HFMonitor(n,de,&e0,moment);
    if(de < sys->eps_converge){ converge = true;  break; }

    e1.copy(e0);
  }
  if(!converge){ std::cerr << "max iteration reached" << std::endl; }


//---------------------------------------
//      Output Results

  /*** print BCS parameters */
  if(pmf.bcs) HFPrintBCS(bcs);

  /*** print calculated energies */
  if(pmf.energy){
    HFPrintMoment(moment);
    HFPrintEnergy(&e0);
  }

  /*** print 3D potential data on the Gauss mesh */
  if(pmf.potential) HFPrintPotential(basis,pot,&fl,&fh);

  /*** expand single particle wavefunctions by spherical harmonics */
  for(int p=0 ; p<=1 ; p++){
    int np = (p == 0) ? sys->getN() : sys->getZ();
    HFEigenvalue(np,&bcs[p],basis,&pot[p],&spe[p]);
    HOSphericalExpansion(basis,&spe[p],matH,&spspec[p]);
  }

  /*** single particle spectra */
  if(pmf.spstate) HFPrintSPState(basis,spe);

  /*** density distribution */
  if(pmf.HFdensity){
    /*** set up mesh for printing density distribution */
    GridData dens[2];
    double dr = 0.25;  double rmax = 15.0;
    double dz = 0.25;  double zmax = 15.0;

    int nr =   rmax/dr + 1;
    int nz = 2*zmax/dz + 1;

    dens[0].memalloc(nr,nz);
    dens[1].memalloc(nr,nz);

    for(int p=0 ; p<=1 ; p++){
      int np = (p == 0) ? sys->getN() : sys->getZ();
      HFEigenvalue(np,&bcs[p],basis,&pot[p],&spe[p]);
      HFPrintDensityCalc(nr,nz,dr,dz,basis,matV,&dens[p]);
    }
    HFPrintDensity(nr,nz,dr,dz,&dens[0],&dens[1]);
  }


  /*** free allocated memory */
  HFDeleteAllocated();
}


/**********************************************************/
/*      Hartree-Fock Iteration                            */
/**********************************************************/
void HFSolve(MFTSystem *sys, BCSParameter *bcs,
             HFEnergy *energy, Basis *basis, Skyrme *skyrme, Moment *moment)
{
  /*** for neutron and proton shells */
  for(int p=0 ; p<=1 ; p++){
    int np = (p == 0) ? sys->getN() : sys->getZ();
    HFEigenvalue(np,&bcs[p],basis,&pot[p],&spe[p]);
    HFLocalDensity(basis,matV,&rho[p],&d2rho[p],&tau[p],&divJ[p],&fl,&fh);
  }

  /*** multipole momemnts, and renormalization of local densities */
  HFMultipoleMoment(basis,moment,rho,&fl,&fh);
  HFRenormalizeLocalDensity(sys->getN(),sys->getZ(),rho,tau,d2rho,divJ,&fl,&fh);

  /*** Coulomb energy */
  HFCoulomb(basis,&d2rho[1],&vcoul,&fl,&fh);

  /*** each energy component, and total energy */
  HFCalcEnergy(sys->getA(),basis,rho,tau,d2rho,divJ,&vcoul,skyrme,bcs,energy,&fl,&fh);

  /*** calculate potential */
  HFPotential(sys->getA(),basis,pot,rho,tau,d2rho,divJ,&vcoul,skyrme,moment,&fl,&fh);
}


/**********************************************************/
/*      Neutron/Proton Separable Part                     */
/**********************************************************/
void HFEigenvalue(int np, BCSParameter *bcs, Basis *basis,
                  OneBodyPotential *pot, SPEnergy *spe)
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

  /*** BCS calculation */
  HFSolveBCS(np,basis,bcs,spe,v2,matV,matH,rho,spd,&fl,&fh);

  /*** density matrix */
  HFDensityMatrix(basis,v2,matV,matH);
}


/**********************************************************/
/*      Remove Q-Constraint at Given Iteration Number     */
/**********************************************************/
inline void HFConstraintControl(int n, Moment *m)
{
  for(int i=0 ; i<5 ; i++){
    if(m->getW(i) > 0.0){
      if(m->getN(i) < n) m->setConstraint(i, 0.0, 0.0, 0);
    }
  }
}


/**********************************************************/
/*      Allocate Large Arrays                             */
/**********************************************************/
void HFAllocateMemory(const int nr, const int nz, Basis *basis)
{
  const int ndim = basis->nv;
  try{
    /*** allocate operator H, Nc blocks, and each contains NxN matrix */
    matH = new Operator [basis->nb];

    /*** allocate single particle density */
    spd  = new GridData [basis->nv];

    /*** for other memory strage */
    membuf = new double [ndim + ndim*(ndim+1)/2];
  }
  catch(std::bad_alloc &e){
    message << "MFT memory allocation error for Vector data " << e.what();
    cohTerminateCode("HFAllocateMemory");
  }

  int icon = 0;
  for(int k=0 ; k<basis->nb ; k++) icon += matH[k].memalloc(basis->ndata[k]);
  icon += spe[0].memalloc(basis->nv);
  icon += spe[1].memalloc(basis->nv);
  if(icon > 0){
    message << "memory allocation error for Operator";
    HFDeleteAllocated();
    cohTerminateCode("HFAllocateMemory");
  }

  /*** allocate arrays defined on the (nz, nr) grid */
  for(int i=0 ; i<2 ; i++){
    icon += rho[i].memalloc(nr,nz);
    icon += d2rho[i].memalloc(nr,nz);
    icon += tau[i].memalloc(nr,nz);
    icon += divJ[i].memalloc(nr,nz);
  }
  for(int i=0 ; i<basis->nv ; i++){
    icon += spd[i].memalloc(nr,nz);
  }
  icon += vcoul.memalloc(nr,nz);
  if(icon > 0){
    message << "memory allocation error for GridData";
    HFDeleteAllocated();
    cohTerminateCode("HFAllocateMemory");
  }

  int pbuf = 0;
  v2   = &membuf[pbuf]; pbuf += ndim;
  matV = &membuf[pbuf]; pbuf += ndim*(ndim+1)/2;
}


/**********************************************************/
/*      Delete Allocated Arrays                           */
/**********************************************************/
void HFDeleteAllocated()
{
  for(int i=0 ; i<2 ; i++){
    rho[i].memfree();
    d2rho[i].memfree();
    tau[i].memfree();
    divJ[i].memfree();
    spe[i].memfree();
    pot[i].memfree();
  }
  fl.memfree();
  fh.memfree();
  
  delete [] spd;
  delete [] matH;
  delete [] membuf;
}



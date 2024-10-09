/******************************************************************************/
/*  cccalc.cpp                                                                */
/*        coupled channels model calculation                                  */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "structur.h"
#include "optical.h"
#include "omcalc.h"
#include "global.h"
#include "output.h"
#include "omoutput.h"
#include "eispack.h"
#include "terminate.h"

#define CC_TOPLEVEL
#include "ewtrans.h"

static Potential        *pot, **pc1, **pc2;
static Wavefunc         *wfn;
static std::complex<double>  *Smat, *Cmat, *work[3], **Vpot;
static double           *Vmat[MAX_LAMBDA];
static int               Jmin = 0, Jmax = 0, Lmax = 0, Jcut = 0;
static NuclearModel      nucmodel = spherical;

static void ccSetCoupledLevels (const double, const double, Pdata *, Direct *);
static void ccSetCouplingPotential (const int, const double, Pdata *, ZAnumber *);
static int  ccTransmissionIndex    (Particle, const int, const int);

#undef DEBUG_CC_MATRIX
#ifdef DEBUG_CC_MATRIX
static void  ccCheckMatrix(const int m, std::complex<double> *);
#endif

#undef PRINT_CC_POTENTIAL


/**********************************************************/
/*      Coupled-Channels Calculation for Entrance         */
/**********************************************************/
int ccCalc
(const double  energy,        // CMS incident energy
 const double  excitation,    // target excitation energy
 Pdata        *proj,          // particle data for projectile
 ZAnumber     *targ,          // target ZA
 const double  mu,            // reduced mass for Targ+Proj system
 Direct       *dir,           // collective state data
 double      **tran,          // transmission coefficients (output)
 CrossSection *crx)           // calculated cross sections (output)
{
  /*** preparation */
  if( ccPreset(energy,excitation,proj,targ,mu,dir) < 0 ) return(-1);

  /*** memory allocation if ang-dist is calculated */
  if(prn.angdist && ctl.entrance) ccAllocateAmplitute(&gCol);

  /*** main JPi loop part */
  ccMain(energy,excitation,proj->pid,targ,tran,crx);

  if(prn.angdist && ctl.entrance) ccFreeAmplitude(&gCol);

  /*** clean up */
  ccCleanUp();

  return(Lmax);
}


/**********************************************************/
/*      Memory Allocation and Save Static Pointers        */
/**********************************************************/
int ccPreset(const double energy, const double excitation, Pdata *proj, ZAnumber *targ, const double mu, Direct *dir)
{
  if(energy <= 0.0) return(-1);
  if( (proj->za.getZ() > 0) && (energy < ECUT_CHARGED) ) return(-1);


//---------------------------------------
//     Memory Allocation

  try{
    /*** allocate arrays and save static pointers */
    pot = new Potential [MAX_DIRECT];
    pc1 = new Potential * [MAX_DIRECT];
    pc2 = new Potential * [MAX_DIRECT];
    for(int n=0 ; n<MAX_DIRECT ; n++){
      pc1[n] = new Potential [MAX_LAMBDA];
      pc2[n] = new Potential [MAX_LAMBDA];
    }

    pot[0].memalloc_radius(MAX_POINTS);

    for(int n=0 ; n<MAX_DIRECT ; n++){
      pot[n].memalloc_potential(MAX_POINTS);
      pot[n].memalloc_spinorbit(MAX_POINTS);
      pot[n].memalloc_coulomb(MAX_POINTS);

      pc1[n][0].memalloc_spinorbit(MAX_POINTS);
      pc2[n][0].memalloc_spinorbit(MAX_POINTS);

      for(int l=0 ; l<MAX_LAMBDA ; l++){
        pc1[n][l].memalloc_potential(MAX_POINTS);
        pc2[n][l].memalloc_potential(MAX_POINTS);
        pc1[n][l].memalloc_coulomb(MAX_POINTS);
        pc2[n][l].memalloc_coulomb(MAX_POINTS);
      }
    }

    /*** here, we do something crafty.
         external function: [number of levels] x [angular momentum]
         internal function: [total channel] x [7 radial points]     */
    wfn = new Wavefunc [MAX_CCMATRIX*MAX_CCMATRIX];
    for(int i=0 ; i<MAX_DIRECT ; i++) wfn[i].memalloc_ext(MAX_L);
    for(int i=0 ; i<MAX_CCMATRIX*MAX_CCMATRIX ; i++)  wfn[i].memalloc_int(7);

    gCdt = new CCdata[MAX_CCMATRIX];
    Smat = new std::complex<double> [MAX_CCMATRIX*MAX_CCMATRIX];
    Cmat = new std::complex<double> [MAX_CCMATRIX*MAX_CCMATRIX];
    Vpot = new std::complex<double> * [MAX_POINTS];
    for(int i=0 ; i<MAX_POINTS ; i++){
      Vpot[i] = new std::complex<double> [MAX_CCMATRIX*(MAX_CCMATRIX+1)/2];
    }
    for(int i=0 ; i<MAX_LAMBDA ; i++){
      Vmat[i] = new double [MAX_CCMATRIX*(MAX_CCMATRIX+1)/2];
    }
    for(int i=0 ; i<3 ; i++){
      work[i] = new std::complex<double> [MAX_CCMATRIX*MAX_CCMATRIX];
    }

    gCol.memalloc(MAX_DIRECT,MAX_LAMBDA);

    gHermiteEigenvalue = new double [MAX_CCMATRIX];
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("ccPreset");
  }


//---------------------------------------
//     Coupled States 

  ccSetCoupledLevels(energy,excitation,proj,dir);


//---------------------------------------
//     Energy and Coulomb Parameter

  int zzprod = targ->getZ() * proj->za.getZ();

  ccSetEnergy(energy,zzprod,mu,&gCol);
  if(prn.system && ctl.entrance){
    outCoupledState(nucmodel,gCol.nlevel,gCol.lev);
    outDeformation(nucmodel,gCol.max_lambda,gCol.beta);
  }


//---------------------------------------
//     Setup Optical Potential Geometry

  ccSetCouplingPotential(zzprod,mu,proj,targ);

  /***  free space wave function */
  Lmax = ccExternalFunction(pot->rad_match,wfn,&gCol);
  if(Lmax == 0){
    message << "coulomb function calculation error";
    cohTerminateCode("ccPreset");
  }

  /***  max J-value and equations */
  Jmin = ( (int)((fabs(gCol.lev[gCol.target_index].spin)+gCol.pspin)*2)%2 != 0 ) ? 1:0;
  Jmax = Jmin + 2*Lmax + (int)(2*fabs(gCol.lev[gCol.target_index].spin));

  /*** check maximum coupled channels */
  int nmax  = ccNumberOfChannels(&gCol);
  if(nmax >= MAX_CCMATRIX){
    message << "number of coupled equations " << nmax << " too large";
    cohTerminateCode("ccPreset");
  }

  return(nmax);
}


/**********************************************************/
/*      Copy Level Data                                   */
/**********************************************************/
void ccSetCoupledLevels(const double energy, const double excitation, Pdata *proj, Direct *dir)
{
  gCol.max_lambda = (MAX_LAMBDA-1)*2;
  gCol.beta       = (dir->type[0] == ccrot) ? dir->defstatic : dir->defdynamic;
  for(int l=0 ; l<MAX_LAMBDA ; l++){
    if(gCol.beta[l] == 0.0){
      gCol.max_lambda = (dir->type[0] == ccrot)  ? l*2 : l;
      break;
    }
  }

  
  gCol.nlevel = 1;
  gCol.lev[0].excitation = dir->lev[0].energy;
  gCol.lev[0].spin       = dir->lev[0].spin; 
  gCol.lev[0].energy     = energy;
  gCol.lev[0].phonon     = 0;

  bool vib = false;
  for(int i=1 ; i<dir->ncc ; i++){
    gCol.lev[i].excitation = dir->lev[i].energy;
    gCol.lev[i].spin       = dir->lev[i].spin;

    if(     dir->type[i] == ccvib1){ gCol.lev[i].phonon = 1; vib = true; }
    else if(dir->type[i] == ccvib2){ gCol.lev[i].phonon = 2; vib = true; }

    gCol.nlevel++;
  }

  /*** determine nuclear excitation model */
  nucmodel = (vib) ? vibration : rotation;

  gCol.target_index = (excitation > 0.0) ? ccFindTargetLevel(excitation,&gCol) : 0;
  if(gCol.target_index < 0){
    message << "target level " << excitation << " not found in the band";
    cohTerminateCode("ccSetCoupledLevels");
  }

  /*** save projectile spin */
  gCol.pspin = 0.5*proj->spin2;
}


/**********************************************************/
/*      Coupling Potential and OM Parameters              */
/**********************************************************/
void ccSetCouplingPotential(const int zzprod, const double mu, Pdata *proj, ZAnumber *targ)
{
  Optical *omp;
  omp = new Optical [gCol.nlevel];

  /*** reuse the same location for Coulomb and spin-orbit */
  for(int n=0 ; n<MAX_DIRECT ; n++){
    if(n != 0) pot[n].r2inv = pot[0].r2inv;
    for(int l=0 ; l<MAX_LAMBDA ; l++){
      pc1[n][l].r2inv = pc2[n][l].r2inv = pot[0].r2inv;
      if(l == 0) continue;
      pc1[n][l].spin_orbit = pc1[n][0].spin_orbit;
      pc2[n][l].spin_orbit = pc2[n][0].spin_orbit;
    }
  }

  /*** store energy-dependent optical potential parameters for all levels */
  for(int i=0 ; i<gCol.nlevel ; i++){

    int j = i;
    if(opt.groundstateomp) j = gCol.target_index;

    omSetOmp(proj->omp,gCol.lev[j].energy,targ->getZ(),targ->getA(),proj->za.getZ(),proj->za.getA(),&omp[i]);

    if(prn.system && ctl.entrance) outOMP(i,&omp[i]);
#ifdef PRINT_CC_POTENTIAL
    outOMPtable(i,&omp[i]);
#endif
  }

  /*** determine matching point, point pot[0] */
  omPotentialForm(zzprod,&omp[0],&gCol.lev[0],&pot[0]);

  for(int n=0 ; n<gCol.nlevel ; n++){
    if(n > 0){
      pot[n].setMatching(pot[0].n_match);  pot[n].rho_match = pot[0].rho_match;
    }
    for(int l=0 ; l<MAX_LAMBDA ; l++){
      pc1[n][l].setMatching(pot[0].n_match);  pc1[n][l].rho_match = pot[0].rho_match;
      pc2[n][l].setMatching(pot[0].n_match);  pc2[n][l].rho_match = pot[0].rho_match;
    }
  }

  /*** radial part optical potentials and nuclear matrix elements */
  if(nucmodel == vibration){
    ccPotentialFormVib(gCol.nlevel,zzprod,omp,mu,pot,pc1,pc2);
    ccMatrixElementVib(&gCol);
  }
  else{
    /*** band-head K2 value for rotational model */
    int k2 = ccBandheadSpin(&gCol);

    ccPotentialFormRot(gCol.nlevel,gCol.max_lambda,gCol.beta,zzprod,omp,mu,pot,pc1);
    ccMatrixElementRot(k2,&gCol);
  }
  
  delete [] omp;
}


/**********************************************************/
/*      Coupled-Channels Calculation for Entrance         */
/**********************************************************/
void ccMain(double energy, double excitation, Particle pid, ZAnumber *targ,
            double **tran, CrossSection *crx)
{
  double cx[MAX_DIRECT];
  std::complex<double> *s0;

  s0 = new std::complex<double> [3*MAX_L];

  /*** clear transmission and S-matrix arrays */
  for(int i=0 ; i<MAX_DIRECT ; i++){
    for(int j=0 ; j<3*MAX_L ; j++) tran[i][j]=0.0;
  }
  for(int j=0 ; j<3*MAX_L ; j++) s0[j] = std::complex<double>(0.0,0.0);


//---------------------------------------
//    Main Calculation

  /*** clear cross section array */
  crx->elastic = crx->reaction = crx->total = 0.0;
  for(int i=0 ; i<gCol.nlevel ; i++) crx->direct[i] = 0.0;

  /***  for all total spin J and parity P   */
  Jcut = Jmax;
  for(int jj=Jmin ; jj<=Jmax ; jj+=2){
    double sfact = (jj+1.0)/ ((2*gCol.pspin+1)*(2*fabs(gCol.lev[gCol.target_index].spin)+1));
    double sreac = crx->reaction;

    for(int p=-1 ; p<=1 ; p+=2){
      int nmax = ccSetChannel(jj,p,Lmax,&gCol,gCdt);

      /*** solving coupled-equations */
      ccCouplingPotential(nucmodel,nmax,jj,Vmat,Vpot,gCdt,&gCol,pot,pc1,pc2);
      ccInternalFunction(nmax,pot->n_match,pot->width,gCdt,wfn,Vpot,work);
      ccSmatrix(nmax,pot->width,gCdt,wfn,work[0],Smat,Cmat);

      /*** count number of channels that belong to the incoming channel */
      int ngs = 0;
      for(int i=0 ; i<nmax ; i++) if(gCdt[i].level == gCol.target_index) ngs++;
      if(ngs == 0) continue;

      /*** for solution of coupled-equations index I */
      for(int i=0 ; i<nmax ; i++){
        if(gCdt[i].level != gCol.target_index) continue;

        /*** for coupled-equations index J, sum |C|^2 */
        double ssum = 0.0;
        for(int j=0 ; j<nmax ; j++){
          if(!gCdt[j].open) continue;

          crx->direct[gCdt[j].level] += 4*sfact*norm(Cmat[i*nmax+j]);
          ssum += norm(Smat[i*nmax+j]);
        }
        /*** cross sections */
        crx->reaction  += sfact*(1-ssum);
        crx->total     += 4*sfact*Cmat[i*nmax+i].imag();
      }

      /*** transmission coefficients including all inverse channels */
      for(int i=0 ; i<nmax ; i++){
        for(int n=0 ; n<gCol.nlevel ; n++) cx[n] = 0.0;
        for(int j=0 ; j<nmax ; j++){
          if(gCdt[j].open) cx[gCdt[i].level] += norm(Smat[i*nmax+j]);
        }

        /*** index for transmission array for L of i-channel */
        int    k = ccTransmissionIndex(pid,gCdt[i].chn.l,gCdt[i].chn.j2);
        double x = (jj+1.0)/(gCdt[i].chn.j2+1.0)/(2*fabs(gCol.lev[gCdt[i].level].spin)+1);
        /*** avoid round-off error */
        if(cx[gCdt[i].level] > 1.0) cx[gCdt[i].level] = 1.0;

        if(gCdt[i].open) tran[gCdt[i].level][k] += x*(1.0-cx[gCdt[i].level]);

        /*** approximated average diagonal S-matrix elements this does not so good at higher L-values */
        if(gCdt[i].level == gCol.target_index){
          s0[k] += x * Smat[i*nmax+i];
        }
      }

      /*** scattering amplitudes */
      if(prn.angdist && ctl.entrance){
        ccScatteringAmplitude(jj,nmax,crx->costh,&gCol,gCdt,Cmat);
      }

      /*** print s matrix elements */
      if(prn.smatrix && ctl.entrance) outCCSmatrix(nmax,jj,p,&gCol,gCdt,Smat);

    }
    if( (sreac > 0.0) && (fabs(crx->reaction/sreac-1.0) < CRIT_LCUTCC) ){
      Jcut = jj;
      break;
    }
  }

  if(prn.smatrix && ctl.entrance) outSmatrix(Lmax,(int)(2*gCol.pspin),s0);


//---------------------------------------
//     Angular Distributions

  if(prn.angdist && ctl.entrance){
    for(int i=0 ; i<gCol.nlevel ; i++){
      for(int j=0 ; j<MAX_ANGDIST ; j++){
        crx->angdist[i][j]
          = NORM_FACT*ccAngularDistribution(i,j,crx->costh,&gCol);
      }
    }
    /*** print out every ANGLE_STEP angle points */
    outAngularDistribution(1,0,gCol.nlevel,ANGLE_STEP,targ);
  }


//---------------------------------------
//      Strength Function if E < 100keV

  double wavesq = gCol.lev[gCol.target_index].wavesq;

  if(prn.system && (energy < 0.1) && ctl.entrance && (pid == neutron)){
    omStrengthFunction(Lmax, targ->getA(), wavesq, energy, &wfn[gCol.target_index], s0);
  }

//---------------------------------------
//     Output Cross Sections

  if((pid != neutron) && (crx->reaction == 0.0)) Lmax=1;

  crx->reaction *= NORM_FACT*PI/wavesq;
  crx->total    *= NORM_FACT*PI/wavesq;
  for(int i=0 ; i<gCol.nlevel ; i++){
    crx->direct[i] *= NORM_FACT*PI/wavesq;
    if(i != gCol.target_index) crx->totaldir += crx->direct[i];
  }

  if(pid == neutron)  crx->elastic = crx->direct[gCol.target_index];

  if(ctl.entrance){
    if(prn.xsection){
      if(pid == neutron)
        outCrossSection(0,crx->total,crx->elastic,crx->reaction);
      else
        outCrossSection(1,crx->total,crx->elastic,crx->reaction);
      outLevelExcite(0,gCol.nlevel, energy, excitation, gCol.lev, crx->direct);
    }
    if(prn.transmission) outTransmission(Lmax,(int)2*gCol.pspin,NORM_FACT*PI/wavesq,tran[gCol.target_index]);
  }

  /*** CC direct changed into negative, so that SigR is not scaled */
  for(int i=0 ; i<gCol.nlevel ; i++) crx->direct[i] *= -1;

  delete [] s0;
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void ccCleanUp(void)
{
  for(int n=0 ; n<MAX_DIRECT ; n++){
    pot[n].memfree();
    for(int l=0 ; l<MAX_LAMBDA ; l++){
      pc1[n][l].memfree();
      pc2[n][l].memfree();
    }
    delete [] pc1[n];
    delete [] pc2[n];
  }
  delete [] pot;

  for(int i=0 ; i<MAX_CCMATRIX*MAX_CCMATRIX ; i++) wfn[i].memfree();
  delete [] wfn;

  delete [] gCdt;
  delete [] Smat;
  delete [] Cmat;
  for(int i=0 ; i<MAX_POINTS ; i++){
    delete [] Vpot[i];
  }
  delete [] Vpot;
  for(int i=0 ; i<MAX_LAMBDA ; i++) delete [] Vmat[i];
  for(int i=0 ; i<3 ; i++) delete [] work[i];

  gCol.memfree();

  delete [] gHermiteEigenvalue;
}


/**********************************************************/
/*      Coupled-Channels Hermitian Matrix                 */
/*      --------                                          */
/*      diagonalize P = 1-SS' matrix for open channels    */ 
/**********************************************************/
int ccHermiteMatrix(const int jj, const int p)
{
  if(jj > Jcut) return(-1);
  
  std::complex<double> *pmat = work[0];
  gHermiteEigenvector   = work[1];
  gChannelSigma         = Vmat[0]; // recycle array
  gSPhase               = Vmat[1]; // recycle

  int m = ccSetChannel(jj,p,Lmax,&gCol,gCdt); // m: total number of channels

  /*** solve schroedinger equation */
  ccCouplingPotential(nucmodel,m,jj,Vmat,Vpot,gCdt,&gCol,pot,pc1,pc2);
  ccInternalFunction(m,pot->n_match,pot->width,gCdt,wfn,Vpot,work);
  ccSmatrix(m,pot->width,gCdt,wfn,work[0],Smat,Cmat);

  /*** Satchler penetration (P) matrix */
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++) pmat[i*m+j] = std::complex<double>(0.0,0.0);
  }
  for(int i=0 ; i<m ; i++){ if(!gCdt[i].open) continue;
    for(int j=0 ; j<m ; j++){ if(!gCdt[j].open) continue;
      int ij = i*m+j;

      std::complex<double> w(0.0,0.0);
      for(int k=0 ; k<m ; k++){
        int ik = i*m+k;
        int jk = j*m+k;
        w += Smat[ik] * conj(Smat[jk]);
      }
      pmat[ij] = ((i==j) ? std::complex<double>(1.0,0.0) : std::complex<double>(0.0,0.0)) - w;
    }
  }

  /*** diagonalize P matrix, unitary matrix stored in gHermiteEigenvector */
  EISPACKHermiteDiag(m,pmat,gHermiteEigenvalue,gHermiteEigenvector);

  /*** calculate and store phase of diag-S */
  for(int i=0 ; i<m ; i++){
    std::complex <double> w(0.0,0.0);
    for(int l=0 ; l<m ; l++){
      std::complex <double>c(0.0,0.0);
      for(int k=0 ; k<m ; k++) c += gHermiteEigenvector[k*m+i]*Smat[k*m+l];
      w += c*gHermiteEigenvector[l*m+i];
    }
    gSPhase[i] = arg(w);
  }

#ifdef DEBUG_CC_MATRIX
  ccCheckMatrix(m,pmat);
#endif

  return(m);
}


/**********************************************************/
/*      Index for transmission array for L                */
/**********************************************************/
int ccTransmissionIndex(Particle p, const int l, const int j2)
{
  int k = 0;

  switch(p){
  case neutron :
  case proton  :
  case triton  :
  case helion  :
    k = l*3 + (int)((1-j2)/2.0) + l;   // [L-1/2][L+1/2][0]
    break;
  case alpha   :
    k = l*3;                           // [L][0][0]
    break;
  case deuteron:
    k = l*3 + (int)((2-j2)/2.0) + l;   // [L-1][0][L+1]
    break;
  default      :
    k = l*3;
    break;
  }

  return(k);
}


/**********************************************************/
/*      Matrix Calculation Debugging                      */
/**********************************************************/
#ifdef DEBUG_CC_MATRIX
void ccCheckMatrix(const int m, std::complex<double> *pmat)
{
  cout.setf(ios::scientific, ios::floatfield);
  cout << setprecision(6);
  cout << "M = " << m << endl;

  cout <<" Smat " << endl;
  double x = 0.0;
  for(int i=0 ; i<m ; i++){
    cout << setw(3) << gCdt[i].level;
    cout << setw(3) << gCdt[i].chn.l;
    cout << setw(3) << gCdt[i].chn.j2;
    for(int j=0 ; j<m ; j++){
      int ij = i*m+j;
      cout << setw(14) << Smat[ij].real() << setw(14) << Smat[ij].imag();
    }
    cout << endl;
  }
  x = 0.0;
  for(int i=0 ; i<m ; i++){
    cout << setw(14) << 1.0-norm(Smat[i*m+i]) << endl;
    x += 1.0-norm(Smat[i*m+i]);
  }
  cout << setw(14) << x <<" sum" << endl;;

  cout <<" Pmat " << endl;
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++){
      int ij = i*m+j;
      cout << setw(14) << pmat[ij].real() << setw(14) << pmat[ij].imag();
    }
    cout << endl;
  }

  cout <<" Eigenvalues " << endl;
  x = 0.0;
  for(int i=0 ; i<m ; i++){
    cout << setw(14) << gHermiteEigenvalue[i] << endl;
    x += gHermiteEigenvalue[i];
  }
  cout << setw(14) << x <<" sum" << endl;;

  cout <<" Umat " << endl;
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++){
      int ij = i*m+j;
      cout << setw(14) << gHermiteEigenvector[ij].real() << setw(14) << gHermiteEigenvector[ij].imag();
    }
    cout << endl;
  }

  cout <<" Unitary " << endl;
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++){
      std::complex<double> w = std::complex<double>(0.0,0.0);
      for(int k=0 ; k<m ; k++){
        int ik = i*m+k;
        int jk = j*m+k;
        w += gHermiteEigenvector[ik]*conj(gHermiteEigenvector[jk]);
      }
      if(fabs(w.real()) < 1e-15) w.real(0.0);
      if(fabs(w.imag()) < 1e-15) w.imag(0.0);
      cout << setw(14) << w.real() << setw(14) << w.imag();
    }
    cout << endl;
  }

  std::complex<double> *wmat = work[2];
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++){
      wmat[i*m+j] = std::complex <double>(0.0,0.0);
      for(int l=0 ; l<m ; l++){
        std::complex <double>c(0.0,0.0);
        for(int k=0 ; k<m ; k++) c += gHermiteEigenvector[k*m+i]*pmat[k*m+l];
        wmat[i*m+j] += c*conj(gHermiteEigenvector[l*m+j]);
      }
    }
  }

  cout <<" UPU^{-1} " << endl;
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++){
      std::complex<double> w = wmat[i*m+j];
      if(fabs(w.real()) < 1e-15) w.real(0.0);
      if(fabs(w.imag()) < 1e-15) w.imag(0.0);
      cout << setw(14) << w.real() << setw(14) << w.imag();
    }
    cout << endl;
  }


  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++){
      wmat[i*m+j] = std::complex <double>(0.0,0.0);
      for(int l=0 ; l<m ; l++){
        std::complex <double>c(0.0,0.0);
        for(int k=0 ; k<m ; k++){
          c += gHermiteEigenvector[k*m+i]*Smat[k*m+l];
        }
        wmat[i*m+j] += c*gHermiteEigenvector[l*m+j];
      }
    }
  }

  cout <<" USU " << endl;
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<m ; j++){
      std::complex<double> w = wmat[i*m+j];
      if(fabs(w.real()) < 1e-15) w.real(0.0);
      if(fabs(w.imag()) < 1e-15) w.imag(0.0);
      cout << setw(14) << w.real() << setw(14) << w.imag();
    }
    cout << endl;
  }
  x = 0.0;
  for(int i=0 ; i<m ; i++){
    cout << setw(14) << 1.0-norm(wmat[i*m+i]) << endl;
    x += 1.0-norm(wmat[i*m+i]);
  }
  cout << setw(14) << x << " sum" << endl;
}
#endif

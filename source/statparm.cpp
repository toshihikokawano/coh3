/******************************************************************************/
/*  statparm.cpp                                                              */
/*        setting up model parameters in statistical model calculation        */
/*        GDR parameter setting                                               */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "statmodel.h"
#include "levden.h"
#include "kcksyst.h"
#include "parameter.h"
#include "terminate.h"

static inline double binenergy(const double d, const int n)
{ return( d*(n+0.5) ); }


/**********************************************************/
/*      Determine Energy Bin Width in Continuum           */
/**********************************************************/
double statBinWidth(const double e)
{
  return(  (  0.0<e && e<=  0.1)* 0.02
          +(  0.1<e && e<=  2.0)* 0.05
          +(  2.0<e && e<= 10.0)* 0.1
          +( 10.0<e && e<= 20.0)* 0.1
          +( 20.0<e && e<= 50.0)* 0.2
          +( 50.0<e && e<=100.0)* 0.5
          +(100.0<e            )* 1.0 );
}


/**********************************************************/
/*     Generate Energy Mesh in the Continuum              */
/*      ----                                              */
/*          generate continuum bins above the maximum     */
/*          level excitation (elmax)                      */
/**********************************************************/
int statSetupEnergyBin(Nucleus *n)
{
  double elmax = n->lev[n->ndisc-1].energy;
  int    klow = 0, khigh = MAX_ENERGY_BIN;

  n->ncont = 0;

  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){

    /*** bin energy is the mid-point of each bin,
         Khigh is the bin number of highest excitation */
    if(binenergy(n->de,k) > n->max_energy){
      khigh = k-1;
      break;
    }
  }
  if(khigh == MAX_ENERGY_BIN){
    message << "continuum energy bins " << khigh << " too large";
    cohTerminateCode("statSetupEnergyBin");
  }

  /*** small difference between Emax and N x dE */
  double e0 = n->max_energy - binenergy(n->de,khigh);
    
  /*** find a bin just above discrete level */
  if(elmax < n->max_energy){
    for(int k=0 ; k<=khigh ; k++){
      if(binenergy(n->de,k)+e0 > elmax){
        klow = k;
        break;
      }
    }
  }else{
    klow = khigh+1;
  }

  n->ntotal = khigh+2;      // number of bins, incl. discrete region
  n->ncont  = khigh+1-klow; // number of bins in continuum only

  /*** special case, when no discrete, it should be ntotal - 1, but set equal to ntotal */
  if(elmax == 0.0) n->ncont = n->ntotal;

  /*** determine the mid-point energies in each bin.
       the sequence starts from the highest toward zero,
       so that n.excitation[0] = n.max_energy */

  for(int k=0 ; k<n->ntotal ; k++){
    n->excitation[k] = e0 + binenergy(n->de,khigh-k);
  }
  if(n->excitation[n->ntotal-1] < 0.0) n->excitation[n->ntotal-1] = 0.0;

  return(n->ncont);
}


/**********************************************************/
/*      Store Level Densities in the Continuum Array      */
/**********************************************************/
void statSetupLevelDensity(Nucleus *n, LevelDensity *ldp)
{
  double fa = parmGetFactor(parmLD  ,n->za); if(fa < 0.0) fa = 1.0;
  double fp = parmGetFactor(parmPAIR,n->za); if(fp < 0.0) fp = 1.0;
  double fs = parmGetFactor(parmSPIN,n->za); if(fs < 0.0) fs = 1.0;
  double fd = parmGetFactor(parmPDST,n->za); if((fd < 0.0) || (fd > 2.0)) fd = 1.0;

  /*** read in level density parameters */
  kckDataRead(&n->za, ldp);
  ldp->a = kckAsymptoticLevelDensity(n->za.getA()) * fa;
  ldp->pairing_energy *= fp;

  /*** determine sigma2 from discrete levels */
  double s0 = ldLevelSpinCutoff(n->ndisc, n->lev);

  /*** from systematics */
  double s1 = sqrt(kckSpinCutoff(n->za.getA()));

  if(s0 == 0.0) ldp->sigma0 = s1;
  else{
    /*** if sigma2 from the levels is very different from the systematics
         we take average of them */
    double rs = s0*s0 /(s1*s1);
    if((rs < 0.5) || (rs > 2.0)) ldp->sigma0 = sqrt((s0*s0 + s1*s1) * 0.5);
    else ldp->sigma0 = s0;
//  ldp->sigma0 = s0;
  }

  /*** re-connect ConstTemp and Fermi Gas, because Nlevel can be changed */
  ldTGconnect((double)n->za.getA(), n->lev[n->ndisc-1].energy, n->ndisc, ldp);

  /*** if no connection, use systematics */
  if(ldp->match_energy == 0.0){
    ldp->temperature = kckTemperature(n->za.getA(),ldp->shell_correct);

    if(fa != 1.0){
      /*** if f!=1 and not connected, the tweak parameter is applied to T, but its inverse */
      ldp->temperature /= fa;
    }
    /*** use ConstTemp at all excitation energies*/
    ldp->E0 = kckE0(n->za.getA(),ldp->pairing_energy,ldp->shell_correct);
    ldTextrapolate((double)n->za.getA(), n->lev[n->ndisc-1].energy, ldp);
  }

  for(int k=0 ; k<n->ntotal ; k++){
    double ex = n->excitation[k];

    /*** level density adjustment, multiplied by a factor */
    double r  = ldLevelDensity(ex,(double)n->za.getA(),ldp);
    double r1 = r * ldParityDistributionEdepend( 1,ex,fd);
    double r2 = r * ldParityDistributionEdepend(-1,ex,fd);

    for(int j=0 ; j<MAX_J ; j++){
      double sd = ldSpinDistribution(j+halfint(n->lev[0].spin),ldp->spin_cutoff,ldp->sigma0,fs);
      n->density[k][j].even = r1*sd;
      n->density[k][j].odd  = r2*sd;
    }
  }
}


/**********************************************************/
/*      Get GDR Parameters                                */
/**********************************************************/
void statSetupGdrParameter(Nucleus *n, GDR *gdr, const double beta)
{
  /*** scan GDR parameters if provided */
  bool gdrinput = false;
  for(int i=0 ; i<MAX_GDR ; i++){
    if(gdr[i].getXL() != "  "){
      gdrinput = true;
      break;
    }
  }

  /*** if no GDR parameters are given */
  if(!gdrinput) statSetupGdrSystematics(n,beta);

  /*** use input parameters */
  else{
    for(int i=0 ; i<MAX_GDR ; i++){
      if(gdr[i].getXL() != "  ") n->gdr[i] = gdr[i];
    }
  }

 /*** copy GDR parameters in the gtrans scope */
  gdrParameterSave(n->gdr,n->za.getA());
}


/**********************************************************/
/*      GDR Parameters from Internal Systematics          */
/**********************************************************/
void statSetupGdrSystematics(Nucleus *n, const double beta)
{
  int m = 0, km1 = 0, ke2 = 0;

  /*** get E1 parameters from systematics */
  if(beta == 0.0){
    gdrE1((double)n->za.getA(),beta,&n->gdr[m++]);
  }else{
    /*** second hump might be ignored if deformation is small */
    gdrE1DoubleHump0((double)n->za.getA(),beta,&n->gdr[m]);  m++;
    gdrE1DoubleHump1((double)n->za.getA(),beta,&n->gdr[m]);  if(n->gdr[m].getXL() == "E1") m++;
  }
  /*** M1 and E2 */
  gdrM1((double)n->za.getA(),&n->gdr[m]);                      km1 = m; m++;
  gdrE2((double)n->za.getZ(),(double)n->za.getA(),&n->gdr[m]); ke2 = m; m++;

  /*** renormalize M1 photo cross section */
  gdrM1norm((double)n->za.getA(),n->ldp.a,n->gdr);

  /*** M2 and E3, if needed */
  gdrM2(&n->gdr[km1],&n->gdr[m++]);
  gdrE3(&n->gdr[ke2],&n->gdr[m++]);
  
  /*** M1 scissors */
  if(beta != 0.0){
    gdrM1scissors((double)n->za.getA(),beta,&n->gdr[m]);
    n->gdr[m].setSigma(parmGetFactor(parmGDRM1) * n->gdr[m].getSigma());
//  n->gdr[m].setEnergy(parmGetFactor(parmGDRM1) * n->gdr[m].getEnergy());
//  n->gdr[m].setWidth(parmGetFactor(parmGDRM1) * n->gdr[m].getWidth());
    m++;
  }
}


/**********************************************************/
/*      Reset GDR Parameters When Given by External Data  */
/**********************************************************/
void statSetupResetGdrParameter(const int m, Nucleus *n)
{
  static std::string XL[MAX_MULTIPOL] = {"E1","M1","E2","M2","E3"};

  for(int c = 0 ; c < MAX_MULTIPOL ; c++){

    if(m > c){
      for(int i=0 ; i<MAX_GDR ; i++){
        if(n->gdr[i].getXL() == "  ") continue;

        /*** when XL is given as external data, set this as "external" */
        if(n->gdr[i].getXL() == XL[c]) n->gdr[i].setProfile(EX);
      }
    }
  }
}


/**********************************************************/
/*      Fission Related Parameter Set Up                  */
/**********************************************************/
void statSetupFissionParameter(Fission *fb)
{
  const double INERTIA_PARAMETER = 0.005;  // h-bar^2/(2I) = 5keV

  /*** tweak fission barriers */
  for(int i=0 ; i<MAX_HUMP ; i++){
    if(fb->barrier[i].height == 0.0) continue;

    /*** if inertia parameter is not given, default of 5 keV */
    if(fb->barrier[i].inertia == 0.0) fb->barrier[i].inertia = INERTIA_PARAMETER;
  }
}

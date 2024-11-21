/******************************************************************************/
/*  setupinit.cpp                                                             */
/*        setting up reaction chain                                           */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>

#include "physicalconstant.h"
#include "structur.h"
#include "coh.h"
#include "setupinit.h"
#include "nucleus.h"
#include "masstable.h"
#include "terminate.h"
#include "mt19937.h"

static void setupInitMT (void);
static bool firstcall = true;

// user-defined uncommon incident particle
static const int    userdefZ    = 0;
static const int    userdefA    = 4;
static const int    userdefJ2   = 0;
static const double userdefMass = 4.0;
static const double userdefMassExcess = 4*ENEUTRON - 0.42;
//static const double userdefMassExcess = 4*ENEUTRON - 1.0;

// static const int    userdefZ    = 0;
// static const int    userdefA    = 6;
// static const int    userdefJ2   = 0;
// static const double userdefMass = 6.0;
// static const double userdefMassExcess = 6*ENEUTRON - 1.0;


/**********************************************************/
/*     Initialize All Data                                */
/**********************************************************/
void setupInitSystem(const int mc, Pdata *pdt)
{
  /*** store particle mass, spin, and identifiers */
  for(int p=0 ; p<mc ; p++){
    pdt[p].pid = unknown;
    pdt[p].omp = 0;
    switch(p){
    case neutron:
      pdt[p].za.setZA(0,1);    pdt[p].pid = neutron;    pdt[p].spin2 = 1;
      pdt[p].mass = MNEUTRON;  pdt[p].mass_excess = ENEUTRON;
      break;
    case proton:
      pdt[p].za.setZA(1,1);    pdt[p].pid = proton;     pdt[p].spin2 = 1;
      pdt[p].mass = MPROTON;   pdt[p].mass_excess = EPROTON;
      break;
    case alpha:
      pdt[p].za.setZA(2,4);    pdt[p].pid = alpha;      pdt[p].spin2 = 0;
      pdt[p].mass = MALPHA;    pdt[p].mass_excess = EALPHA;
      break;
    case deuteron:
      pdt[p].za.setZA(1,2);    pdt[p].pid = deuteron;   pdt[p].spin2 = 2;
      pdt[p].mass = MDEUTERON; pdt[p].mass_excess = EDEUTERON;
      break;
    case triton:
      pdt[p].za.setZA(1,3);    pdt[p].pid = triton;     pdt[p].spin2 = 1;
      pdt[p].mass = MTRITON;   pdt[p].mass_excess = ETRITON;
      break;
    case helion:
      pdt[p].za.setZA(2,3);    pdt[p].pid = helion;     pdt[p].spin2 = 1;
      pdt[p].mass = MHELIUM3;  pdt[p].mass_excess = EHELIUM3;
      break;
    case gammaray:
      pdt[p].za.setZA(0,0);    pdt[p].pid = gammaray;   pdt[p].spin2 = 0;
      pdt[p].mass = 0.0;       pdt[p].mass_excess = 0.0;
      break;
    case userdef:
      pdt[p].za.setZA(userdefZ,userdefA);
      pdt[p].pid         = userdef;
      pdt[p].spin2       = userdefJ2;
      pdt[p].mass        = userdefMass;
      pdt[p].mass_excess = userdefMassExcess;
      break;
    default:
      pdt[p].za.setZA(0,0);    pdt[p].pid = unknown;    pdt[p].spin2 = 0;
      pdt[p].mass = 0.0;       pdt[p].mass_excess = 0.0;
      break;
    }
  }

  /*** initialize random number generator */
  if(firstcall){
    setupInitMT();
    firstcall = false;
  }
}


/**********************************************************/
/*     Setup Scattering Angles                            */
/*     -----                                              */
/*            as default, differential cross sections     */
/*            are calculated at every 1-deg in the CMS    */
/**********************************************************/
void setupScatteringAngle()
{
  double angstep = 1.0;
  for(int i=0 ; i<MAX_ANGDIST ; i++){
    crx.theta[i] = (double)(i+1)*angstep;
    crx.costh[i] = cos(crx.theta[i]*PI/180.0);
  }
}


/**********************************************************/
/*     Setup Reaction Chains, Mass, Binding Energy        */
/**********************************************************/
int setupChain(const int pmax, const double ex, ZAnumber *cZA, Pdata *pdt)
{
  int      n=0, *gen;
  bool     found;
  double   mass, sepa, emax, eps=1e-7;

  gen = new int [MAX_COMPOUND]; // generation of i-th compound

  /*** first compound */
  ncl[n].za          = *cZA;
  ncl[n].max_energy  = ex;
  ncl[n].mass_excess = mass_excess(ncl[n].za.getZ(),ncl[n].za.getA());
  gen[n] = 0;

  n++;

  /*** add all compounds until energy becomes negative */
  for(int i=0 ; i<MAX_COMPOUND ; i++){

    if(n >= MAX_COMPOUND){
      message << "maximum number of compound nuclei reached " << n;
      cohTerminateCode("setupChain");
    }

    if((n == i) && (i != 0)) break;

    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(pdt[j].omp == 0) continue;

      ZAnumber za = ncl[i].za - pdt[j].za;

      /*** check if residual nucleus is bigger than He */
      if( (za.getZ() <= 2) || (za.getA() <= 4) ) continue;

      /*** check if residual nucleus is already included */
      found = false;
      for(int k=0 ; k<n ; k++){
        if(za == ncl[k].za){
          sepa = ncl[k].mass_excess + pdt[j].mass_excess - ncl[i].mass_excess;
          emax = ncl[i].max_energy - sepa;
          if(fabs(1.0 - ncl[k].max_energy/emax) < eps){
            found = true;
            break;
          }
        }
      }

      /*** if this is a new nucleus, add to the NCL array */
      if( !found ){
        bool checkmass = true;
        mass = mass_excess(za.getZ(),za.getA(),&checkmass);
        sepa = mass + pdt[j].mass_excess - ncl[i].mass_excess;
        emax = ncl[i].max_energy - sepa;
        if(checkmass && (emax > 0.0) && (gen[i] < pmax)){
          cohAllocateNucleus();
          ncl[n].za.setZA(za.getZ(),za.getA());
          ncl[n].max_energy  = emax;
          ncl[n].mass_excess = mass;
          gen[n] = gen[i]+1;
          n++;
          if(n >= MAX_COMPOUND) break;
        }
      }
    }
  }

  /*** search for all open channels */
  for(int i=0 ; i<n ; i++){
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      ncl[i].cdt[j].status         = false;
      ncl[i].cdt[j].pid            = (Particle)j;
      ncl[i].cdt[j].next           =  -1;
      ncl[i].cdt[j].binding_energy = 0.0;
      ncl[i].cdt[j].spin2          = pdt[j].spin2;

      ZAnumber za = ncl[i].za - pdt[j].za;

      for(int k=0 ; k<n ; k++){
        if(za == ncl[k].za){
          sepa = ncl[k].mass_excess + pdt[j].mass_excess - ncl[i].mass_excess;
          emax = ncl[i].max_energy - sepa;
          if(fabs(1-ncl[k].max_energy/emax)<eps){
            ncl[i].cdt[j].next           = k;
            ncl[i].cdt[j].binding_energy = sepa;
            if(emax>0.0) ncl[i].cdt[j].status = true;
            break;
          }
        }
      }
    }
    ncl[i].mass = ncl[i].za.getA() + ncl[i].mass_excess / AMUNIT;
  }
/*
  for(int i=0 ; i<n ; i++){
    std::cout << std::setw(3) << ncl[i].za.getZ();
    std::cout << std::setw(3) << ncl[i].za.getA() << std::endl;
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(ncl[i].cdt[j].status){
        std::cout << std::setw(5) << i << std::setw(4) << j;
        std::setprecision(4);
        std::cout << std::setw(12) << ncl[i].cdt[j].binding_energy;
        std::cout << std::setw(12) << pdt[j].mass_excess;
        std::cout << std::setw(12) << ncl[i].max_energy-ncl[i].cdt[j].binding_energy << std::endl;
      }
    }
  }
*/
  /*** count number of particles emitted */
  for(int i=0 ; i<MAX_COMPOUND ; i++){
    gen[i] = 0;
    for(int j=0 ; j<MAX_CHANNEL ; j++) crx.prod[i].par[j] = 0;
    crx.prod[i].xsec = crx.prod[i].fiss = 0.0;
  }


  for(int i=0 ; i<n ; i++){
    for(int j=1 ; j<MAX_CHANNEL ; j++){
      int k = ncl[i].cdt[j].next;
      if(k>0){
        crx.prod[k] = crx.prod[i];
        crx.prod[k].par[j]++;
      }
    }
  }
/*
  for(int i=0 ; i<n; i++){
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      std::cout << std::setw(1) << (int)crx.prod[i].par[j];
    }
    std::cout << std::endl;
  }
*/
 
  delete [] gen;

  message << "max excitation energy " << ex << " total number of CN " << n;
  cohNotice("setupChain");

  return(n);
}


/*************************************************/
/*  Count Number of Actual Compound Nuclei       */
/*************************************************/
int setupCountCompound(const int n)
{
  /*** count number of unique nuclide appears in the chain */
  int c = 0;
  for(int i=0 ; i<n ; i++){
    bool found = false;
    for(int j=0 ; j<i ; j++){
      if(ncl[i].za == ncl[j].za){
        found = true;
        break;
      }
    }
    if(!found) c++;
  }

  if(c >= MAX_NUCLIDE){
    message << "maximum number of unique nuclide reached " << c;
    cohTerminateCode("setupCountCompound");
  }

  message << "number of unique nuclide " << c;
  cohNotice("setupCountCompound");

  return(c);
}


/*************************************************/
/*  Initialize MT Random Number Generator        */
/*************************************************/
void setupInitMT(void)
{
  unsigned long seed = 87545L;

#ifdef RANDOM_SEED_BY_SYSTEM
  random_device seedgen;
  seed = (unsigned long)seedgen();
#endif

  genrand_init(seed);
}

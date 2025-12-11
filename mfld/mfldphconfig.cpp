/******************************************************************************/
/*  mfldphconfig.cpp                                                          */
/*        Create 1-particle 1-hole Configurations                             */
/******************************************************************************/

#include <iostream>
#include <iomanip>

using namespace std;

#include "mfld.h"

static const double hole_cutoff     = 0.05;  // hole exists if v^2 > 0.1
static const double particle_cutoff = 0.95;  // particle exists if u^2 = 1 - v^2 > 0.1


/**********************************************************/
/*      Total Number of 1p-1h Configurations              */
/**********************************************************/
void PHCount(const bool sharp_cutoff, SingleParticle *sp)
{
  sp->h0 = 0;
  sp->p1 = sp->getNmax()-1;

  if(sharp_cutoff){
    /*** scan s-p energy up to Fermi energy */
    for(int i=0 ; i<sp->getNmax() ; i++){
      if(sp->getE(i) > sp->lambda){
        sp->h1 = i - 1;
        sp->p0 = i;
        break;
      }
    }
    /*** reset V2 */
    for(int i=0 ; i<sp->getNmax() ; i++){
      if(0 <= i && i <= sp->h1)  sp->setV2(i,1.0);
      else                       sp->setV2(i,0.0);
    }
  }

  else{
    /*** indices for starting and ending hole states */
    for(int i=0 ; i<sp->getNmax() ; i++){
      if(sp->getV2(i) <= hole_cutoff){ sp->h1 = i - 1; break; }
    }

    /*** indices for starting and ending particle states */
    for(int i=sp->getNmax()-1 ; i>=0 ; i--){
      if(sp->getV2(i) >= particle_cutoff){ sp->p0 = i + 1; break; }
    }
  }
/*
  std::cout << "Hole" << std::endl;
  for(int i=sp->h0 ; i<=sp->h1 ; i++){
    std::cout << std::setw(4) << i;
    std::cout << std::setw(12) << sp->getE(i);
    std::cout << std::setw(12) << sp->getV2(i);
    std::cout << std::setw(3) << sp->getP(i);
    std::cout << std::setw(4) << sp->getM(i) << std::endl;
  }

  std::cout << "Particle" << std::endl;
  for(int i=sp->p0 ; i<=sp->p1 ; i++){
    std::cout << std::setw(4) << i;
    std::cout << std::setw(12) << sp->getE(i);
    std::cout << std::setw(12) << sp->getV2(i);
    std::cout << std::setw(3) << sp->getP(i);
    std::cout << std::setw(4) << sp->getM(i) << std::endl;
  }
*/
}


/**********************************************************/
/*      Calculate Pairing Gap                             */
/**********************************************************/
double PHPairingGap(SingleParticle *sp)
{
  double dz, dn, e1, e2;

  double efz = sp[proton ].lambda;
  double ehz = sp[proton ].getE(sp[proton ].h1);
  double epz = sp[proton ].getE(sp[proton ].p0);

  e1 = abs(efz - ehz);
  e2 = abs(efz - epz);
  dz = (e1 > e2) ? e1 : e2;  
  
  double ehn = sp[neutron].getE(sp[neutron].h1);
  double epn = sp[neutron].getE(sp[neutron].p0);
  double efn = sp[neutron].lambda;

  e1 = abs(efn - ehn);
  e2 = abs(efn - epn);
  dn = (e1 > e2) ? e1 : e2;  

  double delta = (dz > dn) ? dz : dn;

  return(delta);
}

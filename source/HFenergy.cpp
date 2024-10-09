/******************************************************************************/
/*  HFenergy.cpp                                                              */
/*        calculate each energy component                                     */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "HF.h"

/**********************************************************/
/*      Calculate Energy Components                       */
/**********************************************************/
void HFCalcEnergy(const int a, Basis *basis, GridData *rho, GridData *tau, GridData *d2rho, GridData *divJ, GridData *vcoul, Skyrme *skyrme, BCSParameter *bcs, HFEnergy *energy, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  const double hbn = HBARSQ*VLIGHTSQ / (2*MNEUTRON*AMUNIT);
  const double hbz = HBARSQ*VLIGHTSQ / (2*MPROTON *AMUNIT);

  double ek = 0.0, ev = 0.0, es = 0.0, eso = 0.0, ecd = 0.0, ece = 0.0;

  for(int i=0 ; i<fl->nr ; i++){
    for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;

      double w = fl->G[i] * fh->G[j];

      double rn = rho[0].p[i][j];
      double rp = rho[1].p[i][j];  if(rp < 0.0) rp = 0.0;
      double rt = rn + rp;   if(rt < 0.0) rt = 0.0;

      double rn2 = rn * rn;
      double rp2 = rp * rp;
      double rt2 = rt * rt;

      double tn = tau[0].p[i][j];
      double tp = tau[1].p[i][j];
      double tt = tn + tp;

      double drn = d2rho[0].p[i][j];
      double drp = d2rho[1].p[i][j];
      double drt = drn + drp;
      double jn = divJ[0].p[i][j];
      double jp = divJ[1].p[i][j];
      double jt = jn + jp;

      ek  += w*(hbn*tn + hbz*tp);
      ev  += w*( skyrme->b1 * rt2
                +skyrme->b2 * (rn2 + rp2)
                +skyrme->b3 * rt * tt
                +skyrme->b4 * (rn * tn + rp * tp)
                +pow(rt,skyrme->alpha)
                 *( skyrme->b7 * rt2 + skyrme->b8 * (rn2 + rp2)) );

      es  += w*( skyrme->b5 * rt * drt
                +skyrme->b6 * (rn * drn + rp * drp));
      eso += w * skyrme->b9 * (rt * jt + rn * jn + rp * jp);
      ecd += 0.5*w * rp * vcoul->p[i][j];
      ece += w * pow(rp, 4.0/3.0);
    }
  }

  double c1 = PI/(basis->bz * basis->bp2);
  double c2 = 3.0/4.0 * COULOMBSQ * PERMITTIV * pow(3.0/PI, 1.0/3.0);

  energy->kinetic   = c1 * ek * (1.0 - 1.0/(double)a);
  energy->volume    = c1 * ev;
  energy->surface   = c1 * es;
  energy->spinorbit = c1 * eso;
  energy->coulomb   = c1 * (ecd - c2 * ece);
  energy->pairing   = bcs[0].pairing_energy + bcs[1].pairing_energy;

  energy->setTotal();

/*
  std::cout <<"Kinetic Energy     " << std::setw(12) << energy->kinetic << std:endl;
  std::cout <<"Volume Energy      " << std::setw(12) << energy->volume << std:endl;
  std::cout <<"Surface Energy     " << std::setw(12) << energy->surface << std:endl;
  std::cout <<"Spin-Orbit Energy  " << std::setw(12) << energy->spinorbit << std:endl;
  std::cout <<"Coulomb Energy     " << std::setw(12) << energy->coulomb << std:endl;
  std::cout <<"        direct     " << std::setw(12) << c1*ecd << std:endl;
  std::cout <<"        exchange   " << std::setw(12) << -c1*c2*ece << std:endl;
  std::cout <<"Pairing Energy     " << std::setw(12) << energy->pairing << std:endl;
  std::cout <<"Total Energy       " << std::setw(12) << energy->total << std:endl;
*/
}

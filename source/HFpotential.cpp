/******************************************************************************/
/*  HFpotential.cpp                                                           */
/*        Hartree-Fock potential calculation                                  */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "HF.h"
#include "FRDM.h"

#undef DEBUG

static void HFInitialPotentialWS (const int, Basis *, OneBodyPotential *, LaguerrePolynomial *, HermitePolynomial *);


/**********************************************************/
/*      Initial Potential                                 */
/**********************************************************/
void HFInitialPotential(const bool frdm, MFTSystem *sys, Basis *basis, OneBodyPotential *pot, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  /*** initial potential taken from FRDM potential */ 
  if(frdm) FRDMPotential(sys,basis,pot,fl,fh);

  /*** initial potential to be Woods-Saxon */
  else{
    HFInitialPotentialWS(sys->getA(),basis,pot,fl,fh);
  }
}


/**********************************************************/
/*      Initial Potential in the Woods Saxon Form         */
/**********************************************************/
void HFInitialPotentialWS(const int a, Basis *basis, OneBodyPotential *pot, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  const double hbn = HBARSQ*VLIGHTSQ / (2*MNEUTRON*AMUNIT);
  const double hbz = HBARSQ*VLIGHTSQ / (2*MPROTON *AMUNIT);

  double v0, vso, diffuse, r0, rso, r, f;

  const int initpot = 2;

  if(initpot == 1){
    /*** Woods-Saxon, E.Rost, PL 26B, 184 (1968) */
    v0      = -40.6;
    diffuse =   0.7;
    vso     =  31.5;
    r0      = 1.347;
    rso     = 1.28 ;
  }
  else{
    /*** Woods-Saxon, Ross-Mark-Lawson */
    v0      = -42.8;
    diffuse =  0.69;
    vso     =  39.5;
    r0      =   1.3;
    rso     =   1.3;
  }

  double a13 = pow((double)a,1.0/3.0);
  r0  *= a13;
  rso *= a13;

  double cn = HBAR * VLIGHT / (2.0*AMUNIT * MNEUTRON);
  double cz = HBAR * VLIGHT / (2.0*AMUNIT * MPROTON);

  for(int i=0 ; i<fl->nr ; i++){
    for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j == 0) continue;

      /*** radius from the center */
      r = sqrt(fl->x[i]/basis->bp2 + fh->x[j]*fh->x[j]/basis->bz2);

      /*** Woods-Saxon formfactor */
      f = 1.0/(1.0 + exp((r-r0 )/diffuse));

      pot[0].central[i][j] = v0 * f;
      pot[1].central[i][j] = v0 * f;

      /*** spin-orbit potential given not in derivative form */
      f = 1.0/(1.0 + exp((r-rso )/diffuse));

      pot[0].spinorb[i][j] = -vso * cn*cn * v0 * f;
      pot[1].spinorb[i][j] = -vso * cz*cz * v0 * f;

      pot[0].effmass[i][j] = hbn;
      pot[1].effmass[i][j] = hbz;

    }
  }

#ifdef DEBUG
  for(int j=-fh->nz ; j<=fh->nz ; j++){
    if(j == 0) continue;
    for(int i=0 ; i<fl->nr ; i++){
      std::cout << " " << std::setw(11) << sqrt(fl->x[i]/basis->bp2);
      std::cout << " " << std::setw(11) << fh->x[j]/basis->bz;
      std::cout << " " << std::setw(11) << pot[0].central[i][j];
      std::cout << " " << std::setw(11) << pot[0].spinorb[i][j];
      std::cout << " " << std::setw(11) << pot[1].central[i][j];
      std::cout << " " << std::setw(11) << pot[1].spinorb[i][j] << std::endl;
    }
    std::cout << std::endl;
  }
#endif
}


/**********************************************************/
/*      Calculate Potential for Next Iteration            */
/**********************************************************/
void HFPotential(double targa, Basis *basis, OneBodyPotential *pot,
     GridData *rho, GridData *tau, GridData *d2rho, GridData *divJ, GridData *vcoul,
     Skyrme *skyrme, Moment *moment, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  const double hbn = HBARSQ*VLIGHTSQ / (2*MNEUTRON*AMUNIT);
  const double hbz = HBARSQ*VLIGHTSQ / (2*MPROTON *AMUNIT);

  const double e2 = COULOMBSQ * PERMITTIV;
  const double aqn = 1.0;
  double dq[5];

  int iq = 0;
  dq[iq] = moment->getW(iq) * (moment->zcm - moment->getV(iq));         iq++;
  dq[iq] = moment->getW(iq) * (moment->q20 - moment->getV(iq));         iq++;
  dq[iq] = moment->getW(iq) * (moment->q30 - moment->getV(iq)*1000.0);  iq++;
  dq[iq] = moment->getW(iq) * (moment->q40 - moment->getV(iq)*10000.0); iq++;
  dq[iq] = moment->getW(iq) * (moment->qn  - moment->getV(iq));         iq++;

  double w0 = 0.4;
  double w1 = 1.0 - w0;

  for(int i=0 ; i<fl->nr ; i++){
    double r  = fl->x[i] / basis->bp2;
    double r2 = r*r;

    for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;
      double z  = fh->x[j] / basis->bz;
      double z2 = z*z;

      double rn = rho[0].p[i][j];
      double rp = rho[1].p[i][j]; if(rp < 0.0) rp = 0.0;
      double r0 = rn + rp;
      double ra = (r0 <= 0.0) ? 0.0 : pow(r0, skyrme->alpha-1.0);
      double dz = z - moment->zmin;

      /*** constraints on Qs */
      double c = 0.0;
      c += dq[0] * z * PI / (basis->bz * basis->bp2);
      c += dq[1] * (2.0 * z2 - r);
      c += dq[2] * z * (2.0 * z2 - 3.0 * r)/basis->bz;
      c += dq[3] * (0.375 * r2 + z2*z2 - 3.0 * r*z / basis->bz);
      c += dq[4] * exp(-dz*dz/aqn) / (basis->bz*basis->bp2);

      /*** calculate new potentials */
      for(int p=0 ; p<2 ; p++){

        double v = 0.0;
        v += skyrme->b1 * 2.0 * r0;
        v += skyrme->b2 * rho[p].p[i][j] * 2.0;
        v += skyrme->b3 * (tau[0].p[i][j] + tau[1].p[i][j]);
        v += skyrme->b4 * tau[p].p[i][j];
        v += skyrme->b5 * 2.0*(d2rho[0].p[i][j] + d2rho[1].p[i][j]);
        v += skyrme->b6 * d2rho[p].p[i][j] * 2.0;
        v += skyrme->b7 * (2.0 + skyrme->alpha) * ra*r0*r0;
        v += skyrme->b8 * (skyrme->alpha * (rn*rn + rp*rp) + 2.0*r0 * rho[p].p[i][j]) * ra;
        v += skyrme->b9 * (divJ[0].p[i][j] + divJ[1].p[i][j] + divJ[p].p[i][j]);

        v += c;

        if(p == 1) v += vcoul->p[i][j] - e2 * pow((3.0/PI * rp),1.0/3.0);


        pot[p].central[i][j] = w0 * pot[p].central[i][j]
                             + w1 * v;

        pot[p].effmass[i][j] = w0 * pot[p].effmass[i][j]
                             + w1 * (((p == 0) ? hbn : hbz) * (1.0 - 1.0/targa)
                             + skyrme->b3 * r0
                             + skyrme->b4 * rho[p].p[i][j]);

        pot[p].spinorb[i][j] = w0 * pot[p].spinorb[i][j] 
                             - w1 * skyrme->b9 * (r0 + rho[p].p[i][j]);
      }
    }
  }

#ifdef DEBUG
  for(int j=-fh->nz ; j<=fh->nz ; j++){
    if(j == 0) continue;
    for(int i=0 ; i<fl->nr ; i++){
      std::cout << std::setw(12) << fl->x[i] << std::setw(12) << fh->x[j];
      std::cout << std::setw(12) << pot[0].central[i][j];
      std::cout << std::setw(12) << pot[0].effmass[i][j];
      std::cout << std::setw(12) << pot[0].spinorb[i][j] << std::endl;
    }
    std::cout << std::endl;
  }
#endif
}


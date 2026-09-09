/******************************************************************************/
/*  FRDMpotential.cpp                                                         */
/*        FRDM potential calculation                                          */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "FRDM.h"


static double MicroFoldingYukawaPotentialPol (const double, const double, double *, PolarGLPoint *);
static double MicroFoldingCoulombPotentialPol (const double, double *, PolarGLPoint *);

static double MicroAverageAsymmetry(const double, MFTSystem *);
static double MicroAverageDeviation(const double, MFTSystem *, double);


#undef CHECK_SPHERICAL_LIMIT
#undef COULOMB_DIFFUSENESS

/**********************************************************/
/*      FRDM Potential                                    */
/**********************************************************/
void FRDMPotential(MFTSystem *sys, Basis *basis, OneBodyPotential *pot, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  double apot = gYukawaRangeMicro;
  double adif = 0.0;

  /*** h-bar^2/2m */
  const double hbn = HBARSQ*VLIGHTSQ / (2*MNEUTRON*AMUNIT);
  const double hbz = HBARSQ*VLIGHTSQ / (2*MPROTON *AMUNIT);

  /*** scale the shape by the factor of Rpot/R0 */
  PolarGLPoint gp1;

  double c1    = 3.0*gCharge2 / (5.0*gRadius);
  double delta = MicroAverageAsymmetry(c1,sys);
  double eps   = MicroAverageDeviation(c1,sys,delta);
  double rden  = sys->R0 * (1.0 + eps);
  double rpot  = rden + gAdensity - gBdensity / rden;
  double ratio = rpot / sys->R0;

  FRDMShapeNormalization(&gp1,sys);
  gp1.setNormalization(gp1.getNormalization() * ratio);
  FRDMShapeCreateSurface(&gp1,sys);

  /*** potential depths */
  double v0n = gVdepth - gVasymmetry * delta;
  double v0p = gVdepth + gVasymmetry * delta;

  /*** charge density, e x rhoc, corrected by volume change */
  double vrpot = 4.0/3.0 * PI * rpot*rpot*rpot;
  double erho = sys->getZ() * gCharge2 / vrpot;

  /*** spin-orbit potential depths */
  double cn = HBAR * VLIGHT / (2.0*AMUNIT * MNEUTRON);
  double cz = HBAR * VLIGHT / (2.0*AMUNIT * MPROTON);

  double vsn = (gSpinOrbitAn * sys->getA() + gSpinOrbitBn) * cn * cn;
  double vsp = (gSpinOrbitAp * sys->getA() + gSpinOrbitBp) * cz * cz;

  /*** calculate pontentials at each grid point */
  double vr0[3];
  for(int i=0 ; i<fl->nr ; i++){
    double rc = sqrt(fl->x[i]/basis->bp2); // R in cylindrical coordinate

    for(int j=-fh->nz ; j<=fh->nz ; j++){
      if(j == 0) continue;
      double zc = fh->x[j]/basis->bz;      // Z in cylindrical coordinate

      /*** r and theta in the polar coordinate */
      double rp = sqrt(rc*rc + zc*zc);
      double tp = atan(rc/zc); if(tp < 0.0) tp += PI;

      vr0[0] = rp * sin(tp); // X
      vr0[1] = 0.0;          // Y
      vr0[2] = rp * cos(tp); // Z

      double p = MicroFoldingYukawaPotentialPol(apot, adif, vr0, &gp1);
      double c = MicroFoldingCoulombPotentialPol(apot, vr0, &gp1);

      /*** central part */
      pot[0].central[i][j] = -p * v0n;
      pot[1].central[i][j] = -p * v0p;

      /*** spin-orbit part */
      pot[0].spinorb[i][j] = -pot[0].central[i][j] * vsn;
      pot[1].spinorb[i][j] = -pot[1].central[i][j] * vsp;

      /*** add Coulomb */
      pot[1].central[i][j] += c * erho;

      /*** effective mass required for solving Hamiltonian */
      pot[0].effmass[i][j] = hbn;
      pot[1].effmass[i][j] = hbz;
    }
  }

#ifdef CHECK_SPHERICAL_LIMIT
  double dpot = rpot/apot;
  for(int j=-fh->nz ; j<=fh->nz ; j++){
    if(j == 0) continue;
    double zc = fh->x[j]/basis->bz;
    for(int i=0 ; i<fl->nr ; i++){
      double rc = sqrt(fl->x[i]/basis->bp2);
      double rp = sqrt(rc*rc + zc*zc);

      rc = 0.0;

      double p = 0.0, c = 0.0;
      if(rp < rpot){
        p = 1-(1+dpot)*exp(-dpot) * sinh(rp/apot)/(rp/apot);
      }
      else{
        p = (dpot*cosh(dpot) - sinh(dpot))*exp(-rp/apot)/(rp/apot);
      }

      if(rp < rpot){
        c = PI2/3.0*erho * (3.0*rpot*rpot - rp*rp);
      }
      else{
        c = PI4/3.0*erho * rpot*rpot*rpot/rp;
      }
#ifdef COULOMB_DIFFUSENESS
      double d = -PI4 * erho * apot*apot * p;
#else
      double d = 0.0;
#endif
      double vn = -p * v0n;
      double vp = -p * v0p + c + d;

      cout << " " << setw(11) << rc;
      cout << " " << setw(11) << zc;
      cout << " " << setw(11) << pot[0].central[i][j];
      cout << " " << setw(11) << pot[0].spinorb[i][j];
      cout << " " << setw(11) << pot[1].central[i][j];
      cout << " " << setw(11) << pot[1].spinorb[i][j];
      cout << " " << setw(11) << vn;
      cout << " " << setw(11) << vp << endl;
    }
    cout << endl;
  }
#endif
}


/**********************************************************/
/*   Average Bulk Nuclear Asynmmetry, delta-bar           */
/**********************************************************/
double MicroAverageAsymmetry(const double c1, MFTSystem *sys)
{
  double x1 = sys->asym
            + (3.0 * c1 * sys->getZ() * sys->getZ())
             /(8.0 * gSurfaceStiffnessMicro * pow(sys->getA(),5.0/3.0));
  double x2 = 1 
            + (9.0 * gSymmetryEnergyMicro)
             /(4.0 * gSurfaceStiffnessMicro * sys->a13);

  return(x1/x2);
}


/**********************************************************/
/*   Average Relative Deviation in Density, eps-bar       */
/**********************************************************/
double MicroAverageDeviation(const double c1, MFTSystem *sys, double delta)
{
  double x2 = 2.0 * gSurfaceEnergyMicro / sys->a13;
  double x3 = gDensitySymmetryMicro * delta * delta;
  double x4 = c1 * sys->getZ() * sys->getZ() / (sys->a23 * sys->a23);

  return( (-x2+x3+x4)/gCompressibilityMicro );
}


/**********************************************************/
/*   Single Folding Yukawa                                */
/**********************************************************/
double MicroFoldingYukawaPotentialPol(const double a0, const double a1, double *vr0, PolarGLPoint *gp)
{
  double vr1[3],vr2[3],q1[3],q2[3],s[3],d,d2,f;

  /*** in Davies and Nix, PRC 14, 1977 (1976), two equations are given.
       Eq. (3.23) when lambda and a are different
       Eq. (3.26) for lambda = a */

  bool samerange = false;
  if(a0 == a1) samerange = true;


  double z0 = gp->getScaleFactor();
  double v0 = 0.0, v1 = 0.0;
  /*** surface integration over r0 */
  for(int i=0 ; i<gp->ng ; i++){
    double z1 = z0 * gp->w[i];

    for(int j=0 ; j<gp->ng ; j++){
      double z2 = z1 * gp->w[j];

      /*** vector to the point r(i,j) */
      gp->cartesianVector(i, j,vr1);  gp->normalVector(i, j,q1);
      gp->cartesianVector(i,-j,vr2);  gp->normalVector(i,-j,q2);

      /*** s = r0 - r1, d = sqrt |r0-r1|^2 */
      if((d2 = VectorSubtract(vr0,vr1,s)) > 0.0){
        d  = sqrt(d2);
        f = VectorProduct(s,q1) / (d2*d);
        if(samerange){
          v0 += z2 * FuncYukawa(a0,d) * f;
        }else{
          v0 += z2 * FuncDiffusePotential(a0,d) * f;
          v1 += z2 * FuncDiffusePotential(a1,d) * f;
        }
      }
      if((d2 = VectorSubtract(vr0,vr2,s)) > 0.0){
        d  = sqrt(d2);
        f = VectorProduct(s,q2) / (d2*d);
        if(samerange){
          v0 += z2 * FuncYukawa(a0,d) * f;
        }else{
          v0 += z2 * FuncDiffusePotential(a0,d) * f;
          v1 += z2 * FuncDiffusePotential(a1,d) * f;
        }
      }
    }
  }

  if(samerange){
    v0 *= -1.0/(8*PI) * 0.5;
    v1 = 0.0;
  }
  else{
    double x = a1*a1/(a1*a1-a0*a0);
    v0 *= -1.0/(4*PI) * 0.5;
    v1 *= -1.0/(4*PI) * 0.5 * x;
    v1  = v1 - x*v0;
  }
  return(v0 + v1);
}



/**********************************************************/
/*   Single Folding Coulomb                               */
/**********************************************************/
double MicroFoldingCoulombPotentialPol(const double a, double *vr0, PolarGLPoint *gp)
{
  double vr1[3],vr2[3],s[3],q1[3],q2[3],d,d2,f;

  double z0 = gp->getScaleFactor();

  double v = 0.0, dv = 0.0;
  /*** surface integration over r0 */
  for(int i=0 ; i<gp->ng ; i++){
    double z1 = z0 * gp->w[i];

    for(int j=0 ; j<gp->ng ; j++){
      double z2 = z1 * gp->w[j];

      /*** vector to the point r(i,j) */
      gp->cartesianVector(i, j,vr1);  gp->normalVector(i, j,q1);
      gp->cartesianVector(i,-j,vr2);  gp->normalVector(i,-j,q2);

      /*** s = r0 - r1, d = sqrt |r0-r1|^2 */
      if((d2 = VectorSubtract(vr0,vr1,s)) > 0.0){
        d  = sqrt(d2);
        f  = z2 * VectorProduct(s,q1) / d;
        v  += f;
        dv += f * FuncDiffusePotential(a,d) / d2;
      }
      if((d2 = VectorSubtract(vr0,vr2,s)) > 0.0){
        d  = sqrt(d2);
        f  = z2 * VectorProduct(s,q2) / d;
        v  += f;
        dv += f * FuncDiffusePotential(a,d) / d2;
      }
    }
  }

#ifndef COULOMB_DIFFUSENESS
  dv = 0.0;
#endif
  v = (-0.5*v + dv*a*a) * 0.5;

  return(v);
}


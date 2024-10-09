/******************************************************************************/
/*  FRDMmacroshapedepend.cpp                                                  */
/*        Nuclear Shape Dependent Part of Macro Energies                      */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "FRDM.h"

static double sinteg = 0.0, winteg = 0.0, w2integ = 0.0, cinteg = 0.0;

/**********************************************************/
/*   Pre-Calculation of Some Common Quantities            */
/**********************************************************/
void FRDMMacroPrecalc(double **w, PolarGLPoint *gp)
{
  double e1 = 0.0, e2 = 0.0, e3 = 0.0, e4 = 0.0, q[3];

  /*** normalVector returns dS = |V|, and
       it will be r^2 sin(t) for a spherical nucleus */
  double z0 = gp->getScaleFactor();
  for(int i=0 ; i<gp->ng ; i++){
    double z1 = z0 * gp->w[i];
    for(int j=0 ; j<gp->ng ; j++){
      double z2 = z1 * gp->w[j] * gp->normalVector(i,j,q);
      e1 += z2;                      // int dS
      e2 += z2 * w[i][j];            // int W(r) dS
      e3 += z2 * w[i][j] * w[i][j];  // int W^2(r) dS
      e4 += z2 * gp->c[i][j];        // int (1/R1 + 1/R2) dS
    }
  }

  /*** save pre-calculated quantities within this file */
  sinteg  = e1;
  winteg  = e2;
  w2integ = e3;
  cinteg  = e4;
}


/**********************************************************/
/*   W(r) = int 1/|r-r'| {dS . d(r-r')}                   */
/**********************************************************/
double FRDMMacroW(double *vr0, PolarGLPoint *gp)
{
  double vr1[3],vr2[3],q1[3],q2[3],s[3],d2;

  double z0 = gp->getScaleFactor();

  double e = 0.0;
  for(int i1=0 ; i1<gp->ng ; i1++){
    double z1 = z0 * gp->w[i1];

    for(int j1=0 ; j1<gp->ng ; j1++){
      double z2 = z1 * gp->w[j1];

      gp->cartesianVector(i1, j1,vr1); gp->normalVector(i1, j1,q1);
      gp->cartesianVector(i1,-j1,vr2); gp->normalVector(i1,-j1,q2);

      if((d2 = VectorSubtract(vr0,vr1,s)) > 0.0){
        e += z2 *  VectorProduct(s,q1) / sqrt(d2);
      }

      if((d2 = VectorSubtract(vr0,vr2,s)) > 0.0){
        e += z2 *  VectorProduct(s,q2) / sqrt(d2);
      }
    }
  }
  return(-0.25*e);
}


/**********************************************************/
/*   W-bar, by Double Surface Integration                 */
/**********************************************************/
double FRDMMacroWbar(PolarGLPoint *gp, MFTSystem *sys)
{
  double vr0[3],vr1[3],vr2[3],s[3],q0[3],q1[3],q2[3],d2;
  double w = 0.0;

  double z0 = gp->getScaleFactor() * gp->getScaleFactor();
  /*** surface integration over r0 */
  for(int i0=0 ; i0<gp->ng ; i0++){
    double z1 = z0 * gp->w[i0];

    for(int j0=0 ; j0<gp->ng ; j0++){
      double z2 = z1 * gp->w[j0];

      /*** vector to the point r0(i,j), and normal vector q0 */
      gp->cartesianVector(i0,j0,vr0);  gp->normalVector(i0,j0,q0);

      /*** surface integration over r1 */
      for(int i1=0 ; i1<gp->ng ; i1++){
        double z3 = z2 * gp->w[i1];

        for(int j1=0 ; j1<gp->ng ; j1++){
          double z4 = z3 * gp->w[j1];

          /*** vector to the point r1(i,j) and r1(i,-j) */
          gp->cartesianVector(i1, j1,vr1);  gp->normalVector(i1, j1,q1);
          gp->cartesianVector(i1,-j1,vr2);  gp->normalVector(i1,-j1,q2);

           /*** s = r0 - r1, d = sqrt |r0-r1|^2 */
          if((d2 = VectorSubtract(vr0,vr1,s)) > 0.0){
            w += z4 * VectorProduct(s,q0) * VectorProduct(s,q1) / sqrt(d2);
          }
          if((d2 = VectorSubtract(vr0,vr2,s)) > 0.0){
            w += z4 * VectorProduct(s,q0) * VectorProduct(s,q2) / sqrt(d2);
          }
        }
      }
    }
  }

  /*** factor -1/6 comes from conversion from volume to surface integral,
       factor of 1/2 comes from double counting of phi direction */
  w *= - 1.0/6.0 /sys->V0 * 0.5;

  return(w);
}


/**********************************************************/
/*   Surface Energy, Bs                                   */
/**********************************************************/
double FRDMMacroSurfaceEnergy(MFTSystem *sys)
{
  return(sinteg / sys->S0);
}


/**********************************************************/
/*   Neutron Skin Energy, Bv                              */
/**********************************************************/
double FRDMMacroNeutronSkinEnergy(const double w0, MFTSystem *sys)
{
  double e = winteg - sinteg * w0;
  return(-e * 15.0/(sys->S0 * sys->S0));
}


/**********************************************************/
/*   Surface Redistribution Energy, Bw                    */
/**********************************************************/
double FRDMMacroSurfaceRedistributionEnergy(const double w0, MFTSystem *sys)
{
  double e = w2integ - 2.0*w0*winteg + w0*w0*sinteg;
  return(e * 225.0/(sys->S0 * sys->S0 * sys->S0));
}


/**********************************************************/
/*   Curvature Energy, Bk                                 */
/**********************************************************/
double FRDMMacroCurvatureEnergy(MFTSystem *sys)
{
  double e = cinteg;
  return(e/(8.0*PI*sys->R0));
}


/**********************************************************/
/*   Volume Redistribution Energy, Br                     */
/**********************************************************/
double FRDMMacroVolumeRedistributionEnergy(const double w0, PolarGLPoint *gp, MFTSystem *sys)
{
  /*** lower order of GL quadrature is just fine, since Br behaves like 
       Br(r) = A - B r^2 +  Cr^4 */
  const int ng = 4;
  static double glx[ng/2] = {0.861136311594053, 0.339981043584856};
  static double gla[ng/2] = {0.347854845137454,	0.652145154862546};
/*
  const int ng = 16;
  static double glx[ng/2] = {
  0.98940093499164993 , 0.94457502307323258 , 0.86563120238783174 ,
  0.75540440835500303 , 0.61787624440264375 , 0.45801677765722739 ,
  0.28160355077925891 , 0.095012509837637440};
  static double gla[ng/2] = {
  0.027152459411754095, 0.062253523938647893, 0.095158511682492785,
  0.12462897125553387 , 0.14959598881657673 , 0.16915651939500254 ,
  0.18260341504492359 , 0.18945061045506850 };
*/

  double vr0[3], vr1[3], vr2[3], w1, w2, gr1, gr2, f1, f2, r1, r2;

  double e = 0.0;
  double z0 = gp->getScaleFactor();
  for(int i=0 ; i<gp->ng ; i++){
    double z1 = z0 * gp->w[i] * sin(gp->t[i]);

    for(int j=0 ; j<gp->ng ; j++){
      double z2 = z1 * gp->w[j];
      double r0 = gp->cartesianVector(i,j,vr0);

      /*** Gauss-Legendre points in the r-direction */
      double rm = r0 * 0.5;
      double z3 = rm * z2;

      double se = 0.0;
      for(int ig=0 ; ig<ng/2 ; ig++){
        gr1 = rm*(1.0 - glx[ig]);  f1 = gr1/r0;  r1 = gr1*gr1;
        gr2 = rm*(1.0 + glx[ig]);  f2 = gr2/r0;  r2 = gr2*gr2;
        for(int p=0 ; p<3 ; p++){
          vr1[p] = vr0[p] * f1;
          vr2[p] = vr0[p] * f2;
        }
        w1 = FRDMMacroW(vr1,gp) - w0;
        w2 = FRDMMacroW(vr2,gp) - w0;
        se += (r1*w1*w1 + r2*w2*w2) * gla[ig];
      }
      e += z3 * se;
    }
  }
  double c = 1575.0 /(64.0*PI*PI*PI*pow(sys->R0,7.0));

  return(c*e);
}

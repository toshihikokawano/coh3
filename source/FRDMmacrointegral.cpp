/******************************************************************************/
/*  FRDMmacrointegral.cpp                                                     */
/*        Integration Over The Shape For Macro Energy                         */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "FRDM.h"

static double MacroDerivSurfaceEnergy (const double, MFTSystem *);
static double MacroDerivCoulombEnergy (const double, MFTSystem *);
static double MacroDoubleIntegSurfaceEnergy (const double, PolarGLPoint *, MFTSystem *);
static double MacroDoubleIntegCoulombEnergy (const double, PolarGLPoint *, MFTSystem *);

#undef FOLDINGCALCTEST
#ifdef FOLDINGCALCTEST
static void   MacroIntegTest (const double, const double, PolarGLPoint *, MFTSystem *);
#endif


/**********************************************************/
/*   Double Folding to Calculate B1, B2, B3, and B4       */
/**********************************************************/
void FRDMMacroDoubleVolumeIntegral(PolarGLPoint *gp, double *b, MFTSystem *sys)
{
  double arng = gYukawaRange;
  double aden = gYukawaCharge;

  /*** B1: generalized surface (nuclear) energy */
  b[0] = MacroDoubleIntegSurfaceEnergy(arng, gp, sys);
  //MacroIntegTest(arng,aden,gp,sys);

  /*** B3: Coulomb energy */
  b[2] = MacroDoubleIntegCoulombEnergy(aden, gp, sys);

  /*** B2: derivative of surface energy */
  b[1] = MacroDerivSurfaceEnergy(arng, sys);

  /*** B4: derivative of Coulomb energy */
  b[3] = MacroDerivCoulombEnergy(aden, sys);
}


/**********************************************************/
/*   Analytical Solutions for B for Spherical Nucleus     */
/**********************************************************/
void FRDMMacroVolumeIntegralSpherical(double *b, MFTSystem *sys)
{
  double x = sys->R0/gYukawaRange;
  double y = sys->R0/gYukawaCharge;

  b[0] = 1.0 - 3/(x*x) + (1.0+x) * (2.0+3.0/x+3.0/(x*x))*exp(-2.0*x);
  b[1] = 1.0 - (1.0 + 2.0*x + 2.0*x*x) * exp(-2.0*x);
  b[2] = 1.0 - 5/(y*y) * (1.0 - 15/(8.0*y) + 21/(8.0*y*y*y) -3/4.0*(1.0 + 9/(2.0*y) + 7/(y*y) + 7/(2.0*y*y*y))*exp(-2.0*y));
  b[3] = 1.0 + 5*(-3/(y*y) + 15/(2.0*y*y*y) - 63/(4.0*y*y*y*y*y) + 3/4.0*(2/y + 12/(y*y) + 32/(y*y*y) + 42/(y*y*y*y) + 21/(y*y*y*y*y))*exp(-2.0*y));
}


/**********************************************************/
/*   Numerical Derivative of x^2 B1 to Get B2             */
/**********************************************************/
double MacroDerivSurfaceEnergy(const double a, MFTSystem *sys)
{
  double b[4], dx = 0.1, x[4];
  PolarGLPoint gd;

  /*** original shape */
  FRDMShapeNormalization(&gd,sys);
  double r0 = sys->R0;               // save R0
  double f0 = gd.getNormalization(); // shape normalization
  double w0 = r0/f0;                 // omega0 ratio
  double x0 = r0/a;

  x[0] = x0 - 2*dx;
  x[1] = x0 - dx;
  x[2] = x0 + dx;
  x[3] = x0 + 2*dx;

  for(int i=0 ; i<4 ; i++){
    sys->R0 = a*x[i];
    gd.setNormalization(sys->R0/w0);
    FRDMShapeCreateSurface(&gd,sys);
    b[i] = MacroDoubleIntegSurfaceEnergy(a,&gd,sys) * x[i]*x[i];
    //std::cout << x[i] <<" "<< b[i]/(x[i]*x[i]) <<" "<< b[i] << std::endl;
  }

  /*** restore R0 */
  sys->R0 = r0;

  double b2 = derivative(b,dx) / (2.0*x0);
  return(b2);
}


/**********************************************************/
/*   Numerical Derivative of B3/y to Get B4               */
/**********************************************************/
double MacroDerivCoulombEnergy(const double a, MFTSystem *sys)
{
  double b[4], dy = 0.1, y[4];
  PolarGLPoint gd;

  /*** original shape */
  FRDMShapeNormalization(&gd,sys);
  double r0 = sys->R0;
  double f0 = gd.getNormalization();
  double w0 = r0/f0;
  double y0 = r0/a;

  y[0] = y0 - 2*dy;
  y[1] = y0 - dy;
  y[2] = y0 + dy;
  y[3] = y0 + 2*dy;

  for(int i=0 ; i<4 ; i++){
    sys->R0 = a*y[i];
    gd.setNormalization(sys->R0/w0);
    FRDMShapeCreateSurface(&gd,sys);
    b[i] = MacroDoubleIntegCoulombEnergy(a,&gd,sys) / y[i];
  }

  /*** restore R0 */
  sys->R0 = r0;

  double b4 = derivative(b,dy) * (-y0*y0);
  return(b4);
}


/**********************************************************/
/*   Double Folding Yukawa                                */
/**********************************************************/
double MacroDoubleIntegSurfaceEnergy(const double a, PolarGLPoint *gp, MFTSystem *sys)
{
  double vr0[3],vr1[3],vr2[3],s[3],q0[3],q1[3],q2[3],d2;
  double b1 = 0.0;

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
            double d4 = d2*d2;
            double z5 = z4 * VectorProduct(s,q0) * VectorProduct(s,q1);
            b1 += FuncYukawa(a,sqrt(d2)) * z5 / d4;
          }
          if((d2 = VectorSubtract(vr0,vr2,s)) > 0.0){
            double d4 = d2*d2;
            double z5 = z4 * VectorProduct(s,q0) * VectorProduct(s,q2);
            b1 += FuncYukawa(a,sqrt(d2)) * z5 / d4;
          }
        }
      }
    }
  }

  /*** factor of 1/2 accounts for the double counting of phi1 direction */
  b1 *= -1.0/ (sys->R0 * sys->R0 * 8*PI*PI) * 0.5;

  return(b1);
}


/**********************************************************/
/*   Double Folding Coulomb                               */
/**********************************************************/
double MacroDoubleIntegCoulombEnergy(const double a, PolarGLPoint *gp, MFTSystem *sys)
{
  double vr0[3],vr1[3],vr2[3],s[3],q0[3],q1[3],q2[3],d,d2;
  double b3 = 0.0, db = 0.0;

  double z0 = gp->getScaleFactor() * gp->getScaleFactor();
  /*** surface integration over r0 */
  for(int i0=0 ; i0<gp->ng ; i0++){
    double z1 = z0 * gp->w[i0];

    for(int j0=0 ; j0<gp->ng ; j0++){
      double z2 = z1 * gp->w[j0];

      /*** vector to the point r0(i,j) */
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
            double d4 = d2*d2;
            d = sqrt(d2);
            double z5 = z4 * VectorProduct(s,q0) * VectorProduct(s,q1);
            b3 += z5 / d;
            db += FuncDiffuseCoulomb(a,d) * z5 / d4;
          }
          if((d2 = VectorSubtract(vr0,vr2,s)) > 0.0){
            double d4 = d2*d2;
            d = sqrt(d2);
            double z5 = z4 * VectorProduct(s,q0) * VectorProduct(s,q2);
            b3 += z5 / d;
            db += FuncDiffuseCoulomb(a,d) * z5 / d4;
          }
        }
      }
    }
  }

  b3 = - (b3/6.0 - db*a*a*a) * 15.0/(32.0*PI*PI*pow(sys->R0,5.0)) * 0.5;

  return(b3);
}


#ifdef FOLDINGCALCTEST

static double inline gf(const int ng, double r, double *gx, double *x){
  double d  = r*0.5;
  for(int i=0 ; i<ng/2 ; i++){
    x[i]      = d - gx[i]*d;
    x[ng-i-1] = d + gx[i]*d;
  }
  return(d);
}

void MacroIntegTest(const double ar, const double ac, PolarGLPoint *gp, MFTSystem *sys)
{
  const int ng = 32;
  static double glx[ng/2] = {
  0.997263861849481563545 ,  0.9856115115452683354002,  0.9647622555875064307738,
  0.9349060759377396891709,  0.8963211557660521239653,  0.8493676137325699701337,
  0.7944837959679424069631,  0.7321821187402896803874,  0.6630442669302152009751,
  0.5877157572407623290407,  0.5068999089322293900237,  0.4213512761306353453641,
  0.3318686022821276497799,  0.2392873622521370745446,  0.1444719615827964934852,
  0.0483076656877383162348};
  static double gla[ng/2] = {
  0.0070186100094700966004 , 0.0162743947309056706052 , 0.0253920653092620594558 ,
  0.0342738629130214331027 , 0.04283589802222668065688, 0.05099805926237617619616,
  0.0586840934785355471453 , 0.06582222277636184683765, 0.0723457941088485062254 ,
  0.0781938957870703064717 , 0.0833119242269467552222 , 0.087652093004403811143  ,
  0.0911738786957638847129 , 0.093844399080804565639  , 0.0956387200792748594191 ,
  0.0965400885147278005668};


  double wr[ng],gr0[ng],gr1[ng];
  for(int i=0 ; i<ng/2 ; i++) wr[i] = wr[ng-i-1] = gla[i];

  double v0[3],v1a[3],v1b[3];
  double u0[3],u1a[3],u1b[3], s[3], d2, d, x;

  double z0 = gp->getScaleFactor();

  double ss0 = 0.0, sc0 = 0.0;
  for(int i0=0 ; i0<gp->ng ; i0++){
    double z1 = z0 * gp->w[i0] * sin(gp->t[i0]);

    for(int j0=0 ; j0<gp->ng ; j0++){
      double r0  = gp->cartesianVector(i0, j0,v0);
      double z2  = z1 * gp->w[j0] * gf(ng,r0,glx,gr0);

      double srs0 = 0.0, src0 = 0.0;
      for(int ig0=0 ; ig0<ng ; ig0++){
        double f0 = gr0[ig0]/r0;
        for(int p0=0 ; p0<3 ; p0++) u0[p0] = v0[p0] * f0;


        double w0 = gp->getScaleFactor();
        double ss1 = 0.0, sc1 = 0.0;
        for(int i1=0 ; i1<gp->ng ; i1++){
          double w1 = w0 * gp->w[i1] * sin(gp->t[i1]);

          for(int j1=0 ; j1<gp->ng ; j1++){
            double r1 = gp->cartesianVector(i1, j1,v1a);
                        gp->cartesianVector(i1,-j1,v1b);

            double w2 = w1 * gp->w[j1] * gf(ng,r1,glx,gr1);

            double srs1 = 0.0, src1 = 0.0;
            for(int ig1=0 ; ig1<ng ; ig1++){
              double f1 = gr1[ig1]/r1;
              for(int p1=0 ; p1<3 ; p1++){
                u1a[p1] = v1a[p1] * f1;
                u1b[p1] = v1b[p1] * f1;
              }

              double ys = 0.0, yc = 0.0;
              if((d2 = VectorSubtract(u0,u1a,s)) > 0.0){
                d = sqrt(d2);
                x = d/ar;  ys += (2.0 - x)*exp(-x)/x;
                x = d/ac;  yc += (1.0 - (1.0 + 0.5*x)*exp(-x))/d;
              }

              if((d2 = VectorSubtract(u0,u1b,s)) > 0.0){
                d = sqrt(d2);
                x = d/ar;  ys += (2.0 - x)*exp(-x)/x;
                x = d/ac;  yc += (1.0 - (1.0 + 0.5*x)*exp(-x))/d;
              }

              srs1 += ys * gr1[ig1]*gr1[ig1] * wr[ig1];
              src1 += yc * gr1[ig1]*gr1[ig1] * wr[ig1];
            }
            ss1 += w2 * srs1;
            sc1 += w2 * src1;
          }
        }

        srs0 += ss1 * gr0[ig0]*gr0[ig0] * wr[ig0];
        src0 += sc1 * gr0[ig0]*gr0[ig0] * wr[ig0];
      }
      ss0 += z2 * srs0;
      sc0 += z2 * src0;
    }

    std::cout << i0 <<" " << ss0 << " " << sc0 << std::endl;
  }

  double b1 = ss0 / (sys->R0 * sys->R0 * 8*PI*PI * pow(ar,4.0)) * 0.5;
  double b3 = sc0 * 15.0 / (pow(sys->R0,5.0) * 32*PI*PI) * 0.5;

  std::cout << ss0 << " " << b1 << std::endl;
  std::cout << sc0 << " " << b3 << std::endl;
}
#endif

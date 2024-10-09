/******************************************************************************/
/*  FRDMshape.cpp                                                             */
/*        Create 3D Shape in Stretched Polar Coordinate                       */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "FRDM.h"
#include "etc.h"

static void   ShapeNormalVector (const double, const double, const double, double *, MFTSystem *);
static void   ShapeNormalVectorExact (const double, const double, double *, PolarGLPoint *, MFTSystem *);
static double ShapeCurvature (const double, const double, const double, MFTSystem *);
static double ShapeStretchedPolarAngle (const double, const double, const double, const double);
static double ShapeStretchedAzimuthalAngle (const double, const double, const double);
static double ShapeHexadecapole (const double, const double, const double, const int);

static bool exact_normalvector = true;

#undef DEBUG


/**********************************************************/
/*   Nucleus Shape                                        */
/**********************************************************/
double FRDMShape(const double theta, const double phi, MFTSystem *sys)
{
  double x1, x2, x3, x4, y1, y2, y3;

  double eps   = sys->getEps(2);
  double gamma = sys->getGamma();

  x1 = 1.0/3.0       * eps * cos(gamma);
  x2 = sqrt(1.0/3.0) * eps * sin(gamma);
  x3 = cos(theta);
  x4 = cos(2.0*phi);

  double u  = ShapeStretchedPolarAngle(x1,x2,x3,x4);
  double v  = ShapeStretchedAzimuthalAngle(x1,x2,x4);
  double uu = u*u;

  y1 = 1.0 - 2.0/3.0 * eps * cos(gamma + PI2/3.0);
  y2 = 1.0 - 2.0/3.0 * eps * cos(gamma - PI2/3.0);

  double z1 = 1.0/sqrt(y1*y2*(1.0 - 2.0*x1));

  y1 = 1.0 - x1 - 2.0*x1*x1;
  y2 = eps * (cos(gamma) + eps * cos(2.0*gamma)/ 3.0);
  y3 = x2 * (1.0 - 2.0*x1);

  double z2 = sqrt(y1 + y2*uu - y3*(1.0 - uu)*v);

  y1 = 1.0 - x1*(3.0*uu - 1.0);
  y2 = x2*(1.0 - uu)*v;
  y3 = 0.0;
  for(int i=1 ; i<=6 ; i++){
    if( (i == 2) || (i == 4) ) continue;

    eps = sys->getEps(i);
    if(eps != 0.0){
      y3 += eps * legendre(i,u);
    }
  }
  if( (eps = sys->getEps(4)) != 0.0 ){
    y3 += eps * ShapeHexadecapole(u,v,sys->getGamma(),0);
  }

  double z3 = 1.0 / sqrt(y1 + y2 + 2.0*y3);
  double r  = z1 * z2 * z3;

  return(r);
}


/**********************************************************/
/*   Noramlization Factor by Integrating Nucleus Shape    */
/**********************************************************/
void FRDMShapeNormalization(PolarGLPoint *gp, MFTSystem *sys)
{
  double eps   = sys->getEps(2);
  double gamma = sys->getGamma();

  double x1,x2,x3,y1,y2,y3;

  x1 = 1.0/3.0       * eps * cos(gamma);
  x2 = sqrt(PI*8/45.0)*eps * sin(gamma);
  x3 = sqrt(15.0/(PI*8));

  y1 = 1.0 - 2.0/3.0 * eps * cos(gamma + PI2/3.0);
  y2 = 1.0 - 2.0/3.0 * eps * cos(gamma - PI2/3.0);

  double f1 = 1.0/sqrt(y1*y2*(1.0 - 2.0*x1));
  double f2 = 0.0;

  double z0 = gp->getScaleFactor(); // Gauss-Legendre scaling factor

  /*** summation for theta, Ng points */
  for(int i=0 ; i<gp->ng ; i++){
    double z1 = z0 * gp->w[i] * sin(gp->t[i]); // GL weight x sin(theta)
    double u  = cos(gp->t[i]);

    double p1 = sys->getEps(1) * legendre(1,u);
    double p2 =                  legendre(2,u);
    double p3 = sys->getEps(3) * legendre(3,u);
    double p5 = sys->getEps(5) * legendre(5,u);
    double p6 = sys->getEps(6) * legendre(6,u);

    y1 = 1.0 - 2.0*x1*p2 + 2.0*(p1 + p3 + p5 + p6);

    /*** summation for phi, Ng points */
    for(int j=0 ; j<gp->ng ; j++){
      double z2 = z1 * gp->w[j];
      double v  = cos(2.0*gp->p[j]);

      y2 = x2*x3*(1-u*u)*v;
      y3 = 2.0 * sys->getEps(4) * ShapeHexadecapole(u,v,gamma,0);

      f2 += z2*pow(y1 + y2 + y3, -1.5);
    }
  }

  /*** f1*f2/PI4 will be the (omega-ratio)^3
       normalization factor is defined as R0/(omega-ratio) */
  double f = pow(f1*f2/PI4,-1.0/3.0) * sys->R0;
  gp->setNormalization(f);
}


/**********************************************************/
/*   3D points at given (theta,phi)                       */
/**********************************************************/
void FRDMShapeCreateSurface(PolarGLPoint *gp, MFTSystem *sys)
{
  for(int i = 0 ; i < gp->ng ; i ++){
    for(int j = 0 ; j < gp->ng ; j ++){
      gp->r[i][j] = gp->getNormalization() * FRDMShape(gp->t[i], gp->p[j], sys);
      if(exact_normalvector){
        ShapeNormalVectorExact(gp->t[i], gp->p[j], gp->v[i][j], gp, sys);
      }else{
        ShapeNormalVector(gp->t[i], gp->p[j], gp->getNormalization(), gp->v[i][j], sys);
      }
    }
  }
}


/**********************************************************/
/*   Curvatures at 3D points                              */
/**********************************************************/
void FRDMShapeStoreCurvature(PolarGLPoint *gp, MFTSystem *sys)
{
  for(int i = 0 ; i < gp->ng ; i ++){
    for(int j = 0 ; j < gp->ng ; j ++){
      gp->c[i][j] = ShapeCurvature(gp->t[i], gp->p[j], gp->getNormalization(), sys);
    }
  }
}


/**********************************************************/
/*   Normal Line From Grid, Numerical Derivertive         */
/**********************************************************/
void ShapeNormalVector(const double theta, const double phi, const double f, double *vec, MFTSystem *sys)
{
  /*** 5-point derivative = {f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)}/12h */
  double dt = 0.5/180.0 * PI;
  double dp = 0.5/180.0 * PI;

  double t[4],p[4],x[4],y[4],z[4];
  double a[2][3];

  /*** derivertive of theta */
  t[0] = theta + 2*dt;
  t[1] = theta +   dt;
  t[2] = theta -   dt;
  t[3] = theta - 2*dt;

  for(int i=0 ; i<4 ; i++){
    double r = f*FRDMShape(t[i],phi,sys);
    x[i] = r * sin(t[i]) * cos(phi);
    y[i] = r * sin(t[i]) * sin(phi);
    z[i] = r * cos(t[i]);
  }

  a[0][0] = derivative(x,dt);
  a[0][1] = derivative(y,dt);
  a[0][2] = derivative(z,dt);

  /*** derivertive of phi */
  p[0] = phi - 2*dp;
  p[1] = phi -   dp;
  p[2] = phi +   dp;
  p[3] = phi + 2*dp;

  for(int i=0 ; i<4 ; i++){
    double r = f*FRDMShape(theta,p[i],sys);
    x[i] = r * sin(theta) * cos(p[i]);
    y[i] = r * sin(theta) * sin(p[i]);
    z[i] = r * cos(theta);
  }

  a[1][0] = derivative(x,dp);
  a[1][1] = derivative(y,dp);
  a[1][2] = derivative(z,dp);

  /*** normal vector: v = dS/dt x dS/dp */
  /*** for a spherical case, |v| will be r*r*sin(theta) */

  vec[0] = a[0][2]*a[1][1] - a[0][1]*a[1][2];
  vec[1] = a[0][0]*a[1][2] - a[0][2]*a[1][0];
  vec[2] = a[0][1]*a[1][0] - a[0][0]*a[1][1];

#ifdef DEBUG
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4);
  std::cout << std::setw(12) << theta << std::setw(12) << phi << std::setw(12) << f*f*sin(theta);
  std::cout << std::setw(12) << a[0][0] << std::setw(12) << a[0][1] << std::setw(12) << a[0][2];
  std::cout << std::setw(12) << a[1][0] << std::setw(12) << a[1][1] << std::setw(12) << a[1][2];
  std::cout << std::setw(12) << vec[0]  << std::setw(12) << vec[1]  << std::setw(12) << vec[2];
  std::cout << std::setw(12) << sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]) << std::endl;
#endif
}


/**********************************************************/
/*   Normal Line From Grid, Analytical Calculation        */
/**********************************************************/
void ShapeNormalVectorExact(const double theta, const double phi, double *vec, PolarGLPoint *gp, MFTSystem *sys)
{
  double a[2][3];
  double eps   = sys->getEps(2);
  double gamma = sys->getGamma();

  double x1, x2, x3, x4, x5, y1, y2, y3, y4;

  x1 = 1.0/3.0       * eps * cos(gamma);
  x2 = sqrt(1.0/3.0) * eps * sin(gamma);
  x3 = cos(theta);
  x4 = cos(2.0*phi);
  x5 = 1.0 - x3*x3; // sin^2 theta

  /*** variables u and v */
  double u  = ShapeStretchedPolarAngle(x1,x2,x3,x4);
  double v  = ShapeStretchedAzimuthalAngle(x1,x2,x4);
  double uu = u*u;

  /*** constant C */
  y1 = 1.0 - 2.0/3.0 * eps * cos(gamma + PI2/3.0);
  y2 = 1.0 - 2.0/3.0 * eps * cos(gamma - PI2/3.0);
  double c = gp->getNormalization() / sqrt(y1*y2*(1.0 - 2.0*x1));

  /*** function f(u,v) and df/du, df/dv */
  y1 = 1.0 - x1 - 2.0*x1*x1;
  y2 = eps * (cos(gamma) + eps * cos(2.0*gamma)/ 3.0);
  y3 = x2 * (1.0 - 2.0*x1);

  double f  = y1 + y2*uu - y3*(1.0 - uu)*v;
  double fu = 2.0*y2*u + 2.0*y3*u*v;
  double fv = -y3*(1-uu);

  /*** function g(u,v) and dg/du, dg/dv */
  y1 = 1.0 - x1*(3.0*uu - 1.0);
  y2 = x2*(1.0 - uu)*v;
  y3 = 0.0;
  for(int i=1 ; i<=6 ; i++){
    if( (i == 2) || (i == 4) ) continue;
    if(sys->getEps(i) != 0.0){
      y3 += sys->getEps(i) * legendre(i,u);
    }
  }
  if( sys->getEps(4) != 0.0){
    y3 += sys->getEps(4) * ShapeHexadecapole(u,v,gamma,0);
  }
  double g = y1 + y2 + 2.0*y3;

  y1 = -6.0 * x1 * u;
  y2 = -2.0 * x2 * u*v;
  y3 = y4 = 0.0;
  for(int i=1 ; i<=6 ; i++){
    if( (i == 2) || (i == 4) ) continue;
    if(sys->getEps(i) != 0.0){
      y3 += sys->getEps(i) * dlegendre(i,u);
    }
  }
  if( sys->getEps(4) != 0.0){
    y3 += sys->getEps(4) * ShapeHexadecapole(u,v,gamma,1);
    y4  = sys->getEps(4) * ShapeHexadecapole(u,v,gamma,2);
  }
  double gu = y1 + y2 + 2.0*y3;
  double gv = x2*(1-uu) + 2.0*y4;

  /*** dr/du and dr/dv */
  y1 = 2.0*sqrt(f*g)*g;
  double ru = c * (fu*g - f*gu)/y1;
  double rv = c * (fv*g - f*gv)/y1;

  /*** du/d(cos t) and du/d(cos 2p) */
  y1 = sqrt(1.0 - 2.0*x1);
  y2 = 1.0 + x1 + x2*x4;
  y3 = pow(1.0 - x1*(3.0*x3*x3 - 1.0) + x2*x5*x4, -1.5);
  double ut = y1 * y2 * y3;
  double up = -0.5*y1 * x3 * x2 * x5 * y3;

  /*** dv/d(cos 2p) */
  y1 = 1.0 + 2.0*x1 + x1*x1 - x2*x2;
  y2 = 1.0 + x1 + x2*x4;
  double vp = y1/(y2*y2);

  /*** dr/dt and dr/dp */
  double rt = -sin(theta) * ru * ut;
  double rp = -2.0*sin(2.0*phi) * (ru*up + rv*vp);

  /*** dr(x,y,z)/dt */
  double r  = gp->getNormalization() * FRDMShape(theta,phi,sys);

  double st = sin(theta);
  double ct = cos(theta);
  double sp = sin(phi);
  double cp = cos(phi);

  a[0][0] = rt * st * cp + r * ct * cp;
  a[0][1] = rt * st * sp + r * ct * sp;
  a[0][2] = rt * ct      - r * st;

  a[1][0] = rp * st * cp - r * st * sp;
  a[1][1] = rp * st * sp + r * st * cp;
  a[1][2] = rp * ct;


  vec[0] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
  vec[1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
  vec[2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

#ifdef DEBUG
   std::cout.setf(std::ios::scientific, std::ios::floatfield);
   std::cout << std::setprecision(4);
   f = gp->getNormalization();
   std::cout << std::setw(12) << theta << std::setw(12) << phi << std::setw(12) << f*f*st;
  std::cout << std::setw(12) << a[0][0] << std::setw(12) << a[0][1] << std::setw(12) << a[0][2];
  std::cout << std::setw(12) << a[1][0] << std::setw(12) << a[1][1] << std::setw(12) << a[1][2];
  std::cout << std::setw(12) << vec[0]  << std::setw(12) << vec[1]  << std::setw(12) << vec[2];
   std::cout << std::setw(12) << sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]) << std::endl;
#endif
}


/**********************************************************/
/*   Curvature At Each Grid Point                         */
/**********************************************************/
double ShapeCurvature(const double theta, const double phi, const double f, MFTSystem *sys)
{
  double dt = 1.0/180.0 * PI;

  double t[4],p[4],q0,q1,l1,l2,d1,d2,r1,r2;
  double a[3],b[4][3],m[3],c1,c2;

  /*** center point (x,y,z) */
  q0 = f*FRDMShape(theta,phi,sys);
  a[0] = q0 * sin(theta) * cos(phi);
  a[1] = q0 * sin(theta) * sin(phi);
  a[2] = q0 * cos(theta);

  /*** plus/minus 2 delta theta */
  t[0] = theta + 2*dt;
  t[1] = theta +   dt;
  t[2] = theta -   dt;
  t[3] = theta - 2*dt;

  for(int i=0 ; i<4 ; i++){
    q1 = f*FRDMShape(t[i],phi,sys);
    b[i][0] = q1 * sin(t[i]) * cos(phi);
    b[i][1] = q1 * sin(t[i]) * sin(phi);
    b[i][2] = q1 * cos(t[i]);
  }

  /*** curvature = 1/radius in theta direction */
  midpoint(a,b[0],m);
  l1 = distance(a,b[0]);
  d1 = distance(m,b[1]);

  midpoint(a,b[3],m);
  l2 = distance(a,b[3]);
  d2 = distance(m,b[2]);

  r1 = (4*d1+l1)/(8.0*sqrt(d1));
  r2 = (4*d2+l2)/(8.0*sqrt(d2));
  r2 = (r1+r2)*0.5;
  c2 = 1.0/r2;


  /*** curvature in the perpendicular to theta
       a plain that contains a[] as well as the origin
       find 2 points on the X-Y surface */
  m[0] = m[1] = m[2] = 0.0;

  p[0] = phi + PI_2;
  p[1] = phi - PI_2;

  q1 = f*FRDMShape(PI_2,p[0],sys);
  b[0][0] = q1 * sin(p[0]);
  b[0][1] = q1 * cos(p[0]);
  b[0][2] = 0.0;

  q1 = f*FRDMShape(PI_2,p[1],sys);
  b[1][0] = q1 * sin(p[1]);
  b[1][1] = q1 * cos(p[1]);
  b[1][2] = 0.0;

  /*** curvature = 1/radius in phi direction */
  r1 = (4*q0*q0 + distance(b[0],b[1]))/(8.0*q0);
  c1 = 1.0/r1;

  return(c1 + c2);
}


/**********************************************************/
/*   Hexadecapole Deformation and Derivatives             */
/**********************************************************/
double ShapeHexadecapole(const double u, const double v, const double gamma, const int ctl)
{
  double cg = cos(gamma);
  double sg = sin(gamma);

  double a40 = (5.0 * cg * cg + 1.0) / 6.0;
  double a42 = -sqrt(30.0) * sin(2.0 * gamma) / 12.0;
  double a44 = sqrt(70.0) * sg * sg / 12.0;

  double ct = u*u;
  double st = 1.0 - ct;

  double p0 = 0.0, p2 = 0.0, p4 = 0.0;
  if(ctl == 0){
    p0 = legendre(4,u);
    p2 = sqrt( 45.0/( 32.0*PI)) * st*(7.0*ct-1.0) * v;
    p4 = sqrt(315.0/(128.0*PI)) * st*st * (2.0*v*v-1.0);
  }
  /*** dV/du */
  else if(ctl == 1){
    p0 = dlegendre(4,u);
    p2 = sqrt( 45.0/( 32.0*PI)) * (16.0- 28.0*ct) * u * v;
    p4 = sqrt(315.0/(128.0*PI)) * (-4.0)*u*st * (2.0*v*v-1.0);
  }
  /*** dV/dv */
  else{
    p0 = 0.0;
    p2 = sqrt( 45.0/( 32.0*PI)) * st*(7.0*ct-1.0);
    p4 = sqrt(315.0/(128.0*PI)) * st*st * 4.0*v;
  }
  double z = a40*p0 + sqrt(PI4/9.0)*(a42*p2 + a44*p4);

  return(z);
}


/**********************************************************/
/*   Coordinate Variable, u                               */
/**********************************************************/
double ShapeStretchedPolarAngle(const double a, const double b, const double c, const  double d)
{
  double x1 = 1.0 - 2.0 * a;
  double x2 = 1.0 - a * (3.0 * c * c - 1.0);
  double x3 = b * (1.0 - c * c) * d;

  double u  = sqrt(x1/(x2+x3)) * c;

  return(u);
}


/**********************************************************/
/*   Coordinate Variable, v                               */
/**********************************************************/
double ShapeStretchedAzimuthalAngle(const double a, const double b, const double d)
{
  double x1 = (1.0 + a) * d + b;
  double x2 = 1.0 + a + b * d;
  double v  = x1/x2;

  return(v);
}


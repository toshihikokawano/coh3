/******************************************************************************/
/*  bstatsolv.cpp                                                             */
/*        solve Schroedinger equation for bound state                         */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "optical.h"
#include "dsd.h"
#include "coupling.h"
#include "terminate.h"

extern double *fact;

static inline double lagrange(const double h, const double a, const double b, const double c, const double e, const double f, const double g){
  return( ((g-a)/60+0.15*(b-f)+0.75*(e-c))/h );
}

static int    bstateTurningPoint (const int, const int, const double, const double, Potential *);
static void   bstateExternalIntegral (const int, const int, const int, double, double, const double, const double, double *, Potential *);
static int    bstateInternalIntegral (const int, const int, const int, double, const  double, const double, double *, Potential *);
static void   bstateNonlocalCorrection (const double, const double, double *, Potential *);
static void   bstateRenormSpwave (const int, const double, double *);
static double hankel1 (const int, const double);


/**********************************************************/
/*  Search for Potential Depth for Given Binding Energy   */
/**********************************************************/
double bstateFindDepth(const int n, const int l, const int j2, const double k2, const double v0, const double alpha, double *bw, Potential *pot)
{
  int    mx,nod,im=100,jm=100,km=100;
  double dp,q,dq,qt[3],dv,purtb=0.01;

  double beta = sqrt(fabs(k2));
  double h    = pot->width;
  double q1   = hankel1(l,beta* pot->n_match   * h);
  double q2   = hankel1(l,beta*(pot->n_match-1)* h);
  double v1   = v0;

//---------------------------------------
//     Find roughly the number of nodes

  int k=0;
  for(int i=0 ; ; i++){
    if(v1 < 0.0){
      message << "abnormal potential depth " << v1;
      cohTerminateCode("bstateFindDepth");
    }
    mx  = bstateTurningPoint(l,j2,k2,v1/v0,pot);
    nod = bstateInternalIntegral(mx,l,j2,0.0,k2,v1/v0,bw,pot);
    if(nod == n) break;
    v1 += 5.0*(n-nod);
    if(i > im){
      message << "search potential depth failed (outer loop)";
      cohTerminateCode("bstateFindDepth");
    }
  }
  if(v1 == 0.0) v1=1.0;

//---------------------------------------
// Find v0 using numerical derivative

  int j=0;
  for(int i=0 ; ; i++){
    for(j=0 ; j<jm ; j++){
      double v2 =v1*(1.0+purtb);
      mx = bstateTurningPoint(l,j2,k2,v2/v0,pot);
      bstateExternalIntegral(mx,l,j2,q1,q2,k2,v2/v0,bw,pot);
      q  = bw[mx];
      dq = lagrange(h,bw[mx-3],bw[mx-2],bw[mx-1],bw[mx+1],bw[mx+2],bw[mx+3]);
      nod = bstateInternalIntegral(mx,l,j2,q,k2,v2/v0,bw,pot);
      if(nod == n) break;
      if(nod  > n) purtb *= -0.5; 
      else         purtb *=  0.5;
    }
    if(j == jm){
      v1 += 2.0*(n-nod);
      purtb = 0.01;
      k++;
      if(k > km){
        message << "search potential depth failed (inner loop)";
        cohTerminateCode("bstateFindDepth");
      }
      continue;
    }

    dp  = lagrange(h,bw[mx-3],bw[mx-2],bw[mx-1],bw[mx+1],bw[mx+2],bw[mx+3]);
    double w2  = (dp-dq)/q;

    bstateExternalIntegral(mx,l,j2,q1,q2,k2,v1/v0,bw,pot);
    q   = bw[mx]; qt[0]=bw[mx+1]; qt[1]=bw[mx+2]; qt[2]=bw[mx+3];
    dq  = lagrange(h,bw[mx-3],bw[mx-2],bw[mx-1],bw[mx+1],bw[mx+2],bw[mx+3]);
    nod = bstateInternalIntegral(mx,l,j2,q,k2,v1/v0,bw,pot);
    dp  = lagrange(h,bw[mx-3],bw[mx-2],bw[mx-1],bw[mx+1],bw[mx+2],bw[mx+3]);
    double w1  = (dp-dq)/q;

    if(fabs((dv=(w2-(1.0+purtb)*w1)/(w2-w1))-1.0) < 1.0e-05){
      bw[mx+1]=qt[0]; bw[mx+2]=qt[1]; bw[mx+3]=qt[2];
      break;
    }
    if(dv > 1.5) dv=1.1; else if(dv < 0.5) dv=0.9;
    v1 *= dv;

    if(i > im){
      message << "search potential depth failed";
      cohTerminateCode("bstateFindDepth");
    }
  }

  if(alpha > 0.0) bstateNonlocalCorrection(v0, alpha, bw, pot);
  bstateRenormSpwave(pot->n_match,h,bw);
  
  return(v1);
}


/**********************************************************/
/*      Determine Matching Point Between Ext. Int. Waves  */
/**********************************************************/
int bstateTurningPoint(const int l, const int j2, const double k2, const double v0, Potential *pot)
{
  double qr1,qr2;
  int i;

  int m   =  pot->n_match;
  double xl  =  l     *(l     +1.0);
  double so  =  j2/2.0*(j2/2.0+1.0) - xl - 3.0/4.0;
  double r   =  pot->n_match*pot->width;
  qr1= xl/(r*r)-k2-v0*pot->mean_field[m].real() - so*pot->spin_orbit[m].real();
  for(i=m-1 ; i>7 ; i--){
    r   = i*pot->width;
    qr2 = xl/(r*r)-k2-v0*pot->mean_field[i].real() - so*pot->spin_orbit[i].real();
    if((qr2*qr1) < 0.0) break;
    qr1 = qr2;
  }
/*
  if(i == 7 || i > m-6){
    message << "bound state calculation failed " << i;
    cohTerminateCode("bstateTurningPoint");
  }
*/
  return(i+3);
}


/**********************************************************/
/*      Solve External Function, --- Hankel Function      */
/**********************************************************/
void bstateExternalIntegral(const int mx, const int l, const int j2, double pr1, double pr2, const double k2, const double v0, double *bw, Potential *pot)
{
  double r,rr,ar1,ar2,ar3,qr1,qr2,qr3,pr3;

  int    m   =  pot->n_match;
  double d12 =  pot->width*pot->width/12.0;
  double xl  =  l     *(l     +1.0);
  double so  =  j2/2.0*(j2/2.0+1.0)-xl- 3.0/4.0;

  bw[m  ]  = pr1;
  bw[m-1]  = pr2;
  r   = m  *pot->width; rr=r*r;
  qr1 = xl/rr-k2-v0*pot->mean_field[m  ].real() - so*pot->spin_orbit[m  ].real();
  r   =(m-1)*pot->width; rr=r*r;
  qr2 = xl/rr-k2-v0*pot->mean_field[m-2].real() - so*pot->spin_orbit[m-1].real();

  for(int i=m-2 ; i>=mx-3 ; i--){
    r   = i*pot->width; rr=r*r;
    qr3 = xl/rr-k2-v0*pot->mean_field[i].real() - so*pot->spin_orbit[i].real();

    ar1 = 1.0     -d12*qr1;
    ar2 = 2.0+10.0*d12*qr2;
    ar3 = 1.0     -d12*qr3;

    bw[i] = pr3 = (ar2*pr2-ar1*pr1)/ar3;

    qr1 = qr2; qr2 = qr3;
    pr1 = pr2; pr2 = pr3;
  }
  return;
}


/**********************************************************/
/*  Solve Internal Function,  Fox-Goodwin method          */
/**********************************************************/
int bstateInternalIntegral(const int m, const int l, const int j2, double prm, const double k2, const double v0, double *bw, Potential *pot)
{
  double r,rr,ar1,ar2,ar3,qr1,qr2,qr3,pr1,pr2,pr3;

  int    i0  =  (int)(l/3.0+1.0);
  double d12 =  pot->width*pot->width/12.0;
  double xl  =  l     *(l     +1.0);
  double so  =  j2/2.0*(j2/2.0+1.0)-xl- 3.0/4.0;

  bw[0] = pr1 = 0.0;
  bw[1] = pr2 = pow(pot->width,l+1.0);

  int node=0;
  r   = pot->width; rr=r*r;
  qr1 = 0.0;
  qr2 = xl/rr-k2-v0*pot->mean_field[1].real() - so*pot->spin_orbit[1].real();

  for(int i=2 ; i<=m+3 ; i++){
    r   = i*pot->width; rr=r*r;
    qr3 = xl/rr-k2-v0*pot->mean_field[i].real() - so*pot->spin_orbit[i].real();

    ar1 = 1.0     -d12*qr1;
    ar2 = 2.0+10.0*d12*qr2;
    ar3 = 1.0     -d12*qr3;

    bw[i] = pr3= (i<i0) ? pow(r,l+1.0) : (ar2*pr2-ar1*pr1)/ar3;

    qr1 = qr2; qr2 = qr3;
    pr1 = pr2; pr2 = pr3;
    if((pr1*pr2) < 0.0) node++;
  }

  if(prm != 0.0){
    prm = prm/bw[m];
    for(int i=0 ; i<=m+3 ; i++) bw[i] *= prm;
  }
  return(node);
}


/**********************************************************/
/*  Non-Locality Correction                               */
/**********************************************************/
void bstateNonlocalCorrection(const double v0, const double a, double *bw, Potential *pot)
{
   for(int i=0 ; i<=pot->n_match ; i++){
     bw[i] *= 1.0/sqrt(1.0 + a*v0*pot->mean_field[i].real());
   }
}


/**********************************************************/
/*  Re-Normalize Bound Satate Wave Function               */
/**********************************************************/
void bstateRenormSpwave(const int m, const double h, double *bw)
{
   double sum = 0.0;
   for(int i=0 ; i<=m ; i++) sum += bw[i]*bw[i]*h;
   sum = sqrt(1.0/sum);
   if(bw[1] < 0.0) sum *= -1.0;
   for(int i=1 ; i<=m ; i++) bw[i] *= sum;
}


/**********************************************************/
/*  Hankel function                                       */
/**********************************************************/
double hankel1(const int l, const double r)
{
  double x, y, sum;
  int k;

  x=y=1.0 / (2*r);
  sum=1.0;
  for(k=1;k<=l;k++){
    sum += exp(fact[l+k]-fact[k]-fact[l-k]) * x;
    x   *= y;
  }
  return( exp(-r) * sum );
}

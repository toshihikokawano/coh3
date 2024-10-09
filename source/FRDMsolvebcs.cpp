/******************************************************************************/
/*  FRDMsolvebcs.cpp                                                          */
/*        Lipkin-Nogami BCS Model for FRDM                                    */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "FRDM.h"
#include "polysq.h"

static void   LNSolveBCS          (LNParameter *, SPEnergy *);
static void   LNEvaluate          (LNParameter *, SPEnergy *, double *, double *);
static double LNLambda2           (const int, const int, double, double *);
static double LNchisquare         (const int, double *, double *, double *);
static void   LNInteractingOrbit  (const int, LNParameter *, SPEnergy *);
static double LNPairingStrength   (const double, const double, LNParameter *);

static const bool calcmonitor = false;

#undef SOLUTION_CHECK
#ifdef SOLUTION_CHECK
static void   LNSolutionCheck     (LNParameter *, SPEnergy *);
#endif

/**********************************************************/
/*      Lipkin-Nogami BCS model                           */
/**********************************************************/
void FRDMSoleveBCS(MFTSystem *sys, double rhobar, LNParameter *bcs, SPEnergy *spe)
{
  const double ecut = 5.0;    // determine active orbits, starting below Ef - ecut
  const int nv      = spe->ne;

  bcs->cutoff = ecut;

  /*** initial parameters */
  bcs->pairing_eff  = sys->bs * gPairingGapMicro * pow((double)bcs->np,-1/3.0);
  bcs->pairing_mac  = ( (bcs->np % 2) == 0)
                    ? 0.0 : sys->bs * gPairingGap * pow((double)bcs->np,-1/3.0);

  bcs->pairing_gap = sys->bs * gPairingGap * pow((double)bcs->np,-1/3.0);
  bcs->chemical_potential = (spe->energy[bcs->np/2-1]+spe->energy[bcs->np/2])*0.5;
  bcs->lambda2  = 2.0;

  /*** determine single particle states to be included */
  LNInteractingOrbit(nv,bcs,spe);

  /*** pairing interaction strength, G */
  bcs->strength = LNPairingStrength(rhobar,bcs->pairing_eff,bcs);

  /*** determine chemical potential and pairing gap */
#ifdef SOLUTION_CHECK
  LNSolutionCheck(bcs,spe);
#else
  LNSolveBCS(bcs,spe);

  if(calcmonitor){
    std::cout << std::setw(12) << bcs->strength;
    std::cout << std::setw(12) << bcs->chemical_potential;
    std::cout << std::setw(12) << bcs->pairing_gap;
    std::cout << std::setw(12) << bcs->lambda2;
    std::cout << std::endl;
  }
#endif
}


/**********************************************************/
/*      Solve Gap Equation for lambda and delta           */
/*      ----                                              */
/*      Solving the LN gap equation for lambda and delta  */
/*      with the Levenberg-Marquardt method.              */
/*      Lambda2 can be determined by the given G and      */
/*      single-particle energies, when lambda and delta   */
/*      are provided.                                     */
/**********************************************************/
void LNSolveBCS(LNParameter *bcs, SPEnergy *spe)
{
  const double eps = 1e-10, pt = 0.005;
  const int n = 2, m = 2, kmax = 1000;
  double p0[m], p1[m], pa[m], pb[m], pc[m], pd[m],
         d[n], y[n], ya[n], yb[n], yc[n], yd[n], w[n], z[n*m];

  /*** initial parameters */
  p0[0] =  p0[1] = 0.0;

  /*** normalized data and weight */
  d[0] = 1.0;      d[1] = 1.0;
  w[0] = 1.0;      w[1] = 1.0;

  /*** non-linear LS equation by Levenberg-Marquardt method */
  double damping = 100.0;

  LNEvaluate(bcs,spe,p0,y);
  double x0 = LNchisquare(n,y,d,w);

  int k = 0;
  for(k=0 ; k<kmax ; k++){

    /*** monitor iteration for convergence */
    if(calcmonitor){
      std::cout << std::setw(6) << k;
      std::cout << std::setw(12) << y[0] << std::setw(12) << y[1];
      std::cout << std::setw(12) << bcs->chemical_potential + p0[0];
      std::cout << std::setw(12) << bcs->pairing_gap + p0[1];
      std::cout << std::setw(12) << bcs->lambda2;
      std::cout << std::setw(12) << damping << std::setw(12) << x0 << std::endl;
    }

    /*** output difference, data - calc. */
    for(int i=0 ; i<n ; i++) y[i] = d[i] - y[i];

    /*** Jacobian */
    for(int j=0 ; j<m ; j++){
      for(int k=0 ; k<m ; k++) pa[k] = pb[k] = pc[k] = pd[k] = p0[k];

      pa[j] = p0[j] - 2.0*pt;
      pb[j] = p0[j] -     pt;
      pc[j] = p0[j] +     pt;
      pd[j] = p0[j] + 2.0*pt;

      LNEvaluate(bcs,spe,pa,ya);
      LNEvaluate(bcs,spe,pb,yb);
      LNEvaluate(bcs,spe,pc,yc);
      LNEvaluate(bcs,spe,pd,yd);

      /*** 4-point Lagrange derivative */
      for(int i=0 ; i<n ; i++) z[i*m+j] = (ya[i] - 8.0*yb[i] + 8.0*yc[i] - yd[i]) / (pt * 12.0);
    }

    LSQMarquardt(n,m,damping,p0,p1,y,w,z);

    if(p1[1] + bcs->pairing_gap < 0.0) p1[1] = fabs(p1[1] + bcs->pairing_gap) -  + bcs->pairing_gap;

    LNEvaluate(bcs,spe,p1,y);
    double x1 = LNchisquare(n,y,d,w);

    /*** re-adjust damping factor */
    if(x1 < x0){
      if(x1 < eps) break;
      else{
        damping *= 0.5;
        for(int j=0 ; j<m ; j++) p0[j] = p1[j];
      }
    }
    else damping *= 2.0;

    x0 = x1;
  }
  if(k == kmax) std::cerr << "WARNNING : Lipkin-Nogami didnt converge" << std::endl;

  bcs->chemical_potential += p1[0];
  bcs->pairing_gap        += p1[1];//  if(bcs->pairing_gap < 0.0) bcs->pairing_gap *= -1.0;
}


/**********************************************************/
/*      Solve Lambda2 Equation                            */
/**********************************************************/
void LNEvaluate(LNParameter *bcs, SPEnergy *spe, double *p, double *y)
{
  const double eps = 1e-8;
  const int kmax = 10;

  double lambda  = p[0] + bcs->chemical_potential;
  double delta   = p[1] + bcs->pairing_gap;
  double deltasq = delta * delta;

  bcs->lambda2 = LNLambda2(bcs->n1,bcs->n2,bcs->strength,spe->v2);

  double px = 0.0, gx = 0.0;
  for(int k1=0 ; k1<kmax ; k1++){
    px = 0.0;
    gx = 0.0;
    for(int i=bcs->n1 ; i<=bcs->n2 ; i++){

      double e  = spe->energy[i] + (4*bcs->lambda2 - bcs->strength)*spe->v2[i];
      double x1 = e - lambda;
      double x2 = sqrt(x1 * x1 + deltasq);
      spe->v2[i] = (1.0 - x1/x2)*0.5;
      px += spe->v2[i] * 2.0;
      gx += 1.0 / x2;
    }
    px += 2.0*bcs->n1;
    gx  = 2.0/gx;

    double p = LNLambda2(bcs->n1,bcs->n2,bcs->strength,spe->v2);
    double z = fabs(1.0 - bcs->lambda2/p);

    bcs->lambda2 = (p+bcs->lambda2)*0.5;

    if(z < eps) break;
  }

  y[0] = px / bcs->np;
  y[1] = gx / bcs->strength;
}


/**********************************************************/
/*      Calculate Chi-Square                              */
/**********************************************************/
double LNchisquare(const int n, double *y, double *d, double *w)
{
  double x = 0.0;
  for(int i=0 ; i<n ; i++) x += (y[i] - d[i]) * (y[i] - d[i]) / (d[i]*d[i] * w[i]*w[i]);

  return(x / n);
}


/**********************************************************/
/*      lambda2, number fluctuation                       */
/**********************************************************/
double LNLambda2(const int n1, const int n2, double g, double *v2)
{
  double lambda2 = 0.0;

  /*** to avoid round-off error, the products of summations 
       are first expanded, and dropped u^4 v^4 terms in the sum */
  double y1 = 0.0, y2 = 0.0;
  for(int i=n1 ; i<=n2 ; i++){
    double uui = 1.0 - v2[i];
    double vvi = v2[i];
    double ui  = sqrt(uui);
    double vi  = sqrt(v2[i]);

    for(int j=n1 ; j<=n2 ; j++){
      if(i == j) continue;
      double uuj = 1.0 - v2[j];
      double vvj = v2[j];
      double uj  = sqrt(uuj);
      double vj  = sqrt(v2[j]);

      y1 += uui * ui * vi * uj * vvj * vj;
      y2 += uui * vvi * uuj * vvj;
    }
  }

  if(y2 != 0.0) lambda2 = 0.25 * g * y1/y2;

//  if(lambda2 < 0.25*g) lambda2 = 0.25 * g;

  return(lambda2);
}


/**********************************************************/
/*      Active Particles                                  */
/*      -----                                             */
/*      Note on determining N1 and N2                     */
/*                                                        */
/*  Even case: Np = 8           Odd case: Np = 9          */
/*  7 -----                     7 -----                   */
/*  6 -----                     6 -x-x- < N2=5=2Nf-N1     */
/*  5 -x-x- < N2=5 = 2Nf-N1+1   5 -x-x-                   */
/*  4 -x-x-                     4 -o-x- < Nf              */
/*  3 -o-o- < Nf                3 -o-o-                   */
/*  2 -o-o- < N1=2 = 2Nf-N2+1   2 -o-o- < N1=2=2Nf-N2     */
/*  1 -----                     1 -----                   */
/*  0 -----                     0 -----                   */
/**********************************************************/
void LNInteractingOrbit(const int n, LNParameter *bcs, SPEnergy *spe)
{
  bcs->n1 = 0; // index of the lowest level
  bcs->n2 = 0; // highest level
  bcs->nf = 0; // level just above Ef

  /*** determine symmetric energy range, lambda +/- dN */
  bcs->nf = ((bcs->np%2) == 0) ? bcs->np/2 - 1 : bcs->np/2;

  int n2a = 0, n2b = 0;
  for(int i=bcs->nf ; i<n ; i++){
    if(spe->energy[i] > bcs->chemical_potential + bcs->cutoff){ n2a = i; break; }
  }
  for(int i=bcs->nf ; i<n ; i++){
    if(spe->energy[i] > 0.0){ n2b = i-1; break; }
  }
  bcs->n2 = (n2a > n2b) ? n2a : n2b;

  /*** even case, nf is just below Ef */
  if( (bcs->np % 2) == 0 ){
    bcs->n1 = 2*bcs->nf - bcs->n2 + 1;
  }
  /*** odd case, the last particle on the Nf-th level */
  else{
    bcs->n1 = 2*bcs->nf - bcs->n2;
  }

  /*** if too few levels, extend the range */
  if( ((bcs->nf - bcs->n1) < 2) && (bcs->n1 >= 1) ){
    bcs->n1 --;
    bcs->n2 ++;
  }

  /** in case n1 < 0, narrower the region */
  if(bcs->n1 < 0){
    bcs->n2 += bcs->n1;
    bcs->n1 = 0;
  }

  for(int i=0 ; i<n ; i++) spe->v2[i] = (i<=bcs->nf) ? 1.0 : 0.0;
/*
  for(int i=bcs->n1 ; i<=bcs->n2 ; i++){
    if(i == bcs->nf){ cout << i << " > "; }
    else{             cout << i << "   "; }
    cout << spe->energy[i] << " " << spe->v2[i] << endl;
  }
*/
}


/**********************************************************/
/*      Pairing Strength G                                */
/**********************************************************/
double LNPairingStrength(const double rhobar, const double deltag, LNParameter *bcs)
{
  double y1 = (-0.5*bcs->np - 0.5 + bcs->n1)/rhobar;
  double y2 = (-0.5*bcs->np + 0.5 + bcs->n2)/rhobar;
//double y1 = (-0.5*bcs->np +1 + bcs->n1)/rhobar;
//double y2 = (-0.5*bcs->np    + bcs->n2)/rhobar;

  double x1 = log(sqrt(y1 * y1 + deltag * deltag) + y1);
  double x2 = log(sqrt(y2 * y2 + deltag * deltag) + y2);

  double g = 2.0 / ( rhobar * (x2 - x1) );

  return(g);
}


#ifdef SOLUTION_CHECK
/**********************************************************/
/*      For Debugging, Check Equation Solution            */
/**********************************************************/
void LNSolutionCheck(LNParameter *bcs, SPEnergy *spe)
{
  const int n = 2, m = 2;
  double p0[m], d[n], y[n], w[n];

  d[0] = bcs->np;  d[1] = bcs->strength;
  w[0] = 1.0;      w[1] = 1.0;

  double lambda_min = -10.0;
  double lambda_max =  -5.0;
  double delta_min = 0.0;
  double delta_max = 1.2;

  int div = 50;
  double dl = (lambda_max - lambda_min)/div;
  double dd = (delta_max  - delta_min )/div;

  for(int i = 1 ; i<= div ; i++){
    p0[0] = i*dl + lambda_min  - bcs->chemical_potential;
    for(int j = 1 ; j<= div ; j++){
      p0[1] = j*dd + delta_min - bcs->pairing_gap;

      LNEvaluate(bcs,spe,p0,y);
      double x0 = LNchisquare(n,y,d,w);

      std::cout << p0[0] + bcs->chemical_potential <<" " << p0[1] + bcs->pairing_gap << " ";
      std::cout << bcs->lambda2 << " ";
      std::cout << y[0]/d[0] << " " << y[1]/d[1] << " " << x0 << std::endl;
    }
    std::cout << std::endl;
  }
  return;
}
#endif

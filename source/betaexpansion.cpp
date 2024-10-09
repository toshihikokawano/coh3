/******************************************************************************/
/*  betaexpansion.cpp                                                         */
/*        expand FRDM shape in spherical harmonics                            */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "FRDM.h"
#include "global.h"
#include "etc.h"


static double BetaIntegration (const int, const int, PolarGLPoint *);


/**********************************************************/
/*   Spherical Harmonics Beta Expansion                   */
/**********************************************************/
void BetaExpansion(const int lmax, double **alm, PolarGLPoint *gp, MFTSystem *sys)
{
  /*** calculate expansion coefficient a(l,m) for the 3D shape */
  for(int l=0 ; l<=lmax ; l++){
    for(int m=-l ; m<=l ; m++){
      BetaIntegration(l,m,gp);
      alm[l][m] = BetaIntegration(l,m,gp);
    }
  }

  /*** create new 3D shape reconsutructed by the expansion */
  if(pmf.FRDMshape){
    PolarCoordinate q;

    for(int i = 0 ; i < q.nt ; i ++){
      for(int j = 0 ; j < q.np ; j ++) q.r[i][j] = sys->R0;
    }

    for(int l=1 ; l<=lmax ; l++){
      for(int m=-l ; m<=l ; m++){
        for(int i = 0 ; i < q.nt ; i ++){
          for(int j = 0 ; j < q.np ; j ++){
            q.r[i][j] += alm[l][m] * SphericalHarmonics(l,m,q.t[i],q.p[j]);
          }
        }
      }
    }
    FRDMWrite3DMesh(fileFRDMBeta,&q);
  }
}


/**********************************************************/
/*   Calculate Expansion Coefficient                      */
/**********************************************************/
double BetaIntegration(const int l, const int m, PolarGLPoint *gp)
{
  double a = 0.0;
  double z0 = gp->st * gp->sp;

  for(int i=0 ; i<gp->ng ; i++){
    double z1 = z0 * gp->w[i] * sin(gp->t[i]);

    for(int j=0 ; j<gp->ng ; j++){
      double z2 = z1 * gp->w[j];
      double r  = gp->getR(i,j);

      /*** [0,pi] and [0,-pi] are symmetric */
      a += z2 * r * SphericalHarmonics(l,m,gp->t[i], gp->p[j]);
      a += z2 * r * SphericalHarmonics(l,m,gp->t[i],-gp->p[j]);
    }
  }
  return(a);
}


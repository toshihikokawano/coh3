/******************************************************************************/
/*  nilsson.cpp                                                               */
/*        Nilson model (n,l,j2) - Bengtsson                                   */
/******************************************************************************/

#include <iostream>

#include "structur.h"
#include "optical.h"
#include "dsd.h"


/**************************************/
/*   Shell Model Quantum Numbers      */
/**************************************/
class ShellOrbit{
 public:
    int     n     ;
    int     l     ;
    int     j2    ;
    double  energy;
};

static const int N_LEVEL = 45;
static const int Z_LEVEL = 45;

static ShellOrbit nshell[N_LEVEL]={
{0, 0, 1, 1.5000},{0, 1, 3, 2.3800},{0, 1, 1, 2.7400},{0, 2, 5, 3.2900},{1, 0, 1, 3.5000},
{0, 2, 3, 3.8150},{0, 3, 7, 4.1625},{1, 1, 3, 4.5675},{0, 3, 5, 4.7925},{1, 1, 1, 4.8375},
{0, 4, 9, 5.0562},{1, 2, 5, 5.5784},{0, 4, 7, 5.6862},{2, 0, 1, 5.8822},{0, 5,11, 5.9234},
{1, 2, 3, 5.9284},{1, 3, 7, 6.5273},{0, 5, 9, 6.6054},{0, 6,13, 6.8118},{2, 1, 3, 6.9179},
{1, 3, 5, 6.9613},{2, 1, 1, 7.1039},{1, 4, 9, 7.3996},{0, 6,11, 7.6178},{0, 7,15, 7.7275},
{2, 2, 5, 7.8187},{1, 4, 7, 7.9576},{3, 0, 1, 8.0692},{2, 2, 3, 8.1287},{1, 5,11, 8.2706},
{0, 8,17, 8.5526},{0, 7,13, 8.6575},{2, 3, 7, 8.6848},{1, 5, 9, 8.9526},{3, 1, 3, 8.9700},
{2, 3, 5, 9.1188},{3, 1, 1, 9.1560},{1, 6,13, 9.1602},{0, 8,15, 9.6066},{2, 4, 9, 9.6389},
{1, 6,11, 9.9662},{3, 2, 5, 9.9886},{2, 4, 7,10.1969},{4, 0, 1,10.2093},{3, 2, 3,10.2986}};
      	 
static ShellOrbit zshell[Z_LEVEL]={
{0, 0, 1, 1.5000},{0, 1, 3, 2.3800},{0, 1, 1, 2.7400},{0, 2, 5, 3.2900},{1, 0, 1, 3.5000},
{0, 2, 3, 3.8150},{0, 3, 7, 4.1490},{1, 1, 3, 4.5990},{0, 3, 5, 4.7790},{1, 1, 1, 4.8690},
{0, 4, 9, 5.0177},{0, 4, 7, 5.6027},{1, 2, 5, 5.6664},{0, 5,11, 5.8100},{1, 2, 3, 5.9914},
{2, 0, 1, 6.0187},{0, 5, 9, 6.4700},{0, 6,13, 6.6171},{1, 3, 7, 6.6320},{1, 3, 5, 7.0520},
{2, 1, 3, 7.1420},{0, 6,11, 7.3191},{2, 1, 1, 7.3220},{0, 7,15, 7.3395},{1, 4, 9, 7.5448},
{0, 8,17, 8.0247},{1, 4, 7, 8.0308},{0, 7,13, 8.1495},{2, 2, 5, 8.1745},{1, 5,11, 8.4163},
{2, 2, 3, 8.4445},{3, 0, 1, 8.5060},{0, 8,15, 8.9427},{1, 5, 9, 9.0103},{2, 3, 7, 9.1950},
{1, 6,13, 9.2505},{2, 3, 5, 9.5730},{3, 1, 3, 9.6756},{3, 1, 1, 9.8376},{1, 6,11, 9.9525},
{2, 4, 9,10.1782},{2, 4, 7,10.6642},{3, 2, 5,10.8079},{3, 2, 3,11.0779},{4, 0, 1,11.1394}};


/**********************************************************/
/*      Fill Nucleons in Each Shell                       */
/*      If shell is closed, return zero (no pairing)      */
/**********************************************************/
int nilssonFillOrbit(const Particle p, const int n)
{
  int nfermi,nmax;
  ShellOrbit *orb;

  orb  = (p == proton) ? zshell : nshell;
  nmax = (p == proton) ? Z_LEVEL : N_LEVEL;

  int m=0;
  for(nfermi=0 ; nfermi<nmax ; nfermi++){
    m += orb[nfermi].j2+1; if(m >= n) break;
  }
  if(m == n) nfermi++;

  return(nfermi);
}


/**********************************************************/
/*      Nillson Levels                                    */
/**********************************************************/

static const double VN = 41.0;
static const double VP = 41.0;

int nilssonLevel(const Particle p, const int m, const int mf, const unsigned int a, Bstate *bst, double sep, double bmax, double fshift)
{
  double vn,vp;

  if(m < 0) return(-2);

  vn  = VN*pow((double)a,-1.0/3.0);
  vp  = VP*pow((double)a,-1.0/3.0);

  if(p == neutron){
    if(m >= N_LEVEL) return(-1);
    bst->n     = nshell[m].n;
    bst->l     = nshell[m].l;
    bst->j2    = nshell[m].j2;
    bst->bind  = -sep - (nshell[mf].energy - nshell[m].energy)*vn + fshift;
  }
  else{
    if(m >= Z_LEVEL) return(-1);
    bst->n     = zshell[m].n;
    bst->l     = zshell[m].l;
    bst->j2    = zshell[m].j2;
    bst->bind  = -sep - (zshell[mf].energy - zshell[m].energy)*vp;
  }

/* weak bound assumption */
  if(bmax >= 0.0){
    if(bst->bind > 0.0 && bst->bind < bmax){
      bst->bind  = -0.1;
    }
  }

  return(0);
}

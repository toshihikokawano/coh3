/******************************************************************************/
/*  parleg.cpp                                                                */
/*        particle emission angular distribution                              */
/******************************************************************************/

#include <iostream>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "nucleus.h"
#include "global.h"
#include "coupling.h"
#include "etc.h"
#include "polysq.h"

static void specBiedenharn (const int, const int, const double, const int, const int, const int, const int, const int, const int, Transmission *, Nucleus *);

static void specBiedenharnContinuum (const double, const int, const int, const int, const int, const int, const int, Transmission *, Nucleus *, double ***);


/**********************************************************/
/*      Particle Angular Distribution : Level Excite      */
/*      -----                                             */
/*             Legendre Coefficients for                  */
/*             Compound Nuclear Reaction                  */
/**********************************************************/

static int i0 = 0, spin0 = 0, ir = 0, spin1 = 0;

void specLegendreCoefficient(const int targid, const int incid, const int targlev, const double coef, const int pcn, const int jcn, const int lp0, const int jp0, Transmission **td)
{
  Nucleus *n1;
  int c0 = 0;

  /*** set global scope parameters */
  spin0 = ncl[c0].cdt[incid].spin2;                 // projectile spin
  i0    = (int)(2.0*ncl[targid].lev[targlev].spin); // target state spin

  double g = coef*(jcn+1.0)/PI4;

  /*** For all particle emission channels */
  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!ncl[c0].cdt[id].status) continue;

    n1    = &ncl[ncl[c0].cdt[id].next];
    spin1 = ncl[c0].cdt[id].spin2;
    specBiedenharn(incid,targlev,g,ncl[c0].jmax,id,pcn,jcn,lp0,jp0,td[id],n1);
  }
}


/**********************************************************/
/*      Legendre Coefficients in Continuum (optional)     */
/**********************************************************/
void specLegendreCoefficientContinuum(const double coef, const int pcn, const int jcn, const int lp0, const int jp0, Transmission **tc, double ***cleg)
{
  Nucleus *n1;
  int c0 = 0;
  double g = coef*(jcn+1.0)/PI4;

  /*** For all particle emission channels */
  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!ncl[c0].cdt[id].status) continue;

    n1    = &ncl[ncl[c0].cdt[id].next];
    spin1 = ncl[c0].cdt[id].spin2;
    ir    = (int)(2.0*halfint(n1->lev[0].spin));    // residual ground state spin (half or int)
    specBiedenharnContinuum(g,ncl[c0].jmax,id,pcn,jcn,lp0,jp0,tc[id],n1,cleg);
  }
}


/**********************************************************/
/*      Legendre Coefficients for fixed J, l0,j0          */
/*      using the Biedenharn formula                      */
/**********************************************************/
void specBiedenharn(const int incid, const int targlev, const double g, const int jmax, const int id, const int pcn, const int jcn, const int lp0, const int jp0, Transmission *td, Nucleus *n1)
{
  /*** Loop over Legendre coefficient index L */
  /*** all odd-terms are zero, so that we dont calculate */
  for(int l=0 ; l<2*jmax ; l+=4){
    int idx = l/2;

    double z0 = z_coefficient(lp0,jp0,lp0,jp0,spin0,l); if(z0 == 0.0) continue;
    double w0 = racah(jp0,jcn,jp0,jcn,i0,l);            if(w0 == 0.0) continue;
    double c0 = z0*w0*g;

    /*** Loop over discrete state of residual nucleus up to max number */
    int k1max = std::min(n1->ndisc, MAX_ANGDISTLEVELS);

    for(int k1=0 ; k1<k1max ; k1++){
      if(td[k1].lmax <= 0) continue;

      int i1  = (int)(2.0*n1->lev[k1].spin);
      int p1  = n1->lev[k1].parity;
      int phi = parity( spin0-spin1-i0+i1 );

      /*** For exit channel Tlj coupled to JP */
      for(int lp1=0 ; lp1<=2*td[k1].lmax ; lp1+=2){
        if(parity(lp1) != pcn*p1) continue;

        for(int sp1=spin1 ; sp1>=-spin1 ; sp1-=2){
          int jp1 = lp1 + sp1;  if(jp1 < 0) continue;

          if((i1 >= std::abs(jcn-jp1)) && (i1 <= jcn+jp1)){
            double z1 = z_coefficient(lp1,jp1,lp1,jp1,spin1,l); if(z1 == 0.0) continue;
            double w1 = racah(jp1,jcn,jp1,jcn,i1,l);            if(w1 == 0.0) continue;
            double t1 = td[k1].tran[tj_index(lp1,sp1,spin1)];   if(t1 == 0.0) continue;
            double wf = (ctl.fluctuation) ? statMoldauer(lp1,jp1,id,k1,t1,0.0) : 1.0;
            double c1 = z1*w1*t1;

            /*** Check if compound elastic */
            double c2 = 0.0;
            if(ctl.fluctuation && (id == incid) && (k1 == targlev)){
              if((lp0 != lp1) || (jp0 != jp1)){
                double z2 = z_coefficient(lp0,jp0,lp1,jp1,spin0,l);
                double w2 = racah(jcn,jp0,jcn,jp1,i0,l);
                c2 = z2*z2*w2*w2 * g*t1;
              }
            }
            crx.legcoef[id][k1][idx] += phi*(c0*c1 + c2)*wf;
          }
        }
      }
    }
  }
}


/**********************************************************/
/*      Biedenharn formula for Continuum (optional)       */
/**********************************************************/
void specBiedenharnContinuum(const double g, const int jmax, const int id, const int pcn, const int jcn, const int lp0, const int jp0, Transmission *tc, Nucleus *n1, double ***cleg)
{
  /*** Loop over Legendre coefficient index L */
  for(int l=0 ; l<2*jmax ; l+=4){
    int idx = l/2;

    double z0 = z_coefficient(lp0,jp0,lp0,jp0,spin0,l); if(z0 == 0.0) continue;
    double w0 = racah(jp0,jcn,jp0,jcn,i0,l);            if(w0 == 0.0) continue;
    double c0 = z0*w0*g;

    /*** Loop over continuum bins in the residual nucleus */
    for(int k1=0 ; k1<n1->ncont ; k1++){
      if(tc[k1].lmax <= 0) continue;

      for(int i1=ir ; i1<=2*jmax+ir ; i1+=2){ // spin distribution at excitation[k]
        int jdx = (i1-ir)/2;

        for(int p1=-1 ; p1<=1 ; p1+=2){       // parity distribution

          int phi = parity( spin0-spin1-i0+i1 );
          double rho = (p1 == 1) ? n1->density[k1][jdx].even : n1->density[k1][jdx].odd;

          /*** For exit channel Tlj coupled to JP */
          for(int lp1=0 ; lp1<=2*tc[k1].lmax ; lp1+=2){
            if(parity(lp1) != pcn*p1) continue;

            for(int sp1=spin1 ; sp1>=-spin1 ; sp1-=2){
              int jp1 = lp1 + sp1;  if(jp1 < 0) continue;

              if((i1 >= std::abs(jcn-jp1)) && (i1 <= jcn+jp1)){
                double z1 = z_coefficient(lp1,jp1,lp1,jp1,spin1,l); if(z1 == 0.0) continue;
                double w1 = racah(jp1,jcn,jp1,jcn,i1,l);            if(w1 == 0.0) continue;
                double t1 = tc[k1].tran[tj_index(lp1,sp1,spin1)];   if(t1 == 0.0) continue;
                double c1 = z1*w1*t1*rho;
                cleg[id][k1][idx] += phi*c0*c1;
              }
            }
          }
        }
      }
    }
  }
}


/**********************************************************/
/*      Calculate Angular Distribution from Leg. Coef.    */
/*      for Inelastic Channel, and Add Direct             */
/**********************************************************/
int statInelAngularDistribution(const int id, const int jmax, Transmission *td)
{
  Nucleus *n1;

  n1 = &ncl[ncl[0].cdt[id].next];

  int k, kmax = std::min(n1->ndisc, MAX_ANGDISTLEVELS);
  for(k=0 ; k<kmax ; k++){
    if(td[k].lmax <= 0) break;

    /*** if direct is given */
    bool dcx = false;
    if(crx.angdist[k][0] != 0.0) dcx = true;

    for(int i=0 ; i<MAX_ANGDIST ; i++){

      double x = 0.0;
      for(int l=0 ; l<jmax ; l+=2){
        x += crx.legcoef[id][k][l] * legendre(l,crx.costh[i]);
        if(crx.legcoef[id][k][l] == 0.0) break;
      }

      crx.angdist[k][i] += x;
    }

    /*** Legendre coefficients for the differential cross section by LESQ */
    if(dcx){
      for(int j=0 ; j<MAX_J ; j++) crx.legcoef[id][k][j] = 0.0;
      LSQLegendre(true,MAX_ANGDIST,jmax,crx.theta,crx.angdist[k],crx.legcoef[id][k]);
    }
  }

  return(k);
}


/**********************************************************/
/*      Calculate Angular Distribution                    */
/*      for All Other Channels                            */
/**********************************************************/
int statAngularDistribution(const int id, const int jmax, Transmission *td)
{
  Nucleus *n1;

  n1 = &ncl[ncl[0].cdt[id].next];

  int k, kmax = std::min(n1->ndisc, MAX_ANGDISTLEVELS);
  for(k=0 ; k<kmax ; k++){
    if(td[k].lmax <= 0) break;

    for(int i=0 ; i<MAX_ANGDIST ; i++){

      double x = 0.0;
      for(int l=0 ; l<jmax ; l+=2){
        x += crx.legcoef[id][k][l] * legendre(l,crx.costh[i]);
      }

      crx.angdist[k][i] = x;
    }
  }

  return(k);
}

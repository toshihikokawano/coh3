/******************************************************************************/
/*  knockout.cpp                                                              */
/*        alpha particle knockout reaction                                    */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "statmodel.h"
#include "nucleus.h"
#include "parameter.h"


/**********************************************************/
/*      Alpha Knockout in Pre-Equilibrium Process         */
/*      model and parameters taken from C. Kalbach Walker */
/*      Users Manual for PRECO-2006 (2001)                */
/**********************************************************/
void preqAlphaKnockout(System *sys, Pdata *pdt, Transmission **tc, Transmission **td, Spectra *spc)
{
  if((sys->incident.pid != neutron) && (sys->incident.pid != proton)) return;
  if(!ncl[0].cdt[alpha].status) return;


/************************************************/
/*      Parameter Setting                       */
/************************************************/

  /*** Index of nucleus for alpha emission */
  int idin  = sys->inc_id;
  int idout = (int)alpha;
  int idres = ncl[0].cdt[alpha].next;

  /*** Default single particle state density */
  double gz = sys->target.getZ()/13.0;
  double gn = sys->target.getN()/13.0;
  double ga = (sys->incident.pid == neutron) ? gn : gz;
  double gb = sys->target.getA()/208.0;

  /*** Coulomb barrier */
  double a13 = pow(sys->compound.getA(),-1.0/3.0);
  double bca = sys->compound.getZ() * pdt[idin ].za.getZ() * a13;
  double bcb = sys->compound.getZ() * pdt[idout].za.getZ() * a13;

  /*** Average inverse cross section */
  double siga = 0.0;
  for(int k=1 ; k<ncl[idin ].ntotal ; k++){
    if(tc[idin ][k].ecms>bca) siga += tc[idin ][k].sigr * ncl[0].de;
  }
  double sigb = 0.0;
  for(int k=1 ; k<ncl[idres].ntotal ; k++){
    if(tc[idout][k].ecms>bcb) sigb += tc[idout][k].sigr * ncl[0].de;
  }

  /*** Probability for exciting p-h pair */
  double p;
  int    n = sys->target.getN();;
  if( (116 < n) && (n < 126) )       p = 0.02 + 0.06*(126-n)/10.0;
  else if( (126 <= n) && (n < 129) ) p = 0.02 + 0.06*(n-126)/ 3.0;
  else                               p = 0.08;
  double pb = p*sys->target.getZ()*0.5 / (sys->target.getA() - 1.5*p*sys->target.getZ());

  /*** Pauli correction */
  double ak = 1.0/(2.0*ga*ga) - 1.0/(2.0*gb*gb);

  /*** Clear spc array */
  spc->memclear("pe");

/************************************************/
/*      For Alpha Emission                      */
/************************************************/

  double c1,c2,c3;

  c1 = crx.reaction * (pdt[idout].spin2+1.0) /12.0 * pdt[idout].za.getA();
  c2 = pb * ga * gb;
  c3 = (pdt[idin ].spin2+1.0) * pdt[idin ].za.getA() * siga * gb *gb
      *(sys->ex_total+2.0*bca)*(sys->ex_total-bca)*(sys->ex_total-bca)/6.0
      + (pdt[idout].spin2+1.0) * pdt[idout].za.getA() * sigb * ga *gb
      *(sys->ex_total+2.0*bcb)*(sys->ex_total-bcb)*(sys->ex_total-bcb)/6.0;

  if(c3 == 0.0) return;

  c1 *= parmGetFactor(parmKO);

  /*** Calculation for continuum region only */
  for(int k=1 ; k<ncl[idout].ncont ; k++){
    spc->pe[idout][k] = c1* c2/c3*(ncl[idres].excitation[k]-ak)*tc[idout][k].ecms*tc[idout][k].sigr;
/*
    printf("%2d  %10.4f %10.4f %11.4e\n",
            k,tc[idout][k].ecms,ncl[idout].excitation[k], spc->pe[idout][k]);
*/
  }

  preqExcitonPopulation(idout,tc[idout],td[idout],spc->pe[idout]);
}



/******************************************************************************/
/*  stattrans.cpp                                                             */
/*        set-up transmission coefficients                                    */
/******************************************************************************/
#include <iostream>
#include <iomanip>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "nucleus.h"
#include "omcalc.h"
#include "levden.h"
#include "masstable.h"
#include "parameter.h"

static void statStoreGammaTransmission (const double, const double, double *, Nucleus *);

static const double ecms_cutoff = 1e-9; // lowest emission energy 1 meV

#undef DEBUG_TRANS


/**********************************************************/
/*      Store Particle Continuum Transmission into Array  */
/**********************************************************/
void statStoreContinuumTransmission(const int ip, const double ex, Pdata *p, Transmission **tc)
{
  CrossSection cx;

  /*** Clear transmission array */
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[j][k].memclear();
  }

  /*** For all particle decay channels */
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    if (!ncl[ip].cdt[j].status) continue;

    double sp = ncl[ip].cdt[j].binding_energy;
    int    id = ncl[ip].cdt[j].next;
    double mu = ncl[id].mass * p[j].mass / (ncl[id].mass + p[j].mass);
    double f  = parmGetFactor(parmTJ,ncl[ip].za,(Particle)j);

    for(int k=0 ; k<ncl[id].ntotal ; k++){
      tc[j][k].ecms = ex - ncl[id].excitation[k] - sp;

      /*** For the lowest transition, represented by 1/4 mesh point */
      if(k==0){
        tc[j][k].ecms += (ncl[id].excitation[0] - ncl[id].excitation[1])/4.0;
      }
      if(tc[j][k].ecms < ecms_cutoff) continue;
      
      tc[j][k].lmax = omCalc(tc[j][k].ecms, &p[j], &ncl[id].za, mu, tc[j][k].tran, &cx);
      tc[j][k].sigr = cx.reaction * f;

      for(int l=0 ; l<=3*tc[j][k].lmax ; l++) tc[j][k].tran[l] *= f;
#ifdef DEBUG_TRANS
      std::cout << "CONT " << std::setw(3) << j << std::setw(3) << k << std::setw(3) << id;
      std::cout << std::setw(5) << tc[j][k].lmax;
      std::cout.setf(std::ios::fixed, std::ios::floatfield);
      std::cout << std::setprecision(5) << std::setw(10) << tc[j][k].ecms;
      std::cout << std::setprecision(5) << std::setw(10) << ncl[id].excitation[k];
      std::cout.setf(std::ios::scientific, std::ios::floatfield);
      std::cout << std::setprecision(4) << std::setw(11) << tc[j][k].tran[0];
      std::cout << std::setprecision(4) << std::setw(11) << tc[j][k].tran[3];
      std::cout << std::setprecision(4) << std::setw(11) << tc[j][k].tran[4];
      std::cout << std::setprecision(4) << std::setw(11) << tc[j][k].tran[7] << std::endl;
#endif
    }

    /*** Account for quater-size mesh for k=0 */
    for(int l=0 ; l<=3*tc[j][0].lmax ; l++) tc[j][0].tran[l] *= 0.25;
  }
}


/**********************************************************/
/*      Store Particle Discrete Transmission into Array   */
/**********************************************************/
void statStoreDiscreteTransmission(const int ip, const double ex, Pdata *p, Transmission **td)
{
  CrossSection cx;

  /*** Clear transmission array */
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    for(int k=0 ; k<MAX_LEVELS ; k++) td[j][k].memclear();
  }

  /*** For all particle decay channels */
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    if (!ncl[ip].cdt[j].status) continue;

    double sp = ncl[ip].cdt[j].binding_energy;
    int    id = ncl[ip].cdt[j].next;
    double mu = ncl[id].mass * p[j].mass / (ncl[id].mass + p[j].mass);
    double f  = parmGetFactor(parmTJ,ncl[ip].za,(Particle)j);

    /*** Transition to discrete levels */
    for(int k=0 ; k<ncl[id].ndisc ; k++){
      td[j][k].ecms = ex - ncl[id].lev[k].energy - sp;
      if(td[j][k].ecms < ecms_cutoff) continue;

      td[j][k].lmax = omCalc(td[j][k].ecms, &p[j], &ncl[id].za, mu, td[j][k].tran, &cx);
      td[j][k].sigr = cx.reaction * f;

      for(int l=0 ; l<=3*td[j][k].lmax ; l++) td[j][k].tran[l] *= f;
#ifdef DEBUG_TRANS
      std::cout << "DISC " << std::setw(3) << j << std::setw(3) << k << std::setw(3) << id;
      std::cout << std::setw(5) << td[j][k].lmax;
      std::cout.setf(std::ios::fixed, std::ios::floatfield);
      std::cout << std::setprecision(5) << std::setw(10) << td[j][k].ecms;
      std::cout << std::setprecision(5) << std::setw(10) << ncl[id].lev[k].energy;
      std::cout.setf(std::ios::scientific, std::ios::floatfield);
      std::cout << std::setprecision(4) << std::setw(11) << td[j][k].tran[0];
      std::cout << std::setprecision(4) << std::setw(11) << td[j][k].tran[3];
      std::cout << std::setprecision(4) << std::setw(11) << td[j][k].tran[6];
      std::cout << std::setprecision(4) << std::setw(11) << td[j][k].tran[7] << std::endl;
#endif
    }
  }
}


/**********************************************************/
/*  Store Gamma-ray Transmissions into 2-dim Array        */
/**********************************************************/
void statStoreContinuumGammaTransmission(const int k0, double **tg, Nucleus *n)
{
  /*** once copy GDR parameters into the gtrans.cpp scope */
  if(k0 == 0) gdrParameterSave(n->gdr,n->za.getA());

  /*** Clear transmission array */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    for(int i=0 ; i<MAX_MULTIPOL ; i++) tg[k][i] = 0.0;
  }

  /*** Neutron separation energy for temperature dependent Gamma */
  double sn = 0.0;
  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if( n->cdt[id].pid == neutron ){ sn = n->cdt[id].binding_energy; break; }
  }
  if(sn == 0.0){
    bool checkmass = true;
    double mx = mass_excess(n->za.getZ(),n->za.getA()-1,&checkmass);
    sn = (checkmass) ? mx - n->mass_excess + ENEUTRON : 0.0;
  }

  for(int k1=k0+1 ; k1<n->ntotal ; k1++){
    double eg = n->excitation[k0] - n->excitation[k1];
    statStoreGammaTransmission(eg,sn,tg[k1],n);
  }
}


/**********************************************************/
/*  Store Gamma-ray Transmissions into 1-dim Array        */
/**********************************************************/
void statStoreDiscreteGammaTransmission(const double eg, double *tg, Nucleus *n)
{
  /*** Clear transmission array */
  for(int i=0 ; i<MAX_MULTIPOL ; i++) tg[i] = 0.0;

  double sn = 0.0;
  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if( n->cdt[id].pid == neutron ){ sn = n->cdt[id].binding_energy; break; }
  }
  if(sn == 0.0){
    bool checkmass = true;
    double mx = mass_excess(n->za.getZ(),n->za.getA()-1,&checkmass);
    sn = (checkmass) ? mx - n->mass_excess + ENEUTRON : 0.0;
  }
  statStoreGammaTransmission(eg,sn,tg,n);
}


/**********************************************************/
/*  Calculate Gamma-ray Transmissions                     */
/**********************************************************/
void statStoreGammaTransmission(const double eg, const double sn, double *tg, Nucleus *n)
{
  double ex = (eg <= sn) ? sn - eg : 0.0;
  double a  = ldDensityParameter(ex,(double)n->za.getA(),&n->ldp);

  for(int m=0 ; m<MAX_MULTIPOL ; m++){
    if(m == E1) tg[m] = gdrGammaTransmission((GammaMultipolarity)m,eg,a  ,ex );
    else        tg[m] = gdrGammaTransmission((GammaMultipolarity)m,eg,0.0,0.0);
  }

#ifdef DEBUG_TRANS
  std::cout << "GAMMA ";
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout << std::setprecision(5) << std::setw(10) << eg;
  std::cout << std::setprecision(5) << std::setw(10) << ex;
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4) << std::setw(11) << tg[E1]/(PI2*pow(eg,3.0));
  std::cout << std::setprecision(4) << std::setw(11) << tg[M1]/(PI2*pow(eg,3.0));
  std::cout << std::setprecision(4) << std::setw(11) << tg[E2]/(PI2*pow(eg,5.0)) << std::endl;
#endif

}

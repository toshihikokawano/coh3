/******************************************************************************/
/*  MFTmain.cpp                                                               */
/*        mean-field theory calculation                                       */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstring>

#include "structur.h"
#include "coh.h"
#include "HF.h"
#include "FRDM.h"
#include "global.h"
#include "outformat.h"

static void MFTPrintSectionHead (const char *);
static void MFTPrintNucleus (MFTSystem *);
static void MFTPrintExpansion (HFInterface *);


/**********************************************************/
/*     Hartree-Fock BCS Calculation                       */
/**********************************************************/
int mftCalc(MFTparm *parm, HFInterface *spspec, double *beta)
{
  /*** create harmonic oscillator basis */
  Basis basis;

  basis.setN0(parm->n0);
  basis.b = parm->b;
  basis.q = parm->q;

  /*** system parameters */
  MFTSystem mftsys(parm->target.getZ(), parm->target.getA());
  mftsys.setNuclearSize(gRadius);

  /*** copy input parameters */
  for(int i=0 ; i<6 ; i++) mftsys.setEps(i+1,parm->eps[i]);
  mftsys.setGamma(parm->gamma);

  if(mftsys.getEps(1) < 0.0) FRDMReadShape(&mftsys);

  if(pmf.system){
    MFTPrintNucleus(&mftsys);
    if(ctl.frdm && pmf.deformation) FRDMPrintShape(&mftsys);
  }

  /*** FRDM calculation */
  if(ctl.frdm && !ctl.hartreefock) FRDMmain(&mftsys, &basis, beta, spspec);

  /*** Hartree-Fock BCS calculaion */
  if(ctl.hartreefock){

    Skyrme   skyrme;
    Moment   moment;
    BCSParameter bcs[2];

    /*** copy constraints parameters */
    for(int i=0 ; i<5 ; i++){
      if(parm->value[i] > 0.0 && parm->weight[i] > 0.0){
        moment.setConstraint(i,parm->value[i],parm->weight[i],parm->nfree[i]);
      }
    }

    mftsys.setIteration(parm->max_iteration,parm->eps_converge);

    /*** BCS model parameters */
    if(parm->strength[0] != 0.0) bcs[0].strength = parm->strength[0];
    if(parm->strength[1] != 0.0) bcs[1].strength = parm->strength[1];

    bcs[0].pairing = parm->pairing;
    bcs[1].pairing = parm->pairing;

    /*** setup Skyrme parameters */
    SkyrmeParameter(parm->force, &skyrme);

    /*** solve Hartree-Fock equation */
    HartreeFock(&mftsys, &basis, &skyrme, bcs, &moment, spspec);
  }


  /*** print expansion coefficients */
  if(pmf.expansion)  MFTPrintExpansion(spspec);

  return 0;
}


/**********************************************************/
/*      Write Header of Each Section                      */
/**********************************************************/
void MFTPrintSectionHead(const char *sec)
{
  const std::string s = "[ Mean Field Theory Calculation ]";
  int n1 = s.length();
  int n2 = strlen(sec);

  std::string bar1 = "#";

  for(int i=0 ; i<DisplayWidth-n1-1 ; i++) bar1 += ".";
  std::cout << bar1 << s << std::endl;

  std::string bar2 = "#";
  for(int i=0 ; i<DisplayWidth-n2-1 ; i++) bar2 += " ";
  std::cout << bar2 << sec << std::endl;
}


/**********************************************************/
/*      Nucleus Size Parameters                           */
/**********************************************************/
void MFTPrintNucleus(MFTSystem *s)
{
  MFTPrintSectionHead("TARGET NUCLEUS");
  std::cout << cline << "  AtomicNum    MassNum Radius[fm] Surf[fm2]  Vol[fm3]" << std::endl;
  std::cout << blank;
  outVal(11,s->getZ());
  outVal(11,s->getA());
  outVal(s->R0);
  outVal(s->S0);
  outVal(s->V0);
 nl();
}


/**********************************************************/
/*      Print Spherical HO Expansion Coefficients         */
/**********************************************************/
void MFTPrintExpansion(HFInterface *hfsp)
{
  MFTPrintSectionHead("SPHERICAL HO EXPANSION");

  for(int p=0 ; p<=1 ; p++){

    if(p == 0) std::cout << "# neutron shell" << std::endl;
    else       std::cout << "# proton shell" << std::endl;

    /*** for all deformed states */
    for(int i=0 ; i<hfsp[p].n ; i++){

      /*** expansion convergence check */
      double s = 0.0;
      for(int j=0 ; j<hfsp[p].state[i].m ; j++){
        s += hfsp[p].state[i].q[j].c * hfsp[p].state[i].q[j].c;
      }

      std::cout << "# Ebind      V2         U           Ns   K2 Lmax" << std::endl;
      outVal(11,hfsp[p].state[i].energy);
      outVal(11, hfsp[p].state[i].v2);
      outVal(11, sqrt(1.0-hfsp[p].state[i].v2));
      outVal(5,hfsp[p].state[i].m);
      outVal(5,hfsp[p].state[i].k2);
      outVal(5,hfsp[p].state[i].lmax); nl();

      /*** print spherical state and expansion coefficient */
      for(int j=0 ; j<hfsp[p].state[i].m ; j++){
        outVal( 3,hfsp[p].state[i].q[j].n);
        outVal( 2,hfsp[p].state[i].q[j].l);
        outVal( 3,hfsp[p].state[i].q[j].j);
        outVal(11,hfsp[p].state[i].q[j].c);
        if(j%5 == 4) std::cout << std::endl;
      }
      if((hfsp[p].state[i].m)%5 != 0) std::cout << std::endl;
      std::cout << "#    Sum ";
      outVal(s); nl();
    }
    nl();
  }
}



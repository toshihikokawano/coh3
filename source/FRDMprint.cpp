/******************************************************************************/
/*  FRDMprint.cpp                                                             */
/*        Output FRDM Calculation Results                                     */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cstring>

#include "FRDM.h"
#include "outformat.h"

static void FRDMPrintSectionHead(const char *);

/**********************************************************/
/*      Write Header of Each Section                      */
/**********************************************************/
void FRDMPrintSectionHead(const char *sec)
{
  const std::string s = "[ FRDM Calculation ]";
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
/*      FRDM System Parameters                            */
/**********************************************************/
void FRDMPrintShape(MFTSystem *sys){
  FRDMPrintSectionHead("SHAPE PARAMETERS");

  std::cout << cline << "   epsilon1   epsilon2   epsilon3";
  std::cout <<"   epsilon4   epsilon5   epsilon6      gamma" << std::endl;
  std::cout << " DeformParm";
  outVal(11,4,sys->getEps(1));
  outVal(11,4,sys->getEps(2));
  outVal(11,4,sys->getEps(3));
  outVal(11,4,sys->getEps(4));
  outVal(11,4,sys->getEps(5));
  outVal(11,4,sys->getEps(6));
  outVal(11,4,sys->getGamma()); nl();
}


/**********************************************************/
/*      Print Macroscopic Energy                          */
/**********************************************************/
void FRDMPrintMacroQuantity(const double eps, const double delta, const double w0, double *b, const double bs, const double bv, const double bw, const double br)
{
  FRDMPrintSectionHead("MACRO-PART QUANTITY");

  std::cout << cline << " eps-bar    delta-bar  W0/W0(sph)" <<  std::endl;
  std::cout << blank;
  outVal(eps);
  outVal(delta);
  outVal(w0); nl();

  std::cout << cline << " B1         B2         B3         B4        ";
  std::cout << " Bs         Bv         Bw         Br" << std::endl;
  std::cout << blank;
  for(int i=0 ; i<4 ; i++) outVal(b[i]);
  outVal(bs);
  outVal(bv);
  outVal(bw);
  outVal(br); nl();
}


void FRDMPrintMacroEnergy(double *e)
{
  FRDMPrintSectionHead("MACROSCOPIC ENERGIES");

  int i = 0;
  std::cout << " Mass Exess                      ";  outVal(e[i++]); nl();
  std::cout << " Volume Energy                   ";  outVal(e[i++]); nl();
  std::cout << " Surface Energy                  ";  outVal(e[i++]); nl();
  //  std::cout << " Curvature Energy                ";  outVal(e[i++]); nl();
  //  std::cout << " A0 Energy                       ";  outVal(e[i++]); nl();
  i += 2;
  std::cout << " Coulomb Energy                  ";  outVal(e[i++]); nl();
  std::cout << " Volume Redistribution Energy    ";  outVal(e[i++]); nl();
  std::cout << " Coulomb Exchange Energy         ";  outVal(e[i++]); nl();
  std::cout << " Surface Redistribution Energy   ";  outVal(e[i++]); nl();
  std::cout << " Proton Formfactor Energy        ";  outVal(e[i++]); nl();
  std::cout << " Charge Asymmetry Energy         ";  outVal(e[i++]); nl();
  std::cout << " Wigner Term Energy              ";  outVal(e[i++]); nl();
  std::cout << " Pairing Energy                  ";  outVal(e[i++]); nl();
  std::cout << " Bound Electrons                 ";  outVal(e[i++]); nl();
}


/**********************************************************/
/*      Print Spherical Harmonics Expansion               */
/**********************************************************/
void FRDMPrintBeta(const int lmax, double **alm, const double R0)
{
  const double eps = 1e-10;

  FRDMPrintSectionHead("SHAPE BETA EXPANSION");

  std::cout << "#   L   M    a(L,M)      beta(L)   " << std::endl;

  for(int l=0 ; l<=lmax ; l++){
    for(int m=-l ; m<=l ; m++){
      double beta = alm[l][m] / R0;
      if(fabs(beta) > eps){
        std::cout << std::setw(5) << l << std::setw(4) << m;
        std::cout << " "; outVal(alm[l][m]);
        std::cout << " "; outVal(beta);  nl();
      }
    }
  }
}


/**********************************************************/
/*      Print Generated FRDM Potential                    */
/**********************************************************/
void FRDMPrintPotential(Basis *basis, OneBodyPotential *pot, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4);

  FRDMPrintSectionHead("POTENTIAL SHAPE");

  std::cout << "# R[fm]       Z[fm]     ";
  std::cout << "  V0(n)       Vso(n)      V0(p)       Vso(p)"<<std::endl;
  for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;
    for(int i=0 ; i<fl->nr ; i++){
      std::cout << std::setw(12) << sqrt(fl->x[i]/basis->bp2);
      std::cout << std::setw(12) << fh->x[j]/basis->bz;
      std::cout << std::setw(12) << pot[0].central[i][j];
      std::cout << std::setw(12) << pot[0].spinorb[i][j];
      std::cout << std::setw(12) << pot[1].central[i][j];
      std::cout << std::setw(12) << pot[1].spinorb[i][j] << std::endl;
    }
    std::cout << std::endl;
  }
}


/**********************************************************/
/*      Print Single-Particle States                      */
/**********************************************************/
void FRDMPrintSPState(const bool bcsflag, Basis *b, SPEnergy *s, LNParameter *bcs)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4);

  if(bcsflag) FRDMPrintSectionHead("SHIFTED SINGLE PARTICLE SPECTRA");
  else        FRDMPrintSectionHead("SINGLE PARTICLE SPECTRA");

  const double ecut = 10.0;
  const double vcut = 1.0e-10;

  std::cout << "#     neutron shell                     proton shell" << std::endl;
  std::cout << "#     Omg2    Energy[MeV]  v2          Omg2    Energy[MeV]  v2" << std::endl;
  for(int k=0 ; k<b->nv ; k++){

    double en = s[0].energy[k];
    double ez = s[1].energy[k];

    if(bcsflag){
      if( (bcs[0].n1 <= k) && (k<=bcs[0].n2) ){
        en += (4.0*bcs[0].lambda2 - bcs[0].strength)*s[0].v2[k];
      }
      if( (bcs[1].n1 <= k) && (k<=bcs[1].n2) ){
        ez += (4.0*bcs[1].lambda2 - bcs[1].strength)*s[1].v2[k];
      }
    }

    std::cout << std::setw(4) << k;
    std::cout << std::setw(6) << b->omega[s[0].block[k]];
    std::cout << std::setw(3) << (s[0].parity[k] > 0.0 ? " + " : " - ");
    std::cout << std::setw(12) << en;
    std::cout << std::setw(12) << s[0].v2[k];
    std::cout << std::setw(6) << b->omega[s[1].block[k]];
    std::cout << std::setw(3) << (s[1].parity[k] > 0.0 ? " + " : " - ");
    std::cout << std::setw(12) << ez;
    std::cout << std::setw(12) << s[1].v2[k] << std::endl;

    if( (s[0].energy[k] > ecut) && (s[1].energy[k] > ecut) ) break;
    if( (s[0].v2[k] < vcut) && (s[1].v2[k] < vcut) ) break;
  }
}


/**********************************************************/
/*      Print BCS Parameters                              */
/**********************************************************/
void FRDMPrintBCS(LNParameter *bcs)
{
  FRDMPrintSectionHead("BCS PARAMETERS");

  std::cout << "#             Strength     Lambda    Lambda2";
  std::cout << "      Delta    Delta-G  DeltaMacr" << std::endl;
  std::cout << " Neutron   ";
  outVal(11,4,bcs[0].strength);
  outVal(11,4,bcs[0].chemical_potential);
  outVal(11,4,bcs[0].lambda2);
  outVal(11,4,bcs[0].pairing_gap);
  outVal(11,4,bcs[0].pairing_eff);
  outVal(11,4,bcs[0].pairing_mac);
  nl();
  std::cout << " Proton    ";
  outVal(11,4,bcs[1].strength);
  outVal(11,4,bcs[1].chemical_potential);
  outVal(11,4,bcs[1].lambda2);
  outVal(11,4,bcs[1].pairing_gap);
  outVal(11,4,bcs[1].pairing_eff);
  outVal(11,4,bcs[1].pairing_mac);
  nl();
}


/**********************************************************/
/*      Print Shell Correction Parameters                 */
/**********************************************************/
void FRDMPrintSC(MFTSystem *sys, SCParameter *scp)
{
  FRDMPrintSectionHead("SHELL CORRECTION PARAMETERS");

  std::cout << "#               Lambda    Rho-bar sp-density   nucl./13" << std::endl;
  std::cout << " Neutron   ";
  outVal(11,4,scp[0].chemical_potential);
  outVal(11,4,scp[0].rho);
  outVal(11,4,scp[0].rho*2.0);
  outVal(11,4,sys->getN()/13.0);
  nl();
  std::cout << " Proton    ";
  outVal(11,4,scp[1].chemical_potential);
  outVal(11,4,scp[1].rho);
  outVal(11,4,scp[1].rho*2.0);
  outVal(11,4,sys->getZ()/13.0);
  nl();
}


/**********************************************************/
/*      Print Microscopic Energy                          */
/**********************************************************/
void FRDMPrintMicroEnergy(LNParameter *bcs, SCParameter *scp)
{
  FRDMPrintSectionHead("MICROSCOPIC ENERGIES");

  std::cout << "#              Pairing    Average    Ep-<Ep>" << std::endl;
  std::cout << " Neutron   ";
  outVal(11,4,bcs[0].energy_micro);
  outVal(11,4,bcs[0].energy_average);
  outVal(11,4,bcs[0].energy_micro - bcs[0].energy_average); nl();

  std::cout << " Proton    ";
  outVal(11,4,bcs[1].energy_micro);
  outVal(11,4,bcs[1].energy_average);
  outVal(11,4,bcs[1].energy_micro - bcs[1].energy_average); nl();
  nl();

  std::cout << "#             ShellCor    Average    Es-<Es>" << std::endl;
  std::cout << " Neutron   ";
  outVal(11,4,scp[0].energy_micro);
  outVal(11,4,scp[0].energy_average);
  outVal(11,4,scp[0].energy_micro - scp[0].energy_average); nl();

  std::cout << " Proton    ";
  outVal(11,4,scp[1].energy_micro);
  outVal(11,4,scp[1].energy_average);
  outVal(11,4,scp[1].energy_micro - scp[1].energy_average); nl();
}


/**********************************************************/
/*      Print Total Micro / Macro Energies                */
/**********************************************************/
void FRDMPrintTotalEnergy(FRDMEnergy *e)
{
  FRDMPrintSectionHead("TOTAL ENERGY");

  std::cout << "                                    deformed    spherical"; nl();
  std::cout << " macroscopic energy              ";  outVal(e->macro);
  std::cout << "  ";                                 outVal(e->spherical); nl();
  std::cout << " microscopic pair+shell energies ";  outVal(e->micro); nl();
  std::cout << " zero-point energy               ";  outVal(e->zero); nl();
  std::cout << "                                 -----------";  nl();
  std::cout << " TOTAL                           ";  outVal(e->total); nl();
}


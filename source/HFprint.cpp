/******************************************************************************/
/*  HFprint.cpp                                                               */
/*        output calculated results                                           */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstring>

#include "HF.h"
#include "outformat.h"

static void   HFPrintSectionHead (const char *);
static double PRNPointDensity (const double, const double, Basis *, double *);
static double PRNHermite (const int, const  double);
static double PRNLaguerre (const int, const int, const double);


/**********************************************************/
/*      Write Header of Each Section                      */
/**********************************************************/
void HFPrintSectionHead(const char *sec)
{
  const std::string s = "[ Hartree Fock BCS Calculation ]";
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
/*      Print Internal Values                             */
/**********************************************************/
void HFMonitor(const int i, const double d, HFEnergy *e, Moment *m)
{
  std::cout << "#----" << std::endl;
  std::cout << "#" << std::setw(4) << i;
  std::cout << "  zcm fm      Q20 fm^2    Q30 b3/2  ";
  std::cout << "  Q40 b2      QN          Zmin      " << std::endl;
  std::cout << "#    ";
  outVal(12,m->zcm);
  outVal(12,m->q20);
  outVal(12,m->q30 * 1.0e-3);
  outVal(12,m->q40 * 1.0e-4);
  outVal(12,m->qn);
  outVal(12,m->zmin);
  std::cout << std::endl;

  std::cout << "#    ";
  std::cout << "  Kinetic     Volume      Surface   ";
  std::cout << "  SpinOrbit   Coulomb     Total MeV   epsilon" << std::endl;
  std::cout << "#    ";
  outVal(12, e->kinetic);
  outVal(12, e->volume);
  outVal(12, e->surface);
  outVal(12, e->spinorbit);
  outVal(12, e->coulomb);
  outVal(12, e->total);
  outVal(12, d);
  std::cout << std::endl;
}


/**********************************************************/
/*      Print BCS Parameters                              */
/**********************************************************/
void HFPrintBCS(BCSParameter *b)
{
  HFPrintSectionHead("BCS PARAMETERS");

  std::cout << "#               Lambda     Delta      Gap" << std::endl;
  std::cout << " Neutron   ";
  outVal(11,4,b[0].chemical_potential);
  outVal(11,4,b[0].pairing_energy);
  outVal(11,4,b[0].average_gap);
  nl();
  std::cout << " Proton    ";
  outVal(11,4,b[1].chemical_potential);
  outVal(11,4,b[1].pairing_energy);
  outVal(11,4,b[1].average_gap);
  nl();
}


/**********************************************************/
/*      Print Calculated Energies                         */
/**********************************************************/
void HFPrintEnergy(HFEnergy *e)
{
  HFPrintSectionHead("TOTAL ENERGY");

  std::cout << " Nucleon Masses                  ";  outVal(e->nucleon);   nl();
  std::cout << " Kinetic Energy                  ";  outVal(e->kinetic);   nl();
  std::cout << " Volume Energy                   ";  outVal(e->volume);    nl();
  std::cout << " Surface Energy                  ";  outVal(e->surface);   nl();
  std::cout << " Spin-Orbit Energy               ";  outVal(e->spinorbit); nl();
  std::cout << " Coulomb Energy                  ";  outVal(e->coulomb);   nl();
  std::cout << "                                 -----------"; nl();
  std::cout << " TOTAL                           ";  outVal(e->total);  nl();
}


/**********************************************************/
/*      Print Calculated Moments                          */
/**********************************************************/
void HFPrintMoment(Moment *m)
{
  HFPrintSectionHead("NUCLEAR MOMENTS");

  std::cout << cline;
  std::cout << "  zcm[fm]     Q20[fm^2]   Q30[b3/2] ";
  std::cout << "  Q40[b2]     QN          Zmin[fm]" << std::endl;
  std::cout << blank;
  outVal(12,m->zcm);
  outVal(12,m->q20);
  outVal(12,m->q30 * 1.0e-3);
  outVal(12,m->q40 * 1.0e-4);
  outVal(12,m->qn);
  outVal(12,m->zmin);
  std::cout << std::endl;
}


/**********************************************************/
/*      Print Calculated Potential                        */
/**********************************************************/
void HFPrintPotential(Basis *basis, OneBodyPotential *pot, LaguerrePolynomial *fl, HermitePolynomial *fh)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4);

  HFPrintSectionHead("POTENTIAL SHAPE");

  std::cout << "# R[fm]       Z[fm]     ";
  std::cout << "  V0(n)       EffMass(n)  Vso(n)      V0(p)       EffMass(p)  Vso(p)"<<std::endl;
  for(int j=-fh->nz ; j<=fh->nz ; j++){ if(j==0) continue;
    for(int i=0 ; i<fl->nr ; i++){
      std::cout << std::setw(12) << sqrt(fl->x[i]/basis->bp2);
      std::cout << std::setw(12) << fh->x[j]/basis->bz;
      std::cout << std::setw(12) << pot[0].central[i][j];
      std::cout << std::setw(12) << pot[0].effmass[i][j];
      std::cout << std::setw(12) << pot[0].spinorb[i][j];
      std::cout << std::setw(12) << pot[1].central[i][j];
      std::cout << std::setw(12) << pot[1].effmass[i][j];
      std::cout << std::setw(12) << pot[1].spinorb[i][j] << std::endl;
    }
    std::cout << std::endl;
  }
}


/**********************************************************/
/*      Print Single-Particle States                      */
/**********************************************************/
void HFPrintSPState(Basis *b, SPEnergy *s)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4);

  HFPrintSectionHead("SINGLE PARTICLE SPECTRA");

  const double ecut = 10.0;
  const double vcut = 1.0e-10;

  std::cout << "#     neutron shell                 proton shell" << std::endl;
  std::cout << "#     Omg2 Energy[MeV]  v2          Omg2 Energy[MeV]  v2" << std::endl;
  for(int k=0 ; k<b->nv ; k++){
    std::cout << std::setw(4) << k;
    std::cout << std::setw(6) << b->omega[s[0].block[k]];
    std::cout << std::setw(12) << s[0].energy[k];
    std::cout << std::setw(12) << s[0].v2[k];
    std::cout << std::setw(6) << b->omega[s[1].block[k]];
    std::cout << std::setw(12) << s[1].energy[k];
    std::cout << std::setw(12) << s[1].v2[k] << std::endl;

    if( (s[0].energy[k] > ecut) && (s[1].energy[k] > ecut) ) break;
    if( (s[0].v2[k] < vcut) && (s[1].v2[k] < vcut) ) break;
  }
}


/**********************************************************/
/*      Print Nuclear Density                             */
/**********************************************************/
void HFPrintDensityCalc(const int nr, const int nz, const double dr, const double dz, Basis *basis, double *matV, GridData *rho)
{
  for(int i=0 ; i<nr ; i++){
    std::cerr << std::setw(6) << i*dr;
    for(int j=-nz/2 ; j<nz/2 ; j++){
      if(j == 0) continue;
      rho->p[i][j] = PRNPointDensity(dr * i, dz * j, basis, matV);
    }
  }
  std::cerr << std::endl;
}


void HFPrintDensity(int nr, int nz, double dr, double dz, GridData *rhon, GridData *rhop)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4);

  HFPrintSectionHead("DENSITY");

  /*** at Z=0, take average of both adjucent mesh */
  for(int i=0 ; i<nr ; i++){
    rhon->p[i][0] = 0.5*(rhon->p[i][-1] + rhon->p[i][1]);
    rhop->p[i][0] = 0.5*(rhop->p[i][-1] + rhop->p[i][1]);
  }

  for(int i=-nr+1 ; i<nr ; i++){
    int k = std::abs(i);
    for(int j=-nz/2 ; j<nz/2 ; j++){
      std::cout << std::setw(12) << i*dr << std::setw(12) << j*dz;
      std::cout << std::setw(12) << rhon->p[k][j];
      std::cout << std::setw(12) << rhop->p[k][j];
      std::cout << std::setw(12) << rhon->p[k][j] + rhop->p[k][j] << std::endl;
    }
    std::cout << std::endl;
  }
}


/**********************************************************/
/*      rho at Given (R,Z)                                */
/**********************************************************/
double PRNPointDensity(const double r, const double z, Basis *basis, double *matV)
{
  const double rhocut = 1.0e-07;

  double zeta = basis->bz * z;
  double eta  = basis->bp2 * r*r;

  double sum = 0.0;

  for(int kb=0 ; kb<basis->nb ; kb++){
    for(int k1=basis->index[kb] ; k1<basis->index[kb]+basis->ndata[kb] ; k1++){

      int nr1 = basis->state[k1].getNr();
      int nz1 = basis->state[k1].getNz();
      int l1  = basis->state[k1].getL();
      int s1  = basis->state[k1].getS();

      for(int k2=basis->index[kb] ; k2<basis->index[kb]+basis->ndata[kb] ; k2++){

        int nr2 = basis->state[k2].getNr();
        int nz2 = basis->state[k2].getNz();
        int l2  = basis->state[k2].getL();
        int s2  = basis->state[k2].getS();

        if(s1 != s2) continue;
        if(l1 != l2) continue;

        int k0 = (k2<=k1) ? k1*(k1+1)/2 + k2 : k2*(k2+1)/2 + k1;
        double rho12 = matV[k0];
        if( fabs(rho12) < rhocut ) continue;

        double yh = PRNHermite(nz1,zeta) * PRNHermite(nz2,zeta)
          / sqrt(PI * pow(2.0, (double)(nz1+nz2)) * exp(fact[nz1] + fact[nz2]));
        double yl = PRNLaguerre(nr1,l1,eta) * PRNLaguerre(nr2,l2,eta)
          * sqrt(exp(fact[nr1] + fact[nr2] - fact[nr1 + l1] - fact[nr2 + l2]));
 
        sum += rho12 * yh * yl * pow(eta, (double)l1);
      }
    }
  }

  double d = sum * exp(-zeta*zeta - eta) * basis->bz * basis->bp2/PI;

  return(d);
}


double PRNHermite(const int n, const double x)
{
  double fh = 0.0;

  if(n == 0) fh = 1.0;
  else{
    int p = n%2;

    if((x == 0.0) && (p == 0)){
      fh = ((n/2)%2 == 0 ? 1.0 : -1.0) * exp( fact[n] - fact[n/2] );
    }
    else{
      double s = 0.0;
      for(int m=0 ; m<=(n-p)/2 ; m++){
        s += (m%2 == 0 ? 1.0 : -1.0) * pow(2.0*x,(n-2.0*m)) / exp(fact[m] + fact[n-2*m]);
      }
      fh = s * exp(fact[n]);
    }
  }
  return(fh);
}


double PRNLaguerre(const int n, const int k, const double x)
{
  double fl = 0.0;

  if(n == 0) fl = 1.0;
  else{
    double s = 0.0;
    if(x == 0.0){
      s = 1.0 / exp(fact[n] + fact[k]);
    }
    else{
      for(int m=0 ; m<=n ; m++){
        s += pow(-x,(double)m) / exp(fact[m] + fact[n-m] + fact[k+m]);
      }
    }
    fl = s * exp(fact[n+k]);
  }

  return(fl);
}

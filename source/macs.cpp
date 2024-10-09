/******************************************************************************/
/*  macs.cpp                                                                  */
/*        Maxwellian average cross sections                                   */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "coh.h"
#include "nucleus.h"
#include "statmodel.h"
#include "levden.h"
#include "macs.h"


static void   macsPartitionFunction (int, double *, double *);
static double macsAverage (const int, double *, double *, double);
static inline double lowfilter(double x)
{ if(fabs(x) < 1.0e-99) return 0.0; else return(x); }

/**********************************************************/
/*      Store Calculated Cross Sections                   */
/**********************************************************/
void macsStoreCrossSection(const int ie, const int nc)
{
  std::string str;

  for(int i=0 ; i<N_REACTION_RATE ; i++) crx.macs[i][ie] = 0.0;

  for(int i=0 ; i<nc ; i++) {
    std::ostringstream os;
    for(int j=0 ; j<MAX_CHANNEL ; j++) os << std::setw(1) << (int)crx.prod[i].par[j];
    str = os.str();

    if(     str=="0000000"){ crx.macs[0][ie] = crx.prod[i].xsec; } // capture
    else if(str=="0100000"){ double z =crx.prod[i].xsec - crx.compound;
                             crx.macs[1][ie] = (z<1.0e-10) ? 0.0 : z;} // (n,n')
    else if(str=="0010000"){ crx.macs[2][ie] = crx.prod[i].xsec; } // (n,p)
    else if(str=="0001000"){ crx.macs[3][ie] = crx.prod[i].xsec; } // (n,alpha)
    else if(str=="0200000"){ crx.macs[4][ie] = crx.prod[i].xsec; } // (n,2n)
    else if(str=="0300000"){ crx.macs[5][ie] = crx.prod[i].xsec; } // (n,3n)
    else if(str=="0400000"){ crx.macs[6][ie] = crx.prod[i].xsec; } // (n,4n)
    else if(str=="0110000"){ crx.macs[7][ie] = crx.prod[i].xsec; } // (n,np)
    else if(str=="0101000"){ crx.macs[8][ie] = crx.prod[i].xsec; } // (n,na)

    crx.macs[9][ie] += crx.prod[i].fiss;        // total fission
  }

  crx.macs[10][ie] = crx.total;                  // total
  crx.macs[11][ie] = crx.elastic + crx.compound; // elastic
}


/**********************************************************/
/*      Tabulate MACS at Given Temperatures               */
/**********************************************************/
void macsReactionRate(const int tid, const Particle pin, double *t)
{
  double *g;
  g = new double [MAX_MACSTEMP];

  macsPartitionFunction(tid,t,g);


  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(5);
  std::cout << "# Cross Sections" << std::endl;

  if(pin == gammaray){
    std::cout << "# Ecms[MeV]  (g,g')[mb]   (g,n)[mb] ";
    std::cout << "  (g,p)[mb]   (g,a)[mb]  (g,2n)[mb]  (g,np)[mb]  (g,na)[mb]   (g,f)[mb] " << std::endl;

    for(int i=0 ; i<N_MACS_GGRID ; i++){
      std::cout << std::setw(12) << macs_ggrid[i];
      for(int j=0 ; j<N_REACTION_RATE ; j++) std::cout << std::setw(12) << lowfilter(crx.macs[j][i]);
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  else{
    std::cout << "# Ecms[MeV]    (n,g)[mb]  (n,n')[mb]";
    std::cout << "   (n,p)[mb]   (n,a)[mb]  (n,2n)[mb]  (n,3n)[mb]  (n,4n)[mb]  (n,np)[mb]  (n,na)[mb]   (n,f)[mb]";
    std::cout << "   total[mb] elastic[mb]" << std::endl;

    for(int i=0 ; i<N_MACS_EGRID ; i++){
      std::cout << std::setw(12) << macs_egrid[i];
      for(int j=0 ; j<N_REACTION_RATE ; j++) std::cout << std::setw(12) << lowfilter(crx.macs[j][i]);
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }


  std::cout << "# Maxwellian Average Cross Sections" << std::endl;

  if(pin == gammaray){
    std::cout << "#   E[MeV]       T9[K]     PartFunc   (g,g')[mb]   (g,n)[mb]";
    std::cout << "   (g,p)[mb]   (g,a)[mb]  (g,2n)[mb]  (g,np)[mb]  (g,na)[mb]   (g,f)[mb]" << std::endl;
  }
  else{
    std::cout << "#   E[MeV]       T9[K]     PartFunc    (n,g)[mb]  (n,n')[mb]";
    std::cout << "   (n,p)[mb]   (n,a)[mb]  (n,2n)[mb]  (n,3n)[mb]  (n,4n)[mb]  (n,np)[mb]  (n,na)[mb]   (n,f)[mb]";
    std::cout << "   total[mb] elastic[mb]" << std::endl;
  }

  for(int i=0 ; i<MAX_MACSTEMP ; i++){
    if(t[i] == 0.0) break;

    std::cout << std::setw(12) << t[i] << std::setw(12) << t[i] * 1e+6 / BOLTZMANN * 1e-9;
    std::cout << std::setw(12) << g[i];

    for(int j=0 ; j<N_REACTION_RATE ; j++){
      double rate = macsAverage(N_MACS_EGRID,macs_egrid,crx.macs[j],t[i]);
      std::cout << std::setw(12) << lowfilter(rate);
    }
    std::cout << std::endl;
  }

  delete [] g;
}


/**********************************************************/
/*      Partition Function                                */
/**********************************************************/
// M: Normalized (to the gs spin) partition function otherwise set gspin to unity.
// References: (1) Fowler, Caughlan, Zimmermann Annual Reviews 1967 Eq. 11 (not normalized)
//             (2) Bahcall, Fowler ApJ 1970 Eq. 9 (not normalized)
//             (3) Rauscher & Thielemann At. Data Nuc. Data Table 75, 1 (2000)
//                 Eqs. 12-15 (normalized)

void macsPartitionFunction(int tid, double *t, double *g)
{
  const double emax = 10.0;   // highest excitation energy for integration
  const double gcut = 1000.0; // ad hoc maximum value for the partition function
  const double eps  = 1.0e-5; // cut off for spin distribution

  for(int k=0 ; k<MAX_MACSTEMP ; k++) g[k] = 0.0;

  double de = 0.1;
  double elmax = ncl[tid].lev[ncl[tid].ndisc-1].energy;
  double gspin = 2.0*ncl[tid].lev[0].spin + 1.0;

  for(int n=0 ; n<MAX_MACSTEMP ; n++){
    if(t[n] == 0.0) break;

    for(int i=0 ; i<ncl[tid].ndisc ; i++){
      double ex = ncl[tid].lev[i].energy;
      g[n] += (2.0*ncl[tid].lev[i].spin + 1.0)/gspin * exp(-ex/t[n]);
    }

    for(int k=0 ; ; k++){
      double ex = (k+0.5)*de;
      if(ex > emax) break;
      if(ex < elmax) continue;

      double r = ldLevelDensity(ex,(double)ncl[tid].za.getA(),&ncl[tid].ldp);

      for(int j=0 ; j<MAX_J ; j++){
        double sd = ldSpinDistribution(j+halfint(ncl[tid].lev[0].spin),ncl[tid].ldp.spin_cutoff,ncl[tid].ldp.sigma0,1.0);
        if(sd < eps) break;
        double r1 = r*sd*ldParityDistribution( 1);
        double r2 = r*sd*ldParityDistribution(-1);
        double p  = (2.0*j+1.0) * exp(-ex/t[n]) * de;
        g[n] +=  p * r1;
        g[n] +=  p * r2;
      }
      if(g[n] > gcut){
        g[n] = gcut;
        break;
      }
    }
  }
}


/**********************************************************/
/*      Maxwellian Average, Simpson Integration           */
/**********************************************************/
double macsAverage(const int n, double *x, double *y, double t)
{
  const int subdiv = 50;

  double s1 = 0.0, s2 = 0.0;
  for(int i=0 ; i<n-1 ; i++){
    double d = x[i+1] - x[i];
    if(d == 0.0) continue;

    double d1 = d/(double)subdiv;
    double y1 = (y[i+1] - y[i])/(double)subdiv;

    double f = 0.0, p1 = 0.0, p2 = 0.0;
    for(int j=0 ; j<=subdiv ; j++){
      if( (j == 0) || (j == subdiv) ) f = d1    /3.0;
      else if((j%2) == 0)             f = d1*2.0/3.0;
      else                            f = d1*4.0/3.0;

      double e = x[i] + d1* j;
      double u = y[i] + y1* j;
      double w = e * exp(-e/t);

      p1 += f * w * u;
      p2 += f * w;
    }
    s1 += p1;
    s2 += p2;
  }

  double macs = s1/s2 * 2.0/sqrt(PI);

  return(macs);
}

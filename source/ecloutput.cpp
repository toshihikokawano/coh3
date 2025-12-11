/******************************************************************************/
/*  ecloutput.cpp                                                             */
/*        print out calculated exclusive spectra and DDX                      */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>

#include "physicalconstant.h"
#include "structur.h"
#include "nucleus.h"
#include "eclipse.h"
#include "output.h"
#include "outformat.h"


/**********************************************************/
/*      Top Part                                          */
/**********************************************************/
void eclOutHead(const int n, const int c, const double elab)
{
  std::cout << std::setprecision(6) << std::setiosflags(std::ios::scientific);
  std::cout << "# ECLIPSE    "
            << std::setw(13) << elab
            << std::setw(13) << n
            << std::setw(13) << c    << std::endl;
}


/**********************************************************/
/*      Ouput Individual Spectrum                         */
/**********************************************************/
void eclOutSpectra(const int nm, const int km, const int cm, const double de, double **np, EXSpectra *dat)
{
  std::cout << std::setprecision(6) << std::setiosflags(std::ios::scientific);

  for(int n=0 ; n<nm ; n++){

    /*** find non-zero element from the high-side */
    int  km2  = eclFindKmax(km,cm,dat[n].spec);

    if(km2 <= 0) continue;

    std::cout << "# ";
    outZA(&ncl[n].za);
    for(int c=0 ; c<cm ; c++) std::cout << std::setw(13) << particle_name[c];
    std::cout << " Disc.Gamma   Bin.React    Bin.Gamma" <<  std::endl;

    std::cout << "#   Multipl. ";
    for(int c=0 ; c<cm ; c++) std::cout << std::setw(13) << np[n][c];
    nl();

    for(int k=0 ; k<km2 ; k++){
      std::cout << std::setw(13) << k*de;
      for(int c=0 ; c<cm ; c++) std::cout << std::setw(13) << dat[n].spec[c][k];
      std::cout << std::setw(13) << dat[n].glin[k];
      if(n < cm){
        std::cout << std::setw(13) << dat[n].sbin[k];
        std::cout << std::setw(13) << dat[n].gbin[k];
      }
      nl();
    }
    nl();
    nl();
    nl();
  }
}


/**********************************************************/
/*      Ouput Total Spectra                               */
/**********************************************************/
void eclTotalSpectra(const int nm, const int km, const int cm, const double de, EXSpectra *dat)
{
  double *sum;
  sum = new double [cm];

  /*** sum all the partial spectra to the first section in dat */
  for(int n=1 ; n<nm ; n++){
    for(int k=0 ; k<km ; k++){
      for(int c=0 ; c<cm ; c++) dat[0].spec[c][k] += dat[n].spec[c][k];
    }
  }

  /*** find non-zero element from the high-side */
  int  km2  = eclFindKmax(km,cm,dat[0].spec);

  if(km2 > 0){
    for(int c=0 ; c<cm ; c++){
      sum[c] = 0.0;
      for(int k=0 ; k<km2 ; k++) sum[c] += dat[0].spec[c][k];
    }

    std::cout << "# Total Spectra" << std::endl;

    for(int k=0 ; k<km2 ; k++){
      std::cout << std::setw(13) << k*de;
      for(int c=0 ; c<cm ; c++) std::cout << std::setw(13) << dat[0].spec[c][k];
      nl();
    }
    std::cout << "# Sum        ";
    for(int c=0 ; c<cm ; c++) std::cout << std::setw(13) << sum[c];
    nl();
    nl();
    nl();
  }

  delete [] sum;
}


/**********************************************************/
/*      Write Legendre Coefficient, Header Part           */
/**********************************************************/
void eclOutNucleusHead(const int n, const int cm, double *np)
{
  static char p[] = {'g','n','p','a','d','t','h'};
  int j = 0;

  std::cout << "# Nucleus    ";
  std::cout <<" ";

  if(np[0] < 0)  std::cout << "            "; // for inclusive case
  else{
    for(int i=1 ; i<7 ; i++){
      j = (i >= cm) ? 0 : (int)np[i];
      std::cout << std::setw(1) << p[i] << std::setw(1)<< j;
    }
  }
  std::cout << std::setw(13) << n << std::endl;
}


/**********************************************************/
/*      Write Legendre Coefficient, Sub-Header Part       */
/**********************************************************/
void eclOutChannelHead(const int c, const double mp)
{
  int nsec = 0;

  if(c == 0) nsec = (mp > 0.0) ? 2 : 0;
  else       nsec = (mp > 0.0) ? 1 : 0;

  std::cout << "# Particle   ";
  std::cout << std::setw(13) << c
            << std::setw(13) << nsec
            << std::setw(13) << mp << std::endl;
}


/**********************************************************/
/*      Write Legendre Coefficient                        */
/**********************************************************/
void eclOutLegCoeff(const int km, const int nleg, double *ep, double **cleg)
{
  if(km <= 0){
    std::cout << "# Continuum  " << std::setw(13) << 0 << std::setw(13) << 0 << std::endl;
  }
  else{
    int k0 = 1, k1 = km-2, k2 = km;
    /*** cut low side zeros, keep the first two points */
    for(int k=2 ; k<km ; k++){
      bool p = false;
      for(int l=0 ; l<nleg ; l++) if(cleg[l][k] != 0.0) p = true;
      if(p){
        k0 = k-1;
        break;
      }
    }

    /*** high side check */
    for(int k=km-3 ; k>=k0 ; k--){
      bool p = false;
      for(int l=0 ; l<nleg ; l++) if(cleg[l][k] != 0.0) p = true;
      if(p){
        k1 = k+1;
        break;
      }
    }
    k2 = k1 - k0 + 3;

    std::cout << "# Continuum  " << std::setw(13) << k2 << std::setw(13) << nleg << std::endl;

    std::cout << std::setw(13) << ep[0];
    for(int l=0 ; l<nleg ; l++) std::cout << std::setw(13) << cleg[l][0];
    nl();

    for(int k=k0 ; k<=k1 ; k++){
      std::cout << std::setw(13) << ep[k];
      for(int l=0 ; l<nleg ; l++) std::cout << std::setw(13) << cleg[l][k];
      nl();
    }

    std::cout << std::setw(13) << ep[km-1];
    for(int l=0 ; l<nleg ; l++) std::cout << std::setw(13) << cleg[l][km-1];
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Write Discrete Gamma Lines                        */
/**********************************************************/
void eclOutGammaLine(const int ng, const int nleg, double *e, double *x)
{
  std::cout << "# Discrete   " << std::setw(13) << ng << std::setw(13) << nleg << std::endl;
  for(int j=0 ; j<ng ; j++){
    std::cout << std::setw(13) << e[j] << std::setw(13) << x[j]<< std::endl;
  }
  nl();
  nl();
}


/**********************************************************/
/*      Print Out Decay Probability Table                 */
/**********************************************************/
void eclOutPtable(const int n, const int m, int *ktot, int **cidx, double ****ptbl)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);

  for(int c0=0 ; c0<n ; c0++){
    for(int k0=0 ; k0<=ktot[c0] ; k0++){

      std::cout << std::setw(5) << c0 << std::setw(4) << k0 << std::setw(4) << ktot[c0];
      for(int j=0 ; j<m ; j++) std::cout << std::setw(13) <<  cidx[c0][j];
      nl();

      for(int k1=k0 ; k1<=ktot[c0] ; k1++){
        std::cout << std::setw(13) << k1;
        for(int c1=0 ; c1<m ; c1++){
          std::cout << std::setw(13) << ptbl[c0][k0][c1][k1];
        }
        nl();
      }
    }
  }
}


/**********************************************************/
/*      Print Inclusive DDX                               */
/**********************************************************/
void eclOutInclusiveSpectrum(const int c, const int kmax, const int nang, double *xang, int *nleg, const double de, double **cleg, double **fleg, double *espc)
{
  std::ostringstream title;
  title << particle_name[c] << " DOUBLE DIFFERENTIAL CROSS SECTION";
  outSectionHead( (title.str()).c_str() );
  std::cout << cline << std::endl;

  /*** gamma-ray spectrum, isotropic */
  if(c == 0){
    std::cout << "# Emin[MeV]  Emax[MeV]";
    std::cout << "  EDX [mb/MeV]" << std::endl;
  }
  /*** particle spectra */
  else{
    std::cout << "#          Angle [deg]";
    for(int n=0 ; n<nang ; n++) outVal(14,1,xang[n]);
    nl();
    std::cout << "# Emin[MeV]  Emax[MeV]";
    std::cout << " [mb/MeV/sr]" << std::endl;
  }

  for(int k=0 ; k<kmax ; k++){
    if(k == 0){
      outVal(11,4,0.0);
      outVal(11,4,0.5*de);
    }
    else{
      outVal(11,4,((double)k-0.5)*de);
      outVal(11,4,((double)k+0.5)*de);
    }

    if(c == 0){
      outVal(14,lowfilter( espc[k] * cleg[0][k] ));
    }
    else{
      for(int n=0 ; n<nang ; n++){
        double y = 0.0;
        for(int j=0 ; j<nleg[k] ; j++) y += cleg[j][k] * fleg[j][n];
        outVal(14,lowfilter( espc[k] * y ));
      }
    }
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Look For Non-Zero Component                       */
/**********************************************************/
int eclFindKmax(const int km, const int cm, double **dat)
{
  bool kend = false;
  int  km2  = km-1;

  for(int k=km-1 ; k>0 ; k--){
    for(int c=0 ; c<cm ; c++){
      if(dat[c][k] > 0.0){
        kend = true;
        break;
      }
    }
    if(kend){
      km2 = k+2;
      break;
    }
  }
  return( (kend) ? km2 : -1 );
}



/******************************************************************************/
/*  omoutput.cpp                                                              */
/*        data output for optical model calculations                          */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "optical.h"
#include "omoutput.h"
#include "outformat.h"


/**********************************************************/
/*     Output optical potential parameters                */
/**********************************************************/
void outOMP(const int ctl, Optical *omp)
{
  outSectionHead("OPTICAL MODEL POTENTIAL");

  if(ctl >= 0){
    std::cout << "# OMP" << std::setw(3) << ctl << "   ";
    std::cout << "   V         Vs        Wv        Ws        Vso       Wso";
  }
  else  
    std::cout << cline << "   V         Vs        Wv        Ws        Vso       Wso";

  if(omp->rc != 0.0) std::cout << "       Coul" << std::endl;
  else std::cout << std::endl;

  std::cout << "      r[fm]";
  outVal(10,5,omp->r0);   outVal(10,5,omp->r0s);
  outVal(10,5,omp->rv);   outVal(10,5,omp->rs);
  outVal(10,5,omp->rvso); outVal(10,5,omp->rwso); if(omp->rc != 0.0) outVal(10,2,omp->rc);
  std::cout << std::endl;

  std::cout << "      a[fm]";
  outVal(10,5,omp->a0);   outVal(10,5,omp->a0s);
  outVal(10,5,omp->av);   outVal(10,5,omp->as);
  outVal(10,5,omp->avso); outVal(10,5,omp->awso);
  std::cout << std::endl;

  if(ctl < 0){
    std::cout << "  Depth(E0)";
    outVal(10,5,omp->v1);   outVal(10,5,omp->vs1);
    outVal(10,5,omp->wv1);  outVal(10,5,omp->ws1);
    outVal(10,5,omp->vso1); outVal(10,5,omp->wso1);
    std::cout << std::endl;

    std::cout << "  Depth(E1)";
    outVal(10,5,omp->v2);   outVal(10,5,omp->vs2);
    outVal(10,5,omp->wv2);  outVal(10,5,omp->ws2);
    outVal(10,5,omp->vso2); outVal(10,5,omp->wso2);
    std::cout << std::endl;

    std::cout << "  Depth(E2)";
    outVal(10,5,omp->v3);   outVal(10,5,omp->vs3);
    outVal(10,5,omp->wv3);  outVal(10,5,omp->ws3);
    outVal(10,5,omp->vso3); outVal(10,5,omp->wso3);
    std::cout << std::endl;
  }
  std::cout << " Depth[MeV]";

  outVal(10,5,omp->volume.real());     outVal(10,5,omp->surface.real());
  outVal(10,5,omp->volume.imag());     outVal(10,5,omp->surface.imag());
  outVal(10,5,omp->spin_orbit.real()); outVal(10,5,omp->spin_orbit.imag());
  std::cout << std::endl;
}


void outOMPtable(const int ctl, Optical *omp)
{
  if(ctl == 0) outVal(8,3,outRetrieveLabE());
  else std::cout << "        ";

  outVal(8,3,omp->volume.real());     outVal(8,3,omp->r0);   outVal(8,3,omp->a0);
  outVal(8,3,omp->surface.real());    outVal(8,3,omp->r0s);  outVal(8,3,omp->a0s);
  outVal(8,3,omp->volume.imag());     outVal(8,3,omp->rv);   outVal(8,3,omp->av);
  outVal(8,3,omp->surface.imag());    outVal(8,3,omp->rs);   outVal(8,3,omp->as);
  outVal(8,3,omp->spin_orbit.real()); outVal(8,3,omp->rvso); outVal(8,3,omp->avso);
  outVal(8,3,omp->spin_orbit.imag()); outVal(8,3,omp->rwso); outVal(8,3,omp->awso);
  outVal(8,2,omp->rc);
  std::cout << std::endl;
}


/**********************************************************/
/*     Output collective states                           */
/**********************************************************/
void outCoupledState(const NuclearModel m, const int n, LevelData *lev)
{
  if(m == vibration) outSectionHead("VIBRATIONAL STATE DATA");
  else               outSectionHead("ROTATIONAL STATE DATA");
  std::cout << cline << "   Energy       Spin   WaveNum Eout[MeV]   Coulomb   Sig0   ";
  if(m == vibration) std::cout << " Phonon";
  std::cout << std::endl;
  for(int i=0 ; i<n ; i++){
    std::cout << "      Level";
    outVal(10,5,lev[i].excitation);
    outVal(10,1,lev[i].spin);
    outVal(10,5,lev[i].wave_number);
    outVal(10,5,lev[i].energy);
    outVal(10,5,lev[i].coulomb.eta);
    outVal(10,5,lev[i].coulomb.sigma0);
    if(m == vibration){ std::cout << std::setw(5) << lev[i].phonon; }
    std::cout << std::endl;
  }
}


/**********************************************************/
/*     Output deformation parameters                      */
/**********************************************************/
void outDeformation(const NuclearModel m, const int n, double *beta)
{
  outSectionHead("DEFORMATION PARAMETERS");
  std::cout << cline << "    Lambda     Beta" << std::endl;
  if(m == rotation){
    for(int i=0 ; i<n/2 ; i++){
      std::cout << "       Beta";
      outVal(10,(i+1)*2); outVal(10,5,beta[i]);
      std::cout << std::endl;
    }
  }
  else{
    for(int i=0 ; i<n ; i++){
      std::cout << "       Beta";
      outVal(10,i+2); outVal(10,5,beta[i]);
      std::cout << std::endl;
    }
  }
}


/**********************************************************/
/*      Discrete Level Population by Direct               */
/**********************************************************/
void outLevelExcite(const int n0, const int n1, const double ecms, const double ex, LevelData *lev, double *cdir)
{
  outSectionHead("DISCRETE DIRECT CROSS SECTION");
  std::cout << cline << "    Ex[MeV]     Spin    Eout[MeV]  Sigma[mb]" << std::endl;

  for(int i=n0 ; i<n1 ; i++){
    if(cdir[i] < 0.0) continue;
    std::cout << "  Direct" << std::setw(3) << i;
    outVal(11,5,lev[i].excitation);
    outVal( 9,1,lev[i].spin);
    std::cout << "  ";
    outVal(ecms-lev[i].excitation+ex);
    outVal(cdir[i]);
    std::cout << std::endl;
  }
}


/**********************************************************/
/*      Output R-Matrix and Strength Functions            */
/**********************************************************/
void outRmatrix(const int n, std::complex<double> *s, std::complex<double> *r)
{
  outSectionHead("R MATRIX AND STRENGTH FUNCTION");

  std::cout << cline << "   Rr         Ri         R'       PoleStr.   Str. Func."<<std::endl;
  for(int j=0 ; j<=n ; j++){
    std::cout << std::setw(11) << j;
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout << std::setprecision(3);
    std::cout << std::setw(11) << r[j].real() << std::setw(11) << r[j].imag();
    std::cout << std::setw(11) << s[j].real() << std::setw(11) << r[j].imag()/PI;
    std::cout << std::setw(11) << s[j].imag() << std::endl;
  }
}


/**********************************************************/
/*      Output S-Matrix Elements                          */
/**********************************************************/
void outSmatrix(const int n, const int spin2, std::complex<double> *s)
{
  outSectionHead("S MATRIX ELEMENTS");

  std::cout << "#        L ";
  switch(spin2){
  case 0 : std::cout << "  Re(L)"; break;
  case 1 : std::cout << "  Re(L+1/2)  Im(L+1/2)  Re(L-1/2)  Im(L-1/2)"; break;
  case 2 : std::cout << "  Re(L+1)    Im(L+1)    Re(L)      Im(L)      Re(L-1)    Im(L-1)"; break;
  default: break;
  }
  std::cout << std::endl;

  for(int l=0 ; l<=n ; l++){
    std::cout << std::setw(11) << l;
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout << std::setprecision(3);
    for(int j=0 ; j<=2 ; j++){
      if(j <= spin2){
        std::cout << std::setw(11) << s[3*l+j].real();
        std::cout << std::setw(11) << s[3*l+j].imag();
      }
    }
    std::cout << std::endl;
  }
}


/**********************************************************/
/*      Transmission Coefficients                         */
/**********************************************************/
void outTransmission(const int n, const int spin2, const double f, double *tran)
{
  outSectionHead("TRANSMISSION COEFFICIENTS");

  std::cout << "#        L ";
  switch(spin2){
  case 0 : std::cout << "   T(L)                          "; break;
  case 1 : std::cout << "   T(L+1/2)   T(L-1/2)           "; break;
  case 2 : std::cout << "   T(L+1)     T(L)       T(L-1)  "; break;
  default: break;
  }
  std::cout << "Sig Partial" << std::endl;

  double t = 0.0;
  for(int l=0 ; l<=n ; l++){
    std::cout << std::setw(11) << l;
    for(int j=0 ; j<=2 ; j++){
      if(j <= spin2) outVal(tran[3*l+j]);
      else std::cout << blank;
    }

    double s = 0.0;
    if(spin2 == 0)      s =  (2*l+1)*tran[l*3  ];
    else if(spin2 == 1) s =    (l+1)*tran[l*3  ] +   l   *tran[l*3+1];
    else if(spin2 == 2) s = ((2*l+3)*tran[l*3  ] +(2*l+1)*tran[l*3+1] + (2*l-1)*tran[l*3+2])/3.0;
    outVal(s*f);
    t += s*f;
    std::cout << std::endl;

    if(s == 0.0) break;
  }

  std::cout << blank << blank << blank << "        Sum";
  outVal(t);
  std::cout << std::endl;
}


/**********************************************************/
/*      Output S-Matrix Elements for Coupled-Channels     */
/**********************************************************/
void outCCSmatrix(const int nc, const int j, const int p, Collective *col, CCdata *cdt, std::complex<double> *s)
{
  outSectionHead("COUPLED-CHANNELS S MATRIX ELEMENTS");

  std::cout << "#   JPi  Nc";
  outVal(9,1,j/2.0);
  std::cout << std::setw(2) << ((p>0) ? '+' : '-') << std::setw(11) << nc << std::endl;
  std::cout << "#  n   I   l   j" << std::endl;

  for(int i=0 ; i<nc ; i++){

    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout << std::setprecision(1);

    std::cout << std::setw(4) << cdt[i].level;
    std::cout << std::setw(4) << fabs(col->lev[cdt[i].level].spin);
    std::cout << std::setw(4) << cdt[i].chn.l;
    std::cout << std::setw(4) << cdt[i].chn.j2 * 0.5;

    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout << std::setprecision(5);
    for(int j=0 ; j<nc ; j++){
      int ij = i*nc+j;
      std::cout << ' ' << std::setw(12) << s[ij].real() << ' ' << std::setw(12) << s[ij].imag();
    }
    std::cout << std::endl;
  }
}

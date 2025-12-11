/******************************************************************************/
/*  mfldprint.cpp                                                             */
/*        Output Calculated Results                                           */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "mfld.h"

static int findJmax(const int, Density *);

static const std::string blank = "              ";

/**********************************************************/
/*      Print Calculated Level / State Density            */
/**********************************************************/
void printLevelDensity(const bool bothparity, Density *sd)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);

  std::cout << "# Total and Particle-Hole Level Densities" << std::endl;

  std::cout << "# E [MeV]  " << " Total [1/MeV]";
  if(!bothparity) std::cout << blank;
  for(int p=0 ; p<sd->getPsize() ; p++){
    std::cout << std::setw(2) << p+1 << "p" << " - " << std::setw(2) << p+1 << "h     ";
    if(!bothparity) std::cout << blank;
  }
  std::cout << std::endl;
  std::cout << "#          ";
  for(int p=0 ; p<sd->getPsize()+1 ; p++){
    if(bothparity) std::cout << " Even + Odd   ";
    else           std::cout << " Even          Odd          ";
  }
  std::cout << std::endl;

  double z0 = 0.0, z1 = 0.0;
  for(int n=0 ; n<sd->getNsize() ; n++){
    std::cout << std::setw(11) << std::setprecision(4) << sd->binwidth() * n;
    std::cout << std::setprecision(7);

    /*** total density */
    z0 = z1 = 0.0;
    for(int p=0 ; p<sd->getPsize() ; p++){
      for(int j=0 ; j<sd->getJsize() ; j++){
        z0 += sd->r0[p][n][j];
        z1 += sd->r1[p][n][j];
      }
    }
    if(bothparity) std::cout << std::setw(14) << z0 + z1;
    else           std::cout << std::setw(14) << z0 << std::setw(14) << z1;

    /*** for each p-h */
    for(int p=0 ; p<sd->getPsize() ; p++){
      z0 = z1 = 0.0;
      for(int j=0 ; j<sd->getJsize() ; j++){
        z0 += sd->r0[p][n][j];
        z1 += sd->r1[p][n][j];
      }
      if(bothparity) std::cout << std::setw(14) << z0 + z1;
      else           std::cout << std::setw(14) << z0 << std::setw(14) << z1;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;
}


/**********************************************************/
/*      Print Calculated Partial Level Density            */
/**********************************************************/
void printPartialLevelDensity(const bool halfint, const bool bothparity, Density *sd)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);

  std::cout << "# Partial Level Densities" << std::endl;
  std::cout << "# DataSize";
  std::cout << std::setw(5) << sd->getPsize();
  std::cout << std::setw(5) << sd->getNsize();
  std::cout << std::setw(5) << sd->getJsize() << std::endl;

  /*** for each p-h */
  for(int p=0 ; p<sd->getPsize() ; p++){

    /*** find maximum J */
    int jmax = findJmax(p,sd);

    if(bothparity){
      std::cout << "# PH" << std::setw(3) << p+1 << " Jmax" << std::setw(3) << jmax << " Even + Odd" << std::endl;

      /*** energy bins */
      std::cout << "# E [MeV]  ";
      for(int j=0 ; j<jmax ; j++){
        if(halfint) std::cout << std::setw(2) << 2*j+1 << "/2+         ";
        else        std::cout << std::setw(2) << j << "+           ";
      }
      std::cout << std::endl;

      for(int n=0 ; n<sd->getNsize() ; n++){
        std::cout << std::setw(11) << std::setprecision(4) << sd->binwidth() * n;
        std::cout << std::setprecision(7);

        /*** spins */
        for(int j=0 ; j<jmax ; j++) std::cout << std::setw(14) << sd->r0[p][n][j] + sd->r1[p][n][j];
        std::cout << std::endl;
      }
      std::cout << std::endl;
      std::cout << std::endl;
    }
    else{
      std::cout << "# PH" << std::setw(3) << p+1 << " Jmax" << std::setw(3) << jmax << " Even" << std::endl;

      /*** even parities, energy bins */
      std::cout << "# E [MeV]  ";
      for(int j=0 ; j<jmax ; j++){
        if(halfint) std::cout << std::setw(2) << 2*j+1 << "/2+         ";
        else        std::cout << std::setw(2) << j << "+           ";
      }
      std::cout << std::endl;

      for(int n=0 ; n<sd->getNsize() ; n++){
        std::cout << std::setw(11) << std::setprecision(4) << sd->binwidth() * n;
        std::cout << std::setprecision(7);

        /*** spins */
        for(int j=0 ; j<jmax ; j++) std::cout << std::setw(14) << sd->r0[p][n][j];
        std::cout << std::endl;
      }
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << "# PH" << std::setw(3) << p+1 << " Jmax" << std::setw(3) << jmax << " Odd" << std::endl;

      /*** for odd parities */
      std::cout << "# E [MeV]  ";
      for(int j=0 ; j<jmax ; j++){
        if(halfint) std::cout << std::setw(2) << 2*j+1 << "/2-         ";
        else        std::cout << std::setw(2) << j << "-           ";
      }
      std::cout << std::endl;

      for(int n=0 ; n<sd->getNsize() ; n++){
        std::cout << std::setw(11) << std::setprecision(4) << sd->binwidth() * n;
        std::cout << std::setprecision(7);

        for(int j=0 ; j<jmax ; j++) std::cout << std::setw(14) << sd->r1[p][n][j];
        std::cout << std::endl;
      }

      std::cout << std::endl;
      std::cout << std::endl;
    }
  }
}


/**********************************************************/
/*      Better Format for Data Plotting                   */
/**********************************************************/
void printPlottingLevelDensity(const bool bothparity, Density *sd)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);

  /*** for each p-h */
  for(int p=0 ; p<sd->getPsize() ; p++){

    int jmax = findJmax(p,sd);

    std::cout << "# PH" << std::setw(3) << p+1 << std::endl;

    for(int j=0 ; j<jmax ; j++){

      for(int n=0 ; n<sd->getNsize() ; n++){
        std::cout << std::setw(11) << std::setprecision(4) << sd->binwidth() * n;
        std::cout << std::setw(5) << j;

        std::cout << std::setprecision(7);
        if(bothparity) std::cout << std::setw(14) << sd->r0[p][n][j] + sd->r1[p][n][j];
        else           std::cout << std::setw(14) << sd->r0[p][n][j] << std::setw(14) << sd->r1[p][n][j];
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
}


/**********************************************************/
/*      Print JPi-Depend Level Density                    */
/**********************************************************/
void printJDependentLevelDensity(const bool halfint, const bool bothparity, Density *sd)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);

  std::cout << "# JDependent Level Densities" << std::endl;
  std::cout << "# DataSize";
  std::cout << std::setw(5) << sd->getNsize();
  std::cout << std::setw(5) << sd->getJsize() << std::endl;

  int jmax = 0;
  for(int p=0 ; p<sd->getPsize() ; p++){
    int jm = findJmax(p,sd);
    if(jm > jmax) jmax = jm;
  }

  std::cout << "# E [MeV]  ";
  for(int j=0 ; j<jmax ; j++){
    if(halfint) std::cout << std::setw(2) << 2*j+1 << "/2+         ";
    else        std::cout << std::setw(2) << j << "+           ";
    if(!bothparity) std::cout << blank;
  }
  std::cout << std::endl;

  for(int n=0 ; n<sd->getNsize() ; n++){
    std::cout << std::setw(11) << std::setprecision(4) << sd->binwidth() * n;

    std::cout << std::setprecision(7);
    for(int j=0 ; j<sd->getJsize() ; j++){
      double x0 = 0.0, x1 = 0.0;
      for(int p=0 ; p<sd->getPsize() ; p++){
        x0 += sd->r0[p][n][j];
        x1 += sd->r1[p][n][j];
      }
      if(bothparity) std::cout << std::setw(14) << x0 + x1;
      else           std::cout << std::setw(14) << x0 << std::setw(14) << x1;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;
}


/**********************************************************/
/*      Spin Distribution of Level Density                */
/**********************************************************/
void printSpinDistribution(const bool spincutoff, const double ewidth, Density *sd)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);

  double *z0 = new double [sd->getJsize()];
  double *z1 = new double [sd->getJsize()];

  int kmax = (sd->binwidth() * sd->getNsize()) / ewidth + 1;

  std::cout << "# Spin Distributoin" << std::endl;
  std::cout << "# DataSize";
  std::cout << std::setw(5) << sd->getPsize();
  std::cout << std::setw(5) << kmax << std::endl;

  double emin, emax;
  for(int p=0 ; p<sd->getPsize() ; p++){

    int jmax = findJmax(p,sd);

    if(spincutoff){
      std::cout << "# PH" << std::setw(3) << p+1 << std::endl;
      std::cout << "# Emin      Emax       <J>        sigma1     sigma2" << std::endl;
    }
    else{
      std::cout << "# PH" << std::setw(3) << p+1 << " Jmax" << std::setw(3) << jmax << std::endl;
      std::cout << "# Energy     J " << std::endl;
    }

    for(int k=0 ; k<kmax ; k++){
      emin = k * ewidth;
      emax = emin + ewidth;

      for(int j=0 ; j<sd->getJsize() ; j++) z0[j] = z1[j] = 0.0;

      /*** re-binning */
      for(int n=0 ; n<sd->getNsize() ; n++){
        double e0 = sd->binwidth() * n;
        double e1 = sd->binwidth() * (n+1);
        double em = (e0 + e1) / 2.0;

        if( (emin <= em) && (em < emax) ){
          double de = sd->binwidth();
          if( (e0 <= emin) && (emin < e1) ) de = e1 - emin;
          else if( (e0 <= emax) && (emax < e1) ) de = emax - e0;

          de = de / ewidth;

          for(int j=0 ; j<sd->getJsize() ; j++){
            z0[j] += sd->r0[p][n][j] * de;
            z1[j] += sd->r1[p][n][j] * de;
          }
        }
      }

      /*** normalize and average J */
      double x = 0.0, y = 0.0;
      for(int j=0 ; j<jmax ; j++){
        x +=       z0[j] + z1[j];
        y +=  j * (z0[j] + z1[j]);
      }
      if(x > 0.0) y = y/x;

      /*** spin cut-off parameter (not squared) */
      double s1 = 0.796884 * y + 0.414469; // relation between <J> and sigma1, no (2J+1) equation
      double s2 = 0.626089 * y + 0.319161; // relation between <J> and sigma2
 
      if(y == 0.0) s1 = s2 = 0.0;
      if(s2 < 0.0) s2 = 0.0;

      if(spincutoff){
        std::cout << std::setprecision(4);
        std::cout << std::setw(11) << emin;
        std::cout << std::setw(11) << emax;
        std::cout << std::setw(11) << y;
        std::cout << std::setw(11) << s1;
        std::cout << std::setw(11) << s2 << std::endl;
      }
      else{
        double q0 = 0.0, q1 = 0.0;
        for(int j=0 ; j<jmax ; j++){

          std::cout << std::setprecision(4) << std::setw(11) << (emin + emax) / 2.0;
          std::cout << std::setw(3) << j;

          double z = (x == 0.0) ? 0.0 : (z0[j] + z1[j]) / x;
          double r = (s1 == 0.0) ? 0.0 : (j+0.5)/(s1*s1) * exp( -(j+0.5)*(j+0.5)/(2.0*s1*s1) );

          q0 += z;
          q1 += r;

          std::cout << std::setprecision(4);
          std::cout << std::setw(11) << z;
          std::cout << std::setw(11) << r;
          std::cout << std::endl;
        }

        std::cout << "# Sum/Average ";
        std::cout << std::setw(11) << q0 << std::setw(11) << q1 << std::endl;
        std::cout << std::endl;
      }
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }

  delete [] z0;
  delete [] z1;
}


/**********************************************************/
/*      Average Number of Particle-Hole                   */
/**********************************************************/
void printAverageConfiguration(Density *sd)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);

  std::cout << "# Average Particle-Hole Configuration" << std::endl;

  std::cout << "# E [MeV]  " << "     Average N" << std::endl;

  for(int n=0 ; n<sd->getNsize() ; n++){
    std::cout << std::setw(11) << std::setprecision(4) << sd->binwidth() * n;
    std::cout << std::setprecision(7);

    /*** total density */
    double zt = 0.0;
    for(int p=0 ; p<sd->getPsize() ; p++){
      for(int j=0 ; j<sd->getJsize() ; j++) zt += sd->r0[p][n][j] + sd->r1[p][n][j];
    }

    /*** average weighted by p+h */
    double zp = 0.0;
    if(zt > 0.0){
      for(int p=0 ; p<sd->getPsize() ; p++){
        for(int j=0 ; j<sd->getJsize() ; j++){
          zp += (sd->r0[p][n][j] + sd->r1[p][n][j]) * (p+1) * 2;
        }
      }
      zp /= zt;
    }

    std::cout << std::setw(14) << zp << std::endl;
  }
}


/**********************************************************/
/*      Find Maximum J value                              */
/**********************************************************/
int findJmax(const int p, Density *sd)
{
  int jmax = 0;
  for(int n=0 ; n<sd->getNsize() ; n++){
    int jm = 0;
    for(jm = sd->getJsize()-1 ; jm >= 0 ; jm--){
      if((sd->r0[p][n][jm] + sd->r1[p][n][jm]) != 0.0) break;
    }
    if(jmax < jm+1) jmax = jm+1;
  }
  jmax ++;
  if(jmax > sd->getJsize()) jmax = sd->getJsize();

  return jmax;
}


/**********************************************************/
/*      Monitoring of Number of States                    */
/**********************************************************/
void printMState(std::string c, const int q, const int p, MCount *mc)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(2);

  int mmax = mc->getMsize(); if(mmax > 26) mmax = 26;
  int nmax = mc->getNsize();
  
  std::cout << "# " << std::setw(2) << c << std::setw(6) << p << std::endl;
  std::cout << "     ";
  for(int m=q ; m<=mmax ; m+=2) std::cout << std::setw(10) << m;
  std::cout << std::endl;

  for(int i=0 ; i<nmax ; i++){
    std::cout << std::setw(5) << i;
    for(int m=q ; m<=mmax ; m+=2){
      unsigned long mstate = mc->m0[i][m] + mc->m1[i][m];
      double x = mstate / 10.0;
      if(x == 0.0) std::cout << "          ";
      else if(x > 999909999.0) std::cout << " xxxxxxxxx";
      else std::cout << std::setw(10) << (int)x;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;
}




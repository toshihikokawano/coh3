/******************************************************************************/
/*  eclenergy.cpp                                                             */
/*        check energy balance                                                */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "physicalconstant.h"
#include "structur.h"
#include "nucleus.h"
#include "eclipse.h"
#include "outformat.h"


/**********************************************************/
/*      Check Energy Balance and Fix Gamma Multiplicity   */
/**********************************************************/
void eclCheckEnergyBalance(const int cm, double *np, int *kmax, double **epar, const int ngam, double *gx, double *gy, double ***cleg, EXSpectra *dat, Nucleus *n)
{
  double *espc = new double [cm];
  double *edsc = new double [cm];
  double *egam = new double [cm];

  /*** available energy for this reaction */
  double emax = dat->getEm();

  /*** integrate continuum spectra, sum(cleg[0] + gy) are normalized to unity */
  for(int c=0 ; c<cm ; c++){
    espc[c] = 0.0;
    if( (c > 0) && (np[c] == 0.0) ) continue;
    for(int k=0 ; k<kmax[c]-1 ; k++){
      double e0 = epar[c][k  ];
      double e1 = epar[c][k+1];
      double s0 = cleg[c][0][k  ];
      double s1 = cleg[c][0][k+1];
      espc[c] += (s0*(e1+2.0*e0) + s1*(e0+2.0*e1)) * (e1-e0) / 6.0;
    }

    espc[c] *= np[c];
  }

  /*** sum discrete gammas */
  for(int k=0 ; k<ngam ; k++) espc[0] += gx[k] * gy[k] * np[0];

  double reac = crx.total - crx.elastic;

  /*** discrete particle emission */
  for(int c=0 ; c<cm ; c++){
    edsc[c] = egam[c] = 0.0;

    if( (c > 0) && (np[c] == 0.0) ) continue;
    if(dat->isBinary()){
      for(int k=0 ; k<n->ndisc ; k++){
        double p = crx.levexcite[c][k] / reac;
        edsc[c] += p * (emax -  n->lev[k].energy);  // energy of discrete particle transition
        egam[c] += p * n->lev[k].energy;            // energy emitted by gammas
      }
    }
  }

  double esum = 0.0;

  std::cout << "# Total " << std::setw(14) << emax << std::endl;
  std::cout << "#         cont          disc          gamma" << std::endl;
  for(int c=0 ; c<cm ; c++){
    std::cout << "#      ";
    std::cout << std::setw( 1) << p_name[c];
    std::cout << std::setw(14) << espc[c];
    std::cout << std::setw(14) << edsc[c];
    std::cout << std::setw(14) << egam[c] << std::endl;
    esum += espc[c] + edsc[c] + egam[c];
  }

  std::cout << "# ESum  " << std::setw(14) << esum << std::endl;

  delete [] espc;
  delete [] edsc;
  delete [] egam;
}


/**********************************************************/
/*      Reconstruct Total Energy Spectrum                 */
/**********************************************************/
void eclTotalEnergySpectrum(const int cm, const int km, const double sigma, double *np, EXSpectra *dat, Nucleus *n)
{
  double **spec = new double * [cm];
  for(int c=0 ; c<cm ; c++){
    spec[c] = new double [km];
    for(int k=0 ; k<km ; k++) spec[c][k] = 0.0;
  }

  /*** copy and normalize continumm spectrum */
  for(int c=0 ; c<cm ; c++){
    double sum = 0.0;
    for(int k=0 ; k<km ; k++){
      if(c == 0) sum += (dat->spec[c][k] + dat->glin[k] + dat->gbin[k]) * n->de;
      else       sum += (dat->spec[c][k] + dat->sbin[k]) * n->de;
    }
    if(sum > 0.0){
      double f = sigma * np[c] / sum;
      for(int k=0 ; k<km ; k++){
        if(c == 0) spec[c][k] = f * (dat->spec[c][k] + dat->glin[k] + dat->gbin[k]);
        else       spec[c][k] = f * (dat->spec[c][k] + dat->sbin[k]);
      }
    }
  }


  std::cout << "# Energy      ";
  for(int c=0 ; c<cm ; c++){
    std::cout << std::setw(14) << p_name[c];
  }
  std::cout << std::endl;

  int km2 = km;
  for(km2=km ; km2>=0 ; km2--){
    bool allzero = true;
    for(int c=0 ; c<cm ; c++){
      if(spec[c][km2] != 0.0) allzero = false;
    }
    if(!allzero) break;
  }
  km2++;

  for(int k=0 ; k<=km2 ; k++){
    std::cout << std::setw(14) << k * n->de;
    for(int c=0 ; c<cm ; c++){
      std::cout << std::setw(14) << spec[c][k];
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;

  for(int c=0 ; c<cm ; c++){
    delete [] spec[c];
  }
  delete [] spec;
}

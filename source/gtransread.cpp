/******************************************************************************/
/*  gtrans.cpp                                                                */
/*        Gamma-ray transmission coefficients                                 */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "datatable.h"
#include "terminate.h"

// local variables for reading photo-absorption cross section from an external file
const std::string photoabsfile = "CoHPhotoAbsorption.dat";
static DataTable pabs;


/**********************************************************/
/*      Read Photo-Absorption Table                       */
/**********************************************************/
void gdrAbsorptionDataRead()
{
  /*** set max size, 1000 energy points and MAX_multipolarities */
  pabs.setLimit(1000,MAX_MULTIPOL);

  if(dataTableRead(photoabsfile,&pabs) < 0){
    message << "photo-absorption data file " << photoabsfile << " not found";
    cohTerminateCode("photoAbsorptionRead");
  }

  if((pabs.getXsize() == 0) || (pabs.getYsize() == 0)){
    message << "no photo-absorption cross section given in " << photoabsfile;
    cohTerminateCode("photoAbsorptionRead");
  }
//dataTableWrite(&pabs);

  /*** convert cross section data into strength functions */
  double c0 = 1.0 / (HBARSQ * VLIGHTSQ * PI * PI) / NORM_FACT;

  for(int i=0 ; i<pabs.getXsize() ; i++){
    double e2 = pabs.xdata[i] * pabs.xdata[i];

    for(int j=0 ; j<pabs.getYsize() ; j++){
      int l = j/2 + 1; // L
      int n = 2*l + 1; // 2L + 1
      double el = pow(pabs.xdata[i],(double)n); // E^{2L+1}
      pabs.ydata[i][j] *= c0 * e2 / (double)n / el; // T = sigma E^2 const / (2L+1) / E^{2L+1}
    }
  }
//dataTableWrite(&pabs);
}


/**********************************************************/
/*      Interpolate Gamma-ray Transmission Coefficients   */
/**********************************************************/
double gdrGammaTransmissionInterpolate(const GammaMultipolarity gm, const double eg)
{
  double tg = 0.0;

  /*** interpolate table */
  int col = pabs.getYsize();

  dataTableInterpolate(eg,&pabs);

  switch(gm){
  case E1:
    if(col >= 1) tg = pabs.q[0];
    break;
  case M1:
    if(col >= 2) tg = pabs.q[1];
    break;
  case E2:
    if(col >= 3) tg = pabs.q[2];
    break;
  case M2:
    if(col >= 4) tg = pabs.q[3];
    break;
  default:
    break;
  }

  return(tg);
}


/**********************************************************/
/*      Highest Multipolarity Given in File               */
/**********************************************************/
int gdrAbsorptionDataYsize()
{  return( pabs.getYsize() ); }



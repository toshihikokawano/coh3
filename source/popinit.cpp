/******************************************************************************/
/*  popinit.cpp                                                               */
/*        initial population in the first compound                            */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "structur.h"
#include "statmodel.h"
#include "nucleus.h"

#undef DEBUG_POPCHECK

/**********************************************************/
/*      Store initial population at the given excitation  */
/**********************************************************/
double statStoreInitPopulation(const int spin0, const int targlev, const int targid, const double coef, Transmission *tin)
{
  int ig    = (int)(2.0*halfint(ncl[0].lev[0].spin));   // half-int ofset for CN
  int it    = (int)(2.0*ncl[targid].lev[targlev].spin); // target state spin
  int pt    = ncl[targid].lev[targlev].parity;          // target state parity
  double g  = coef /( (it+1.0)*(spin0+1.0) );           // spin factor without (2J+1)

  /*** Loop over CN J and Parity*/
  int jdxmax = 0;
  for(int j0=ig ; j0<=ncl[0].jmax*2+ig ; j0+=2){
    int   jdx = (j0-ig)/2;
    if(jdx > jdxmax) jdxmax = jdx;
    for(int p0=-1 ; p0<=1 ; p0+=2){

      /*** For entrance channel Tlj coupled to JP*/
      for(int lp0=0 ; lp0<=tin->lmax*2 ; lp0+=2){
        if(parity(lp0) != p0*pt) continue;

        for(int sp0=spin0 ; sp0>=-spin0 ; sp0-=2){
          int jp0 = lp0 + sp0;
          if(jp0 < 0) continue;

          if( std::abs(jp0-it) > j0 || j0 > (jp0+it) ) continue;

          double tlj  = tin->tran[tj_index(lp0,sp0,spin0)];
          if(tlj == 0.0) continue;

          if(p0 == 1) ncl[0].pop[0][jdx].even += (j0+1.0)*g*tlj;
          else        ncl[0].pop[0][jdx].odd  += (j0+1.0)*g*tlj;
        }
      }
    }
  }

  /*** calculate sum */
  double s = 0.0;
  for(int j=0 ; j<=jdxmax ; j++) s += ncl[0].pop[0][j].even + ncl[0].pop[0][j].odd;

#ifdef DEBUG_POPCHECK
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4);

  for(int j=0 ; j<20 ; j++) {
    std::cout << std::setw(10) << j;
    std::cout << std::setw(12) << ncl[0].pop[0][j].even;
    std::cout << std::setw(12) << ncl[0].pop[0][j].odd;
    std::cout << std::endl;
  }
  std::cout << std::setw(12) << s << std::endl;
#endif

  return(s);
}


/**********************************************************/
/*      Store initial population for photo-reaction       */
/**********************************************************/
double statStoreInitPopulationPhoton(const int targlev, const int targid, const double coef, Transmission *tin)
{
  int ig    = (int)(2.0*halfint(ncl[0].lev[0].spin));   // half-int ofset for CN
  int it    = (int)(2.0*ncl[targid].lev[targlev].spin); // target state spin
  int pt    = ncl[targid].lev[targlev].parity;          // target state parity
  double g  = coef /(it+1.0) / 2.0;                     // spin factor without (2J+1)
                                                        // factor of 1/2, PR 126, 671 (1961)

  /*** Loop over CN J and Parity*/
  int jdxmax = 0;
  for(int j0=ig ; j0<=ncl[0].jmax*2+ig ; j0+=2){
    int   jdx = (j0-ig)/2;
    if(jdx > jdxmax) jdxmax = jdx;
    for(int p0=-1 ; p0<=1 ; p0+=2){

      /*** For entrance channel Tlj coupled to JP, skip L=0 */
      for(int lp0=2 ; lp0<=tin->lmax*2 ; lp0+=2){

        if( std::abs(lp0-it)>j0 || j0>(lp0+it) ) continue;

        int pl = parity(lp0);
        for(int pg=-1 ; pg <=1 ; pg+=2){ // E: 1, M: -1

          if(pl*pg != p0*pt) continue;

          /*** Tj for photo reaction contains L=0(E,M,*), L=1(E,M,*), L=2(E,M,*) ... */
          int index = lp0/2*3 + ( (pg == -1) ? 1 : 0 );
          double tlj  = tin->tran[index];
          if(tlj == 0.0) continue;

          if(p0 == 1) ncl[0].pop[0][jdx].even += (j0+1.0)*g*tlj;
          else        ncl[0].pop[0][jdx].odd  += (j0+1.0)*g*tlj;
        }
      }
    }
  }

  /*** calculate sum */
  double s = 0.0;
  for(int j=0 ; j<=jdxmax ; j++) s += ncl[0].pop[0][j].even + ncl[0].pop[0][j].odd;

#ifdef DEBUG_POPCHECK
  std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4);
  for(int j=0 ; j<20 ; j++) {
    std::cout << std::setw(10) << j;
    std::cout << std::setw(12) << ncl[0].pop[0][j].even;
    std::cout << std::setw(12) << ncl[0].pop[0][j].odd;
    std::cout << std::endl;
  }
  std::cout << std::setw(12) << s << std::endl;
#endif

  return(s);
}


/**********************************************************/
/*      Adjust Initial Population to Match Reaction CS    */
/**********************************************************/
void statAdjustInitPopulation(const double sigr)
{
  if(sigr == 0.0) return;

  double f = 0.0;
  for(int j=0 ; j<=ncl[0].jmax ; j++) f += ncl[0].pop[0][j].even + ncl[0].pop[0][j].odd;
  if(f == 0.0) return;

  f = sigr/f;
  for(int j=0 ; j<=ncl[0].jmax ; j++){
    ncl[0].pop[0][j].even *= f;
    ncl[0].pop[0][j].odd  *= f;
  }
  f=0.0;
  for(int j=0;j<=ncl[0].jmax;j++)  f += ncl[0].pop[0][j].even + ncl[0].pop[0][j].odd;
} 

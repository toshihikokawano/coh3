/******************************************************************************/
/*  outprep.cpp                                                               */
/*        some pre-calculations for data output                               */
/******************************************************************************/

#include <iostream>

#include "structur.h"
#include "nucleus.h"
#include "output.h"

static CumulativeResidualProduct res;


/**********************************************************/
/*      Prepare and Print Cumulative Long-Lived States    */
/**********************************************************/
void outPrepTotalResidual(const int n, const int m)
{
  outPrepAllocateResidual(m);
  outPrepCumulativeResidual(n,1.0);
  outPrepPrintResidual();
  outPrepPrintIsomericRatio();
  outPrepFreeResidual();
}


/**********************************************************/
/*      Sum Produced Long-Lived States                    */
/**********************************************************/
void outPrepCumulativeResidual(const int n, const double fraction)
{
  if(!res.isAlloc()) return;

  /*** for all residual nuclides */
  for(int i=0 ; i<n ; i++){

    bool found = res.exist(ncl[i].za);

    /*** found a new residual */
    if(!found){
      /*** ground state production */
      int meta = 0;
      res.setData(ncl[i].za,meta,
                  ncl[i].lev[meta].energy,
                  ncl[i].lev[meta].spin,
                  ncl[i].lev[meta].parity,
                  ncl[i].lev[meta].halflife,
                  ncl[i].lpop[meta] * fraction);


      /*** scan discrete states */
      meta = 1;
      for(int k=1 ; k<ncl[i].ndisc ; k++){
        if(ncl[i].lev[k].halflife >= thalfmin){
          res.setData(ncl[i].za,meta,
                      ncl[i].lev[k].energy,
                      ncl[i].lev[k].spin,
                      ncl[i].lev[k].parity,
                      ncl[i].lev[k].halflife,
                      ncl[i].lpop[k] * fraction);
          meta++;
        }
      }
    }

    /*** if data already exist */
    else{
      for(int k=0 ; k<ncl[i].ndisc ; k++){
        /*** find the index of the level, and add production cross section */
        int idx = res.getIndex(ncl[i].za,ncl[i].lev[k].energy);
        if(idx >= 0) res.addProduction(idx,ncl[i].lpop[k] * fraction);
      }
    }
  }
}


/**********************************************************/
/*      Surrogate Call of Output Routines                 */
/**********************************************************/
void outPrepPrintResidual()
{
  res.setRange();
  outTotalResidual(&res);
}

void outPrepPrintIsomericRatio()
{
  res.setRange();
  outIsomericRatio(&res);
}


/**********************************************************/
/*      Utilities to Interface Local Object               */
/**********************************************************/
CumulativeResidualProduct *outPrepObjectPointer()
{
  return(&res);
}


/**********************************************************/
/*      Memory Allocate / Free                            */
/**********************************************************/
void outPrepAllocateResidual(const int m)
{
  res.memalloc(m * MAX_METASTABLE);
  res.clear();
}

void outPrepFreeResidual()
{
  res.memfree();
}


/******************************************************************************/
/*  HOsorteigenvalue.cpp                                                      */
/*        quick sort of eigenvalue                                            */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "HObasis.h"


/**********************************************************/
/*      Quick Sort of Eigenvalues                         */
/*      -------------------------                         */
/*      idx : original index before sorting               */
/*      ptr : block number that includes the eigenvalue   */
/*      ord : original order inside the block             */
/**********************************************************/
void HOSortEigenvalue(Basis *b, double *spe, int *idx, int *ptr, int *ord)
{
  for(int i=0 ; i<b->nv ; i++) idx[i] = i;

  for(int i=0 ; i<b->nv ; i++){
    int k = i;
    double x = spe[i];
    for(int j=i ; j<b->nv ; j++){
      if(x >= spe[j]){
        k = j;
        x = spe[j];
      }
    }
    int l  = idx[k];
    idx[k] = idx[i];
    idx[i] = l;
    spe[k] = spe[i];
    spe[i] = x;
  }

  for(int i=0 ; i<b->nv ; i++){
    for(int k=0 ; k<b->nb ; k++){

      /*** find original Omega block and its oder in the block */
      for(int j=0 ; j<b->ndata[k] ; j++){
        if(j+b->index[k] == idx[i]){
          ptr[i] = k;
          ord[i] = j;
          break;
        }
      }
    }
  }
/*
  for(int i=0 ; i<b->nv ; i++){
    cout << setw(5) << i << setw(5) << idx[i];
    cout << setw(5) << ptr[i] << setw(5) << ord[i];
    cout << setw(12) << spe[i] << setw(5) << b->omega[ptr[i]] << endl;
  }
*/
}


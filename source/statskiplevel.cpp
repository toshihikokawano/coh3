/******************************************************************************/
/*  statskiplevel.cpp                                                         */
/*        manipulate discrete levels                                          */
/******************************************************************************/

#include <iostream>

#include "structur.h"
#include "statmodel.h"


/**********************************************************/
/*      Remove Discrete Level in RIPL Data                */
/*      --------                                          */
/*      if the first column of the level data is '#',     */
/*      lev[i].flag = 1 : this level will be skipped      */
/*                                                        */
/**********************************************************/
int     statSkipDiscreteLevel(NuclearStructure *d)
{
  int skipid[MAX_LEVELS];
  
  /*** count special flags if given */
  int nskip = 0;
  for(int i=0 ; i<d->nlevel ; i++){
    if(d->lev[i].flag == 1){
      d->lev[i].ngamma = 0;
      skipid[nskip++] = i;
    }
  }
  if(nskip == 0) return(nskip);


  /*** see if there are gamma-lines decaying to skipped levels */
  for(int i=d->nlevel-1 ; i>0 ; i--){

    /*** each decay gamma */
    int    k = 0;
    double s = 0.0;
    for(int j=0 ; j<d->lev[i].ngamma ; j++){

      /*** compare with the skipped level */
      bool p = false;
      for(int n=0 ; n<nskip ; n++){
        if(d->lev[i].fstate[j] == skipid[n]){ p = true; break; }
      }
      if(!p){
        s += d->lev[i].branch[j];
        k ++;
      }
      
    }

    /*** when there was only one gamma-ray, force decaying to g.s. */
    if(k == 0){
      d->lev[i].ngamma    = 1;
      d->lev[i].fstate[0] = 0;
      d->lev[i].branch[0] = 1.0;
    }

    /*** if there are more than one gamma-ray but the sum is 0.0, re-distribute */
    else if(s == 0.0){
      s = 1.0 / k;
      for(int j=0 ; j<d->lev[i].ngamma ; j++){
        for(int n=0 ; n<nskip ; n++){
          if(d->lev[i].fstate[j] == skipid[n]) d->lev[i].branch[j] = 0.0;
          else d->lev[i].branch[j] = s;
        }
      }
    }

    /*** rescale the branching ratios */
    else{
      s = 1.0/s;
      for(int j=0 ; j<d->lev[i].ngamma ; j++){
        d->lev[i].branch[j] *= s;

        for(int n=0 ; n<nskip ; n++){
          if(d->lev[i].fstate[j] == skipid[n]) d->lev[i].branch[j] = 0.0;
        }
      }
    }
  }
  
  return(nskip);
}

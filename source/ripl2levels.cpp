/******************************************************************************/
/*  ripl2levels.cpp                                                           */
/*        read in discrete level data from file                               */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <cstring>
#include <cmath>

#include "structur.h"
#include "ripl2levels.h"
#include "terminate.h"

static void riplSelectCandidate (std::string, double *, int *);
static int  riplFixTransition (int *, double *, double *, const int, Level *);
static int  riplNormalizeBranch (int, double *, double *);
static int  riplSetNmax (const MaxLevelCtl, const int, const int, const int);
static int  riplFixNmax (int, Level *);
static int  riplGSonly (const unsigned int, const unsigned int, Level *);

static std::string datadir = DATADIR;


/**********************************************************/
/*      Read Discrete Level Data from RIPL Database       */
/**********************************************************/
int riplReadDiscreteLevels(ZAnumber *za, Level *lev, MaxLevelCtl ml)
{
  std::ostringstream os;
  std::ifstream      fp;
  std::string   str,file,d;
  int           nlevel = 0;
  int           nf[MAX_GAMMA_BRANCH];
  double        pe[MAX_GAMMA_BRANCH],ic[MAX_GAMMA_BRANCH];

  os << std::setw(4) << std::setfill('0') << za->getZ();
  file = os.str();
  file[0] = 'z';
  str = file + ".dat";

  /*** try current directry first */
  fp.open(&str[0]);

  /*** then local levels directory */
  if(!fp){
    str = LEVELDIRECTORY + file + ".dat";
    fp.open(&str[0]);
  }

  /*** then system data area */
  if(!fp){
    str = datadir + "/" + LEVELDIRECTORY + file + ".dat";
    fp.open(&str[0]);
  }

  if(!fp){
    message << "discrete level data file " << str << " not found, use ground state only";
    cohNotice("riplReadDiscreteLevels");

    nlevel = riplGSonly(za->getZ(), za->getA(), lev);
    return(nlevel);
  }

  for(int i=0 ; i<MAX_LEVELS ; i++) lev[i].energy = lev[i].spin = 0.0;

  bool found=false;
  int  a=0, z=0, nol=0, nc=0, nmax=0;
  while(getline(fp,str)){

    /*** search for Z and A entry in the file */
    d=str.substr( 5, 5);  a    = atoi(&d[0]);
    d=str.substr(10, 5);  z    = atoi(&d[0]);
    d=str.substr(15, 5);  nol  = atoi(&d[0]);
    d=str.substr(20, 5);//nog  = atoi(&d[0]);
    d=str.substr(25, 5);  nmax = atoi(&d[0]);
    d=str.substr(30, 5);  nc   = atoi(&d[0]);
    ZAnumber za1(z,a);

    if((za1.getZ() == za->getZ()) && (za1.getA() == za->getA())) found = true;
    /*
    if(nmax >= MAX_LEVELS){
      message << "discrete levels for Z " << za1.getZ() << " - A " << za1.getA() << " too many (" << nol << ")";
      fp.close();
      cohTerminateCode("riplReadDiscreteLevels");
    }
    */
    /*** set Nmax by the value of MaxLevelCtl */
    nlevel = riplSetNmax(ml, nol, nmax, nc);
    if(nlevel >= MAX_LEVELS) nlevel = MAX_LEVELS - 1;

    /*** for all discrete levels */
    double elev=0.0, s=0.0, thlf=0.0;
    int    p=0, ng=0, f=0;
    for(int i=0 ; i<nol ; i++){
      getline(fp,str);
      d=str.substr( 4,10);  elev = atof(&d[0]);
      d=str.substr(15, 5);  s    = atof(&d[0]);
      d=str.substr(20, 3);  p    = atoi(&d[0]);
      d=str.substr(24,10);  thlf = atof(&d[0]);
      d=str.substr(34, 3);  ng   = atoi(&d[0]);

      f = (str[0] == '#' ) ? 1 : 0;

      /*** for gamma-ray branches */
      for(int j=0 ; j<ng ; j++){
        getline(fp,str);
        if(j < MAX_GAMMA_BRANCH){
          d=str.substr(39, 4);  nf[j] = atoi(&d[0]) -1;
          d=str.substr(66,10);  pe[j] = atof(&d[0]);
          d=str.substr(77,10);  ic[j] = atof(&d[0]);
        }
      }
      if(ng >= MAX_GAMMA_BRANCH) ng = MAX_GAMMA_BRANCH - 1;

      /*** if found, copy data into lev array */
      if(found && i<MAX_LEVELS - 1){
        if( (s < 0.0) && (i == 0) ) riplSelectCandidate(str, &s, &p);

        lev[i].energy   = elev;
        lev[i].spin     = s;
        lev[i].halflife = thlf;
        lev[i].parity   = p;
        lev[i].flag     = f;

        if(i >= 1){
          if(ng == 0) ng = riplFixTransition(nf,pe,ic,i,lev);
          ng = riplNormalizeBranch(ng,pe,ic);
        }

        lev[i].ngamma   = ng;

        for(int j=0 ; j<ng ; j++){
          lev[i].fstate[j] = nf[j];
          lev[i].branch[j] = pe[j];
          lev[i].gratio[j] = 1.0/(1.0 + ic[j]); // Ig/(Ig + Ie) = 1/(1+ICC)
        }
      }
    }
    if(found) break;
  }
  fp.close();


  if(!found){
    if(ml != reassign){
      message << "discrete level data for Z " << za->getZ() << " - A " << za->getA() << " not found";
      cohTerminateCode("riplReadDiscreteLevels");
    }
    else{
      nlevel = riplGSonly(za->getZ(), za->getA(), lev);
    }
  }
  else{
    if( (lev[0].spin < 0.0) || (lev[0].parity == 0) ){
      if(ml == reassign) nlevel = riplGSonly(za->getZ(), za->getA(), lev);
      else{
        message << "discrete level data for Z " << za->getZ() << " - A " << za->getA() << " found but spin/parity not assigned";
        cohTerminateCode("riplReadDiscreteLevels");
      }
    }
  }

  if(ml == extended){
    nlevel = riplFixNmax(nlevel,lev);
    if(nlevel == 0){
      message << "discrete level data for Z " << za->getZ() << " - A " << za->getA() << " not complete";
      cohTerminateCode("riplReadDiscreteLevels");
    }
  }

  return(nlevel);
}


/**********************************************************/
/*      Get Spin and Parity from Candidates               */
/**********************************************************/
void riplSelectCandidate(std::string str, double *s, int *p)
{
  std::string csrc = str.substr(46,18);
  char   cdst[] = "                  ";

  int i0 = csrc.find("(");
  int i1 = csrc.find_first_of(",)");

  if(i0 == (signed int)std::string::npos) return;

  for(int i=i0+1 ; i<i1 ; i++) cdst[i-i0-1] = csrc[i];
  cdst[i1-i0-1] = '\0';

  *s = (strchr(cdst,'/') != nullptr) ?  0.5 : 1.0;
  *p = (strchr(cdst,'-') != nullptr) ? -1   : 1;

  i1 = strlen(cdst);
  for(int i=0 ; i<i1 ; i++){
    if( (cdst[i] == '+') || (cdst[i] == '-') || (cdst[i] == '/') ){
      cdst[i] = '\0';
      break;
    }
  }
  *s *= atof(cdst);
}


/**********************************************************/
/*      Fix Transition If No Decay                        */
/**********************************************************/
int riplFixTransition(int nf[], double pe[], double ic[], const int k, Level *lev)
{
  int n = 0;

  /*** if this is the first excited state */
  if(k == 1){
    nf[0] = 0;
    pe[0] = 1.0;
    ic[0] = 0.0;
    n = 1;
    return(n);
  }

  /*** first, look for states to which decay is possible by E1 */
  int    p = lev[k].parity;
  double s = lev[k].spin;

  int dm = 100;
  for(int i=0 ; i<k ; i++){
    if( (s == 0.0) && (lev[i].spin == 0.0) ) continue;
    int dj = (int)fabs(lev[i].spin - s);
    if(dj < dm) dm = dj;
    /*** E1, different parity and |J0 - J1| <= 1.0 */
    if( (lev[i].parity != p) && (dj <= 1) ){
      nf[n] = i;
      pe[n] = 1.0;
      ic[n] = 0.0;
      n++;
      if(n == MAX_GAMMA_BRANCH-1) break;
    }
  }
  /*** if no transition still,
       find the closest spin state with the lowest energy */
  if(n == 0){
    for(int i=0 ; i<k ; i++){
      if( (s == 0.0) && (lev[i].spin == 0.0) ) continue;
      int dj = (int)fabs(lev[i].spin - s);
      if(dj == dm){
        nf[0] = i;
        pe[0] = 1.0;
        ic[0] = 0.0;
        n = 1;
        break;
      }
    }
  }

  return(n);
}


/**********************************************************/
/*      Renormalize Transition Probability                */
/**********************************************************/
int riplNormalizeBranch(int ng, double pe[], double ic[])
{
  double s = 0.0;

  for(int i=0 ; i<ng ; i++) s += pe[i];

  if(s == 0.0){
    if(ng > 1){
      s = 1.0/(double)ng;
      for(int i=0 ; i<ng ; i++){
        pe[i] = s;
        ic[i] = 0.0;
      }
    }
    else{
      pe[0] = 1.0;  // decay 100% to the first line in the branches
      ic[0] = 0.0;
      ng = 1;
    }
  }else{
    s = 1.0/s;
    for(int i=0 ; i<ng ; i++) pe[i] *= s;
  }

  return(ng);
}


/**********************************************************/
/*      Determine Nmax from Input                         */
/**********************************************************/
int riplSetNmax(const MaxLevelCtl ml, const int nol, const int nmax, const int nc)
{
  int nlevel = 0;

  switch(ml){
  case normal:   // Nc is used for the highest level
    nlevel = nc;
    break;
  case extended: // Nmax is used, but upto unknown spin/parity state
  case reassign: // Nmax is used, and unknown level spin/parity will be re-assigned later
    nlevel = (nmax > nc) ? nmax : nc;
    break;
  case all:      // include all levels given
    nlevel = nol;
    break;
  default:
    nlevel = nc;
    break;
  }

  return(nlevel);
}


/**********************************************************/
/*      Renormalize Transition Probability                */
/**********************************************************/
int riplFixNmax(int nl, Level *lev)
{
  for(int i=0 ; i<nl ; i++){
    /*** spin or parity not assigned */
    if( (lev[i].spin < 0.0) || (lev[i].parity == 0) ){
      nl = i;
      break;
    }
  }

  return(nl);
}


/**********************************************************/
/*      Assume Ground State Spin                          */
/**********************************************************/
int riplGSonly(const unsigned int z, const unsigned int a, Level *lev)
{
  unsigned int n = a -z;

  /*** even-even */
  if(((z%2) == 0) && ((n%2) == 0) ){
    lev[0].spin   = 0.0;
    lev[0].parity = 1;
  }
  /*** odd-odd, 1+ assumed */
  else if(((z%2) != 0) && ((n%2) != 0) ){
    lev[0].spin   = 1.0;
    lev[0].parity = 1;
  }
  /*** even-odd, 1/2+ assumed */
  else{
    lev[0].spin   = 0.5;
    lev[0].parity = 1;
  }

  return(1);
}

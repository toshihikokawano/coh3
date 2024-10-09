/******************************************************************************/
/*  ensdf.cpp                                                                 */
/*        read ENSDF data file, extract beta branching ratio                  */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include "beoh.h"
#include "dir.h"

static std::string datadir = DATADIR;
#define N_RECORD 80

/***********************************************************/
/*      Read Beta Decay Branching Ratio from ENSDF         */
/***********************************************************/
int beohENSDFRead(ZAnumber *za, Beta *beta)
{
  std::string        str,file,line;
  double             ener = 0.0, bran = 0.0, norm = 1.0;
  std::ostringstream os;
  std::ifstream      fp;

  os << std::setw(3) << std::setfill('0') << za->getZ() << std::setw(3) << std::setfill('0') << za->getA();
  file = os.str();
  str = "ENSDF" + file + ".dat";

  /*** try current directry first */
  fp.open(&str[0]);

  /*** then system data area */
  if(!fp){
    str = datadir + "/" + ENSDFDIRECTORY + str;
    fp.open(&str[0]);
  }

  if(!fp){
    beta->nstate=0;
    return(beta->nstate);
  }


  for(int p=0 ; p<N_RECORD ; p++) line += ' ';
  line[N_RECORD] = '\0';

  int    n = 0;
  while(!fp.eof()){
    getline(fp,str);

    for(int p=0 ; p<N_RECORD ; line[p++] = ' ')
    strncpy(&line[0], &str[0], str.length());

    char   c0 = tolower(line[5]);
    char   c1 = tolower(line[6]);
    char   c2 = tolower(line[7]);

    /*** skip comment and history */
    if((c1 == 'c') || (c2 == 'h')) continue;

    /*** normalization record */
    if(c2 == 'n' && c0 == ' ' && c1 == ' '){
      norm = atof((line.substr(41,8)).c_str());
      if(norm == 0.0) norm = 1.0;

      continue;
    }
    /*** excitation energy, in MeV */
    else if(c2 == 'l' && c0 == ' ' && c1 == ' '){
      ener = atof((line.substr(9,10)).c_str()) * 0.001;
      continue;
    }
    /*** beta decay, assume its level is previously found */
    else if(c2 == 'b' && c0 == ' ' && c1 == ' '){
      bran = atof((line.substr(21,8)).c_str()) * 0.01;
      beta->br[n++].setVal(ener,bran);

    }
    if(n >= MAX_BETA_BRANCH){
      n--;
      break;
    }
  }
  fp.close();

  beta->nstate = n;

  /*** re-normalize the branching ratio */
  double sum = 0.0;
  if(norm != 0.0 || norm != 1.0){
    for(int i=0 ; i<beta->nstate ; i++){
      beta->br[i].setR(beta->br[i].getR() * norm);
      sum += beta->br[i].getR();
    }
  }

  /*** if branching ratios are not given, set zero all */
  if(sum==0.0 && beta->nstate>0){
    beta->nstate = 0;
//  for(int i=0 ; i<beta->nstate ; i++) beta->br[i].setR(1.0/beta->nstate);
  }

  /*** in some cases, sum of branching ratio exceed 1.0, because
       the ratios are sometimes given as maximum, like <10%.
       we renormalize the all strength, if this is the case */
  if(sum>1.0){
    for(int i=0 ; i<beta->nstate ; i++) beta->br[i].scaleR(1.0/sum);
  }
/*
  cout << "# ENSDF Sum " << sum << "  " << beta->nstate << endl;
  for(int i=0 ; i<beta->nstate ; i++){
    cout << "#EN "
         << setprecision(4) << setiosflags(ios::scientific) << setw(11)
         << beta->br[i].getE() << setw(11) << beta->br[i].getR() << endl;
  }
*/
  return(beta->nstate);
}


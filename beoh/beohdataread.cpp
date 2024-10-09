/******************************************************************************/
/*  dataread.cpp                                                              */
/*        functions to read input parameters                                  */
/******************************************************************************/

#include <iostream>

#include "beoh.h"
#include "optical.h"
#include "parameter.h"
#include "elements.h"
#include "terminate.h"
#include "global.h"
#include "masstable.h"
#include "dir.h"

static const unsigned int W_HEAD = 10;  // width of keyword

static char *   readExtractData      (const int);
static int      readFgetOneline      (char *);
static Particle readParticleIdentify (const int);

static inline int  readElementToZ    (char *);
static inline int  readGetFbrIndex   (int, unsigned int, unsigned int, Fission *);
#ifndef BeoH
static inline void readParSet1       (const parameterType);
#endif
static inline void readParSet2       (const parameterType);
static inline void readParSet3       (const parameterType);
static inline void readParSet4       (const parameterType);

static std::string head;
static char   inputdata[256];

/**********************************************************/
/*      Read Headline                                     */
/**********************************************************/
static const int NWORD_SEC = 8;
static const std::string table_sec[NWORD_SEC]=
  {"BEGIN    :","DATA     :","BETA     :","HFS      :",
   "FFRAG    :","MCF      :","CFY      :",
   "END      :"};

int readHead(char *s)
{
  if( readFgetOneline(s) < 0 ) return(-1);

  for(int j=0 ; j<NWORD_SEC ; j++) if(head==table_sec[j]) return(j);
  message << "unknown key-word [" << head << "] in HEAD";
  cohTerminateCode("readHead");
  return(-1);
}


/**********************************************************/
/*      Read System Parameter                             */
/**********************************************************/
int readSystem(char *s, System *sys, Pdata *pdt, FFragData *fdt, BData *bdt)
{
  while(1){
    if( readFgetOneline(s) < 0 ) break;

    if(head == "ENDDATA  :"){
      break;

    }else if(head == "precursor:" || head == "compound :" || head == "target   :"){
      int z = readElementToZ(readExtractData(0));
      int a = atoi(          readExtractData(1));

      if(head == "precursor:"){
        sys->target.setZA(z,a);
        bdt->setMode(betadecay);
      }
      else if(head == "compound :"){
        sys->compound.setZA(z,a);
        bdt->setMode(statdecay);
      }
      else{
        sys->target.setZA(z,a);
        bdt->setMode(fissiondecay);
      }

    }else if(head == "exenergy :"){
      bdt->exmean  = atof(readExtractData(0));
      bdt->exwidth = atof(readExtractData(1));

    }else if(head == "qvalue   :" || head == "energy   :"){
      bdt->energy  = atof(readExtractData(0));

    }else if(head == "incident :"){
      sys->incident.pid = readParticleIdentify(0);

    }else if(head == "isomer   :"){
      sys->target_level = atoi(readExtractData(0));

    }else if(head == "initspin :"){
      bdt->initspin = atof(readExtractData(0));
      readExtractData(1);
      if(inputdata[0] == '+') bdt->initparity = 1;
      else if(inputdata[0] == '-') bdt->initparity = -1;
      else bdt->initparity = -2;

    }else if(head == "spinfact :"){
      bdt->spinfactor = atof(readExtractData(0));

    }else if(head == "bin      :"){
      sys->energy_bin_in = atof(readExtractData(0));

    }else if(head == "omp_n    :"){
      readExtractData(0);
      pdt[neutron ].omp = find_omp(inputdata);

    }else if(head == "omp_p    :"){
      if(MAX_CHANNEL <= 2){
        message << "proton not handled as the max channel is " << MAX_CHANNEL;
        cohTerminateCode("readSystem");
      }
      readExtractData(0);
      pdt[proton  ].omp = find_omp(inputdata);

    }else if(head == "omp_a    :"){
      if(MAX_CHANNEL <= 3){
        message << "helium not handled as the max channel is " << MAX_CHANNEL;
        cohTerminateCode("readSystem");
      }
      readExtractData(0);
      pdt[alpha   ].omp = find_omp(inputdata);

    }else if(head == "omp_d    :"){
      if(MAX_CHANNEL <= 4 ){
        message << "deuteron not handled as the max channel is " << MAX_CHANNEL;
        cohTerminateCode("readSystem");
      }
      readExtractData(0);
      pdt[deuteron].omp = find_omp(inputdata);

    }else if(head == "omp_t    :"){
      if(MAX_CHANNEL <= 5){
        message << "triton not handled as the max channel is " << MAX_CHANNEL;
        cohTerminateCode("readSystem");
      }
      readExtractData(0);
      pdt[triton  ].omp = find_omp(inputdata);

    }else if(head == "omp_h    :"){
      if(MAX_CHANNEL <= 6){
        message << "helium3 not handled as the max channel is " << MAX_CHANNEL;
        cohTerminateCode("readSystem");
      }
      readExtractData(0);
      pdt[helion  ].omp = find_omp(inputdata);

    }else if(head == "deform   :"){
      sys->beta2 = atof(readExtractData(0));

    }else if(head == "partmax  :"){
      sys->max_particle = atoi(readExtractData(0));

    }else if(head == "tweakOV  :"){
      readParSet2(parmOV);

    }else if(head == "tweakOW  :"){
      readParSet2(parmOW);

    }else if(head == "tweakORV :"){
      readParSet2(parmORV);

    }else if(head == "tweakORW :"){
      readParSet2(parmORW);

    }else if(head == "tweakOAV :"){
      readParSet2(parmOAV);

    }else if(head == "tweakOAW :"){
      readParSet2(parmOAW);

    }else if(head == "tweakOC  :"){
      readParSet2(parmOC);

    }else if(head == "massfile :"){
      massReadFile(readExtractData(0));

    }else if(head == "select_z :"){
      int z = readElementToZ(readExtractData(0));
      fdt->selectZ = z;

    }else if(head == "select_a :"){
      fdt->selectA = atoi(readExtractData(0));

    }else if(head == "option   :"){
      std::string item = readExtractData(0);
      if(item == "internalconversion"){
        opt.internalconversion = true;
      }
      else if(item == "finegammaspectrum"){
        opt.finegammaspectrum = true;
      }
      else{
        message << "no such option [" << item << "] in DATA";
        cohTerminateCode("readSystem");
      }

    }else{
      message << "unknown key-word [" << head << "] in DATA";
      cohTerminateCode("readSystem");
    }
  }

  /*** set spontaneous fission flag when no incident is given */
  fdt->spontaneous = (sys->incident.pid == unknown) ? true : false;

  return(0);
}


/**********************************************************/
/*      Read Beta Parameters                              */
/**********************************************************/
int readBeta(char *s, Beta *b)
{
  double d1,d2;

  b->nstate = 0;
  while(1){
    if( readFgetOneline(s) < 0 ) break;

    if(head == "ENDBETA  :"){
      return(b->nstate);

    }else if(head == "gtstate  :"){
      if(b->nstate < MAX_BETA_BRANCH){
        d1 = atof(readExtractData(0));
        d2 = atof(readExtractData(1));
        b->br[b->nstate++].setVal(d1,d2);
      }

    }else if(head == "width    :"){
      b->width = atof(readExtractData(0));

    }else{
      message << "unknown key-word [" << head << "] in BETA";
      cohTerminateCode("readBeat");
    }
  }

  return(b->nstate);
}


/**********************************************************/
/*      Read Statistical Model Parameters                 */
/**********************************************************/
int readStatModel(char *s, GDR *gdr, Fission *fbr)
{
  unsigned int z,a;
  int i = -1, i2 = 0, j = 0, m;

  while(1){
    if( readFgetOneline(s) < 0 ) break;

    if(head == "ENDHFS   :"){
      return(0);

    }else if(head == "gdr      :"){
      if(j < MAX_GDR){
        readExtractData(0);
        if(tolower(inputdata[0]) != 'e' && tolower(inputdata[0]) != 'm'){
          message << "gamma-ray emission type " << inputdata[0] << " not defined";
          cohTerminateCode("readStatModel");
        }
        if(inputdata[1] != '1' && inputdata[1] != '2' && inputdata[1] != '3'){
          message << "gamma-ray multipolarity " << inputdata[1] << " should be 1, 2, or 3";
          cohTerminateCode("readStatModel");
        }

        gdr[j].setXL( inputdata );
        gdr[j].setEnergy ( atof(readExtractData(1)) );
        gdr[j].setWidth  ( atof(readExtractData(2)) );
        gdr[j].setSigma  ( atof(readExtractData(3)) );
        j++;
      }

    }else if(head == "barrier  :"){
      z = readElementToZ(readExtractData(0));
      a = atoi(          readExtractData(1));
      m = atoi(          readExtractData(2));

      if( (0 <= m) && (m < MAX_HUMP) && (i < MAX_FISS_CHANCE) ){
        i2 = readGetFbrIndex(i,z,a,fbr);
        if((i2 < 0) && (i < MAX_FISS_CHANCE-1)) i++;
        fbr[i].za.setZA(z,a);
        fbr[i].barrier[m].height       = atof(readExtractData(3));
        fbr[i].barrier[m].curvature    = atof(readExtractData(4));
        fbr[i].barrier[m].inertia      = atof(readExtractData(5));
      }

    }else if(head == "kband    :"){
      z = readElementToZ(readExtractData(0));
      a = atoi(          readExtractData(1));
      m = atoi(          readExtractData(2));
      if( (0 <= m) && (m < MAX_HUMP) ){
        i2 = readGetFbrIndex(i,z,a,fbr);
        if(i2 >= 0){
          if(fbr[i2].barrier[m].nband < MAX_KBAND){
            int k = fbr[i2].barrier[m].nband;
            fbr[i2].barrier[m].kband[k].excitation = atof(readExtractData(3));
            fbr[i2].barrier[m].kband[k].k2         = (int)(atof(readExtractData(4))*2.0);
            fbr[i2].barrier[m].kband[k].parity     = atoi(readExtractData(5));
            fbr[i2].barrier[m].nband++;
          }
        }
        else{
          message << "barrier parameter required before kband";
          cohTerminateCode("readStatModel");
        }
      }

    }else if(head == "levelcut :"){
      readParSet3(parmLCUT);

    }else if(head == "tweakTJ  :"){
      readParSet4(parmTJ);

    }else if(head == "tweakLD  :"){
      readParSet3(parmLD);

    }else if(head == "tweakPAIR:"){
      readParSet3(parmPAIR);

    }else if(head == "tweakFL1 :"){
      readParSet3(parmFL1);

    }else if(head == "tweakFL2 :"){
      readParSet3(parmFL2);

    }else if(head == "tweakFL3 :"){
      readParSet3(parmFL3);

    }else{
      message << "unknown key-word [" << head << "] in HFS";
      cohTerminateCode("readStatModel");
    }
  }
  return(0);
}


/**********************************************************/
/*      Read Prompt Fission Model Parameters              */
/**********************************************************/
int readFissionFragment(char *s, FFragData *fdt, BData *bdt)
{
  while(1){
    if( readFgetOneline(s) < 0 ) break;

    if(head == "ENDFFRAG :"){
      return(0);

    }else if(head == "massreso :"){
      fdt->massresolution = atof(readExtractData(0));

    }else if(head == "maxtemp  :"){
      fdt->maxtemperature = atof(readExtractData(0));

    }else if(head == "ycutoff  :"){
      fdt->ycutoff = atof(readExtractData(0));

    }else if(head == "mcffile  :"){
      fdt->mcffile = readExtractData(0);

    }else if(head == "chicalc  :"){
      bdt->setMode(fissionspec);

    }else if(head == "YAfile   :"){
      fdt->yafile = readExtractData(0);

    }else if(head == "YTKEfile :"){
      fdt->ytkefile = readExtractData(0);

    }else if(head == "RTAfile  :"){
      fdt->rtafile = readExtractData(0);

    }else{
      message << "unknown key-word [" << head << "] in FFRAG";
      cohTerminateCode("readFissionFragment");
    }
  }

  return(0);
}


int readMultiChanceFission(char *s, MultiChanceFissionData *mc)
{
  while(1){
    if( readFgetOneline(s) < 0 ) break;

    if(head == "ENDMCF   :"){
      return(0);

    }else if(head == "fraction :"){
      mc->fraction = atof(readExtractData(0));

    }else if(head == "spinfact :"){
      mc->spinfactor  = atof(readExtractData(0)); // given in the main section too
      mc->spinfactor1 = atof(readExtractData(1));

    }else if(head == "rt       :"){
      mc->rt  = atof(readExtractData(0));
      mc->rt1 = atof(readExtractData(1));

    }else if(head == "tke      :"){
      mc->tke  = atof(readExtractData(0));
      mc->tke1 = atof(readExtractData(1));

    }else if(head == "exfis    :"){
      mc->exfis = atof(readExtractData(0));

    }else if(head == "eprefis  :"){
      mc->eprefis = atof(readExtractData(0));

    }else if(head == "gauss1   :"){
      mc->GaussFract[0] = atof(readExtractData(0));
      mc->GaussSigma[0] = atof(readExtractData(1));
      mc->GaussDelta[0] = atof(readExtractData(2));
      mc->GaussFract[4] = atof(readExtractData(3));
      mc->GaussSigma[4] = atof(readExtractData(4));
      mc->GaussDelta[4] = atof(readExtractData(5));

      mc->ffydata = "input";

    }else if(head == "gauss2   :"){
      mc->GaussFract[1] = atof(readExtractData(0));
      mc->GaussSigma[1] = atof(readExtractData(1));
      mc->GaussDelta[1] = atof(readExtractData(2));
      mc->GaussFract[5] = atof(readExtractData(3));
      mc->GaussSigma[5] = atof(readExtractData(4));
      mc->GaussDelta[5] = atof(readExtractData(5));

    }else if(head == "gauss3   :"){
      mc->GaussFract[2] = atof(readExtractData(0));
      mc->GaussSigma[2] = atof(readExtractData(1));
      mc->GaussDelta[2] = atof(readExtractData(2));
      mc->GaussFract[6] = atof(readExtractData(3));
      mc->GaussSigma[6] = atof(readExtractData(4));
      mc->GaussDelta[6] = atof(readExtractData(5));

    }else if(head == "gausssym :"){
      mc->GaussFract[3] = atof(readExtractData(0));
      mc->GaussSigma[3] = atof(readExtractData(1));
      mc->GaussFract[7] = atof(readExtractData(2));
      mc->GaussSigma[7] = atof(readExtractData(3));

    }else if(head == "finetune :"){
      int    a = atoi(readExtractData(0));
      double f = atof(readExtractData(1));
      if(!mc->setFTune(a,f)){
        message << "too many fine tune mass points";
        cohTerminateCode("readMultiChanceFission");
      }

    }else if(head == "ffydata  :"){
      mc->ffydata = readExtractData(0);

    }else if(head == "zpfactor :"){
      mc->ZpFactor[0] = atof(readExtractData(0));
      mc->ZpFactor[1] = atof(readExtractData(1));
      mc->ZpFactor[2] = atof(readExtractData(2));
      mc->ZpFactor[3] = atof(readExtractData(3));

    }else{
      message << "unknown key-word [" << head << "] in MCF";
      cohTerminateCode("readMultiChanceFission");
    }
  }

  return(0);
}


/**********************************************************/
/*      Read Cumulative Yield Data                        */
/**********************************************************/
int readFissionYield(char *s, FFragData *f)
{
  while(1){
    if( readFgetOneline(s) < 0 ) break;

    if(head == "ENDCFY   :"){
      return(0);
    }else if(head == "longlived:"){
      f->maxhalflife = atof(readExtractData(0));

    }else if(head == "decayfile:"){
      f->branchdatafile = readExtractData(0);

    }else{
      message << "unknown key-word [" << head << "] in CFY";
      cohTerminateCode("readFissionYield");
    }
  }

  return(0);
}


/**********************************************************/
/*      Skip comment-line and read 1-line                 */
/**********************************************************/
static std::string line;

int readFgetOneline(char *s)
{
  /*** read one line, skip comment or blank line */
  while(1){
    std::cin.getline(s,255);
    if(std::cin.eof() != 0) return(-1);
    if(s[0] == '#') continue;
    if(strlen(s) <= 1) continue;
    break;
  }

  /*** copy data line to str, extract line head */
  line = (std::string)s;
  head = line.substr(0, W_HEAD);
  if(head.size() < W_HEAD){
    message << "key-word [" << head << "] too short";
    cohTerminateCode("readFgetOneline");
  }

  return(1);
}


/**********************************************************/
/*      Extract Data from Line                            */
/**********************************************************/
char * readExtractData(const int n)
{
  inputdata[0] = '\0';

  if( strlen(line.c_str())>W_HEAD ){
    std::string work = line;
    char *str = &work[W_HEAD];
    char *tok = nullptr;
    const char *delim = " ";

    tok = strtok(str,delim);

    if(tok!=nullptr){
      for(int i=1 ; i<=n ; i++){
        tok = strtok(nullptr,delim);
        if(tok==nullptr) break;
      }
    }
    if(tok!=nullptr) strncpy(&inputdata[0],tok,sizeof(inputdata)-1);
  }

  return(&inputdata[0]);
}


/**********************************************************/
/*      Convert Element Name into Z Number                */
/**********************************************************/
int readElementToZ(char *d)
{
  int z = element_getZ(d);

  if(z<=0){
    message << "element " << d << " not found";
    cohTerminateCode("readElementToZ");
  }

  return(z);
}


/**********************************************************/
/*      Convert Particle Index into ID                    */
/**********************************************************/
Particle readParticleIdentify(const int n)
{
  readExtractData(n);
  char i = tolower(inputdata[0]);

  Particle p = unknown;
  switch(i){
  case 'g' : p = gammaray; break;
  case 'n' : p = neutron ; break;
  case 'p' : p = proton  ; break;
  case 'a' : p = alpha   ; break;
  case 'd' : p = deuteron; break;
  case 't' : p = triton  ; break;
  case 'h' : p = helion  ; break;
  case 'f' : p = fission ; break;
  default  : break;
  }
  if(p == unknown){
    std::ostringstream os;
    os << "unknown particle " << i;
    cohTerminateCode(os.str());
  }
  return(p);
}


/**********************************************************/
/*      Get Nucleus Index from Z and A                    */
/**********************************************************/
inline int readGetFbrIndex(int n, unsigned int z, unsigned int a, Fission *fbr)
{
  ZAnumber id(z,a);
  int idx = -1;
  for(int i=0 ; i<=n ; i++){
    if( fbr[i].za == id ) idx = i;
  }
  return(idx);
}


/**********************************************************/
/*      Adjustable Parameter Setting                      */
/**********************************************************/
#ifndef BeoH
inline void readParSet1(const parameterType pt)
{
  double   f = atof(readExtractData(0));
  double   w = atof(readExtractData(1));
  if(!adj.addparm(pt,f,w)){
    message << "too many adjustable parameters";
    cohTerminateCode("readParSet1");
  }
}
#endif
inline void readParSet2(const parameterType pt)
{
  Particle p = readParticleIdentify(0);
  double   f = atof(readExtractData(1));
  double   w = atof(readExtractData(2));
  if(!adj.addparm(pt,p,f,w)){
    message << "too many adjustable parameters";
    cohTerminateCode("readParSet2");
  }
}

inline void readParSet3(const parameterType pt)
{
  unsigned int z,a;
  z = readElementToZ(readExtractData(0));
  a = atoi(          readExtractData(1));
  ZAnumber id(z,a);
  double   f = atof(readExtractData(2));
  double   w = atof(readExtractData(3));
  if(!adj.addparm(pt,id,f,w)){
    message << "too many adjustable parameters";
    cohTerminateCode("readParSet3");
  }
}

inline void readParSet4(const parameterType pt)
{ 
  unsigned int z,a;
  z = readElementToZ(readExtractData(0));
  a = atoi(          readExtractData(1));
  ZAnumber id(z,a);
  Particle p = readParticleIdentify(2);
  double   f = atof(readExtractData(3));
  double   w = atof(readExtractData(4));
  if(!adj.addparm(pt,id,p,f,w)){
    message << "too many adjustable parameters";
    cohTerminateCode("readParSet4");
  }
}


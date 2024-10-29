/******************************************************************************/
/*  dataread.cpp                                                              */
/*        functions to read input parameters                                  */
/******************************************************************************/

#include <iostream>

#include "structur.h"
#include "coh.h"
#include "optical.h"
#include "parameter.h"
#include "elements.h"
#include "terminate.h"
#include "global.h"
#include "masstable.h"

static const unsigned int W_HEAD = 10;  // width of keyword

static char *   readExtractData (const int);
static int      readFgetOneline (char *);
static Particle readParticleIdentify (const int);

static inline int  readElementToZ (char *);
static inline DirType readDirect (char *);
static inline int  readGetFbrIndex (const int, const unsigned int, const unsigned int, Fission *);
static inline void readParSet1 (const parameterType);
static inline void readParSet2 (const parameterType);
static inline void readParSet3 (const parameterType);
static inline void readParSet4 (const parameterType);

static std::string head;
static char inputdata[256];

/**********************************************************/
/*      Read Headline                                     */
/**********************************************************/
static const int NWORD_SEC = 9;
static const std::string table_sec[NWORD_SEC]=
  {"BEGIN    :","DATA     :","DWBA     :","DSD      :","PREEQ    :",
   "HFS      :","FNS      :","MFT      :",
   "END      :"};

int readHead(char *s)
{
  if( readFgetOneline(s) < 0 ) return(-1);
  for(int j=0 ; j<NWORD_SEC ; j++) if(head == table_sec[j]) return( j);
  message << "unknown key-word [" << head << "] in HEAD";
  cohTerminateCode("readHead");
  return(-1);
}


/**********************************************************/
/*      Read System Parameter                             */
/**********************************************************/
int readSystem(char *s, System *sys, Pdata *pdt, Direct *dir, double *ein, double *mkt)
{
  int klev = 0;
  int kein = 0;
  int kmkt = 0;
  unsigned char inputcheck = 0x00; // 1: target, 2: energy, 4: incident

  for(int i=0 ; i<MAX_LAMBDA ; i++) dir->defdynamic[i] = dir->defstatic[i] = 0.0;
  for(int i=0 ; i<MAX_DIRECT ; i++) dir->lev[i].spin = dir->lev[i].energy = 0.0;

  while(1){
    if( readFgetOneline(s) < 0 ){
      message << "data structure error in DATA section";
      cohTerminateCode("readSystem");
    }

    if(head == "ENDDATA  :"){
      break;

    }else if(head == "target   :"){
      int z = readElementToZ(readExtractData(0));
      int a = atoi(          readExtractData(1));
      sys->target.setZA(z,a);
      inputcheck = inputcheck | 0x01;

    }else if(head == "incident :"){
      Particle p = readParticleIdentify(0);
      if(     p == gammaray) sys->incident.za.setZA(0,0);
      else if(p == neutron ) sys->incident.za.setZA(0,1);
      else if(p == proton  ) sys->incident.za.setZA(1,1);
      else if(p == alpha   ) sys->incident.za.setZA(2,4);
      else if(p == deuteron) sys->incident.za.setZA(1,2);
      else if(p == triton  ) sys->incident.za.setZA(1,3);
      else if(p == helion  ) sys->incident.za.setZA(2,3);
      else if(p == userdef ) sys->incident.za.setZA(3,0);
      else                   sys->incident.za.setZA(0,1); // neutron assumed;
      inputcheck = inputcheck | 0x04;

    }else if(head == "energy   :"){
      ein[kein++] = atof(readExtractData(0));
      if(kein >= MAX_EINCIDENT){
        message << "too many incident energies " << kein + 1;
        cohTerminateCode("readSystem");
      }
      inputcheck = inputcheck | 0x02;

    }else if(head == "macs     :" || head == "t9       :"){
      double x =  atof(readExtractData(0));
      mkt[kmkt++] = (head == "t9       :") ? x * 1e+3 * BOLTZMANN : x;
      ctl.macs = true;
      if(kmkt >= MAX_MACSTEMP){
        message << "too many Maxwellian temperatures " << kmkt + 1;
        cohTerminateCode("readSystem");
      }
      inputcheck = inputcheck | 0x02;

    }else if(head == "excite   :"){
      sys->excitation = atof(readExtractData(0));

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

    }else if(head == "omp_x    :"){
      if(MAX_CHANNEL <= 6){
        message << "user-defined particle not handled as the max channel is " << MAX_CHANNEL;
        cohTerminateCode("readSystem");
      }
      readExtractData(0);
      pdt[userdef].omp = find_omp(inputdata);

    }else if(head == "level    :"){
      if(klev == MAX_DIRECT){
        message << "too many coupled levels " << klev + 1;
        cohTerminateCode("readSystem");
      }
      dir->lev[klev].energy = atof(readExtractData(0));
      dir->lev[klev].spin   = atof(readExtractData(1));
      dir->type[klev++]     = readDirect(readExtractData(2));

    }else if(head == "deform   :" || head == "beta_rot :"){
      for(int l=0 ; l<MAX_LAMBDA ; l++){
        dir->defstatic[l] = atof(readExtractData(l));
        if(dir->defstatic[l] == 0.0) break;
      }
      if(dir->defstatic[0] > 0.0) sys->beta2 = dir->defstatic[0];

    }else if(head == "beta_vib :"){
      for(int l=0 ; l<MAX_LAMBDA ; l++){
        dir->defdynamic[l] = atof(readExtractData(l));
        if(dir->defdynamic[l] == 0.0) break;
      }

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

    }else if(head == "egrid    :"){
      readParSet1(parmESET);
      inputcheck = inputcheck | 0x02;

    }else if(head == "massfile :"){
      massReadFile(readExtractData(0));

    }else if(head == "option   :"){
      std::string item = readExtractData(0);
      if(item == "ewtransformation"){
        opt.ewtransformation = true;
      }
      else if(item == "internalconversion"){
        opt.internalconversion = true;
      }
      else if(item == "readdensity"){
        opt.readdensity = true;
      }
      else if(item == "readphotoabsorption"){
        opt.readphotoabsorption = true;
      }
      else if(item == "groundstateomp"){
        opt.groundstateomp = true;
      }
      else if(item == "chargediscrete"){
        opt.chargediscrete = true;
      }
      else if(item == "finegammaspectrum"){
        opt.finegammaspectrum = true;
      }
      else if(item == "continuumangdist"){
        opt.continuumangdist = true;
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

  /*** minimun input check, at least target is required */
  if(!(inputcheck & 0x01)){
    message << "target nucleus not given";
    cohTerminateCode("readSystem");
  }

  /*** exclude target only case, which is for mean-field calculation */
  if(inputcheck != 0x01){
    if(!(inputcheck & 0x02)){
      message << "calculation energy not given";
      cohTerminateCode("readSystem");
    }
    if(!(inputcheck & 0x04)){
      message << "incident particle not given";
      cohTerminateCode("readSystem");
    }
  }

  if(klev > 0) message << "number of direct levels " << klev << " ";
  if(kein > 0) message << "number of incident energies " << kein << " ";
  if(kmkt > 0) message << "number of Maxwellian temperature " << kmkt;
  if(message.str() != "") cohNotice("readSystem");

  return(klev);
}


/**********************************************************/
/*      Read DWBA  Parameter                              */
/**********************************************************/
int readDwba(char *s, Direct *dir)
{
  int klev = dir->ncc;

  for(int i=klev ; i<MAX_DIRECT ; i++) dir->lev[i].energy = 0.0;

  while(1){
    if( readFgetOneline(s) < 0 ){
      message << "data structure error in DWBA section";
      cohTerminateCode("readDwba");
    }

    if(head == "ENDDWBA  :"){
      break;

    }else if(head == "level    :"){
      if(klev == MAX_DIRECT){
        message << "too many DWBA levels " << klev + 1;
        cohTerminateCode("readDWBA");
      }
      dir->lev[klev].energy = atof(readExtractData(0));
      dir->lev[klev].spin   = fabs(atof(readExtractData(1)));
      dir->defdwba[klev]    = atof(readExtractData(2));
      dir->type[klev]       = dwba;
      klev++;

    }else if(head == "nonloc   :"){
      dir->nonloc = atof(readExtractData(0));

    }else{
      message << "unknown key-word [" << head << "] in DWBA";
      cohTerminateCode("readDwba");
    }
  }

  if(klev > dir->ncc) message << "number of DWBA levels " << klev - dir->ncc;
  if(message.str() != "") cohNotice("readDWBA");

  return(klev);
}


/**********************************************************/
/*      Read DSD  Parameter                               */
/**********************************************************/
int readDsd(char *s, Dcapt *dsd)
{
  while(1){
    if( readFgetOneline(s) < 0 ){
      message << "data structure error in DSD section";
      cohTerminateCode("readDsd");
    }
    
    if(head == "ENDDSD   :"){
      break;

    }else if(head == "interact :"){
      dsd->v1     =  atof(readExtractData(0));
      dsd->w1     =  atof(readExtractData(1));

    }else if(head == "nonloc   :"){
      dsd->nonloc =  atof(readExtractData(0));

    }else if(head == "boundmax :"){
      dsd->bmax   =  atof(readExtractData(0));

    }else if(head == "efshift  :"){
      dsd->fshift =  atof(readExtractData(0));

    }else if(head == "tweakDSDV:"){
      readParSet1(parmDSDV);

    }else{
      message << "unknown key-word [" << head << "] in DSD";
      cohTerminateCode("readDsd");
    }
  }

  return(0);
}


/**********************************************************/
/*      Read Statistical Model Parameter                  */
/**********************************************************/
int readStatModel(char *s, GDR *gdr, Fission *fbr)
{
  unsigned int z,a;
  int m;
  int ngdr =  0; // count number of GDR
  int nfis = -1; // count number of fissioning nuclei

  while(1){
    if( readFgetOneline(s) < 0 ){
      message << "data structure error in HFS section";
      cohTerminateCode("readStatModel");
    }

    if(head == "ENDHFS   :"){
      break;

    }else if(head == "photo    :"){
      readParSet1(parmGSTR);

    }else if(head == "timecut  :"){
      readParSet1(parmISOM);

    }else if(head == "d0search :"){
      readParSet1(parmD0SR);

    }else if(head == "gdr      :"){
      if(ngdr == MAX_GDR){
        message << "too many GDRs " << ngdr + 1;
        cohTerminateCode("readStatModel");
      }

      readExtractData(0);
      if(tolower(inputdata[0]) != 'e' && tolower(inputdata[0]) != 'm'){
        message << "gamma-ray emission type " << inputdata[0] << " not defined";
        cohTerminateCode("readStatModel");
      }
      if(inputdata[1] != '1' && inputdata[1] != '2' && inputdata[1] != '3'){
        message << "gamma-ray multipolarity " << inputdata[1] << " should be 1, 2, or 3";
        cohTerminateCode("readStatModel");
      }

      gdr[ngdr].setXL( inputdata );
      gdr[ngdr].setEnergy ( atof(readExtractData(1)) );
      gdr[ngdr].setWidth  ( atof(readExtractData(2)) );
      gdr[ngdr].setSigma  ( atof(readExtractData(3)) );
      std::string p = (std::string)readExtractData(4);
      if(p == "SL" || p == "") gdr[ngdr].setProfile(SL);
      else if(p == "GL")       gdr[ngdr].setProfile(GL);
      else if(p == "ML")       gdr[ngdr].setProfile(ML);
      else{
        message << "gamma-ray profile " << p << " not defined";
        cohTerminateCode("readStatModel");
      }

      /*** for backward compatibility */
      if((ngdr <= 1) && (gdr[ngdr].getXL() == "E1") && (p == "")) gdr[ngdr].setProfile(GL);

      ngdr++;

    }else if(head == "barrier  :"){
      z = readElementToZ(readExtractData(0));
      a = atoi(          readExtractData(1));
      m = atoi(          readExtractData(2));

      if( (0 <= m) && (m < MAX_HUMP) && (nfis < MAX_FISS_CHANCE) ){
        int i2 = readGetFbrIndex(nfis,z,a,fbr);
        if((i2 < 0) && (nfis < MAX_FISS_CHANCE-1)) nfis++;

        fbr[nfis].za.setZA(z,a);
        fbr[nfis].barrier[m].height       = atof(readExtractData(3));
        fbr[nfis].barrier[m].curvature    = atof(readExtractData(4));
        fbr[nfis].barrier[m].inertia      = atof(readExtractData(5));
        fbr[nfis].barrier[m].elmax        = atof(readExtractData(6));
      }

    }else if(head == "kband    :"){
      z = readElementToZ(readExtractData(0));
      a = atoi(          readExtractData(1));
      m = atoi(          readExtractData(2));
      if( (0 <= m) && (m < MAX_HUMP) ){
        int i2 = readGetFbrIndex(nfis,z,a,fbr);
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

    }else if(head == "potwell  :"){
      z = readElementToZ(readExtractData(0));
      a = atoi(          readExtractData(1));
      m = atoi(          readExtractData(2));

      if( (0 <= m) && (m < MAX_HUMP-1) && (nfis < MAX_FISS_CHANCE) ){
        int i3 = readGetFbrIndex(nfis,z,a,fbr);
        if(i3 >= 0){
          fbr[i3].potwell[m].height       = atof(readExtractData(3));
          fbr[i3].potwell[m].curvature    = atof(readExtractData(4));
          fbr[i3].potwell[m].absorb       = atof(readExtractData(5));
        }
        else{
          message << "barrier parameter required before potwell";
          cohTerminateCode("readStatModel");
        }
      }

    }else if(head == "levcomp  :" || head == "levcomp+ :"){
      z = readElementToZ(readExtractData(0));
      a = atoi(          readExtractData(1));

      if( nfis < MAX_FISS_CHANCE ){
        int i3 = readGetFbrIndex(nfis,z,a,fbr);
        if(i3 >= 0){
          fbr[i3].za.setZA(z,a);
          fbr[i3].fisenhance.compfact1p = atof(readExtractData(2));
          fbr[i3].fisenhance.compfact2p = atof(readExtractData(3));
        }
      }

    }else if(head == "levcomp- :"){
      z = readElementToZ(readExtractData(0));
      a = atoi(          readExtractData(1));

      if( nfis < MAX_FISS_CHANCE ){
        int i3 = readGetFbrIndex(nfis,z,a,fbr);
        if(i3 >= 0){
          fbr[i3].za.setZA(z,a);
          fbr[i3].fisenhance.compfact1n = atof(readExtractData(2));
          fbr[i3].fisenhance.compfact2n = atof(readExtractData(3));
        }
      }

    }else if(head == "fenhance :"){
      z = readElementToZ(readExtractData(0));
      a = atoi(          readExtractData(1));

      if( nfis < MAX_FISS_CHANCE ){
        int i3 = readGetFbrIndex(nfis,z,a,fbr);
        if(i3 >= 0){
          fbr[i3].za.setZA(z,a);
          fbr[i3].fisenhance.energy = atof(readExtractData(2));
          fbr[i3].fisenhance.width  = atof(readExtractData(3));
          fbr[i3].fisenhance.peak   = atof(readExtractData(4));
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

    }else if(head == "tweakSPIN:"){
      readParSet3(parmSPIN);

    }else if(head == "tweakPDST:"){
      readParSet3(parmPDST);

    }else if(head == "tweakFL1 :"){
      readParSet3(parmFL1);

    }else if(head == "tweakFL2 :"){
      readParSet3(parmFL2);

    }else if(head == "tweakFL3 :"){
      readParSet3(parmFL3);

    }else if(head == "tweakM1  :"){
      readParSet1(parmGDRM1);

    }else if(head == "broaden  :"){
      readParSet1(parmBROD);

    }else{
      message << "unknown key-word [" << head << "] in HFS";
      cohTerminateCode("readStatModel");
    }
  }

  if(ngdr > 0) message << "number of GDRs " << ngdr << " ";
  if(nfis >= 0) message << "number of fissioning nuclei " << nfis + 1;
  if(message.str() != "") cohNotice("readStatModel");

  return(0);
}


/**********************************************************/
/*      Read Exciton Model Parameter                      */
/**********************************************************/
int readExciton(char *s)
{
  while(1){
    if( readFgetOneline(s) < 0 ){
      message << "data structure error in PREEQ section";
      cohTerminateCode("readExcition");
    }

    if(head == "ENDPREEQ :"){
      break;

    }else if(head == "tweakM2  :"){
      readParSet1(parmM2);

    }else if(head == "tweakM2R :"){
      readParSet1(parmM2R);

    }else if(head == "tweakSD  :"){
      readParSet2(parmSD);

    }else if(head == "tweakKO  :"){
      readParSet1(parmKO);

    }else if(head == "tweakCOLL:"){
      readParSet1(parmCOLL);

    }else{
      message << "unknown key-word [" << head << "] in PREEQ";
      cohTerminateCode("readExciton");
    }
  }

  return(0);
}


/**********************************************************/
/*      Read Fission Neutron Spectrum Parameter           */
/**********************************************************/
int readFns(char *s, FNSpec *fns)
{
  int kf = 0, ke = 0, kd = 0, kr = 0, kn = 0, ks = 0, ka = 0;

  for(int i=0 ; i<MAX_FISS_CHANCE ; i++){
    fns->fc[i].tke        = 0.0;
    fns->fc[i].txe        = 0.0;
    fns->fc[i].etotal     = 0.0;
    fns->fc[i].egamma     = 0.0;
    fns->fc[i].meansep    = 0.0;
    fns->fc[i].lf.a       = 0.0;
    fns->fc[i].hf.a       = 0.0;
    fns->fc[i].rt         = 1.0;
    fns->fc[i].nuratio    = 1.0;
    fns->fc[i].tps        = 1.0;
    fns->fc[i].anisotropy = 0.0;
  }

  while(1){
    if( readFgetOneline(s) < 0 ){
      message << "data structure error in FNS section";
      cohTerminateCode("readFns");
    }

    if(head == "ENDFNS   :"){
      break;

    }else if(head == "fragment :"){
      if(kf == MAX_FISS_CHANCE){
        message << "too many fission fragment pairs " << kf + 1;
        cohTerminateCode("readFns");
      }
      int zl = readElementToZ(readExtractData(0));
      int al = atoi(          readExtractData(1));
      int zh = readElementToZ(readExtractData(2));
      int ah = atoi(          readExtractData(3));
      fns->fc[kf].lf.za.setZA(zl,al);
      fns->fc[kf].hf.za.setZA(zh,ah);
      kf ++;

    }else if(head == "energies :"){
      if(ke == MAX_FISS_CHANCE){
        message << "too many fission energies " << ke + 1;
        cohTerminateCode("readFns");
      }
      fns->fc[ke].etotal  =  atof(readExtractData(0));
      fns->fc[ke].tke0    =  atof(readExtractData(1));
      fns->fc[ke].tke1    =  atof(readExtractData(2));
      fns->fc[ke].egamma  =  atof(readExtractData(3));
      fns->fc[ke].meansep =  atof(readExtractData(4));
      ke ++;

    }else if(head == "density  :"){
      if(kd == MAX_FISS_CHANCE){
        message << "too many fission level densities " << kd + 1;
        cohTerminateCode("readFns");
      }
      fns->fc[kd].lf.a    =  atof(readExtractData(0));
      fns->fc[kd].hf.a    =  atof(readExtractData(1));
      fns->fc[kd].ac      =  atof(readExtractData(2));
      kd ++;

    }else if(head == "rt       :"){
      if(kr == MAX_FISS_CHANCE){
        message << "too many Rt parameters " << kr + 1;
        cohTerminateCode("readFns");
      }
      fns->fc[kr].rt = atof(readExtractData(0));
      kr ++;

    }else if(head == "nuratio  :"){
      if(kn == MAX_FISS_CHANCE){
        message << "too many nu-ratio parameters " << kn + 1;
        cohTerminateCode("readFns");
      }
      fns->fc[kn].nuratio = atof(readExtractData(0));
      kn ++;

    }else if(head == "anisotrop:"){
      if(ka == MAX_FISS_CHANCE){
        message << "too many anisotropic parameters " << ka + 1;
        cohTerminateCode("readFns");
      }
      fns->fc[ka].anisotropy =  atof(readExtractData(0));
      ka ++;

    }else if(head == "tpdf_s   :"){
      if(ks == MAX_FISS_CHANCE){
        message << "too many temperature distribution parameters " << ks + 1;
        cohTerminateCode("readFns");
      }
      fns->fc[ks].tps =  atof(readExtractData(0));
      ks ++;

    }else if(head == "maxwell  :"){
      fns->maxwell =  atof(readExtractData(0));

    }else if(head == "omp      :"){
      readExtractData(0);
      fns->omindex = find_omp(inputdata);

    }else{
      message << "unknown key-word [" << head << "] in FNS";
      cohTerminateCode("readFns");
    }
  }

  if(kf > 0) message << "number of multi-chance fission for spectrum " << kf;
  message << " energy " << ke << " density " << kd << " Rt " << kr << " nu-ratio " << kn << " aniso " << ka << " temp dist " << ks;
  if(message.str() != "") cohNotice("readFns");

  return(0);
}


/**********************************************************/
/*      Read Mean Field Theory Parameter                  */
/**********************************************************/
int readMeanfield(char *s, MFTparm *m)
{
  while(1){
    if( readFgetOneline(s) < 0 ){
      message << "data structure error in MFT section";
      cohTerminateCode("readMeanfield");
    }

    if(head == "ENDMFT   :"){
      break;

    }else if(head == "force    :"){
      m->force = readExtractData(0);
      if((m->force != "SIII") && (m->force != "SLy4") && (m->force != "SkM*")){
        message << "unknown Skyrme force [" << m->force << "] in MFT";
        cohTerminateCode("readMeanField");
      }
      ctl.hartreefock = true;

    }else if(head == "pairing  :"){
      m->pairing = readExtractData(0);
      if( (m->pairing != "delta") && (m->pairing != "seniority") ){
        message << "unknown pairing model [" << m->pairing << "] in MFT";
        cohTerminateCode("readMeanField");
      }

    }else if(head == "converge :"){
      m->max_iteration = atoi(readExtractData(0));
      m->eps_converge  = atof(readExtractData(1));

    }else if(head == "strength :"){
      m->strength[0] = atof(readExtractData(0));
      m->strength[1] = atof(readExtractData(1));

    }else if(head == "basis    :"){
      m->n0 = atoi(readExtractData(0));
      m->b  = atof(readExtractData(1));
      m->q  = atof(readExtractData(2));

    }else if(head == "zcm      :"){
      m->value[0]  = atof(readExtractData(0));
      m->weight[0] = atof(readExtractData(1));
      m->nfree[0]  = atoi(readExtractData(2));

    }else if(head == "Q20      :"){
      m->value[1]  = atof(readExtractData(0));
      m->weight[1] = atof(readExtractData(1));
      m->nfree[1]  = atoi(readExtractData(2));

    }else if(head == "Q30      :"){
      m->value[2]  = atof(readExtractData(0));
      m->weight[2] = atof(readExtractData(1));
      m->nfree[2]  = atoi(readExtractData(2));

    }else if(head == "Q40      :"){
      m->value[3]  = atof(readExtractData(0));
      m->weight[3] = atof(readExtractData(1));
      m->nfree[3]  = atoi(readExtractData(2));

    }else if(head == "epsilon  :"){
      m->eps[0]    = atof(readExtractData(0));
      m->eps[1]    = atof(readExtractData(1));
      m->eps[2]    = atof(readExtractData(2));
      m->eps[3]    = atof(readExtractData(3));
      m->eps[4]    = atof(readExtractData(4));
      m->eps[5]    = atof(readExtractData(5));

      ctl.frdm = true;

    }else if(head == "gamma    :"){
      m->gamma     = atof(readExtractData(0));

    }else if(head == "print    :"){
      std::string item = readExtractData(0);
      if(item == "off"){
        pmf.off();
      }
      else if(item == "potential"){
        pmf.potential = true;
      }
      else if(item == "spstate"){
        pmf.spstate = true;
      }
      else if(item == "shiftedstate"){
        pmf.shiftedstate = true;
      }
      else if(item == "bcs"){
        pmf.bcs = true;
      }
      else if(item == "energy"){
        pmf.energy = true;
      }
      else if(item == "expansion"){
        pmf.expansion = true;
      }
      else if(item == "deformation"){
        pmf.deformation = true;
      }
      else if(item == "HFdensity"){
        pmf.HFdensity = true;
      }
      else if(item == "HFmonitor"){
        pmf.HFmonitor = true;
      }
      else if(item == "FRDMzero"){
        pmf.FRDMzero = true;
      }
      else if(item == "FRDMshape"){
        pmf.FRDMshape = true;
      }
      else if(item == "FRDMfileout"){
        pmf.FRDMfileout = true;
      }
      else{
        message << "no such print option [" << item << "] in MFT";
        cohTerminateCode("readMeanField");
      }

    }else{
      message << "unknown key-word [" << head << "] in MFT";
      cohTerminateCode("readMeanField");
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
    cohTerminateCode("readFgetOneLine");
  }

  return(1);
}


/**********************************************************/
/*      Extract Data from Line                            */
/**********************************************************/
char * readExtractData(const int n)
{
  inputdata[0] = '\0';

  if( strlen(line.c_str()) > W_HEAD ){
    std::string work = line;
    char *str = &work[W_HEAD];
    char *tok = nullptr;
    const char *delim = " ";

    tok = strtok(str,delim);

    if(tok != nullptr){
      for(int i=1 ; i<=n ; i++){
        tok = strtok(nullptr,delim);
        if(tok == nullptr) break;
      }
    }
    if(tok != nullptr) strncpy(&inputdata[0],tok,sizeof(inputdata)-1);
  }

  return(&inputdata[0]);
}


/**********************************************************/
/*      Convert Element Name into Z Number                */
/**********************************************************/
int readElementToZ(char *d)
{
  int z = element_getZ(d);

  if(z <= 0){
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
  case 'x' : p = userdef ; break;
  case 'f' : p = fission ; break;
  default  : break;
  }

  if(p == unknown){
    message << "unknown particle " << i;
    cohTerminateCode("readParticleIdentify");
  }

  return(p);
}


/**********************************************************/
/*      Convert Direct Reaction Model into DirType        */
/**********************************************************/
DirType readDirect(char *d)
{
  DirType k = ccrot;
  std::string type = (std::string)d;
  if(     type == "gs")  k = ccvibg;
  if(     type == "1ph") k = ccvib1;
  else if(type == "2ph") k = ccvib2;

  return(k);
}


/**********************************************************/
/*      Get Nucleus Index from Z and A                    */
/**********************************************************/
inline int readGetFbrIndex(const int n, const unsigned int z, const unsigned int a, Fission *fbr)
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
inline void readParSet1(const parameterType pt)
{
  double   f = atof(readExtractData(0));
  double   w = atof(readExtractData(1));
  if(!adj.addparm(pt,f,w)){
    message << "too many adjustable parameters";
    cohTerminateCode("readParSet1");
  }
}

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


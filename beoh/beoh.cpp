/******************************************************************************/
/**                                                                          **/
/**   BeoH   :  Beta-Decay by CoH                                            **/
/**                                            Version 1.2 Naiad    (2017)   **/
/**                                                              T. Kawano   **/
/**   History                                                                **/
/**   1.0  2014 Apr. (Triton)    Imported from CoH3.3.2 Titania              **/
/**             Jul.             Combined with CoH3.4.0 Oberon               **/
/**   1.1  2017 Apr. (Nereid)    Simple Statistical Decay Mode               **/
/**   1.2  2017 Sep. (Naiad)     Prompt Fission Calculation Mode             **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <unistd.h>

#define COH_TOPLEVEL

#include "beoh.h"
#include "beohoutput.h"
#include "global.h"
#include "parameter.h"
#include "coupling.h"
#include "elements.h"
#include "terminate.h"

int main (int, char *[]);
static void cohOutputOptions (const unsigned int);
static void cohHelp (void);
static void cohAllocateMemory (void);
static void cohDeleteAllocated (void);

/**********************************************************/
/*      Global Parameters                                 */
/**********************************************************/

//---------------------------------------
//      Extern variables defined in nucleus.h

CrossSection      crx;               // cross section data
GammaProduction   gml;               // gamma-ray line data 
Nucleus          *ncl;               // nucleus data
NuclearStructure *nst;               // discrete level data
Fission          *fbr;               // fission barrier data

std::string       version = "ver.1.2 (2017 Sep)";


//---------------------------------------
//      Defined in global.h

GlobalControl     ctl;               // booleans to control calculation
GlobalOption      opt;               // booleans to activate optional features
GlobalPrint       prn;               // booleans to control output

//---------------------------------------
//      Extern variables defined in parameter.h

Adjustable        adj;               // adjustable parameters

//---------------------------------------
//      Defined in terminate.h

std::ostringstream message;          // error or warning message string

//---------------------------------------
//      Extern variables defined in coupling.h

double           *fact;              // factorials for coupling constants


/**********************************************************/
/*      Local Scope Variables inside beoh.cpp             */
/**********************************************************/
static int        allocated_ncl = 0; // total number of nucleus
static int        allocated_nst = 0; // total number of nuclear structure data
static bool       verbflag      = false; // verbose output flag


/**********************************************************/
/*      BeoH Starts Here                                  */
/**********************************************************/
int main(int argc, char *argv[])
{

//---------------------------------------
//      Command Line Options

  /*** default output options */
  unsigned int  propt = PRN_SYSTEM | PRN_XSECTION | PRN_SPECTRA;


  /*** excitation energy, Z, and A numbers from command line */
  double        energy   = 0.0; // nuclear excitation energy
  double        exwidth  = 0.0; // distribution of excitation energy
  double        spinfact = 1.0; // multiplication factor for initial spin distribution (sigma)
  double        iniJ     = 0.0; // initial state spin, J
  double        targE    = 0.0; // target kinetic energy, when moving
  double        beta2    = 0.0; // deformation parameter
  int           targZ    = 0;   // target Z number
  int           targA    = 0;   // target A number
  int           iniP     = 0;   // initial state parity, -1 or 1
  unsigned long nsim     = 0;   // number of Monte Carlo sampling
  std::string   elem     = "";  // element name

  int           x;
  while((x = getopt(argc,argv,"p:e:d:a:z:j:J:f:k:b:s:mvh")) != -1){
    switch(x){
    case 'p': propt      = atoi(optarg);                  break;
    case 'e': energy     = atof(optarg);                  break;
    case 'd': exwidth    = atof(optarg);                  break;
    case 'a': targA      = atoi(optarg);                  break;
    case 'z': elem       = optarg;
              targZ = element_getZ(&elem[0]);
              if(targZ == 0) {
                std::cerr << "ERROR     :unknown target name " << elem << std::endl;
                exit(-1);
              }                                           break;
    case 'j': iniJ       = abs(atof(optarg));
              iniP       = (strstr(optarg,"-")) ? -1 : 1; break;
    case 'J': iniJ       = abs(atof(optarg));
              iniP       = -2;                            break;
    case 'f': spinfact   = atof(optarg);                  break;
    case 'k': targE      = atof(optarg);                  break;
    case 'b': beta2      = atof(optarg);                  break;
    case 's': nsim       = atoi(optarg);                  break;
    case 'v': verbflag   = true;                          break;
    case 'h': cohHelp();                                  break;
    case ':': std::cerr << "ERROR     :need a value for option " << x << std::endl;
              cohHelp();                                  break;
    case '?': std::cerr << "ERROR     :unknown option " << argv[argc-1] << std::endl;
              cohHelp();                                  break;
    }
  }

  /*** set print-out global variable */
  cohOutputOptions(propt);

  /*** memory allocation */
  cohAllocateMemory();

  /*** main part */
  if((targZ > 0) && (targA > 0) && (energy > 0.0)){
    beohCGMCompatibleMode(targZ, targA, iniJ, iniP, energy, exwidth, spinfact, targE, beta2, nsim);
  }
  else{
    beohMainLoop(nsim);
  }

  /***  delete allocated memory */
  cohDeleteAllocated();

  return 0;
}


/**********************************************************/
/*      Setting Output Options                            */
/**********************************************************/
void  cohOutputOptions(const unsigned int p)
{
  /*** check and activate each bit */
  prn.system       = p & PRN_SYSTEM;
  prn.xsection     = p & PRN_XSECTION;
  prn.spectra      = p & PRN_SPECTRA;
  prn.gcascade     = p & PRN_GCASCADE;
  prn.mchistory    = p & PRN_MCHISTORY;
}


/**********************************************************/
/*      Help                                              */
/**********************************************************/
void cohHelp()
{
  outBanner();

  std::cerr << std::endl;
  std::cerr << "beoh -e Excitation -a targA -z targZ [-p -j -f -b]" << std::endl;
  std::cerr << "    -p N        :  N is sum of these numbers" << std::endl;
  std::cerr << "                   1 : system parameters" << std::endl;
  std::cerr << "                   2 : ground state production probabilities" << std::endl;
  std::cerr << "                   4 : particle emission energy spectra" << std::endl;
  std::cerr << "    -j Spin_and_Parity" << std::endl;
  std::cerr << "                   spin J=3/2 -> 3.5 as input" << std::endl;
  std::cerr << "                   parity is given by its sign, +3.5, -0, etc" << std::endl;
  std::cerr << "    -f SpinFactor" << std::endl;
  std::cerr << "                   scaling factor to determine the initial spin distribution" << std::endl;
  std::cerr << "    -b Deformation_Parameter" << std::endl;
  std::cerr << "                   deformation parameter, beta2, if deformed" << std::endl;
  exit(0);
}


/**********************************************************/
/*      Allocate Memory                                   */
/**********************************************************/
void cohAllocateMemory()
{
  try{
    /*** compound nucleus, all other nuclei are allocated dynamically */
    ncl = new Nucleus [MAX_COMPOUND];
    ncl[allocated_ncl++].memalloc(MAX_LEVELS,MAX_ENERGY_BIN);

    /*** discrete level data */
    nst = new NuclearStructure [MAX_NUCLIDE];
    nst[allocated_nst++].memalloc(MAX_LEVELS);

    /*** fission barrier data */
    fbr = new Fission [MAX_FISS_CHANCE];
    for(int i=0 ; i<MAX_FISS_CHANCE ; i++){
      fbr[i].memalloc();
    }

    /*** residual nucleus production cross section */
    crx.prod = new ParticleCount [MAX_COMPOUND];

    /*** particle emission spectra */
    crx.spectalloc(MAX_CHANNEL+2,MAX_ENERGY_BIN);

    /*** gamma-ray line array, need larger number for HF3D */
    gml.memalloc(MAX_GAMMALINE*2,CUT_GAMMALINE);

    /*** factorial */
    factorial_allocate();

    /*** parameter adjustment */
    adj.memalloc(MAX_ADJUSTABLES);
  }
  catch(std::bad_alloc &e){
    std::cerr << "ERROR     :memory allocation error " << e.what() << std::endl;
    exit(-1);
  }
}


/**********************************************************/
/*      Instantiation Nucleus Object                      */
/**********************************************************/
void cohAllocateNucleus()
{
  try{
    ncl[allocated_ncl++].memalloc(MAX_LEVELS,MAX_ENERGY_BIN);
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what() << " for nucleus " << allocated_ncl;
    cohTerminateCode("cohAllocateNucleus");
  }
}


/**********************************************************/
/*      Instantiation DiscreteLevel Object                */
/**********************************************************/
void cohAllocateNuclearStructure()
{
  try{
    nst[allocated_nst++].memalloc(MAX_LEVELS);
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what() << " for nuclear structure data " << allocated_nst;
    cohTerminateCode("cohAllocateNuclearStructure");
  }

  if(allocated_nst == MAX_NUCLIDE){
    message << "number of unique nuclides exceed the maximum " << allocated_nst;
    cohTerminateCode("cohAllocateNuclearStructure");
  }
}


/**********************************************************/
/*      Reset Object Pointers                             */
/**********************************************************/
void cohResetAllocationPointer()
{
  /*** Before starting a new calculation,
       move the pointer to the first object. */
  allocated_nst = 1;
  allocated_ncl = 1;
}


/**********************************************************/
/*      Reset Nucleus Object                              */
/**********************************************************/
void cohResetNucleus()
{
  delete [] ncl;
  allocated_ncl = 0;

  /*** re-allocate nucleus up to MAX_COMPOUND */
  ncl = new Nucleus [MAX_COMPOUND];
  ncl[allocated_ncl++].memalloc(MAX_LEVELS,MAX_ENERGY_BIN);
}


/**********************************************************/
/*      Reset DiscreteLevel Object                        */
/**********************************************************/
void cohResetNuclearStructure()
{
  delete [] nst;
  allocated_nst = 0;

  /*** re-allocate discrete level data up to MAX_NUCLIDE */
  nst = new NuclearStructure [MAX_NUCLIDE];
  nst[allocated_nst++].memalloc(MAX_LEVELS);
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void cohDeleteAllocated()
{
  delete [] ncl;
  delete [] nst;
  delete [] fbr;
  delete [] crx.prod;
  crx.spectfree();
  factorial_delete();
}


/**********************************************************/
/*     Warning Message                                    */
/**********************************************************/
void cohWarningMessage(std::string module)
{
  std::cerr << "WARNING   : ";
  if(module != "") std::cerr << "[" << module <<"] ";
  std::cerr << message.str() << std::endl;
  message.str("");
}

void cohNotice(std::string module)
{
  if(module == "NOTE"){
    std::cerr << " (._.) " << message.str() << std::endl;
  }
  else{
    if(verbflag) std::cerr  << " (@_@) [" << module << "] " << message.str() << std::endl;
  }
  message.str("");
}


/**********************************************************/
/*     Emergency Stop                                     */
/**********************************************************/
int cohTerminateCode(std::string module)
{
  /*** Release global storage */
  cohDeleteAllocated();

  /*** Exit code */
  std::cerr << "ERROR     : ";
  if(module != "") std::cerr << "[" << module << "] ";
  std::cerr << message.str() << std::endl;
  exit(-1);
}


/******************************************************************************/
/**                                                                          **/
/**   MFLD   :  Mean-Field Level Density by CoH                              **/
/**                                                            Version 0.1   **/
/**                                                              T. Kawano   **/
/**   History                                                                **/
/**          0.1  2024 Nov. Initial Verision                                 **/
/**                                                                          **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <random>

#include "global.h"
#include "mfld.h"
#include "FRDM.h"
#include "elements.h"
#include "terminate.h"
#include "mt19937.h"

#define COH_TOPLEVEL

int main (int, char **);
static void mfldSystemDefault (SystemData *);
static void mfldSetupParameters (SystemData *, PrintData *, MFTSystem *, ZAnumber *, const bool);
static void mfldHelp ();
static void cohCheckRange (ZAnumber *);
static void cohAllocateMemory (SystemData *);
static void cohDeleteAllocated (void);
static void setupInitMT (void);


//---------------------------------------
//      Global Parameters

GlobalPrintMFT    pmf;                  // print option for FRDM
double          *fact;                  // factorial

//---------------------------------------
//      Defined in terminate.h

std::ostringstream   message;           // error or warning message string

//---------------------------------------
//      Local Scope Variables

static bool          verbflag = false;  // verbose output flag
static Density       sd;                // combinatorial level densities
static HFInterface   hfsp[2];           // extracted data from FRDM for other calculations


/**********************************************************/
/*      Combinatorial Level Density Calculation           */
/**********************************************************/
int main(int argc, char *argv[])
{

//----------------------------------------------------------
//      Command Line Options

  int           targZ    = 0;     // target Z number
  int           targA    = 0;     // target A number
  bool          skipflag = false; // skip data read, and use all default
  std::string   elem     = "";    // element name

  int x;
  while((x = getopt(argc,argv,"a:z:dvh")) != -1){
    switch(x){
    case 'a': targA      = atoi(optarg);                  break;
    case 'z': elem       = optarg;
              targZ = element_getZ(&elem[0]);
              if(targZ == 0) {
                std::cerr << "ERROR     :unknown target name " << elem << std::endl;
                exit(-1);
              }                                           break;
    case 'd': skipflag   = true;                          break;
    case 'v': verbflag   = true;                          break;
    case 'h': mfldHelp();                                 break;
    case ':': std::cerr << "ERROR     :need a value for option " << x << std::endl;
              mfldHelp();                                 break;
    case '?': std::cerr << "ERROR     :unknown option " << argv[argc-1] << std::endl;
              mfldHelp();                                 break;
    }
  }

  /*** model parameters and size of data arrays */
  SystemData sys;
  PrintData  prn;
  MFTSystem  mft;


//----------------------------------------------------------
//  Parameter Initialization and Data Read

  /*** read input and initialize system parameters */
  ZAnumber za(targZ,targA);
  mfldSystemDefault(&sys);
  mfldSetupParameters(&sys,&prn,&mft,&za,skipflag);

  /*** allocate arrays */
  cohAllocateMemory(&sys);


//----------------------------------------------------------
//  Single-Particle Energy by Independent Particle Model


  /*** include unbound states if spmax > 0 */
  if(sys.spmax > 0.0) {
    for(int p=0 ; p<=1 ; p++) hfsp[p].emax = sys.spmax;
  }

  /*** calculate single particle energy spectra */
  FRDMCalc(&mft,hfsp);


//----------------------------------------------------------
//  Combinatorial Level Density from Single-Particle States

  SDUnperturbed(&sys,hfsp,&sd);


//----------------------------------------------------------
//  Print Calculatedl Level Densities

  bool halfint = ((sys.target.getA() % 2) != 0) ? true : false;

  /*** total observed level density */
  if(prn.TotalDensity)      printLevelDensity(prn.BothParity,&sd);
  if(prn.PartialDensity)    printPartialLevelDensity(halfint,prn.BothParity,&sd);
  if(prn.JDependentDensity) printJDependentLevelDensity(halfint,prn.BothParity,&sd);
  if(prn.Plotting)          printPlottingLevelDensity(prn.BothParity,&sd);
  if(prn.SpinDistribution)  printSpinDistribution(false,sys.eave,&sd);
  if(prn.SpinCutOff)        printSpinDistribution(true,sys.eave,&sd);
  if(prn.AverageConfig)     printAverageConfiguration(&sd);


//----------------------------------------------------------
//  Clean Up

  /*** free allocated memory */
  cohDeleteAllocated();

  return 0;
}


/**********************************************************/
/*      Set Default Sysmtem Parameters                    */
/**********************************************************/
void mfldSystemDefault(SystemData *dat)
{
  /*** set default calculation space size */
  dat->lsize = MAX_SPSTATE;             // max number single-particle levels in N or Z shell
  dat->nsize = MAX_ENERGY_BIN;          // max number of excitation energy bins
  dat->psize = MAX_PARTICLE;            // max N for N-particle N-hole configuration
  dat->msize = MAX_M;                   // max M-value in M-scheme level density calculation
  dat->jsize = MAX_J;                   // max J-value for level density

  /*** other parameters, to be overwritten by inputs */
  dat->phmax = dat->psize;
  dat->de    = ENERGY_BIN;
  dat->emax  = dat->de * dat->nsize;

  /*** initialize random number seed */
  setupInitMT();
}


/**********************************************************/
/*      Set Up Sysmtem Parameters                         */
/**********************************************************/
void mfldSetupParameters(SystemData *dat, PrintData *prn, MFTSystem *mft, ZAnumber *t, const bool skipflag)
{
  /*** default FRDM shape, read database by setting -1 */
  double eps[6];
  eps[0] = -1;
  for(int i=1 ; i<6 ; i++) eps[i] = 0.0;

  pmf.off();
  /*** command-line only mode */
  if(skipflag){
    prn->TotalDensity = true;
  }
  /*** read input data */
  else{
    char buf[256];
    readSystem(buf,dat,prn,eps);
  }

  /*** overwrite if Z and A are given */
  if(t->getZ() > 0) dat->target.setZA(t->getZ(),t->getA());
  cohCheckRange(&dat->target);

  /*** check Emax, if less than dE x Nsize */
  double ex = dat->de * dat->nsize; // max energy before data read
  if(dat->emax != ex){
    if(dat->emax > ex) dat->emax = dat->de * dat->nsize;
    else               dat->nsize = (int)(dat->emax / dat->de) + 1;
  }

  sd.de = dat->de;

  /*** if ph-max is given by input */
  if(dat->phmax != dat->psize) dat->psize = dat->phmax;

  /*** nuclear shape parameters */
  mft->init(dat->target.getZ(), dat->target.getA());
  mft->setNuclearSize(gRadius);

  for(int i=0 ; i<6 ; i++) mft->setEps(i+1,eps[i]);
  if(mft->getEps(1) < 0.0) FRDMReadShape(mft);

  message << "highest energy:" << dat->emax << "  bin-width:" << dat->de << "  number of energy bin:" << dat->nsize;
  cohNotice("mfldSetupParameters");

  message << "max number of ph configuration:" << dat->psize << "  max J:" << dat->jsize << "  max single-particle states:" << dat->lsize;
  cohNotice("mfldSetupParameters");
}


/**********************************************************/
/*      Check Z/A Range                                   */
/**********************************************************/
void cohCheckRange(ZAnumber *t)
{
  const unsigned int zmin =    10;
  const unsigned int zmax =   118;
  const unsigned int amin =    20;
  const unsigned int amax =   300;

  bool rangechk = false;
  if((t->getZ() < zmin) || (t->getZ() > zmax)){
    message << "target Z-number " << t->getZ() << " out of range";
    rangechk = true;
  }
  else if((t->getA() < amin) || (t->getA() > amax)){
    message << "target A-number " << t->getA() << " out of range";
    rangechk = true;
  }

  if(rangechk) cohTerminateCode("cohCheckRagne");
}


/**********************************************************/
/*      Help                                              */
/**********************************************************/
void mfldHelp()
{
  std::cerr << "mfld < input_file" << std::endl;
  std::cerr << "mfld -a targA -z targZ < input_file" << std::endl;
  std::cerr << "       Z and A in input_file will be replaced by command-line" << std::endl;
  std::cerr << "mfld -a targA -z targZ -d" << std::endl;
  std::cerr << "       no input_file mode" << std::endl;
  exit(0);
}


/**********************************************************/
/*      Allocate Memory                                   */
/**********************************************************/
void cohAllocateMemory(SystemData *dat)
{
  try{
    /*** level densities */
    sd.memalloc(dat->psize, dat->nsize, dat->jsize);

    hfsp[neutron].memalloc();
    hfsp[proton ].memalloc();
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error" << e.what();
    cohTerminateCode("cohAllocateMemory");
  }

  sd.memclear();
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void cohDeleteAllocated()
{
  sd.memfree();

  hfsp[neutron].memfree();
  hfsp[proton ].memfree();
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


/*************************************************/
/*  Initialize MT Random Number Generator        */
/*************************************************/
#define RANDOM_SEED_BY_SYSTEM

void setupInitMT(void)
{
  unsigned long seed = 87545L;

#ifdef RANDOM_SEED_BY_SYSTEM
  std::random_device seedgen;
  seed = (unsigned long)seedgen();
#endif

  genrand_init(seed);
}

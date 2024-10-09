// global variables, data output and calculation flow control

/**************************************/
/*      Output Control                */
/**************************************/

/*** standard output */
#define PRN_SYSTEM       0x0001    //    1
#define PRN_XSECTION     0x0002    //    2
#define PRN_SPECTRA      0x0004    //    4
#define PRN_LEVELEXCITE  0x0008    //    8
#define PRN_GCASCADE     0x0010    //   16
#define PRN_POPULATION   0x0020    //   32
#define PRN_TRANSMISSION 0x0040    //   64
#define PRN_ANGDIST      0x0080    //  128
#define PRN_PREEQPARM    0x0100    //  256
#define PRN_SMATRIX      0x0200    //  512
#define PRN_MCHISTORY    0x0400    // 1024
#define PRN_FISENERGY    0x0800    // 2048

/*** non-standard output */
#define PEX_GBRANCH      0x0001    //    1
#define PEX_CUMULEVEL    0x0002    //    2
#define PEX_DENSITY      0x0004    //    4
#define PEX_FISDENSITY   0x0008    //    8
#define PEX_DDX          0x0010    //   16
#define PEX_PTABLE       0x0020    //   32
#define PEX_GAMMALINE    0x0040    //   64
#define PEX_DECAYWIDTH   0x0080    //  128
#define PEX_GAMMASTF     0x0100    //  256


class GlobalPrint{
 public:
  bool system;                     // print out system parameters
  bool xsection;                   // print cross sections
  bool spectra;                    // print particle energy spectra
  bool levelexcite;                // print discrete level population
  bool gcascade;                   // print gamma-ray cascading
  bool population;                 // print internal populations
  bool transmission;               // print transmission coefficients
  bool angdist;                    // print angular distributions
  bool preeqparm;                  // print pre-compound parameters
  bool smatrix;                    // print scattering matrix elements
  bool mchistory;                  // print Monte Carlo history
  bool fisenergy;                  // print fission energy table

  GlobalPrint(){
    system       = false;
    xsection     = false;
    spectra      = false;
    levelexcite  = false;
    gcascade     = false;
    population   = false;
    transmission = false;
    angdist      = false;
    preeqparm    = false;
    smatrix      = false;
    mchistory    = false;
    fisenergy    = false;
  }
};


class GlobalPrintExt{
 public:
  bool gbranch;                    // print gamma-ray branching ratios
  bool cumulevel;                  // print cumulative levels
  bool density;                    // print level density
  bool fisdensity;                 // print fission level density
  bool ddx;                        // print DDX data for ENDF6 formatting
  bool ptable;                     // print decay probability matrix
  bool gammaline;                  // print all line gammas
  bool decaywidth;                 // print compound nucleus channel decay width
  bool gammastf;                   // print gamma-ray strength function

  GlobalPrintExt(){
    gbranch      = false;
    cumulevel    = false;
    density      = false;
    fisdensity   = false;
    ddx          = false;
    ptable       = false;
    gammaline    = false;
    decaywidth   = false;
    gammastf     = false;
  }
};


class GlobalPrintMFT{
 public:
  bool system;                     // print system parameters
  bool potential;                  // print potential on 3D mesh
  bool spstate;                    // print single-particle spectrum
  bool shiftedstate;               // print BCS shifted single-particle states
  bool bcs;                        // print BCS parameters
  bool energy;                     // print calculated energies
  bool expansion;                  // print spherical expansion coefficients
  bool deformation;                // print deformation parameters
  bool HFmonitor;                  // print HF iteration monitor
  bool HFdensity;                  // print nuclear density on 3D mesh
  bool FRDMbeta;                   // print spherical harmonics expansion
  bool FRDMzero;                   // calculate zero-point energy
  bool FRDMshape;                  // output nuclear shape on file
  bool FRDMfileout;                // output internal parameters on file


  GlobalPrintMFT(){
    system       = true;
    potential    = false;
    spstate      = false;
    shiftedstate = false;
    bcs          = true;
    energy       = true;
    expansion    = false;
    deformation  = true;
    HFmonitor    = true;
    HFdensity    = false;
    FRDMbeta     = true;
    FRDMzero     = false;
    FRDMshape    = false;
    FRDMfileout  = false;
  }

  void off(){
    system       = false;
    potential    = false;
    spstate      = false;
    shiftedstate = false;
    bcs          = false;
    energy       = false;
    expansion    = false;
    deformation  = false;
    HFmonitor    = false;
    HFdensity    = false;
    FRDMbeta     = false;
    FRDMzero     = false;
    FRDMshape    = false;
    FRDMfileout  = false;
  }
};


/**************************************/
/*      Calculation Flow Control      */
/**************************************/

class GlobalControl{
 public:
  bool entrance;                   // this is entrance channel calculation
  bool deformed;                   // target is deformed
  bool exciton;                    // include exciton model
  bool dwba;                       // include DWBA calculation
  bool dsd;                        // include DSD model
  bool fluctuation;                // include width fluctuation
  bool ewtransform;                // Engelbrecht-Weidenmueller transformation
  bool statmodel;                  // include statistical model
  bool fission;                    // perform fission calculation
  bool montecarlo;                 // do Monte Carlo simulation 
  bool fileout;                    // write cross sections on a file
  bool superlazy;                  // super-lazy mode for quick calculation
  bool exclusive;                  // exclusive particle spectra calculation
  bool ddx;                        // double differential spectra
  bool fns;                        // fission neutron spectra
  bool macs;                       // maxwellian average cross section
  bool meanfield;                  // mean-field theory
  bool hartreefock;                // Hartree-Fock calculation
  bool frdm;                       // finite-range droplet model
  bool expansion;                  // HF/FRDM s.p. state expansion

  GlobalControl(){
    entrance       = false;
    deformed       = false;
    exciton        = false;
    dwba           = false;
    dsd            = false;
    fluctuation    = false;
    ewtransform    = false;
    statmodel      = false;
    fission        = false;
    montecarlo     = false;
    fileout        = false;
    superlazy      = false;
    exclusive      = false;
    ddx            = false;
    fns            = false;
    macs           = false;
    meanfield      = false;
    hartreefock    = false;
    frdm           = false;
    expansion      = false;
  }
};


/**************************************/
/*      Calculation Options           */
/**************************************/
class GlobalOption{
 public:
  bool ewtransformation;       // Engelbrecht-Weidenmueller transformation
  bool internalconversion;     // discrete gamma-ray internal conversion
  bool readdensity;            // read level density table in the current directory
  bool readphotoabsorption;    // read tabulated photo-absorption in the current directory
  bool groundstateomp;         // use G.S. omp for excitated states
  bool chargediscrete;         // explicit charged particle discrete transitions
  bool finegammaspectrum;      // gamma line mapped onto a finer energy grid
  bool continuumangdist;       // angular distributions in continuum

  GlobalOption(){
    ewtransformation    = false;
    internalconversion  = false;
    readdensity         = false;
    readphotoabsorption = false;
    groundstateomp      = false;
    chargediscrete      = false;
    finegammaspectrum   = false;
    continuumangdist    = false;
  }
};


#ifndef COH_TOPLEVEL
extern GlobalPrint    prn;
extern GlobalPrintExt pex;
extern GlobalPrintMFT pmf;
extern GlobalControl  ctl;
extern GlobalOption   opt;
#endif

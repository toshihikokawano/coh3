// array size and some default constants

/****************************/
/*      ARRAY SIZE          */
/****************************/

const int MAX_ENERGY_BIN    =  3000;  // maximum energy bins
const int MAX_CHANNEL       =     7;  // maximum decay channel

#ifdef BeoH
const int MAX_COMPOUND      =   200;  // maximum number of compound nucleus
const int MAX_NUCLIDE       =    50;  // number of unique nuclides in chain
const int MAX_DIRECT        =     1;  // number of direct reaction
const int MAX_FISS_CHANCE   =    20;  // total number of multi-chance fission
const int MAX_GAMMALINE     = 30000;  // maximum number of gamma lines printed
#else
const int MAX_COMPOUND      =   500;  // maximum number of compound nucleus
const int MAX_NUCLIDE       =   150;  // number of unique nuclides in chain
const int MAX_DIRECT        =    20;  // number of direct reaction
const int MAX_FISS_CHANCE   =    20;  // total number of multi-chance fission
const int MAX_GAMMALINE     =  2000;  // maximum number of gamma lines printed
#endif

const int MAX_J             =    60;  // maximum J-value  = MAX_L
const int MAX_LEVELS        =   300;  // maximum discrete levels
const int MAX_GAMMA_BRANCH  =    50;  // maximum gamma-ray branches
const int MAX_MULTIPOL      =     5;  // multipolarity, E1, M1, E2, M2, E3
const int MAX_GDR           =     8;  // maximum GDRs and pygmy resonance
const int MAX_LAMBDA        =     5;  // lambda_max = (MAX_LAMBDA-1)*2
const int MAX_ANGDIST       =   180;  // maximum angular distribution points
const int MAX_ANGDISTLEVELS =    50;  // maximum levels angular dist. given
const int MAX_HUMP          =     3;  // maximum number of fission humps
const int MAX_KBAND         =     5;  // maximum number of transition states
const int MAX_EINCIDENT     =   100;  // maximum number of incident energies
const int MAX_MACSTEMP      =    50;  // maximum number of MACS temperature


/****************************/
/*      DEFAULT PARAMETER   */
/****************************/

const double ENERGY_BIN     =   0.1;  // default energy bin in the continuum
const double NORM_FACT      =  10.0;  // conversion factor to [mb]
const int    MAX_PARTICLE   =    10;  // number of emitted particles in chain
const int    ANGLE_STEP     =     5;  // angle step for angular dist. printing
const double FINE_GRID      =  0.01;  // finer grid for gamma-ray spec, 10 keV
const double CUT_GAMMALINE  =  1e-6;  // cut-off for discrete gamma-lines for printing

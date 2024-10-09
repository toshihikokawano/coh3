// functions for exclusive spectrum calculations

static const int NLEG =    5; // max order of Legendre expansion
static const int NDDX =   61; // calculate angular points
static const int NGLN = 1001; // maximum number of discrete gamma lines
static const int NCNT = 1000; // maximum number of continuum energy grids
static const int NPRG =  100; // maximum number of primary gamma lines
static const int NLEV =  100; // maximum number of discrete levels


/****************************/
/*    Discrete Levels Data  */
/****************************/
class LevPop{
 public:
  double     gamma;    // level population by gamma cascade
  double     particle; // population by binary reaction
  LevPop(){
    gamma    = 0.0;
    particle = 0.0;
  }
};

/****************************/
/*   Spectrum Data          */
/****************************/
class EXSpectra{
 private:
  double   emax; // maximum excitation energy
  double   elev; // higheset discrete level energy
  int      ncon; // number of continuum bins
  int      nlev; // number of discrete levels
  int      cmax; // maximum decay channels
  int      kmax; // maximum energy spectrum bins
  bool     gcon; // discrete population includes both gamma and particle transisions
 public:
  double **spec; // particles and gamma spectra
  double  *glin; // gamma line spectra given in continuum
  double  *sbin; // particle binary reaction spectra in continuum
  double  *pfis; // fission probability
  LevPop  *lpop; // level population by particle or gamma

  EXSpectra(){
    emax = 0.0;
    elev = 0.0;
    ncon = 0;
    cmax = 0;
    kmax = 0;
    nlev = 0;
    gcon = false;
  }

  ~EXSpectra(){
    for(int c=0 ; c<cmax ; c++) delete [] spec[c];
    delete [] spec;
    delete [] glin;
    delete [] sbin;
    delete [] pfis;
    delete [] lpop;
  }

  void memalloc(int c, int k){
    cmax = c;
    kmax = k;

    /*** energy spectra */
    spec = new double * [cmax];
    for(int c=0 ; c<cmax ; c++){
      spec[c] = new double [kmax];
    }
    glin = new double [kmax];
    sbin = new double [kmax];
    pfis = new double [kmax];

    /*** discrete levels */
    lpop = new LevPop [NLEV];
 }

  void set(int nc, int nl, double em, double el){
    ncon = nc;
    nlev = nl;
    emax = em;
    elev = el;
  }
  void setGcon(bool p){
    gcon = p;
  }
  double getEm(){
    return emax;
  }
  double getEl(){
    return elev;
  }
  int    getNc(){
    return ncon;
  }
  int    getNl(){
    return nlev;
  }
  bool   getGcon(){
    return gcon;
  }
  int    getCmax(){
    return cmax;
  }
  int    getKmax(){
    return kmax;
  }
};


/**************************************/
/*      eclipse.cpp                   */
/**************************************/
void    eclCalc (System *, double **, const unsigned long);
void    eclAllocateMemory (const int, const int);
void    eclDeleteAllocated (const int, const int);
void    eclGenerateDecayTable (const int, const int, const int, double **);
void    eclLevelPopulation (const int, double **, double **);
void    eclRetrievePrefissionSpectrum (const int, double *);
void    eclRetrieveFissionProbability (const int, double *);


/**************************************/
/*      eclspectra.cpp                */
/**************************************/
void    eclSpectra (const bool, const bool, const int, int *, int **, double ****, EXSpectra *);


/**************************************/
/*      ecloutput.cpp                 */
/**************************************/
void    eclOutHead (const int, const int, const double);
void    eclOutSpectra (const int, const int, const int, const double, double **, EXSpectra *);
void    eclTotalSpectra (const int, const int, const int, const double, EXSpectra *);
void    eclOutNucleusHead (const int, const int, double **);
void    eclOutChannelHead (const int, const double);
void    eclOutLegCoeff (const int, const int, double *, double **);
void    eclOutGammaLine (const int, const int, double *, double *);
void    eclOutPtable (const int, const int, int *, int **, double ****);
int     eclFindKmax (const int, const int, double **);


/**************************************/
/*      eclddx.cpp                    */
/**************************************/
void    eclDDX (const int, const int, const int, System *, double **, double **, EXSpectra *);


/**************************************/
/*      eclmc.cpp                     */
/**************************************/
void    eclMC (const int, const int, const int, int *, int **, double ****, double **, const unsigned long, double, EXSpectra *);


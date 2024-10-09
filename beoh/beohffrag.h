static const int MAX_FINE_TUNE = 20;

/****************************/
/*   Multi-Chance Fission   */
/****************************/
class MultiChanceFissionData{
 private:
  int    mftune;              // max number of fine tune nuclides
  int    nftune;              // number of fine tune nuclides
  bool   arrayalloc;          // memory allocation flag
 public:
  double fraction;            // fraction for multi-chance fission
  double spinfactor;          // spin distribution scaling for individual chance-fission
  double rta[80];             // RT(A) fitted at thermal
  double spinfactor1;         // energy-dependent term of spin scaling factor
  double rt;                  // RT for energy sorting mechanism
  double rt1;                 // energy-dependent term of RT
  double tke;                 // total kinetic energy, TKE
  double tke1;                // energy-dependent term of TKE
  double exfis;               // average excitaion energy causing fission for chance-fission
  double eprefis;             // energy of prefission neutron
  double GaussSigma[8];       // Gaussian widths for fragment yield
  double GaussDelta[8];       // A-shift of Gaussian distribution
  double GaussFract[8];       // fraction for each Gaussians
  double yap[18];             // parameters for the energy-dependent Y(A) model
  double tkep[4];             // parameters for TKE(Einc)
  double tkea[4]; 	      // parameters for TKE(A)
  double stkea[3];            // parameters for sigma_TKE(A)
  double ZpFactor[4];         // even-odd factor Fz and Fn multiplicatin factor
  std::string ffydata;        // model to be used for Y(Z,A,TKE)
  int    *ftuneMass;          // mass numbers to be modified
  double *ftuneFactor;        // tuning factor

  MultiChanceFissionData(){
    mftune      = MAX_FINE_TUNE;
    ftuneMass   = new int [mftune];
    ftuneFactor = new double [mftune];
    arrayalloc  = true;

    init();
  }

  ~MultiChanceFissionData(){
    if(arrayalloc){
      delete [] ftuneMass;
      delete [] ftuneFactor;
      arrayalloc = false;
    }
  }

  void init(){
    nftune     = 0;
    fraction   = 1.0;
    spinfactor = 1.0;
    spinfactor1= 0.0;
    rt         = 1.0;
    rt1        = 0.0;          
    tke        = 0.0;
    tke1       = 0.0;
    exfis      = 0.0;
    eprefis    = 0.0;
    ffydata    = "internal";
    for(int i=0 ; i<8 ; i++) GaussSigma[i] = GaussDelta[i] = GaussFract[i] = 0.0;
    for(int i=0 ; i<80 ; i++) rta[i] = 0.0;
    for(int i=0 ; i<18 ; i++) yap[i] = 0.0;
    for(int i=0 ; i<4 ; i++) tkep[i] = 0.0;
    for(int i=0 ; i<4 ; i++) tkea[i] = 0.0;
    for(int i=0 ; i<3 ; i++) stkea[i] = 0.0;
    for(int i=0 ; i<4 ; i++) ZpFactor[i] = 0.0;
  }

  bool setFTune(int a, double f){
    bool flag;
    if(nftune >= mftune) flag = false;
    else{
      ftuneMass[nftune] = a;
      ftuneFactor[nftune] = f;
      nftune ++;
      flag = true;
    }
    return flag;
  }

  int getNFTune(){ return nftune; }

};


/****************************/
/*   Fission Decay Data     */
/****************************/
class FFragData{
 private:
  int maxcf;                  // maximum number of chance-fission
  int ncf;                    // number of chance-fission data currently assigned
  bool   arrayalloc;          // memory allocation flag
 public:
  bool   spontaneous;         // flag for spontaneous fission
  double massresolution;      // mass resolution for output
  double maxtemperature;      // Maxwellian temperature for output
  double ycutoff;             // minimun Y(Z,A,TKE) to be included
  double maxhalflife;         // cut-off half-life for long-lived isotopes (year)
  std::string mcffile;        // multi-chance fission data file
  std::string branchdatafile; // beta-decay branching ratio data file
  std::string yafile;	      // user defined Y(A) model parameters
  std::string ytkefile;       // user defined Y(TKE|A) model parameter
  std::string rtafile; 	      // user defined R_T(A) distribuiton
  int    selectZ;             // select Z number of FF to be calculated
  int    selectA;             // select A number

  MultiChanceFissionData *mc; // fission-chance dependent data

  FFragData(){
    maxcf = MAX_FISS_CHANCE;
    mc = new MultiChanceFissionData [maxcf];
    arrayalloc = true;
    init();
  }

  ~FFragData(){
    if(arrayalloc){
      delete [] mc;
      arrayalloc = false;
    }
  }

  void init(){
    ncf            = 0;
    massresolution = 0.0;
    maxtemperature = 1.42;
    ycutoff        = 0.0;
    maxhalflife    = 1000.0;
    mcffile        = "";
    branchdatafile = "";
    yafile         = "";
    ytkefile	   = "";
    rtafile 	   = "";
    selectZ        = 0;
    selectA        = 0;

    if(arrayalloc){
      for(int i=0 ; i<maxcf ; i++) mc[i].init();
    }
  }

  void inclFissionChance(){
    if(ncf <= maxcf) ncf ++;
  }

  int getFissionChance(){ return ncf; }
  int getMaxFissionChance(){ return maxcf; }

  void resetFissionChance(int n){ if(n <= maxcf) ncf = n; }
};


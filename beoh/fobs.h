/****************************/
/*   Output Quantities      */
/****************************/
class FragmentObservable{
 private:
  int       lsize;            // number of total particle emission
  int       nsize;            // number of energy bins
  int       ksize;            // TKE bins
  bool      allocated;
 public:
  double    multiplicity[2];       // neutron and photon multiplicities
  double    etotal[2];             // total energies
  double    eaverage[2];           // average energies
  double    *spectrum[2];          // CMS energy spectra
  double    *speclab;              // LAB energy spectra
  double    *Pn;                   // P(nu)
  double    *yTKE;                 // sum of yield at each TKE bin
  double    *nTKE;                 // nu(TKE)

  FragmentObservable(){
    allocated = false;
    lsize     = 0;
    nsize     = 0;
    ksize     = 0;
    clear();
  }

  ~FragmentObservable(){
    memfree();
  }

  void clear(){
    for(int p=0 ; p<2 ; p++) multiplicity[p] =  eaverage[p] = etotal[p] = 0.0;
    if(allocated){
      for(int i=0 ; i<nsize ; i++) speclab[i] = 0.0;
      for(int p=0 ; p<2 ; p++){
        for(int i=0 ; i<nsize ; i++) spectrum[p][i] = 0.0;
      }
      for(int i=0 ; i<lsize ; i++) Pn[i] = 0.0;
      for(int i=0 ; i<ksize ; i++) nTKE[i] = yTKE[i] = 0.0;
    }
  }

  void memalloc(int l, int n, int k){
    if(!allocated){
      lsize = l;
      nsize = n;
      ksize = k;

      speclab = new double[nsize];
      for(int p=0 ; p<2 ; p++){
        spectrum[p] = new double[nsize];
        for(int i=0 ; i<nsize ; i++) spectrum[p][i] = 0.0;
      }
      Pn  = new double [lsize];
      yTKE = new double [ksize];
      nTKE = new double [ksize];

      clear();

      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      for(int p=0 ; p<2 ; p++){
        delete [] spectrum[p];
      }
      delete [] speclab;
      delete [] Pn;
      delete [] yTKE;
      delete [] nTKE;

      allocated = false;
    }
  }
};


class FissionObservable{
 private:
  int       lsize;            // number of total particle emission
  int       msize;            // number of mass
  int       nsize;            // number of energy bin
  int       ksize;            // TKE bins
  bool      allocated;
 public:
  FragmentObservable L;
  FragmentObservable H;
  double    multiplicity[2];       // neutron and photon multiplicities
  double    etotal[2];             // total energies
  double    eaverage[2];           // average energies
  double    *spectrum[2];          // average CMS energy spectra
  double    *chi;                  // prompt fission neutron spectrum
  double    *multiplicitydist[2];  // multiplicities as function of A
  double    *eaveragedist[2];      // average energies as function of A
  double    *preMassYield;         // Y(A) before neutron emission
  double    *postMassYield;        // Y(A) after neutron emission
  double    *javeragedist;         // distribution of average spin <J>
  double    *Pn;                   // neutron multiplicity distribution
  double    *yTKE;                 // sum of yield at each TKE bin
  double    *nTKE;                 // nu(TKE)
  double    eprefis;               // prefission neutron total energy
  double    nprefis;               // prefission neutron multiplicity

  FissionObservable(){
    allocated = false;
    lsize     = 0;
    msize     = 0;
    ksize     = 0;
  }

  ~FissionObservable(){
    memfree();
  }

  void memalloc(int l, int m, int n){
    if(!allocated){
      lsize = l;
      msize = m;
      nsize = n;
      ksize = 300;  // fixed bin size for TKE, every 1 MeV

      L.memalloc(lsize,nsize,ksize);
      H.memalloc(lsize,nsize,ksize);

      for(int p=0 ; p<2 ; p++) multiplicity[p] =  eaverage[p] = etotal[p] = 0.0;

      for(int p=0 ; p<2 ; p++){
        spectrum[p] = new double[nsize];
        for(int i=0 ; i<nsize ; i++) spectrum[p][i] = 0.0;
      }
      chi = new double [nsize];
      for(int i=0 ; i<nsize ; i++) chi[i] = 0.0;

      preMassYield  = new double [msize];
      postMassYield = new double [msize];
      javeragedist  = new double [msize];
      for(int i=0 ; i<msize ; i++) preMassYield[i] =  postMassYield[i] = javeragedist[i] = 0.0;


      for(int i=0 ; i<2 ; i++){
        multiplicitydist[i] = new double [msize];
        eaveragedist[i] = new double [msize];
        for(int j=0 ; j<msize ; j++) multiplicitydist[i][j] = eaveragedist[i][j] = 0.0;
      }

      Pn  = new double [lsize];
      for(int i=0 ; i<lsize ; i++) Pn[i] = 0.0;

      yTKE = new double [ksize];
      nTKE = new double [ksize];
      for(int i=0 ; i<ksize ; i++) yTKE[i] = nTKE[i] = 0.0;

      eprefis = 0.0;
      nprefis = 0.0;

      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      L.memfree();
      H.memfree();
      for(int p=0 ; p<2 ; p++){
        delete [] multiplicitydist[p];
        delete [] eaveragedist[p];
        delete [] spectrum[p];
      }
      delete [] chi;
      delete [] preMassYield;
      delete [] postMassYield;
      delete [] Pn;
      delete [] yTKE;
      delete [] nTKE;

      allocated = false;
    }
  }

  int getLsize(){ return lsize; }
  int getMsize(){ return msize; }
  int getNsize(){ return nsize; }
  int getKsize(){ return ksize; }
};



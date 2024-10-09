#ifndef __CONSTANT_H__
#define __CONSTANT_H__
#include "constant.h"
#endif

#ifndef __GDR_H__
#define __GDR_H__
#include "gdr.h"
#endif

#ifndef __HFINTERFACE_H__
#define __HFINTERFACE_H__
#include "HFinterface.h"
#endif

//------------------------------------------------------------------------------
//      Enum

enum particle{
  gammaray = 0, neutron = 1, proton  = 2, alpha   = 3, deuteron= 4,
  triton   = 5, helion  = 6, userdef = 7, fission = 8, unknown = 9};
typedef enum particle Particle;

enum dirtype{
  none = 0, dwba = 1, ccrot = 2, ccvibg = 3, ccvib1 = 3, ccvib2 = 4};
typedef enum dirtype DirType;


//------------------------------------------------------------------------------
//     Class

/**********************************************************/
/*   Z and A numbers of Nucleus                           */
/**********************************************************/
class ZAnumber{ 
 private:
    unsigned int Z;
    unsigned int A;
 public:
    ZAnumber(){
      Z = 0;
      A = 0;
    }
    ZAnumber(int z, int a){
      Z = z;
      A = a;
    }
    void setZA(int z, int a){
      Z = z;
      A = a;
    }
    unsigned int getZ(){ return (Z); }
    unsigned int getA(){ return (A); }
    unsigned int getN(){ return (A-Z); }
    ZAnumber operator+(ZAnumber x){
      ZAnumber y;
      y.Z = Z + x.Z;
      y.A = A + x.A;
      return y;
    }
    ZAnumber operator-(ZAnumber x){
      ZAnumber y;
      y.Z = Z - x.Z;
      y.A = A - x.A;
      return y;
    }
    bool operator==(ZAnumber x){
      bool z = false;
      if( (Z == x.Z) && (A == x.A) ) z = true;
      return z;
    }
    bool operator!=(ZAnumber x){
      bool z = false;
      if( (Z != x.Z) || (A != x.A) ) z = true;
      return z;
    }
};


/**********************************************************/
/*   Array for Even and Odd Parities                      */
/**********************************************************/
class Parity{
 public:
    double even;
    double odd;

    Parity(){
      even = 0.0;
      odd  = 0.0;
    }

    Parity(double a, double b){
      even = a;
      odd  = b;
    }
};


/**********************************************************/
/*   Particle Data, spin, mss, optical potential ID       */
/**********************************************************/
class Pdata{
 public:
    Particle  pid;            // particle identifier
    ZAnumber  za;             // particle mass and atomic number
    int       omp;            // optical potential index
    int       spin2;          // intrinsic spin doubled
    double    mass;           // exact particle mass
    double    mass_excess;    // particle mass excess

    Pdata(){
      pid         = unknown;
      za.setZA(0,0);
      omp         = 0;
      spin2       = 0;
      mass        = 0.0;
      mass_excess = 0.0;
    }
};


/**********************************************************/
/*    System Data                                         */
/**********************************************************/
class System{
 public:
    ZAnumber   target;        // target Z and A numbers
    ZAnumber   compound;      // compound nucleus Z and A
    Pdata      incident;      // incident particle
    int        target_id;     // index of target nucleus
    int        inc_id;        // index of incident channel
    int        target_level;  // level index of excited target
    int        max_compound;  // number of compound nucleus
    int        max_particle;  // max number of patricle to be emitted
    int        uniq_compound; // number of unique compounds
    int        max_channel;   // actual maximum number of channes
    double     input_energy;  // energy in input file
    double     cms_energy;    // incident CMS energy
    double     lab_energy;    // incident LAB energy
    double     excitation;    // target excitation energy
    double     ex_total;      // total excitation energy
    double     energy_bin;    // default energy bin width
    double     energy_bin_in; // input energy bin width
    double     reduced_mass;  // reduced mass for entrance channel
    double     wave_number;   // wave number
    double     beta2;         // ground state deformation
    int       *bandidx;       // collective band index

    System(){
      bandidx = new int [MAX_LEVELS];
      init();
    }
    ~System(){
      delete [] bandidx;
    }

    void init(){
      target.setZA(0,0);
      compound.setZA(0,0);
      incident.za.setZA(0,0);
      incident.pid   = unknown;
      target_id      = 0;
      inc_id         = 0;
      target_level   = 0;
      max_compound   = 0;
      max_particle   = MAX_PARTICLE;
      uniq_compound  = 0;
      max_channel    = 1;
      input_energy   = 0.0;
      cms_energy     = 0.0;
      lab_energy     = 0.0;
      excitation     = 0.0;
      ex_total       = 0.0;
      energy_bin     = 0.0;
      energy_bin_in  = 0.0;
      reduced_mass   = 0.0;
      wave_number    = 0.0;
      beta2          = 0.0;
      for(int i=0 ; i<MAX_LEVELS ; i++) bandidx[i] = -1;
    }
};


/**********************************************************/
/*    Discrete Levels, energy, spin, parity               */
/**********************************************************/
class Level{
 public:
    int        flag;          // special flag not defined in RIPL
    double     energy;        // level excitation energy
    double     spin;          // level spin (not doubled)
    int        parity;        // parity
    int        ngamma;        // number of Gamma-rays
    int       *fstate;        // final state index for g-decay
    double    *branch;        // brancing ratio
    double     halflife;      // half-life of the state
    double    *gratio;        // net gamma-ray considering ICC

    Level(){
      flag     = 0;
      energy   = 0.0;
      spin     = 0.0;
      parity   = 0;
      ngamma   = 0;
      halflife = 0.0;
      fstate   = nullptr;
      branch   = nullptr;
      gratio   = nullptr;
    }
};


/**********************************************************/
/*   Channel Depend Data for Compound Reaction            */
/**********************************************************/
class Channel{
 public:
    Particle pid;             // particle identifier
    int      spin2;           // particle spin x2
    int      next;            // index of residual nucleus
    int      lmax;            // maximal angular momentum
    double   energy;          // energy of the emitted particle
    double   binding_energy;  // binding energy of the particle
    bool     status;          // channel status  open / close

    Channel(){
      init();
    }

    void init(){
      pid            = unknown;
      spin2          = 0;
      next           = 0;
      lmax           = 0;
      energy         = 0.0;
      binding_energy = 0.0;
      status         = false;
    }
};


/**********************************************************/
/*    Direct Reaction Data                                */
/**********************************************************/
class Direct{
 public:
    int      nlev;                      // number of total direct levels
    int      ncc;                       // number of coupled levels
    DirType  type[MAX_DIRECT];          // DWBA, coupled-channels
    double   nonloc;                    // non-locality correction
    double   defstatic[MAX_LAMBDA];     // static deformation parameters
    double   defdynamic[MAX_LAMBDA];    // dynamical deformation parameters
    double   defdwba[MAX_DIRECT];       // deformation for DWBA
    Level    lev[MAX_DIRECT];           // discrete levels

    Direct(){
      init();
    }
    void init(){
      nlev           = 0;
      ncc            = 0;
      nonloc         = 0.0;
      for(int i=0 ; i<MAX_LAMBDA ; i++){
        defstatic[i]  = 0.0;
        defdynamic[i] = 0.0;
      }
      for(int i=0 ; i<MAX_DIRECT ; i++){
        defdwba[i] = 0.0;
        type[i] = none;
      }
    }
};


/**********************************************************/
/*    Level Density Parameters                            */
/**********************************************************/
class LevelDensity{
 public:
    double     match_energy;  // matching energy
    double     shell_correct; // shell correction energy
    double     temperature;   // T parameter
    double     E0;            // eergy shift
    double     a;             // a parameter
    double     pairing_energy;// pairing energy
    double     spin_cutoff;   // spin cut-off factor
    double     sigma0;        // spin cut-off factor from levels

    LevelDensity(){
      match_energy   = 0.0;
      shell_correct  = 0.0;
      temperature    = 0.0;
      E0             = 0.0;
      a              = 0.0;
      pairing_energy = 0.0;
      spin_cutoff    = 0.0;
      sigma0         = 0.0;
    }
};


/**********************************************************/
/*   Direct/Semidirect Capture Parameters                 */
/**********************************************************/
class Dcapt{
 public:
    double       v1;          // real isovector potential 
    double       w1;          // imaginary isovector potential
    double       nonloc;      // non-locality correction
    double       bmax;        // max of unbound level energy
    double       fshift;      // shift Fermi enerrgy
    GDR         *gdr;         // GDR parameters given in input
    bool         hf;          // use Hartree-Fock results
    HFInterface *hfsp;        // Hartree-Fock s.p. orbits

    Dcapt(){
      init();
      gdr    = new GDR [MAX_GDR];
    }
    ~Dcapt(){
      delete [] gdr;
    }
    void init(){
      v1     = 0.0;
      w1     = 0.0;
      nonloc = 0.0;
      bmax   = 1.5;
      fshift = 0.0;
      hf     = false;
    }
};


/**********************************************************/
/*   Fission Barrier Parameters                           */
/**********************************************************/
class KBand{
 public:
    int      parity;          // transition state parity
    int      k2;              // band K (doubled)
    double   excitation;      // excitation energy

    KBand(){
      parity     = 0;
      k2         = 0;
      excitation = 0.0;
    }
};

class Barrier{
 public:
    double     height;        // fission barrier height
    double     curvature;     // barrier curvature
    double     inertia;       // inertia parameter, hbar^2/2I
    double     elmax;         // highest discrete state
    int        nband;         // number of transition band-head
    KBand     *kband;         // transition states band-head
    Parity   **density;       // level density on top of barrier

    Barrier(){
      clear();
      kband     = nullptr;
      density   = nullptr;
    }

    void clear(){
      height    = 0.0;
      curvature = 0.0;
      inertia   = 0.0;
      elmax     = 0.0;
      nband     = 0;
    }
};

class PotWell{
 public:
    double     height;        // fission barrier height
    double     curvature;     // barrier curvature
    double     absorb;        // absorption potential depth

    PotWell(){
      clear();
    }

    void clear(){
      height    = 0.0;
      curvature = 0.0;
      absorb    = 0.0; 
    }
};

class FisEnhance{
 public:
    double     energy;        // enhancement of fission
    double     width;
    double     peak;
    double     compfact1p;     // level compression parameter 1
    double     compfact2p;     // level compression parameter 2
    double     compfact1n;     // parameter 1 for negative parity
    double     compfact2n;     // parameter 2 for negative parity

    FisEnhance(){
      energy = 0.0;
      width  = 0.0;
      peak   = 0.0;
      compfact1p = 0.0;
      compfact2p = 0.0;
      compfact1n = 0.0;
      compfact2n = 0.0;
    }
};

class Fission{
 private:
    bool       arrayalloc;    // flag for memory allocation
 public:
    ZAnumber   za;            // Z and A numbers
    int        hump;          // number of fission barriers
    Barrier    *barrier;      // fission barriers
    PotWell    *potwell;      // class-II,III wells
    FisEnhance fisenhance;    // fission enhance (Lorentz shape)

    Fission(){
      za.setZA(0,0);
      hump       = 0;
      barrier    = nullptr;
      potwell    = nullptr;
      arrayalloc = false;
    }
    void memalloc(){
      if(!arrayalloc){
        barrier = new Barrier [MAX_HUMP];
        potwell = new PotWell [MAX_HUMP - 1];
        for(int j=0 ; j<MAX_HUMP ; j++){
          barrier[j].kband   = new KBand [MAX_KBAND];
          barrier[j].density = new Parity *[MAX_ENERGY_BIN];
          for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
            barrier[j].density[k] = new Parity[MAX_J];
          }
        }
        arrayalloc = true;
      }
    }
    ~Fission(){
      if(arrayalloc){
        for(int j=0 ; j<MAX_HUMP ; j++){
          for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
            delete [] barrier[j].density[k];
          }
          delete [] barrier[j].density;
          delete [] barrier[j].kband;
        }
        delete [] barrier;
        delete [] potwell;
        arrayalloc = false;
      }
    }
    void init(){
      za.setZA(0,0);
      for(int j=0 ; j<MAX_HUMP   ; j++) barrier[j].clear();
      for(int j=0 ; j<MAX_HUMP-1 ; j++) potwell[j].clear();
    }
};


/**********************************************************/
/*    Nucleus Data                                        */
/**********************************************************/
class Nucleus{
 private:
    int        ndiscmax;      // max number of discrete levels
    int        ncontmax;      // max number of contonuum bins
    bool       arrayalloc;    // flag for memory allocation
 public:
    ZAnumber   za;            // Z and A numbers
    double     max_energy;    // initial excitation energy
    double     mass;          // exact mass
    double     mass_excess;   // mass excess
    double     de;            // energy bin width in continuum
    LevelDensity ldp;         // level density parameter
    Level     *lev;           // discrete level data
    double    *excitation;    // excitation energy
    Parity   **density;       // level density
    Parity   **pop;           // population
    double    *lpop;          // leve population
    Channel   *cdt;           // channel data
    Fission    *fission;      // fission barriers
    double     *popfis;       // population decreased by fission
    GDR        *gdr;          // GDR parameters
    int        jmax;          // cut-off of angular momentum
    int        ndisc;         // number discrete levels
    int        ncont;         // number of continuum enerby bins
    int        ntotal;        // total number of enerby bins
    bool       fissile;       // flag for fissile nuclide

    Nucleus(){
      za.setZA(0,0);
      max_energy     = 0.0;
      mass           = 0.0;
      mass_excess    = 0.0;
      de             = 0.0;
      jmax           = 0  ;
      ndisc          = 0  ;
      ncont          = 0  ;
      ntotal         = 0  ;
      fission        = nullptr;
      fissile        = false;

      ndiscmax       = 0;
      ncontmax       = 0;
      arrayalloc     = false;
    }

    void memalloc(int n, int m){
      if(!arrayalloc){
        ndiscmax = n;
        ncontmax = m;
        density        = new Parity * [ncontmax];
        pop            = new Parity * [ncontmax];
        popfis         = new double [ncontmax];
        for(int j=0 ; j<ncontmax ; j++){
          density[j]     = new Parity [MAX_J];
          pop[j]         = new Parity [MAX_J];
        }
        lpop           = new double [ndiscmax];
        excitation     = new double [ncontmax];
        cdt            = new Channel [MAX_CHANNEL];
        gdr            = new GDR [MAX_GDR];
        arrayalloc = true;
      }
    }

    void memclear(){
      if(arrayalloc){
        Parity p(0.0,0.0);
        for(int k=0 ; k<ncontmax ; k++){
          excitation[k] = 0.0;
          popfis[k] = 0.0;
          for(int j=0 ; j<MAX_J ; j++){
            pop[k][j]     = p;
            density[k][j] = p;
          }
        }
        for(int k=0 ; k<ndiscmax ; k++) lpop[k] = 0.0;
        for(int k=0 ; k<MAX_GDR ; k++) gdr[k].clear();
      }
    }

    ~Nucleus(){
      if(arrayalloc){
        for(int j=0 ; j<ncontmax ; j++){
          delete [] density[j];
          delete [] pop[j];
        }
        delete [] density;
        delete [] pop;
        delete [] popfis;

        delete [] lpop;
        delete [] excitation;

        delete [] cdt;
        delete [] gdr;
        arrayalloc = false;
      }
    }
};


/**********************************************************/
/*    Nuclear Structure Data, Discrete Levels             */
/**********************************************************/
class NuclearStructure{
 private:
    int        nlevmax;       // max number of levels allocated
    bool       arrayalloc;    // flag for memory allocation
 public:
    ZAnumber   za;            // Z and A numbers
    Level     *lev;           // discrete level data
    int        nlevel;        // number discrete levels

    NuclearStructure(){
      za.setZA(0,0);
      nlevel     = 0;
      nlevmax    = 0;
      arrayalloc = false;
    }
    void memalloc(int n){
      if(!arrayalloc){
        nlevmax = n;
        lev = new Level [nlevmax];
        for(int j=0 ; j<nlevmax ; j++){
          lev[j].fstate  = new int[MAX_GAMMA_BRANCH];
          lev[j].branch  = new double[MAX_GAMMA_BRANCH];
          lev[j].gratio  = new double[MAX_GAMMA_BRANCH];
        }
        arrayalloc = true;
      }
    }
    ~NuclearStructure(){
      if(arrayalloc){
        for(int j=0 ; j<nlevmax ; j++){
          delete [] lev[j].fstate;
          delete [] lev[j].branch;
          delete [] lev[j].gratio;
        }
        delete [] lev;
        arrayalloc = false;
      }
    }
};


/**********************************************************/
/*   Calculated Cross Sections                            */
/**********************************************************/
class ParticleCount{
 public:
    unsigned char par[7];     // number of particles emitted
    double        xsec;       // production cross section
    double        fiss;       // fission cross section

    ParticleCount(){
      init();
    }
    void init(){
      for(int i=0 ; i<7 ; i++) par[i]=0;
      xsec = 0.0;
      fiss = 0.0;
    }
};

class CrossSection{
 private:
    bool    arrayalloc;       // flag for memory allocation
    bool    arrayastro;       // memory for MACS
    bool    arrayspect;       // memory for energy spectra
    int     n_reacrate;       // number of reaction rates
    int     n_spectra;        // number of different spectra
    int     n_specbin;        // number of bins for spectra
 public:
    double  total;            // total cross section
    double  elastic;          // elastic cross section
    double  reaction;         // non-elastic cross section
    double  compound;         // compound elastic cross section
    double  dsd;              // DSD capture cross section
    double  preeq;            // total preequilibrium emission
    double  totaldir;         // sum of direct inelastic
    double  *direct;          // direct inelastic cross sections
    double  **levexcite;      // discrete level excitation
    double  *theta;           // scattering angles
    double  *costh;           // cosine of scattering angles
    double  **angdist;        // scattering angular distribution
    double  ***legcoef;       // Legendre coeff. for CN ang. dist
    double  **spectra;        // total particle energy spectra
    double  **transmission;   // T(lj) for in-coming particle
    double  **macs;           // save cross sections for MACS
    ParticleCount *prod;      // individual reaction cross section

    CrossSection(){
      init();
      arrayalloc = false;
      arrayastro = false;
      arrayspect = false;
      n_reacrate = 0;
      n_spectra  = 0;
      n_specbin  = 0;
    }

    ~CrossSection(){
      if(arrayalloc){
        for(int i=0 ; i<MAX_DIRECT ; i++){
          delete [] transmission[i];
        }
        delete [] transmission;

        delete [] direct;
        delete [] prod;
        delete [] theta;
        delete [] costh;
        for(int j=0 ; j<MAX_ANGDISTLEVELS ; j++){
          delete [] angdist[j];
        }
        delete [] angdist;

        for(int i=0 ; i<MAX_CHANNEL ; i++){
          for(int j=0 ; j<MAX_ANGDISTLEVELS ; j++){
            delete [] legcoef[i][j];
          }
          delete [] levexcite[i];
          delete [] legcoef[i];
        }
        delete [] levexcite;
        delete [] legcoef;
        arrayalloc = false;
      }

      astrofree();
      spectfree();
    }

    void init(){
      total     = 0.0;
      elastic   = 0.0;
      reaction  = 0.0;
      compound  = 0.0;
      dsd       = 0.0;
      preeq     = 0.0;
      totaldir  = 0.0;
    }

    void memalloc(){
      if(!arrayalloc){

        transmission = new double * [MAX_DIRECT];
        for(int i=0 ; i<MAX_DIRECT ; i++){
          transmission[i] = new double [MAX_J*3];
        }

        direct    = new double [MAX_DIRECT];
        prod = new ParticleCount [MAX_COMPOUND];

        theta   = new double [MAX_ANGDIST];
        costh   = new double [MAX_ANGDIST];
        angdist = new double * [MAX_ANGDISTLEVELS];
        for(int j=0 ; j<MAX_ANGDISTLEVELS ; j++){
          angdist[j] = new double [MAX_ANGDIST];
        }

        levexcite = new double * [MAX_CHANNEL];
        legcoef   = new double ** [MAX_CHANNEL];
        for(int i=0 ; i<MAX_CHANNEL ; i++){
          levexcite[i] = new double [MAX_LEVELS];
          legcoef[i]   = new double * [MAX_ANGDISTLEVELS];
          for(int j=0 ; j<MAX_ANGDISTLEVELS ; j++){
            legcoef[i][j] = new double [MAX_J];
          }
        }
        arrayalloc = true;
      }
    }

    void astroalloc(int nrate, int ngrid){
      if(!arrayastro){
        macs = new double * [nrate];
        for(int i=0 ; i<nrate ; i++){
          macs[i] = new double [ngrid];
        }
        n_reacrate = nrate;
        arrayastro = true;
      }
    }

    void spectalloc(int nc, int ne){
      if(!arrayspect){
        spectra = new double * [nc];
        for(int j=0 ; j<nc ; j++){
          spectra[j] = new double [ne];
        }
        n_spectra = nc;
        n_specbin = ne;
        arrayspect = true;
      }
    }

    void specclear(){
      if(arrayspect){
        for(int i=0 ; i<n_spectra ; i++){
          for(int k=0 ; k<n_specbin ; k++) spectra[i][k] = 0.0;
        }
      }
    }

    void astrofree(){
      if(arrayastro){
        for(int i=0 ; i<n_reacrate ; i++){
          delete [] macs[i];
        }
        delete [] macs;
        n_reacrate = 0;
        arrayalloc = false;
      }
    }

    void spectfree(){
      if(arrayspect){
        for(int i=0 ; i<n_spectra ; i++){
          delete [] spectra[i];
        }
        delete [] spectra;
        n_spectra = 0;
        arrayspect = false;
      }
    }

    int getNspectra(){return n_spectra;}
    int getNspecbin(){return n_specbin;}
};


/**********************************************************/
/*   Optional Gamma-Ray Production Cross Sections         */
/**********************************************************/
class GammaLine{
 public:
    ZAnumber   za;            // origin Z,A
    double  energy;           // gamma-ray energy
    double  production;       // production cross section

    GammaLine(){
      init();
    }

    void init(){
      za.setZA(0,0);
      energy     = 0.0;
      production = 0.0;
    }
};


class GammaProduction{
 private:
    int     ngmax;            // maximum number of gamma-rays
    int     ncurr;            // current number of gamma-rays
    bool    arrayalloc;       // flag for memory allocation
    double  cutoff;           // cut off for production
 public:
    GammaLine *line;          // line gammas

    GammaProduction(){
      ngmax = 0;
      ncurr = 0;
      cutoff = 0.0;
      arrayalloc = false;
    }

    ~GammaProduction(){
      memfree();
    }

    void memalloc(int n, double c){
      if(!arrayalloc){
        ngmax = n;
        line  = new GammaLine [ngmax];
        arrayalloc = true;
        init();
      }
      cutoff = c;
    }

    void memfree(){
      if(arrayalloc){
        delete [] line;
        arrayalloc = false;
      }
    }

    void init(){
      if(arrayalloc){
        ncurr = 0;
        for(int i=0 ; i<ngmax ; i++) line[i].init();
      }
    }

    int push(ZAnumber za, double e, double p){
      for(int i=0 ; i<ncurr ; i++){
        if( (line[i].energy == e) && (line[i].za == za) ){
          line[i].production += p;
          return ncurr;
        }
      }

      if( (ncurr < ngmax) && (p >= cutoff) ){
        line[ncurr].za = za;
        line[ncurr].energy = e;
        line[ncurr].production = p;
        ncurr++;
      }
      return ncurr;
    }

    int getN(){ return ncurr; }
};


/**********************************************************/
/*   Prompt Fission Neutron Spectrum Data                 */
/**********************************************************/
class FFragment{
 public:
    ZAnumber     za;          // fragment Z and A numbers
    LevelDensity ldp;         // fragment level density
    double       a;           // actual a-parameter
    double       tmax;        // maximum temperature

    FFragment(){
      init();
    }

    void init(){
      za.setZA(0,0);
      a          = 0.0;
      tmax       = 0.0;
    }
};


class FChance{
 public:
    FFragment  lf;            // light fission fragment
    FFragment  hf;            // heavy fission fragment
    double     ac;            // CN level density parameter
    double     tke;           // total kinetic energy
    double     tke0;          // TKE at zero energy
    double     tke1;          // energy dependence of TKE tke = tke0 + tke1 En
    double     txe;           // total excitation energy
    double     tmax;          // CN maximum temperature
    double     etotal;        // total energy release
    double     egamma;        // average total gamma energy
    double     exfiss;        // average excitation energy causing fission
    double     rt;            // anisothermal parameter, RT
    double     meansep;       // average Sn
    double     nuratio;       // nuL / nuH for average
    double     nubar;         // nu at each chance fission
    double     ecms[5];       // average CMS neutron energy
    double     elab;          // average LAB neutron energy
    double     tps;           // parameter for temperature pdf
    double     anisotropy;    // anisotropy parameter

    FChance(){
      init();
    }

    void init(){
      lf.init();
      hf.init();
      ac         = 0.0;
      tke        = 0.0;
      tke0       = 0.0;
      tke1       = 0.0;
      txe        = 0.0;
      tmax       = 0.0;
      etotal     = 0.0;
      egamma     = 0.0;
      exfiss     = 0.0;
      rt         = 1.0;
      meansep    = 0.0;
      nuratio    = 1.0;
      nubar      = 0.0;
      elab       = 0.0;
      tps        = 1.0;
      anisotropy = 0.0;
      for(int i=0 ; i<5 ; i++){
        ecms[i]  = 0.0;
      }
    }
};


class FNSpec{
 private:
    int        nf;            // max number of fission chances
 public:
    FChance   *fc;            // chance fission-dependent data
    int        omindex;       // optical model potential index
    double     maxwell;       // Maxwellian temperature for ratio calculation

    FNSpec(const int n){
      nf = n;
      fc = new FChance [nf];
      init();
    }

    ~FNSpec(){
      delete [] fc;
    }

    void init(){
      omindex    = 0;
      maxwell    = 0.0;
      for(int i=0 ; i<nf ; i++) fc[i].init();
    }

    int maxchance(){ return nf; }
};


/**********************************************************/
/*      Mean-Field Calculation Printarameters             */
/**********************************************************/
class MFTparm{
 public:
  ZAnumber target;
  int      n0;                // number of Harmonic Oscillator basis
  double   b;                 // b = sqrt(hbar/m w)
  double   q;                 // omega ratio

  /*** parameters for Hartree-Fock */
  double   value[5];
  double   weight[5];
  int      nfree[5];
  ZAnumber nucleus;
  std::string force;
  std::string pairing;
  double   strength[2];
  int      max_iteration;
  double   eps_converge;

  /*** parameters for FRDM */
  double   eps[6];
  double   gamma;

  MFTparm(){
    init();
  }

  void init(){
    b             = 0.5;
    q             = 1.0;
    n0            = 14;

    force         = "";
    pairing       = "";
    strength[0]   = 0.0;
    strength[1]   = 0.0;
    max_iteration = 200;
    eps_converge  = 1.0e-05;

    for(int i=0 ; i<5 ; i++){
      value[i]  = 0.0;
      weight[i] = 0.0;
      nfree[i]  = 0;
    }

    for(int i=0 ; i<6 ; i++) eps[i] = 0.0;
    gamma = 0.0;
  }
};

/******************************************************************************/
/*  fileout.cpp                                                               */
/*        write calculated results on files                                   */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "structur.h"
#include "nucleus.h"
#include "fileout.h"
#include "global.h"
#include "metastable.h"
#include "ripl2levels.h"
#include "terminate.h"

static const int N_PART_KIND       =  7;
static const int N_OUTPUT_LEVELS   = 42;

#undef LONGOUTPUT
#ifdef LONGOUTPUT
static const int COLUMN_WIDTH      = 21;
static const int DECIMAL_WIDTH     = 14;
#else
static const int COLUMN_WIDTH      = 13;
static const int DECIMAL_WIDTH     =  6;
#endif

static const int N_OUTPUT_SIGMA = 25;
  static const std::string reaction_name[N_OUTPUT_SIGMA] = {
  " (x,-)   ",
  " (x,n)   "," (x,p)   "," (x,A)   "," (x,d)   "," (x,t)   "," (x,h)   ",
  " (x,2n)  "," (x,np)  "," (x,nA)  "," (x,nd)  "," (x,nt)  "," (x,nh)  ",
  " (x,2p)  "," (x,pA)  "," (x,pd)  "," (x,2A)  "," (x,Ad)  "
  " (x,3n)  "," (x,2np) "," (x,2nA) "," (x,npA) "," (x,n2p) ",
  " (x,4n)  "," (x,5n)  "};

static const std::string particle_number[N_OUTPUT_SIGMA] = {
  "000000",
  "100000", "010000", "001000", "000100", "000010", "000001",
  "200000", "110000", "101000", "100100", "100010", "100001",
  "020000", "011000", "010100", "002000", "001100",
  "300000", "210000", "201000", "111000", "120000",
  "400000", "500000"};


static void    fioWriteCrossSectionMain        (const bool, const std::string, const int, const double);
static void    fioWriteCrossSectionParticle    (const bool, const std::string, const int, const double);
static void    fioWriteCrossSectionRadioactive (const bool, const std::string, const int, const double);
static void    fioWriteCrossSectionFission     (const bool, const std::string, const int, const double);
static void    fioWriteLevelExcite             (const bool, const std::string, double *, const double, const int, Nucleus *);
static void    fioWriteAngularDistribution     (const std::string, double *, const double);
static void    fioWriteLegendreCoefficient     (const std::string, double *, const double, const int);

static inline void   fioSetLevelEnergy         (double *, Nucleus *);
static inline bool   fioCheckExist             (std::string);
static inline std::string fioNparticle              (unsigned char *);


/*** lower limit of output values, enforce zero */
static const double eps = 1.0e-99;


/**********************************************************/
/*     Output Cross Sections and Angular Distributions    */
/**********************************************************/
void cohFileOut(int nc, int tid, double elab)
// nc    : number of compound nucleus
// tid   : index for the target
// elab  : incident LAB energy 
{
  std::string fname;
  std::ostringstream fn;
  double elev[N_OUTPUT_LEVELS];

  /*** print out angular distributions for elastic and inelastic scattering */
  if(prn.angdist){

    fioSetLevelEnergy(elev,&ncl[tid]);

    fn.str("");
    fn << filename_angular_distribution << "." << file_extension;
    fioWriteAngularDistribution(fn.str(), elev, elab);
  }

  /*** print out cross sections */
  else{
    /*** main cross sections, total, elastic, capture, inelasitc, fission  */
    fn.str("");
    fn << filename_cross_section << "." << file_extension;
    fioWriteCrossSectionMain(fioCheckExist(fn.str()), fn.str(), nc, elab);

    /*** all reaction cross sections */
    fn.str("");
    fn << filename_particle_production << "." << file_extension;
    fioWriteCrossSectionParticle(fioCheckExist(fn.str()), fn.str(), nc, elab);

    /*** radioactive isotopes */
    fn.str("");
    fn << filename_radioactive_production << "." << file_extension;
    fioWriteCrossSectionRadioactive(fioCheckExist(fn.str()), fn.str(), nc, elab);

    /*** multi-chance fission cross sections */
    if(ctl.fission){
      fn.str("");
      fn << filename_fission << "." << file_extension;
      fioWriteCrossSectionFission(fioCheckExist(fn.str()), fn.str(), nc, elab);
    }
  }

  /*** for binary reactions, scattering to discrete levels */
  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;

    int p = ncl[0].cdt[id].next;
    if(ncl[p].ndisc == 0) continue;

    fioSetLevelEnergy(elev,&ncl[p]);

    fn.str("");
    if(prn.angdist){
      fn << filename_legendre_coefficient << std::setw(1) << id << "." << file_extension;
      fioWriteLegendreCoefficient(fn.str(), elev, elab, id);
    }
    else{
      fn << filename_level_excite << std::setw(1) << id << "." << file_extension;
      fioWriteLevelExcite(fioCheckExist(fn.str()), fn.str(), elev, elab, id, &ncl[p]);
    }
  }
}


/**********************************************************/
/*     Store Discrete Level Energies in Array             */
/**********************************************************/
void fioSetLevelEnergy(double *elev, Nucleus *n)
{
  for(int k=0 ; k<N_OUTPUT_LEVELS ; k++) elev[k] = 0.0;

  int m = N_OUTPUT_LEVELS-1;
  if(n->ndisc < m) m = n->ndisc;

  for(int k=1 ; k<m ; k++) elev[k] = n->lev[k].energy;
}


/**********************************************************/
/*     Check If File Already Exists in the Current Dir    */
/**********************************************************/
bool fioCheckExist(std::string fname)
{
  std::ifstream fin;
  bool fexist = false;

  fin.open(&fname[0],std::ios::in);
  if(!fin) fexist = false;
  else     fexist = true;
  fin.close();

  return(fexist);
}


/**********************************************************/
/*     Convert ParticleCount array to std::string              */
/**********************************************************/
std::string fioNparticle(unsigned char *par)
{
  std::ostringstream os;
  for(int j=1 ; j<N_PART_KIND ; j++){
    if(j >= MAX_CHANNEL) os << std::setw(1) << 0;
    else                 os << std::setw(1) << (int)par[j];
  }
  return os.str();
}


/**********************************************************/
/*     Output Selected Cross Sections on File             */
/**********************************************************/
void fioWriteCrossSectionMain(const bool fexist, const std::string fname, const int nc, const double e)
{
  const int N_OUTPUT_CROSSSEC = 6;
  static const std::string rname[N_OUTPUT_CROSSSEC] = {
  " Total       "," Elastic     "," CNFormation "," Capture     "," Inelastic   "," Fission     "};

  std::ofstream fp;
  fp.open(fname.c_str(),std::ios::out | std::ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSectionMain");
  }

  double cx[N_OUTPUT_CROSSSEC];

  for(int i=0 ; i<N_OUTPUT_CROSSSEC ; i++) cx[i]=0.0;

  cx[0] = crx.total;                    // total cross section
  cx[1] = crx.elastic + crx.compound;   // elastic scattering cross section
  cx[2] = crx.reaction;                 // compond formation cross section

  /*** capture, inelastic, and fission */
  for(int i=0 ; i<nc ; i++) {

    std::string str = fioNparticle(crx.prod[i].par);
    if(     str == "000000"){
      cx[3] = crx.prod[i].xsec;
    }
    else if(str == "100000"){
      double z =crx.prod[i].xsec - crx.compound;
      cx[4] = (z < 1.0e-10) ? 0.0 : z;
    }
    cx[5] += crx.prod[i].fiss;
  }

  /*** print header if file not exist */
  if(!fexist){
    fp << "#" << std::setw(COLUMN_WIDTH-1) << " Reaction    ";
    for(int i=0 ; i<N_OUTPUT_CROSSSEC ; i++) fp << std::setw(COLUMN_WIDTH) << rname[i];
    fp << std::endl;
  }

  /*** print cross sections */
  fp << std::setprecision(DECIMAL_WIDTH);
  fp << std::setiosflags(std::ios::scientific);
  fp << std::setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_CROSSSEC ; i++){
    fp << std::setw(COLUMN_WIDTH) << ( (cx[i] < eps) ? 0.0 : cx[i] );
  }
  fp << std::endl;
  fp.close();
}


/**********************************************************/
/*     Output Particle Production Cross Sections          */
/**********************************************************/
void fioWriteCrossSectionParticle(const bool fexist, const std::string fname, const int nc, const double e)
{
  std::ofstream fp;
  fp.open(fname.c_str(),std::ios::out | std::ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSectionParticle");
  }

  double cx[N_OUTPUT_SIGMA];

  for(int k=0 ; k<N_OUTPUT_SIGMA ; k++) cx[k] = 0.0;

  for(int i=0 ; i<nc ; i++) {
    /*** count number of emitted particls */
    std::string str = fioNparticle(crx.prod[i].par);
    
    /*** look into the output table, if exists */
    for(int k=1 ; k<N_OUTPUT_SIGMA ; k++){
      if(str == particle_number[k]){ cx[k] = crx.prod[i].xsec; break; }
    }
  }

  /*** print header if file not exist */
  if(!fexist){
    fp << "#" << std::setw(COLUMN_WIDTH-1) << " Reaction   ";
    for(int i=0 ; i<N_OUTPUT_SIGMA-1 ; i++) fp << std::setw(COLUMN_WIDTH) << reaction_name[i+1];
    fp << std::endl;
  }

  /*** print cross sections */
  fp << std::setprecision(DECIMAL_WIDTH);
  fp << std::setiosflags(std::ios::scientific);
  fp << std::setw(COLUMN_WIDTH) << e;
  for(int k=1 ; k<N_OUTPUT_SIGMA ; k++){
    fp << std::setw(COLUMN_WIDTH) << ( (cx[k] < eps) ? 0.0 : cx[k] );
  }
  fp << std::endl;
  fp.close();
}


/**********************************************************/
/*     Output Radioactive Nuclide Production              */
/**********************************************************/
void fioWriteCrossSectionRadioactive(const bool fexist, const std::string fname, const int nc, const double e)
{
  const int metamax = 3; // number of meta-stable states to be printed, including g.s.

  std::ofstream fp;
  fp.open(fname.c_str(),std::ios::out | std::ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSectionRadioactive");
  }

  /*** initially allocate all the isotopes to make a complete list */
  Radioisotope *isotope = new Radioisotope [N_OUTPUT_SIGMA];
  for(int i=0 ; i<N_OUTPUT_SIGMA ; i++) isotope[i].memalloc(metamax);


  NuclearStructure ns;
  ns.memalloc(MAX_LEVELS);

  /*** first, search meta-stable states in RIPL */
  for(int i=0 ; i<N_OUTPUT_SIGMA ; i++){
    for(int m=0 ; m<isotope[i].getmaxlevel() ; m++) isotope[i].lev[m].init();

    /*** calculate residual ZA */
    ns.za = ncl[0].za;

    for(int j=0 ; j<(int)particle_number[i].length() ; j++){
      if(particle_number[i][j] != '0'){
        ZAnumber zp(0,0);
        switch(j){
        case 0: zp.setZA(0,1); break;
        case 1: zp.setZA(1,1); break;
        case 2: zp.setZA(2,4); break;
        case 3: zp.setZA(1,2); break;
        case 4: zp.setZA(1,3); break;
        case 5: zp.setZA(2,3); break;
        default: break;
        }
        for(int p=0 ; p<(int)(particle_number[i][j] - '0') ; p++){ ns.za = ns.za - zp; }
      }
    }

    /*** read RIPL */
    int nlev = riplReadDiscreteLevels(&ns.za,ns.lev,reassign);

    /*** scan discrete levels to see if radioactive */
    bool active = false;
    for(int k=0 ; k<nlev ; k++){
      if(ns.lev[k].halflife > thalfmin){ active = true; break; }
    }

    if(!active){
      isotope[i].radioactive = false; // flag to skip this isotope to be printed
      continue;
    }

    isotope[i].radioactive = true;
    isotope[i].za = ns.za;

    int m = 0;
    /*** always include ground state */
    isotope[i].lev[m].nx = 0;
    isotope[i].lev[m].energy = ns.lev[0].energy;
    isotope[i].lev[m].halflife = ns.lev[0].halflife;
    m ++; // move to meta-stable state

    /*** look at excited levels */
    for(int k=1 ; k<nlev ; k++){
      if(ns.lev[k].halflife > thalfmin){
        isotope[i].lev[m].nx = k;                        // order of this state
        isotope[i].lev[m].energy = ns.lev[k].energy;     // excitation energy of the state
        isotope[i].lev[m].halflife = ns.lev[k].halflife; // half-life of the state
        m ++; // increment meta-stable state counter
      }
      if(m >= isotope[i].getmaxlevel()) break;
    }

    /*** for each of residuals */
    for(int n=0 ; n<nc ; n++){

      /*** convert emitted particles to std::string */
      std::string str = fioNparticle(crx.prod[n].par);

      /*** this residual was calculated */
      if(str == particle_number[i]){

        for(int m=0 ; m<isotope[i].getmaxlevel() ; m++){
          isotope[i].lev[m].production = 0.0;
          for(int k=0 ; k<ncl[n].ndisc ; k++){
            if(k == isotope[i].lev[m].nx){
              isotope[i].lev[m].production = ncl[n].lpop[k];
              break;
            }
          }
        }

        /*** count number of meta-stable states */
        int nm = 0;
        for(int m=1 ; m<isotope[i].getmaxlevel() ; m++) if(isotope[i].lev[m].halflife != 0.0) nm ++;

        /*** when there is 1 metastable state, subtract this from the g.s. level production */
        if(nm == 1){
          isotope[i].lev[0].production -= isotope[i].lev[1].production;
        }
        /*** when there are 2, calculate m2 decay fraction to m1 state */
        else if(nm == 2){
          double *lpop = new double [ncl[n].ndisc];
          for(int k=0 ; k<ncl[n].ndisc ; k++) lpop[k] = 0.0;
          lpop[ isotope[i].lev[2].nx ] = 1.0;

          for(int i0=isotope[i].lev[2].nx ; i0>0 ; i0--){
            for(int j=0 ; j<ncl[n].lev[i0].ngamma ; j++){
              lpop[ ncl[n].lev[i0].fstate[j] ] += lpop[i0] * ncl[n].lev[i0].branch[j];
            }
          }
          double frac = lpop[ isotope[i].lev[1].nx ];

          isotope[i].lev[1].production -= isotope[i].lev[2].production * frac;
          if(isotope[i].lev[1].production < 0.0) isotope[i].lev[1].production = 0.0;

          isotope[i].lev[0].production -= isotope[i].lev[1].production + isotope[i].lev[2].production;
          delete [] lpop;
        }
        if(isotope[i].lev[0].production < 0.0) isotope[i].lev[0].production = 0.0;
      }
    }
  }

  /*** print header if file not exist */
  if(!fexist){

    fp << std::setprecision(DECIMAL_WIDTH);
    fp << std::setiosflags(std::ios::scientific);
    fp << "#" << std::setw(COLUMN_WIDTH-1) << " Nuclide Ex ";
    for(int i=0 ; i<N_OUTPUT_SIGMA ; i++){
      if(!isotope[i].radioactive) continue;
      int id = (signed int)isotope[i].za.getZ() * 1000 + (signed int)isotope[i].za.getA();
      fp << std::setw(COLUMN_WIDTH) << id;
      for(int m=1 ; m<isotope[i].getmaxlevel() ; m++){
        if(isotope[i].lev[m].nx < 0) fp << std::setw(COLUMN_WIDTH) << "----";
        else fp << std::setw(COLUMN_WIDTH) << isotope[i].lev[m].energy;
      }
    }
    fp << std::endl;

    fp << "#" << std::setw(COLUMN_WIDTH-1) << " Particle L ";
    for(int i=0 ; i<N_OUTPUT_SIGMA ; i++){
      if(!isotope[i].radioactive) continue;
      fp << std::setw(COLUMN_WIDTH) << particle_number[i];
      for(int m=1 ; m<isotope[i].getmaxlevel() ; m++){
        if(isotope[i].lev[m].nx < 0) fp << std::setw(COLUMN_WIDTH) << "----";
        else fp << std::setw(COLUMN_WIDTH) << isotope[i].lev[m].nx;
      }
    }
    fp << std::endl;

    fp << "#" << std::setw(COLUMN_WIDTH-1) << " HalfLife   ";
    for(int i=0 ; i<N_OUTPUT_SIGMA ; i++){
      if(!isotope[i].radioactive) continue;
      for(int m=0 ; m<isotope[i].getmaxlevel() ; m++){
        if(isotope[i].lev[m].halflife < 0.0) fp << std::setw(COLUMN_WIDTH) << 0.0;
        else fp << std::setw(COLUMN_WIDTH) << isotope[i].lev[m].halflife;
      }
    }
    fp << std::endl;
  }

  /*** print cross sections */
  fp << std::setprecision(DECIMAL_WIDTH);
  fp << std::setiosflags(std::ios::scientific);
  fp << std::setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_SIGMA ; i++){
    if(!isotope[i].radioactive) continue;
    for(int m=0 ; m<isotope[i].getmaxlevel() ; m++){
      fp << std::setw(COLUMN_WIDTH) << ( (isotope[i].lev[m].production < eps) ? 0.0 : isotope[i].lev[m].production );
    }
  }
  fp << std::endl;
  fp.close();

  for(int i=0 ; i<N_OUTPUT_SIGMA ; i++) isotope[i].memfree();
  delete [] isotope;
}


/**********************************************************/
/*     Output Fission Cross Sections on File              */
/**********************************************************/
void fioWriteCrossSectionFission(const bool fexist, const std::string fname, const int nc, const double e)
{
  const int N_OUTPUT_FISSION = 5;
  static const std::string rname[N_OUTPUT_FISSION] = {
  " (x,0nf)     "," (x,1nf)     "," (x,2nf)     "," (x,3nf)     "," (x,4nf)     "};

  std::ofstream fp;
  fp.open(fname.c_str(),std::ios::out | std::ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteCrossSectionFission");
  }

  double cx[N_OUTPUT_FISSION];
  
  for(int i=0 ; i<N_OUTPUT_FISSION ; i++) cx[i] = 0.0;
  for(int i=0 ; i<nc ; i++){
    if(i >= N_OUTPUT_FISSION) break;
    cx[i] = crx.prod[i].fiss;
  }

  /*** print header if file not exist */
  if(!fexist){
    fp << "#" << std::setw(COLUMN_WIDTH-1) << " Fission    ";
    for(int i=0 ; i<N_OUTPUT_FISSION ; i++) fp << std::setw(COLUMN_WIDTH) << rname[i];
    fp << std::endl;
  }

  /*** print cross sections */
  fp << std::setprecision(DECIMAL_WIDTH);
  fp << std::setiosflags(std::ios::scientific);
  fp << std::setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_FISSION ; i++){
    fp << std::setw(COLUMN_WIDTH) << ( (cx[i] < eps) ? 0.0 : cx[i] );
  }
  fp << std::endl;
  fp.close();
}


/**********************************************************/
/*     Output Level Excitation Cross Sections on File     */
/**********************************************************/
void fioWriteLevelExcite(const bool fexist, const std::string fname, double *y, const double e, const int id, Nucleus *n)
{
  std::ofstream fp;
  fp.open(fname.c_str(),std::ios::out | std::ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteInelastic");
  }

  double x[N_OUTPUT_LEVELS];

  for(int k=0 ; k<N_OUTPUT_LEVELS ; k++) x[k] = 0.0;
  
  /*** sum all level inelastic up to min(N_OUTPUT_LEVELS-1,ndisc) */
  int m = N_OUTPUT_LEVELS-1;
  if(n->ndisc < m) m = n->ndisc;

  double sum = 0.0;
  for(int k=0 ; k<m ; k++){
    sum += (x[k] = crx.levexcite[id][k]);
  }

  /*** the last element contains continuum inelastic cross section */
  x[N_OUTPUT_LEVELS-1] = n->lpop[0] - sum;
  y[N_OUTPUT_LEVELS-1] = y[m-1] + 0.001;

  /*** not so exact for the continuum, to avoid round off error    */
  if(x[N_OUTPUT_LEVELS-1] < 0.0) x[N_OUTPUT_LEVELS-1] = 0.0;

  /*** if the c/s to the highest level is zero, no continuum can be assumed */
  if(x[m-1] == 0.0) x[N_OUTPUT_LEVELS-1] = 0.0;

  fp << std::setprecision(DECIMAL_WIDTH);
  fp << std::setiosflags(std::ios::scientific);
  if(!fexist){
    fp << "#" << std::setw(COLUMN_WIDTH-1) << "            ";
    for(int i=0 ; i<N_OUTPUT_LEVELS ; i++) fp << std::setw(COLUMN_WIDTH) << y[i];
    fp << std::endl;
  }
  fp << std::setw(COLUMN_WIDTH) << e;
  for(int i=0 ; i<N_OUTPUT_LEVELS ; i++) fp << std::setw(COLUMN_WIDTH) << x[i];
  fp << std::endl;
  fp.close();
}


/**********************************************************/
/*     Output Scattering Angular Distributions on File    */
/**********************************************************/
void fioWriteAngularDistribution(std::string fname, double *y, const double e)
{
  std::ofstream fp;
  fp.open(fname.c_str(),std::ios::out | std::ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteAngularDistribution");
  }
  
  fp << std::setprecision(DECIMAL_WIDTH-1);
  fp << std::setiosflags(std::ios::scientific);
  fp << "#" << std::setw(COLUMN_WIDTH-1) << e << std::setw(COLUMN_WIDTH) << "   Angdist   ";

  fp << std::setprecision(DECIMAL_WIDTH);
  fp << std::setiosflags(std::ios::scientific);
  for(int k=1 ; k<N_OUTPUT_LEVELS-1 ; k++) fp << std::setw(COLUMN_WIDTH) << y[k];
  fp << std::endl;

  for(int i=0 ; i<MAX_ANGDIST ; i++){
    fp << std::setw(COLUMN_WIDTH) << crx.theta[i];
    for(int j=0 ; j<N_OUTPUT_LEVELS-1 ; j++) fp << std::setw(COLUMN_WIDTH) << crx.angdist[j][i];
    fp << std::endl;
  }

  fp.close();
}


/**********************************************************/
/*     Output Legendre Coefficients on File               */
/**********************************************************/
void fioWriteLegendreCoefficient(std::string fname, double *y, const double e, const int id)
{
  std::ofstream fp;
  fp.open(fname.c_str(),std::ios::out | std::ios::app);
  if(!fp){
    message << "data file " << fname << " output error";
    cohTerminateCode("fioWriteLegendreCoefficient");
  }

  fp << std::setprecision(DECIMAL_WIDTH-1);
  fp << std::setiosflags(std::ios::scientific);
  fp << "#" << std::setw(COLUMN_WIDTH-1) << e << std::setw(COLUMN_WIDTH) << "   Legcoef   ";

  fp << std::setprecision(DECIMAL_WIDTH);
  fp << std::setiosflags(std::ios::scientific);
  for(int k=1 ; k<N_OUTPUT_LEVELS-1 ; k++) fp << std::setw(COLUMN_WIDTH) << y[k];
  fp << std::endl;

  fp << std::setprecision(DECIMAL_WIDTH-1);
  for(int i=0 ; i<MAX_J ; i++){
    fp << std::setw(COLUMN_WIDTH) << i;
    for(int j=0 ; j<N_OUTPUT_LEVELS-1 ; j++) fp << std::setw(COLUMN_WIDTH) <<  crx.legcoef[id][j][i];
    fp << std::endl;
  }

  fp.close();
}

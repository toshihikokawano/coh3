// functions for result output in a table format, and default file names

#define filename_cross_section          "CoHCrossSection"
#define filename_particle_production    "CoHParticleProduction"
#define filename_radioactive_production "CoHRadioactiveProduction"
#define filename_fission                "CoHFission"
#define filename_level_excite           "CoHLevelExcite"
#define filename_angular_distribution   "CoHAngularDistribution"
#define filename_legendre_coefficient   "CoHLegendreCoefficient"
#define file_extension                  "dat"


class UnstableLevel{
 public:
  int            nx;        // level number for this meta-stable state
  double     energy;        // level excitation energy
  double   halflife;        // half-life of the state
  double production;        // production cross section

  UnstableLevel(){
    init();
  }

  void init(){
    nx = -1;
    energy = 0.0;
    halflife = 0.0;
    production = 0.0;
  }
};


class Radioisotope{
 private:
  int       maxlevel;        // max number of levels allocated
  bool     allocated;
 public:
  ZAnumber        za;        // Z and A number of this nuclide
  UnstableLevel *lev;        // radio-activity
  bool   radioactive;        // flag if this is radioactive

  Radioisotope(){
    maxlevel = 0;
    allocated = false;
    za = ZAnumber(0,0);
    radioactive = false;
  }

  ~Radioisotope(){
    memfree();
  }

  void memalloc(const int m){
    if(!allocated){
      maxlevel = m;    
      lev = new UnstableLevel [maxlevel];
      for(int i=0 ; i<maxlevel ; i++) lev[i].init();
      allocated = true;
    }
  }

  void memfree(){
    if(allocated){
      delete [] lev;
      allocated = false;
    }
  }

  int getmaxlevel(void){ return maxlevel; }
};




/**************************************/
/*      fileout.cpp                   */
/**************************************/
void    cohFileOut (int, int, double);

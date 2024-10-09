// class for parameter adjustment

const int MAX_ADJUSTABLES   =  100 ;  /* maximum num of adjustable parameters  */

enum parameterType{
     null,
     parmOV,      // real potential depth
     parmOW,      // imaginary potential depth
     parmORV,     // real potential radius
     parmORW,     // imaginary potential radius
     parmOAV,     // real potential diffuseness
     parmOAW,     // imaginary potential diffuseness
     parmOC,      // Coulomb radius
     parmM2,      // M2 (scattering) parameter in pre-equilibrium
     parmM2R,     // M2 (Exchange) parameter
     parmSD,      // single particle state density parameter
     parmKO,      // alpha-particle knock out
     parmTJ,      // transmission coefficients from a specific compound
     parmLD,      // level density parameters
     parmSPIN,    // spin-cutoff parameter in level density formula
     parmPAIR,    // pairing energy in level density formula
     parmPDST,    // parity distribution
     parmLCUT,    // number of discrete levels included
     parmDSDV,    // DSD coupling real strength
     parmFL1,     // first fission barrier level density enhancement factor
     parmFL2,     // second fission barrier level density enhancement factor
     parmFL3,     // third fission barrier level density enhancement factor
     parmGSTR,    // gamma-ray strength function
     parmGDRM1,   // M1 strength function
     parmISOM,    // threshold time for determining isomeric states
     parmD0SR,    // search level density parameter for given D0 value
     parmESET,    // built-in incident energy grid
     parmBROD,    // photon spectrum Gaussian broadening
     parmCOLL     // collective ehnancement in pre-equilibrium level density
};
typedef enum parameterType ParmType;

/*** parameter names for printing out */
static std::string pname[] = 
    {"    ","OV  ","OW  ","ORV ","ORW ","OAV ","OAW ","OC  ","M2  ","M2R ",
     "SD  ","KO  ","TJ  ","LD  ","SPIN","PAIR","PDST","LCUT","DSDV",
     "FL1 ","FL2 ","FL3 ","GSTR","GDRM","ISOM","D0SR","ESET","BROD","COLL"};

/****************************/
/*   Adjustable Parameters  */
/****************************/
class Parameter{
 private:
  double       factor;
  double       width;
 public:
  ParmType     type;
  Particle     particle;
  ZAnumber     za;
  Parameter(){
    init();
  }
  void init(){
    type     = null;
    za.setZA(0,0);
    particle = unknown;
    factor   = 1.0;
    width    = 0.0;
  }
  void setfactor(double f){
    factor = f;
  }
  void setwidth(double w){
    width  = w;
  }
  double getfactor(){
    return(factor);
  }
  double getwidth(){
    return(width);
  }
};


/****************************/
/*   Adjustable Object      */
/****************************/
class Adjustable{
 private:
  int          np;             /* total number of given parameters */
  int          npmax;          /* max  number of given parameters  */
  bool         arrayalloc;     /* flag for memory allocation       */
 public:
  Parameter    *parm;          /* array of individual parameter    */

  Adjustable(){
    np = 0;
    arrayalloc = false;
  }

  void memalloc(int n){
    if(!arrayalloc){
      parm = new Parameter [n];
      npmax = n;
      arrayalloc = true;
      for(int i=0 ; i<n ; i++) parm[i].init();
    }
  }

  ~Adjustable(){
    if(arrayalloc){
      delete [] parm;
      np = 0;
      arrayalloc = false;
    }
  }

  int nparm(){ return(np); }


  bool addparm(ParmType t, double f, double w){
    int k = getindex(t);
    if(k >= 0){
      parm[k].setfactor(f);
      parm[k].setwidth(w);
    }
    else{
      parm[np].type = t;
      parm[np].setfactor(f);
      parm[np].setwidth(w);
      np++;
      if(np == npmax) return false;
    }
    return true;
  }

  bool addparm(ParmType t, Particle p, double f, double w){
    int k = getindex(t,p);
    if(k >= 0){
      parm[k].setfactor(f);
      parm[k].setwidth(w);
    }
    else{
      parm[np].type = t;
      parm[np].particle = p;
      parm[np].setfactor(f);
      parm[np].setwidth(w);
      np++;
      if(np == npmax) return false;
    }
    return true;
  }

  bool addparm(ParmType t, ZAnumber za, double f, double w){
    int k = getindex(t,za);
    if(k >= 0){
      parm[k].setfactor(f);
      parm[k].setwidth(w);
    }
    else{
      parm[np].type = t;
      parm[np].za = za;
      parm[np].setfactor(f);
      parm[np].setwidth(w);
      np++;
      if(np == npmax) return false;
    }
    return true;
  }

  bool addparm(ParmType t, ZAnumber za, Particle p, double f, double w){
    int k = getindex(t,za,p);
    if(k >= 0){
      parm[k].setfactor(f);
      parm[k].setwidth(w);
    }
    else{
      parm[np].type = t;
      parm[np].za = za;
      parm[np].particle = p;
      parm[np].setfactor(f);
      parm[np].setwidth(w);
      np++;
      if(np == npmax) return false;
    }
    return true;
  }

  int getindex(ParmType t){
    int k = -1;
    for(int i=0 ; i < np ; i++){
      if(parm[i].type == t){ k = i; break; }
    }
    return(k);
  }

  int getindex(ParmType t, Particle p){
    int k = -1;
    for(int i=0 ; i < np ; i++){
      if(parm[i].type == t && parm[i].particle == p){ k = i; break; }
    }
    return(k);
  }

  int getindex(ParmType t, ZAnumber za){
    int k = -1;
    for(int i=0 ; i < np ; i++){
      if( (parm[i].type == t) && (parm[i].za == za) ){ k = i; break; }
    }
    return(k);
  }

  int getindex(ParmType t, ZAnumber za, Particle p){
    int k = -1;
    for(int i=0 ; i < np ; i++){
      if( (parm[i].type == t) && (parm[i].za == za) && (parm[i].particle == p) ){ k = i; break; }
    }
    return(k);
  }
};


#ifndef COH_TOPLEVEL
extern Adjustable adj;
#endif


/****************************/
/*   parameter.cpp          */
/****************************/

void   parmCheckFactor (void);
int    parmCheckParameter  (const ParmType);
int    parmCheckParameter  (const ParmType, const Particle);
int    parmCheckParameter  (const ParmType, const ZAnumber);
int    parmCheckParameter  (const ParmType, const ZAnumber, const Particle);
double parmGetFactor       (const ParmType);
double parmGetFactor       (const ParmType, const Particle);
double parmGetFactor       (const ParmType, const ZAnumber);
double parmGetFactor       (const ParmType, const ZAnumber, const Particle);
double parmGetValue        (const ParmType);
int    parmGetValue        (const ParmType, const ZAnumber);

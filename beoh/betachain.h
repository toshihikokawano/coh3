/****************************/
/*   Defaults               */
/****************************/
/*** constants for beta-decay chain calculation */
static const int MAX_BETA_BRANCH = 1000; // maximum beta branches
static const int MAX_DECAYMODE   =    4; // maximum number of decay mode


/****************************/
/*   Beta-Decay Chain Calc. */
/****************************/
inline unsigned int tomza(int m, int z, int a)
{ return (  ((char)m << 16) | ((char)z << 8) | (unsigned int)a  ); }
/*
 MZA alignment
    Z        | A            | M
    01234567 | 0123456789ab | 0123
 */

class DecayMode{
 private:
  int          ndm;                    // number of decay modes
  double       branch[MAX_DECAYMODE];  // branching ratios
  int          next[MAX_DECAYMODE];    // nucide ID after decay
 public:

  DecayMode(){
    ndm = 0;
    for(int i=0 ; i<MAX_DECAYMODE ; i++){
      branch[i]  = 0.0;
      next[i]    = 0;
    }
  }

  void setNDM(int i){
    if(i<MAX_DECAYMODE){
      ndm = i;
    }
  }

  void setData(int i, double b, int k){
    if(i < ndm){
      branch[i]  = b;
      next[i]    = k;
    }
  }

  int getNDM(void){
    return ndm;
  }

  double getBranch(int i){
    if(i < ndm){
      return branch[i];
    }
    else return 0.0;
  }

  int getNext(int i){
    if(i < ndm){
      return next[i];
    }
    else return 0.0;
  }
};


class Isotope{
 private:
  unsigned int MZA;    // ID of this nuclide
  unsigned int M;      // meta-state
  unsigned int Z;      // atomic number
  unsigned int A;      // mass number
  double       spin;   // spin
  int          parity; // parity (+/-) 1
  bool         stable;

 public:
  double    lambda;
  double    initial;
  double    yield;
  double    std;
  DecayMode mode;

  Isotope(){
    init();
  }

  Isotope(int m, int z, int a){
    Z   = z;
    A   = a;
    M   = m;
    MZA = tomza(m,z,a);
    stable   = true;
    yield    = 0.0;
    initial  = 0.0;
    std      = 0.0;
    lambda   = 0.0;
  }

  void init(){
    MZA = 0;
    Z   = 0;
    A   = 0;
    M   = 0;
    spin   = 0.0;
    parity = 0;
    stable   = true;
    yield    = 0.0;
    initial  = 0.0;
    std      = 0.0;
    lambda   = 0.0;
  }

  unsigned int setZA(int m, int z, int a){
    M = m;
    Z = z;
    A = a;
    return(tomza(m,z,a));
  }

  void setLambda(double t){
    if(t > 0.0){
      stable   = false;
      lambda   = t;
    }
    else{
      stable   = true;
      lambda   = 0.0;
    }
  }

  void setJP(double j2, int p){
    spin   = j2;
    parity = (p >= 0) ? 1 : -1;
  }

  unsigned int getZ()  { return Z; }
  unsigned int getA()  { return A; }
  unsigned int getM()  { return M; }
  unsigned int getMZA(){ return(tomza(M,Z,A)); }
  double getJ() { return spin; }
  int getP() { return parity; }
  bool isstable(){ if(stable) return true; else return false; }
  double getLambda()   { return(lambda); }
  double activity()    { return(yield * lambda); }
};


/****************************/
/*   Beta Decay Branching   */
/****************************/
class Branch{
 private:
  double energy;
  double fraction;
 public:
  Branch(){
    energy   = 0.0;
    fraction = 0.0;
  }
  double getE(void){ return(energy); }
  double getR(void){ return(fraction); }
  void setE(double e){ energy   = e; }
  void setR(double b){ fraction = b; }
  void scaleR(double x){ fraction *= x; }
  void setVal(double e, double b){ energy = e; fraction = b; }
};


/****************************/
/*   Beta Decay Parameter   */
/****************************/
class Beta{
 public:
    int     nstate;           // nubmer of final states
    Branch  *br;              // branching ratio
    double  width;            // Gaussian broadening width in MeV

    Beta(){
      nstate   = 0;
      br       = new Branch [MAX_BETA_BRANCH];
      width    = 0.0;
    }

    ~Beta(){
      delete [] br;
    }
};



const int MAX_CASCADE    =    100 ;  /* maximum number of gamma-rays per cascade */
const int MAX_SIMULATION =  10000 ;  /* number of Monte Carlo Simulations */

/****************************/
/*    Monte Carlo Index     */
/****************************/
class MCPoint{
 private:
  int    n;    // nucleus index
  int    p;    // CN state parity
  int    j;    // J-array index
  int    k;    // excitation index
  double e;    // excitation energy
  bool   l;    // true if discrete level
 public:
  MCPoint(){
    n  = 0;
    p  = 0;
    j  = 0;
    k  = 0;
    e  = 0.0;
    l  = false;
  }

  MCPoint(const int n0, const int p0, const int j0, const int k0, const double e0){
    set(n0,p0,j0,k0,e0);
  }

  void set(const int n0, const int p0, const int j0, const int k0, const double e0){
    n = n0;
    p = p0;
    j = j0;
    k = k0;
    e = e0;
  }

  void setE(const double e0){ e = e0; }

  void setDiscrete(){ l = true; }
  void setContinuum(){ l = false; }

  int    getNucleus (){return n;};
  int    getEk      (){return k;};
  int    getParity  (){return p;};
  int    getSpin    (){return j;};
  double getEnergy  (){return e;};
  bool   isDiscrete (){return l;};
};


/****************************/
/*    Particle History      */
/****************************/
class MCHistory{
 private:
  int     id;  // particle identifier
  double  ep;  // emission energy
  MCPoint mp;  // MC index for intermediate state
 public:
  MCHistory(){
    id = 0;
    ep = 0.0;
    mp.set(0,0,0,0,0.0);
  }

  MCHistory(const int id0, const double ep0, MCPoint mp0){
    id = id0;
    ep = ep0;
    mp = mp0;
  }

  void set(const int id0, const double ep0, MCPoint mp0){
    id = id0;
    ep = ep0;
    mp = mp0;
  }

  int     getIndex()  {return id;};
  double  getEnergy() {return ep;};
  MCPoint getMCPoint(){return mp;};
};


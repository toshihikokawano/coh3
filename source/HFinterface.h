const int HF_MAX_EXPANSION =  40;
const int HF_MAX_SPSTATE   = 100;


/**********************************************************/
/*      Single Particle Energies and Spherical Expansion  */
/**********************************************************/
class SphericalExpansion{
 public:
  int n;
  int l;
  int j;
  double c;

  SphericalExpansion(){
    n = 0 ; l = 0 ; j=0 ;
    c = 0.0;
  }
  void setQval(int k1, int k2, int k3){
    n = k1; l = k2;  j = k3;
  }
  void setCval(double x){
    c = x;
  }
};


class HFState{
 public:
  SphericalExpansion *q;
  double energy;
  double v2;
  int    parity;
  int    k2;
  int    m;
  int    lmax;

  HFState(){
    energy = 0.0;
    v2     = 0.0;
    parity = 0;
    k2     = 0;
    m      = 0;
    lmax   = 0;
  }

  void memalloc(){
    q = new SphericalExpansion [HF_MAX_EXPANSION];
  }

  void memfree(){
    delete [] q;
  }

  void set(int n1, int n2, int n3, double c){
    if(m < HF_MAX_EXPANSION){
      q[m].setQval(n1,n2,n3);
      q[m].setCval(c);
      m++;
    }
  }
};


class HFInterface{
 private:
  bool allocated;
 public:
  int      n;
  double   b;
  double   lambda;
  double   delta;
  HFState *state;

  HFInterface(){
    n      = 0;
    b      = 0.0;
    lambda = 0.0;
    delta  = 0.0;
    allocated = false;
  }

  ~HFInterface(){
    memfree();
  }

  void memalloc(){
    if(!allocated){
      state = new HFState [HF_MAX_SPSTATE];
      for(int i=0 ; i<HF_MAX_SPSTATE ; i++){
        state[i].memalloc();
      }
      allocated = true;
    }
  }
  void memfree(){
    if(allocated){
      for(int i=0 ; i<HF_MAX_SPSTATE ; i++){
        state[i].memfree();
      }
      delete [] state;
      allocated = false;
    }
  }

};


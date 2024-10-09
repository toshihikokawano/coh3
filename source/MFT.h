#include <cmath>

/**********************************************************/
/*      Nucleus Parameters                                */
/**********************************************************/
class MFTSystem{
 private:
  double eps[6];
  double gamma;
  int a;
  int z;
  int n;

 public:
  int      max_iteration; // max number of HF iteration
  double   eps_converge;  // HF iteration convergence criterion

  double   R0;    // radius
  double   S0;    // surface area
  double   V0;    // volume
  double   S1;    // deformed surface area
  double   a13;   // A^{1/3}
  double   a23;   // A^{2/3}
  double   asym;  // (Z-N)/A
  double   bs;    // Bs parameter in FRDM macro
  bool     pol;   // polar coordinate

  MFTSystem(){
    for(int i=0 ; i<6 ; i++) eps[i] = 0.0;
    gamma = 0.0;
  }

  MFTSystem(int z0, int a0){
    for(int i=0 ; i<6 ; i++) eps[i] = 0.0;
    gamma = 0.0;

    init(z0,a0);
  }

  void init(int z0, int a0){
    z    = z0;
    a    = a0;
    n    = a - z;

    max_iteration = 200;
    eps_converge  = 1.0e-05;

    asym = (n - z)/(double)a;
    a13  = pow((double)a,1.0/3.0);
    a23  = a13*a13;
    R0   = 0.0;
    S0   = 0.0;
    S1   = 0.0;
    V0   = 0.0;
    bs   = 1.0;
    pol  = true;
  }

  void setIteration(int n, double e){
    max_iteration = n;
    eps_converge  = e;
  }

  void setNuclearSize(double r0){
    R0   = r0 * a13;
    S0   = 4.0 * PI *R0*R0;
    V0   = 4.0/3.0 * PI * R0*R0*R0;
  }

  void setEps(int i, double x){
    if(1 <= i && i <= 6) eps[i-1] = x;
  }

  void setGamma(double x){
    gamma = x;
  }

  double getEps(int i){
    if(i <= 0 || i >= 7) return(0.0);
    return(eps[i-1]);
  }

  double getGamma(){
    return(gamma);
  }

  int    getZ() { return(z);  }
  int    getA() { return(a);  }
  int    getN() { return(a-z);  }
};


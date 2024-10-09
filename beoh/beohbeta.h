
static const int MAX_BETA_BRANCH = 1000; // maximum beta branches

/****************************/
/*   Beta Decay Branching   */
/****************************/
class Branch{
 private:
  double energy;
  double ratio;
 public:
  Branch(){
    energy = 0.0;
    ratio  = 0.0;
  }
  double getE(void){
    return(energy);
  }
  double getR(void){
    return(ratio);
  }
  void setE(double e){
    energy = e;
  }
  void setR(double b){
    ratio  = b;
  }
  void scaleR(double x){
    ratio  *= x;
  }
  void setVal(double e, double b){
    energy = e;
    ratio  = b;
  }
};


/****************************/
/*   Beta Decay Parameter   */
/****************************/
class Beta{
 public:
    int     nstate           ;     /* nubmer of final states            */
    Branch  *br              ;     /* branching ratio                   */

    Beta(){
      nstate   = 0;
      br       = new Branch [MAX_BETA_BRANCH];
    }
    ~Beta(){
      delete [] br;
    }
};

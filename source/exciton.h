#define MAX_PREEQ  6

/**************************************/
/*      Exciton Configuration         */
/**************************************/
class Exconf{
 public:
    int zp;     // proton particle number
    int zh;     // proton hole number
    int np;     // neutron particle number
    int nh;     // neutron hole number

    Exconf(){ zp = zh = np = nh = 0; }
    Exconf(int p1, int h1, int p2, int h2){ set(p1,h1,p2,h2); }

    void set(int p1, int h1, int p2, int h2){
      zp = p1;
      zh = h1;
      np = p2;
      nh = h2;
    }

    int sum(){ return zp + zh + np + nh; }
};


/**************************************/
/*      Singple Particle Density      */
/**************************************/
class SPdens{
 public:
    double     gz;            // state density parameter for proton
    double     gn;            // state density parameter for neutron
    double     pairing;       // pairing energy
    double     well_depth;    // V, finite potential correction

    SPdens(){
      gz           = 0.0;
      gn           = 0.0;
      pairing      = 0.0;
      well_depth   = 0.0;
    }
};


/**************************************/
/*      Preequilibrium Systyem Data   */
/**************************************/
class Preeq{
 public:
    double     ex_total;      // maximum excitation energy
    double     omega_total;   // total state density for p-h
    double     m2zz;          // averaged matrix elements for p-p
    double     m2nn;          // averaged matrix elements for n-n
    double     m2zn;          // averaged matrix elements for p-n
    double     m2nz;          // averaged matrix elements for n-p
    SPdens     *spd;          // pointer to s-p density array

    Preeq(){
      ex_total    = 0.0;
      omega_total = 0.0;
      m2zz        = 0.0;
      m2nn        = 0.0;
      m2zn        = 0.0;
      m2nz        = 0.0;
    }
};


/**************************************/
/*      exciparm.cpp                  */
/**************************************/

double  preqSingleStateDensity (const int, const int);
double  preqPairingDelta (const int, const int);
double  preqPairingCorrection (const int, const double, SPdens *);
double  preqPauliCorrection (const int, const int, const int, const int, const double, const double);
double  preqMSquare (const int, const int, Preeq *);
double  preqFiniteWell (const int, const int, const double, double);
double  preqPotentialDepth (const double, const int, const int);
double  preqStateDensity (const double, SPdens *, Exconf *);
double  preqCollectiveEnhancement (const double, Exconf *);


/**************************************/
/*      excilambda.cpp                */
/**************************************/

double  preqLambdaPlusZ (Preeq *, Exconf *);
double  preqLambdaPlusN (Preeq *, Exconf *);
double  preqLambdaZeroZ (Preeq *, Exconf *);
double  preqLambdaZeroN (Preeq *, Exconf *);


/**************************************/
/*      excipop.cpp                   */
/**************************************/

void    preqPOccupation (double **,double **,double **, double **,double **,double **,double **);

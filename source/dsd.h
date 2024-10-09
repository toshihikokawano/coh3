// functions for DSD model, and constants

#define MAX_POINTS_BW   400   // maximum radial points for bound wave
#define MAX_SPLEVEL      50   // maximum number of sp states
#define VIB_FORMFACTOR    2   // index of coupling formfactor


/****************************/
/*   Bound State Data       */
/****************************/
struct boundstate{
  public:
    int     n;                // bound state principle q.num. n
    int     l;                // bound state q.num. l
    int     j2;               // bound state q.num. j*2
    int     k2;               // band K*2 (or Omega)
    double  bind;             // bound state binding energy
    double  w;                // BCS theory U factor or weight
    double  *wave;            // bound wave func.
};
typedef struct boundstate Bstate;

struct phstate{
  public:
    Bstate  Z;                // proton shell
    Bstate  N;                // neutron shell
};
typedef struct phstate PHstate;


/**************************************/
/*      BCS.CPP                       */
/**************************************/
double  bcsEnergyGap (const Particle, const int, const int);
double  bcsVfactor (const double, const double, const double);
double  bcsGapEquation (const int, const double, const double, Bstate *);


/**************************************/
/*      NILSSON.CPP                   */
/**************************************/
int     nilssonFillOrbit (const Particle, const int);
int     nilssonLevel (const Particle, const int, const int, const unsigned int, Bstate *, double, double, double);


/**************************************/
/*      BSTATWAV.CPP                  */
/**************************************/
int     bstateWavefunc (Pdata *, ZAnumber *, const double, const double, const double, const double, const double, Optical *, Potential *, Bstate *);


/**************************************/
/*      BSTATSOLVE.CPP                */
/**************************************/
double  bstateFindDepth (const int, const int, const int, const double, const double, const double, double *, Potential *);


/**************************************/
/*      DSDOVERLAP.C                  */
/**************************************/
void    dsdFormfactor (const int, Optical *, Potential *, const double, const double, std::complex<double> *, std::complex<double> *);
std::complex<double> dsdRadialIntegral (const int, const double, std::complex<double> *, double *, std::complex<double> *);

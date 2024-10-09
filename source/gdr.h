/**************************************/
/*      Gamma-ray Multipolarity       */
/**************************************/
typedef enum {EX = 0, SL = 1, GL = 2, ML = 3} GammaProfile;
typedef enum {E1 = 0, M1 = 1, E2 = 2, M2 = 3, E3 = 4, M3 = 5, E4 = 6, M4 = 7} GammaMultipolarity;


/****************************/
/*      G(D)R Parameters    */
/****************************/
class GDR{
 private:
    GammaProfile prof;   // Lorentzian shape
    std::string  XL;     // E or M, and multipolarity
    double       energy; // GDR energy
    double       width;  // GDR width
    double       sigma;  // GDR peak-cross section
 public:
    GDR(){
      clear();
    }

    void setGDR(std::string em, double e, double w, double s, GammaProfile p){
      prof    = p;
      XL      = em;
      energy  = e;
      width   = w;
      sigma   = s;
    }

    void clear(){
      prof    = SL;
      XL      = "  ";
      energy  = 0.0;
      width   = 0.0;
      sigma   = 0.0;
    }

    GammaProfile getProfile() {return prof;}
    std::string getXL () {return XL;}
    double getEnergy  () {return energy;}
    double getWidth   () {return width ;}
    double getSigma   () {return sigma ;}
    char   getEM      () {char em = (XL[0]=='E') ? 'E' : 'M'; return(em);}
    int    getL       () {return( (int)(XL[1])-'0');}

    void   setProfile(GammaProfile p) {prof = p;}
    void   setXL     (std::string em) {XL     = em;}
    void   setEnergy (double e ) {energy = e;}
    void   setWidth  (double w ) {width  = w;}
    void   setSigma  (double s ) {sigma  = s;}
};


/**************************************/
/*      gdr.cpp                       */
/**************************************/
void    gdrE1 (const double, const double, GDR *);
void    gdrE1DoubleHump0 (const double, const double, GDR *);
void    gdrE1DoubleHump1 (const double, const double, GDR *);
void    gdrM1 (const double, GDR *);
void    gdrM1scissors (const double, const double, GDR *);
void    gdrE2 (const double, const double, GDR *);
void    gdrM2 (GDR *, GDR *);
void    gdrE3 (GDR *, GDR *);
void    gdrM3 (GDR *, GDR *);
void    gdrE4 (GDR *, GDR *);

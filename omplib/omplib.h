#include "optical.h"
#include "masstable.h"

/*                                 Energy  At   Zt   Ai   Zi           */
unsigned int OMPWilmoreHodgson(            int,      int, int, Optical *);
unsigned int OMPBecchettiGreenlees(        int, int, int, int, Optical *);
unsigned int OMPRapaport(          double, int, int, int, int, Optical *);
unsigned int OMPWalterGuss(        double, int, int, int, int, Optical *);
unsigned int OMPChapelHill89(      double, int, int, int, int, Optical *);
unsigned int OMPKoningDelaroche(   double, int, int, int, int, Optical *);
unsigned int OMPSoukhovitskii(     double, int, int, int,      Optical *);
unsigned int OMPKunieda(           double, int, int, int, int, Optical *);
unsigned int OMPPerey(                     int, int, int, int, Optical *);
unsigned int OMPMadlandSchwandt(   double, int, int, int, int, Optical *);
unsigned int OMPMadlandYoung(      double, int, int,           Optical *);

unsigned int OMPLemos(                               int, int, Optical *);
unsigned int OMPNolte(                     int, int, int, int, Optical *);
unsigned int OMPAvrigeanu(         double, int, int, int, int, Optical *);
unsigned int OMPAvrigeanu2009(     double, int, int, int, int, Optical *);
unsigned int OMPAvrigeanu2014(     double, int, int, int, int, Optical *);
unsigned int OMPTALYS_alpha(       double, int, int, int, int, Optical *);
unsigned int OMPBojowald(          double, int, int, int, int, Optical *);
unsigned int OMPAnHaixia(                  int, int, int, int, Optical *);
unsigned int OMPHanYinlu(                  int, int, int, int, Optical *);
unsigned int OMPBoundState(                               int, Optical *);
unsigned int OMPtest(                                     int, Optical *);

unsigned int OMPspline(            double, int, int, int, int, Optical *);


/*                                 Energy  At   Zt   Ai   Zi           */
unsigned int OMPscratch(           double, int, int, int, int, Optical *);
unsigned int OMPuserdef(           double, int, int, int, int, Optical *);
unsigned int OMPSmithA120(         double, int, int,           Optical *);
unsigned int OMPSmithA136(         double, int, int,           Optical *);
unsigned int OMPSmithA415(                 int, int,           Optical *);
unsigned int OMPFlap22(            double, int, int,           Optical *);
unsigned int OMPYoung_Am(          double, int, int,           Optical *);
unsigned int OMPYoung_Pu(          double, int, int,           Optical *);
unsigned int OMPYoung_Re(          double, int, int,           Optical *);
unsigned int OMPModSoukhovitskii(  double, int, int, int,      Optical *);
unsigned int OMPSoukhovitskii2005( double, int, int, int,      Optical *);
unsigned int OMPDave1p(                    int, int,           Optical *);
unsigned int OMPWLH1(              double, int, int, int, int, Optical *);
unsigned int OMPWLH2(              double, int, int, int, int, Optical *);
unsigned int OMPWLH3(              double, int, int, int, int, Optical *);
unsigned int OMPWLH4(              double, int, int, int, int, Optical *);
unsigned int OMPWLH5(              double, int, int, int, int, Optical *);



double OMPdeltaSurface (const double, const double, const double, const double, const double);
double OMPdeltaVolume (const double, const double, const double, const double);
double OMPasymmetricVolume (const double, const double, const double, const double, const double);

static inline double OMPFermiEnergy(const int at, const int zt, const int zi)
{
  double  sn = 0.0, bn = 0.0;
  bool proton = (zi == 1) ? true : false;

  if(proton){
    sn = mass_excess(zt-1,at-1) + EPROTON  - mass_excess(zt  ,at  );
    bn = mass_excess(zt  ,at)   + EPROTON  - mass_excess(zt+1,at+1);
  }
  else{
    sn = mass_excess(zt  ,at-1) + ENEUTRON - mass_excess(zt  ,at  );
    bn = mass_excess(zt  ,at)   + ENEUTRON - mass_excess(zt  ,at+1);
  }

  return(-(sn+bn)*0.5);
}

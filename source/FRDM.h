#ifndef __PHYSICALCONSTANT_H__
#define __PHYSICALCONSTANT_H__
#include "physicalconstant.h"
#endif

#ifndef __HFINTERFACE_H__
#define __HFINTERFACE_H__
#include "HFinterface.h"
#endif

#ifndef __HOBASIS_H__
#define __HOBASIS_H__
#include "HObasis.h"
#endif

#ifndef __DIR_H__
#define __DIR_H_
#include "dir.h"
#endif

#ifndef __MFT_H__
#define __MFT_H__
#include "MFT.h"
#endif

#ifndef __FRDMCOORDINATE_H__
#define __FRDMCOORDINATE_H__
#include "FRDMcoordinate.h"
#endif


/**********************************************************/
/*      Fixed File Names for Outputs                      */
/**********************************************************/
static const std::string fileFRDMShape     = "CoHFRDMShape.dat";
static const std::string fileFRDMBeta      = "CoHFRDMBeta.dat";
static const std::string fileFRDMParameter = "CoHFRDMParameter.dat";
static const std::string fileFRDMEnergy    = "CoHFRDMEnergy.dat";

static const std::string fileFRLDMShape    = "CoHFRLDMShape.dat";
static const std::string fileFRLDMEnergy   = "CoHFRLDMEnergy.dat";


/**********************************************************/
/*      Macro Part Parameters                             */
/**********************************************************/

const double gCurvature        =   41.0;  // Harmonic oscillator curvature parameter
const double gRadius           =   1.16;
const double gCharge2          =   1.4399764; // square of electronic charge [MeV fm]
const double gYukawaRange      =   0.68;
const double gYukawaCharge     =   0.70;
const double gSurfaceStiffness =  29.21;  // Q: effective surface-stiffness constant
const double gSymmetryEnergy   =  32.73;  // J: symmetry energy constant
const double gCompExp          =  60.0;   // C: pre-exponential compressibility-term constant
const double gCompExpRange     =   0.831; // gamma: compressibility range constant


const double gVolumeEnergy     =  16.247; // a1: volume energy constant
const double gSurfaceEnergy    =  22.92;  // a2: surface energy constant
const double gCurvatureEnergy  =   0.0;   // a3: curvature energy constant
const double gDensitySymmetry  =   0.0;   // L: density-symmetry constant
const double gCompressibility  = 240.0;   // K: nuclear compressibility constant

const double gProtonRadius     =   0.80;  // rp: proton root-mean-square radius

const double gChargeAsymmetry  =   0.436; // ca: charge asymmetry constant
const double gWigner           =  30.0;   // W: Wigner constant
const double gPairingGap       =   4.80;  // rmac: average pairing gap constant
const double gNPInteraction    =   6.6;   // h: neutron proton interaction constant
const double gEBinding         =   1.433e-5; // electronic-binding constnat


/**********************************************************/
/*      Micro Part Parameters                             */
/**********************************************************/

const double gYukawaRangeMicro =   0.80;
const double gAdensity         =   0.82;  // constant of Fermi potential radius
const double gBdensity         =   0.56;  // constant of Fermi potential radius
const double gVdepth           =   52.5;  // central potential depth
const double gVasymmetry       =   48.7;  // asymmetry term for potential depth
const double gSpinOrbitAn      = 0.01875; // spin-orbit coeff. kn
const double gSpinOrbitBn      =   31.5;  // spin-orbit coeff. ln
const double gSpinOrbitAp      = 0.02500; // spin-orbit coeff. kp
const double gSpinOrbitBp      =   28.0;  // spin-orbit coeff. lp
const double gPairingGapMicro  =   3.20;  // rmic: average pairing gap constant

const double gSurfaceEnergyMicro    =  22.0; // a2: surface energy constant
const double gDensitySymmetryMicro  =  99.0; // L: density-symmetry constant
const double gCompressibilityMicro  = 300.0; // K: nuclear compressibility constant


const double gSurfaceStiffnessMicro =  25.0; // Q: effective surface-stiffness constant
const double gSymmetryEnergyMicro   =  35.0; // J: symmetry energy constant


/**********************************************************/
/*      FRLDM Macro Part Parameters                       */
/**********************************************************/

const double fVolumeEnergy     =  16.00126; // av: volume energy constant
const double fVolumeAsymmetry  =   1.92240; // kv: volume asymmetry constant
const double fSurfaceEnergy    =  21.18466; // as: surface energy constant
const double fSurfaceAsymmetry =   2.345  ; // ks: surface asymmetry constant
const double fAZeroEnergy      =   2.615  ; // a0: A0 constant
const double fChargeAsymmetry  =   0.10289; // ca: charge asymmetry constant



/**********************************************************/
/*      Lipkin-Nogami BCS Parameters                      */
/**********************************************************/
class LNParameter{
 public:
  int    n1;   // lower orbit to be included
  int    n2;   // upper orbit
  int    nf;   // index for Fermi level
  int    np;   // total particle
  double strength;           // G
  double pairing_gap;        // Delta
  double pairing_eff;        // Delta_G
  double pairing_mac;        // Delta-bar in macroscopic model
  double chemical_potential; // lambda
  double lambda2;            // lambda2
  double cutoff;             // highest s.p. energy
  double energy_micro;       // microscopic energy
  double energy_average;     // average energy

  LNParameter(){
    n1 = 0;
    n2 = 0;
    nf = 0;
    np = 0;
    strength = 0.0;
    pairing_gap = 0.0;
    pairing_eff = 0.0;
    pairing_mac = 0.0;
    chemical_potential = 0.0;
    lambda2 = 0.0;
    cutoff = 0.0;
    energy_micro = 0.0;
    energy_average = 0.0;
  }
};


/**********************************************************/
/*      Strutinsky Shell Correction Parameters            */
/**********************************************************/
class SCParameter{
 public:
  double chemical_potential; // lambda
  double rho;                // average level density at lambda
  double energy_micro;       // microscopic energy
  double energy_average;     // average energy

  SCParameter(){
    chemical_potential = 0.0;
    rho = 0.0;
    energy_micro = 0.0;
    energy_average = 0.0;
  }
};


/**********************************************************/
/*      Energies                                          */
/**********************************************************/
class FRDMEnergy{
 public:
  double macro;
  double micro;
  double spherical;
  double zero;
  double shellcorr;
  double pairing;
  double total;
  double liquiddrop;

  FRDMEnergy(){
    macro      = 0.0;
    micro      = 0.0;
    spherical  = 0.0;
    zero       = 0.0;
    shellcorr  = 0.0;
    pairing    = 0.0;
    total      = 0.0;
    liquiddrop = 0.0;
  }

  void setTotal(){
    total = (liquiddrop != 0.0) ? liquiddrop + micro : macro + micro + zero;
  }
};


/**********************************************************/
/*      Function Definitions                              */
/**********************************************************/
static double inline FuncYukawa(double a, double d)
{
  double f = 0.0;
  if(a == 0.0 || d == 0.0) return(f);

  double x = d/a;
  f = 2.0 - (x*x + 2.0*x + 2.0) * exp(-x);
  return(f);
}


static double inline FuncDiffuseCoulomb(double a, double d)
{
  double f = 0.0;
  if(a == 0.0 || d == 0.0) return(f);

  double x = d/a;
  f = 2.0*x - 5.0 + (0.5*x*x + 3.0*x + 5.0) * exp(-x);
  return(f);
}


static double inline FuncDiffusePotential(double a, double d)
{
  double f = 0.0;
  if(a == 0.0 || d == 0.0) return(f);

  double x = d/a;
  f = 1.0 - (1.0 + x) * exp(-x);
  return(f);
}


static inline double derivative(double *x, double d){
  double df = (x[0] - 8.0*x[1] + 8.0*x[2] - x[3])/(12.0*d);
  return(df);
}

static inline double distance(double *p, double *q){
  double d = 0.0;
  for(int i=0 ; i<3 ; i++) d += (p[i] - q[i]) * (p[i] - q[i]);
  return(d);
}

static inline void midpoint(double *p, double *q, double *m){
  for(int i=0 ; i<3 ; i++) m[i] = (p[i] + q[i])*0.5;
}

static inline double VectorProduct(double *a, double *b){
  double x = 0.0;
  for(int i=0 ; i<3 ; i++) x += a[i]*b[i];
  return(x);
}

static inline double VectorSubtract(double *a, double *b, double *c){
  double x = 0.0;
  for(int i=0 ; i<3 ; i++){
    c[i] = a[i] - b[i];
    x += c[i]*c[i];
  }
  return(x);
}



/*** FRDMmain.cpp */
void   FRDMmain (MFTSystem *, Basis *, double *, HFInterface *);

/*** FRDMshape.cpp */
double FRDMShape (const double, const double, MFTSystem *);
void   FRDMShapeNormalization (PolarGLPoint *, MFTSystem *);
void   FRDMShapeCreateSurface (PolarGLPoint *, MFTSystem *);
void   FRDMShapeStoreCurvature (PolarGLPoint *, MFTSystem *);

/*** betaexpand.cpp */
void   BetaExpansion (const int, double **, PolarGLPoint *, MFTSystem *);

/*** FRDMmacroenergy.cpp */
double FRDMMacroEnergy (PolarGLPoint *, MFTSystem *);
double FRDMMacroEnergySpherical (MFTSystem *);

/*** FRDMmacrointegral.cpp */
void   FRDMMacroDoubleVolumeIntegral (PolarGLPoint *, double *, MFTSystem *);
void   FRDMMacroVolumeIntegralSpherical (double *, MFTSystem *);

/*** FRDMmacroshapedependent.cpp */
void   FRDMMacroPrecalc (double **, PolarGLPoint *);
double FRDMMacroW (double *, PolarGLPoint *);
double FRDMMacroWbar(PolarGLPoint *, MFTSystem *);
double FRDMMacroSurfaceEnergy (MFTSystem *);
double FRDMMacroNeutronSkinEnergy (const double, MFTSystem *);
double FRDMMacroSurfaceRedistributionEnergy (const double, MFTSystem *);
double FRDMMacroCurvatureEnergy (MFTSystem *);
double FRDMMacroVolumeRedistributionEnergy (const double, PolarGLPoint *, MFTSystem *);

/*** FRDMmicroenergy.cpp */
double FRDMMicroEnergy (MFTSystem *, Basis *, HFInterface *);

/*** FRDMpotential.cpp */
void FRDMPotential(MFTSystem *, Basis *, OneBodyPotential *, LaguerrePolynomial *, HermitePolynomial *);

/*** FRDMsolvebcs.cpp */
void FRDMSoleveBCS (MFTSystem *, double, LNParameter *, SPEnergy *);

/*** FRDMpairing.cpp */
double FRDMPairingEnergy (LNParameter *, SPEnergy *);
double FRDMAveragePairingEnergy (const double, LNParameter *);

/*** FRDMshell.cpp */
double FRDMShellEnergy (MFTSystem *, const int, SPEnergy *, SCParameter *);
double FRDMSumSingleParticleEnergy (const int, SPEnergy *);

/*** FRDMzeroenergy.cpp */
double FRDMZeropointEnergy (PolarGLPoint *, MFTSystem *, Basis *, HFInterface *, FRDMEnergy *);

/*** FRDMprint.cpp */
void FRDMPrintShape (MFTSystem *);
void FRDMPrintMacroQuantity (const double , const double, const double, double *, const double, const double, const double, const double);
void FRDMPrintMacroEnergy (double *);
void FRDMPrintBeta (const int, double **, const double);
void FRDMPrintPotential (Basis *, OneBodyPotential *, LaguerrePolynomial *, HermitePolynomial *);
void FRDMPrintSPState (const bool, Basis *, SPEnergy *, LNParameter *);
void FRDMPrintBCS (LNParameter *);
void FRDMPrintSC (MFTSystem *, SCParameter *);
void FRDMPrintMicroEnergy (LNParameter *, SCParameter *);
void FRDMPrintTotalEnergy (FRDMEnergy *);

/*** FRDMfileout.cpp */
void FRDMWrite3DMesh (std::string, PolarCoordinate *);
void FRDMWriteParameters (const int, const int, LNParameter *, SCParameter *);
void FRDMWriteEnergy (std::string, MFTSystem *, FRDMEnergy *);

/*** FRDMreadshape.cpp */
void FRDMReadShape (MFTSystem *);

/*** FRLDMshape.cpp */
double FRLDMShape(double, QSurface *);
double FRLDMShapeDerivative(double, QSurface *);
int    FRLDMShapeCreateSurface (CylindricalGLPoint *, QSurface *, MFTSystem *);
int    FRLDMShapeCreateSpherical (CylindricalGLPoint *, QSurface *, MFTSystem *);
int    FRLDMCreateShape (QSurface *);
int    FRLDMShapeParameter (QSurface *);
void   FRLDMNeckLocation(QSurface *);
double FRLDMMultipoleMoment(const int, const int, const int, CylindricalGLPoint *);

void   FRLDMPrintShape(QSurface *);
void   FRLDMPrintQuadSurface(QSurface *);
void   FRLDMMonitorParameter(QSurface *);

/*** FRDMmacroenergy.cpp */
double FRLDMMacroEnergy(CylindricalGLPoint *, MFTSystem *);
double FRLDMMacroEnergySpherical (MFTSystem *);
double FRLDMMacroEnergyLiquidDrop (double, CylindricalGLPoint *, MFTSystem *);

/*** FRDMmacrointegral.cpp */
void   FRLDMMacroDoubleVolumeIntegral (double *, CylindricalGLPoint *, MFTSystem *);

/*** FRDMprint.cpp */
void   FRLDMPrintShape (MFTSystem *);
void   FRLDMPrintMacroQuantity (double *, double);
void   FRDMPrintMacroQuantity(const double eps, const double delta, const double w0, double *b, const double bs, const double bv, const double bw, const double br);
void   FRLDMPrintMacroEnergy (double *);
void   FRLDMPrintMacroEnergyLiquidDrop (double *);
void   FRLDMPrintTotalEnergy (FRDMEnergy *);

/*** FRDMfileout.cpp */
void   FRLDMWrite3DMesh (std::string, CylindricalCoordinate *);
void   FRLDMWriteEnergy(std::string, const int, const int, double, FRDMEnergy *, QSurface *, CylindricalGLPoint *);


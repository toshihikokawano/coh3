/******************************************************************************/
/*  FRDMmacroenergy.cpp                                                       */
/*        Macro Energy Part in FRDM                                           */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "FRDM.h"
#include "global.h"

static double EnergyMassExcess (MFTSystem *);
static double AverageAsymmetry (const double, const double, const double, const double, MFTSystem *);
static double AverageDeviation (const double, const double, const double, const double, MFTSystem *);
static double EnergyVolume (const double, const double, MFTSystem *);
static double EnergySurface (const double, const double, const double, MFTSystem *);
static double EnergyCurvature (const double, MFTSystem *);
static double EnergyA0 (void);
static double EnergyCoulomb (const double, const double, MFTSystem *);
static double EnergyVolumeRedistribution (const double, const double, MFTSystem *);
static double EnergyCoulombExchange (const double, MFTSystem *);
static double EnergySurfaceRedistribution (const double, const double, const double, const double, MFTSystem *);
static double EnergyProtonFormfactor (MFTSystem *);
static double EnergyChargeAsymmetry (MFTSystem *);
static double EnergyWigner (MFTSystem *);
static double EnergyPairing (const double, MFTSystem *);
static double EnergyBoundElectron (MFTSystem *);

static inline void FRDMSetMacroConstant(double *);

/**********************************************************/
/*   Macroscopic Energy Emac(Z,N,shape)                   */
/**********************************************************/
double FRDMMacroEnergy(PolarGLPoint *gp, MFTSystem *sys)
{
  /*** constants and each energy component */
  double c[4], b[4], e[15], **w;

  /*** set constants, c1, c2, c4, and c5 */
  FRDMSetMacroConstant(c);

  /*** nuclear energy, Coulomb energy, B1, B2, B3, and B4 */
  FRDMMacroDoubleVolumeIntegral(gp, b, sys);

  /*** pre-calculate W(r) at the surface */
  w = new double * [gp->ng];
  double vr[3];
  for(int i=0 ; i<gp->ng ; i++){
    w[i] = new double [gp->ng];
    for(int j=0 ; j<gp->ng ; j++){
      gp->cartesianVector(i,j,vr);
      w[i][j] = FRDMMacroW(vr,gp);
    }
  }

  FRDMMacroPrecalc(w,gp);

  /*** W-bar at the surface */
  double w0 = FRDMMacroWbar(gp,sys);

  double bs = FRDMMacroSurfaceEnergy(sys);
  double bv = FRDMMacroNeutronSkinEnergy(w0,sys);
  double bw = FRDMMacroSurfaceRedistributionEnergy(w0,sys);
  double bk = FRDMMacroCurvatureEnergy(sys);
  double br = FRDMMacroVolumeRedistributionEnergy(w0,gp,sys);


  /*** save Bs for microscopic model */
  sys->bs = bs;

  /*** each energy component */
  double delta = AverageAsymmetry(c[0],b[0],bv,bs,sys);
  double eps   = AverageDeviation(c[0],b[1],b[3],delta,sys);

  int i=0;
  e[i++] = EnergyMassExcess(sys);
  e[i++] = EnergyVolume(delta,eps,sys);
  e[i++] = EnergySurface(b[0],bs,delta,sys);
  e[i++] = EnergyCurvature(bk,sys);
  e[i++] = EnergyA0();
  e[i++] = EnergyCoulomb(c[0],b[2],sys);
  e[i++] = EnergyVolumeRedistribution(c[1],br,sys);
  e[i++] = EnergyCoulombExchange(c[2],sys);
  e[i++] = EnergySurfaceRedistribution(c[3],b[0],bs,bw,sys);
  e[i++] = EnergyProtonFormfactor(sys);
  e[i++] = EnergyChargeAsymmetry(sys);
  e[i++] = EnergyWigner(sys);
  e[i++] = EnergyPairing(bs,sys);
  e[i++] = EnergyBoundElectron(sys);

  double emac = 0.0;
  for(i=0 ; i<14 ; i++) emac += e[i];

  /*** print calculated quantities */
  if(pmf.energy){
    FRDMPrintMacroQuantity(eps,delta,w0/(1.6*PI*sys->R0*sys->R0),b,bs,bv,bw,br);
    FRDMPrintMacroEnergy(e);
  }

  for(int i=0 ; i<gp->ng ; i++){
    delete [] w[i];
  }
  delete [] w;

  return(emac);
}


/**********************************************************/
/*   Macroscopic Energy for Spherical Nucleus             */
/**********************************************************/
double FRDMMacroEnergySpherical(MFTSystem *sys)
{
  double c[4], b[4], e[15];

  FRDMSetMacroConstant(c);
  FRDMMacroVolumeIntegralSpherical(b, sys);

  double bs = 1.0;
  double bv = 1.0;
  double bw = 1.0;
  double bk = 1.0;
  double br = 1.0;

  double delta = AverageAsymmetry(c[0],b[0],bv,bs,sys);
  double eps   = AverageDeviation(c[0],b[1],b[3],delta,sys);

  int i=0;
  e[i++] = EnergyMassExcess(sys);
  e[i++] = EnergyVolume(delta,eps,sys);
  e[i++] = EnergySurface(b[0],bs,delta,sys);
  e[i++] = EnergyCurvature(bk,sys);
  e[i++] = EnergyA0();
  e[i++] = EnergyCoulomb(c[0],b[2],sys);
  e[i++] = EnergyVolumeRedistribution(c[1],br,sys);
  e[i++] = EnergyCoulombExchange(c[2],sys);
  e[i++] = EnergySurfaceRedistribution(c[3],b[0],bs,bw,sys);
  e[i++] = EnergyProtonFormfactor(sys);
  e[i++] = EnergyChargeAsymmetry(sys);
  e[i++] = EnergyWigner(sys);
  e[i++] = EnergyPairing(bs,sys);
  e[i++] = EnergyBoundElectron(sys);

  double emac = 0.0;
  for(i=0 ; i<14 ; i++) emac += e[i];

  return(emac);
}


/**********************************************************/
/*   Set Some Energy Related Constants                    */
/**********************************************************/
void FRDMSetMacroConstant(double *c)
{
  c[0] = 3.0*gCharge2 / (5.0*gRadius);
  c[1] = (1.0/gSymmetryEnergy + 18.0/gCompressibility)*c[0]*c[0]/336.0;
  c[2] = 5.0/4.0* pow(3.0/PI2, 2.0/3.0) * c[0];
  c[3] = c[0]*c[0]/(64.0*gSurfaceStiffness);
}


/**********************************************************/
/*   Average Bulk Nuclear Asynmmetry, delta-bar           */
/**********************************************************/
double AverageAsymmetry(const double c1, const double b1, const double bv, const double bs, MFTSystem *sys)
{
  double x1 = sys->asym
            + ( 3.0 * c1 * sys->getZ() * bv * bs)
             /(16.0 * gSurfaceStiffness * sys->a23 * b1);
  double x2 = 1 
            + (9.0 * gSymmetryEnergy * bs * bs)
             /(4.0 * gSurfaceStiffness * sys->a13 * b1);

  return(x1/x2);
}


/**********************************************************/
/*   Average Relative Deviation in Density, eps-bar       */
/**********************************************************/
double AverageDeviation(const double c1, const double b2, const double b4, const double delta, MFTSystem *sys)
{
  double x1 = gCompExp * exp(-gCompExpRange * sys->a13);
  double x2 = 2.0 * gSurfaceEnergy * b2 / sys->a13;
  double x3 = gDensitySymmetry * delta * delta;
  double x4 = c1 * sys->getZ() * sys->getZ() * b4 / (sys->a23 * sys->a23);

  return( (x1-x2+x3+x4)/gCompressibility );
}


/**********************************************************/
/*   Mass Excess of Z and N                               */
/**********************************************************/
double EnergyMassExcess(MFTSystem *sys)
{
  return( ENEUTRON * sys->getN() + EPROTON * sys->getZ() );
}


/**********************************************************/
/*   Volume Energy                                        */
/**********************************************************/
double EnergyVolume(const double delta, const double eps, MFTSystem *sys)
{
  double x1 = gVolumeEnergy;
  double x2 = gSymmetryEnergy * delta * delta;
  double x3 = 0.5 * gCompressibility * eps * eps;
  double e  = (-x1 + x2 - x3) * sys->getA();

  return(e);
}


/**********************************************************/
/*   Surface Energy                                       */
/**********************************************************/
double EnergySurface(const double b1, const double bs, const double delta, MFTSystem *sys)
{
  double x1 = gSurfaceEnergy * b1;
  double x2 = 9.0/4.0 * gSymmetryEnergy * gSymmetryEnergy * delta * delta * bs * bs / (gSurfaceStiffness * b1);
  double e  = (x1 + x2) * sys->a23;

  return(e);
}


/**********************************************************/
/*   Curvature Energy                                     */
/**********************************************************/
double EnergyCurvature(const double bk, MFTSystem *sys)
{
  double e = gCurvatureEnergy * bk * sys->a13;
  return(e);
}


/**********************************************************/
/*   A0 Energy                                            */
/**********************************************************/
double EnergyA0(void)
{
  double e = 0.0;
  return(e);
}


/**********************************************************/
/*   Coulomb Energy                                       */
/**********************************************************/
double EnergyCoulomb(const double c1, const double b3, MFTSystem *sys)
{
  double e = c1 * sys->getZ() * sys->getZ() * b3 / sys->a13;
  return(e);
}


/**********************************************************/
/*   Volume Redistribution Energy                         */
/**********************************************************/
double EnergyVolumeRedistribution(const double c2, const double br, MFTSystem *sys)
{
  double e = - c2 * sys->getZ() * sys->getZ() * sys->a13 * br;
  return(e);
}


/**********************************************************/
/*   Coulomb Exchange Correction                          */
/**********************************************************/
double EnergyCoulombExchange(const double c4, MFTSystem *sys)
{
  double e = - c4 * pow((double)sys->getZ(),4.0/3.0) / sys->a13;
  return(e);
}


/**********************************************************/
/*   Surface Redistribution Energy                        */
/**********************************************************/
double EnergySurfaceRedistribution(const double c5, const double b1, const double bs, const double bw, MFTSystem *sys)
{
  double e = - c5 * sys->getZ() * sys->getZ() * bw * bs / b1;
  return(e);
}


/**********************************************************/
/*   Proton-formfactor Correction to Coulomb              */
/**********************************************************/
double EnergyProtonFormfactor(MFTSystem *sys)
{
  const double c = -1.0/8.0 * 145.0/48.0;
  double f0 = c * gProtonRadius * gProtonRadius * gCharge2
              /(gRadius * gRadius * gRadius);
  double e = f0 * sys->getZ() * sys->getZ() / (double)sys->getA();
  return(e);
}


/**********************************************************/
/*   Charge Asymmetry Energy                              */
/**********************************************************/
double EnergyChargeAsymmetry(MFTSystem *sys)
{
  double e = - gChargeAsymmetry * (sys->getN() - sys->getZ());
  return(e);
}


/**********************************************************/
/*   Wigner Energy                                        */
/**********************************************************/
double EnergyWigner(MFTSystem *sys)
{
  double x = fabs(sys->asym);

  if((sys->getZ() == sys->getN()) && ((sys->getZ()%2) == 1)){
    x += 1.0/sys->getA();
  }
  double e = gWigner * x;
  return(e);
}


/**********************************************************/
/*   Pairing Energy                                       */
/**********************************************************/
double EnergyPairing(const double bs, MFTSystem *sys)
{
  double dn = gPairingGap * bs * pow((double)sys->getN(),-1.0/3.0);
  double dp = gPairingGap * bs * pow((double)sys->getZ(),-1.0/3.0);
  double delta = gNPInteraction / (bs * sys->a23);

  bool oddZ = (sys->getZ()%2 == 1) ? true : false;
  bool oddN = (sys->getN()%2 == 1) ? true : false;

  double e = 0.0;
  if(oddZ && oddN){       e = dp + dn - delta; }
  else if(oddZ && !oddN){ e = dp; }
  else if(!oddZ && oddN){ e = dn; }
  return(e);
}


/**********************************************************/
/*   Boundary Electron Energy                             */
/**********************************************************/
double EnergyBoundElectron(MFTSystem *sys)
{
  double e = - gEBinding * pow((double)sys->getZ(),2.39);
  return(e);
}

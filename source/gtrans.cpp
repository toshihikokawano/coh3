/******************************************************************************/
/*  gtrans.cpp                                                                */
/*        Gamma-ray transmission coefficients                                 */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "statmodel.h"
#include "parameter.h"
#include "global.h"

static double  gdrStandardLorentzian (const int, const double, GDR *);
static double  gdrGeneralizedLorentzian (const int, const double, GDR *, const double, const double);
static double  gdrModifiedLorentzian (const int, const double, GDR *, const double, const double);
static inline void    gdrDropletGammaSave (const unsigned int);
static inline double  gdrProfileLorentzian (const double, const double, const double);
static inline double  gdrTempDependGamma (const double, const double, const double, const double);
static inline double  gdrTempDependGamma2 (const double, const double, const double, const double);
static inline double  gdrGammagSystematics (const int);

static double gStrengthFunction = 1.0;
static double gDropletGamma = 0.0;
static const double gStrNormFactor = 1.0 / (HBARSQ * VLIGHTSQ * PI * PI) / NORM_FACT;

static GDR *gdr[MAX_GDR];
static int ngdr = 0;


/**********************************************************/
/*      Normalization Factor for M1 Transmission          */
/**********************************************************/
void gdrM1norm(const double a, const double ldpa, GDR *g)
{
  double fe1 = 0.0, fm1 = 0.0;

  /*** scan all E1 GDR */
  for(int i=0 ; i<MAX_GDR ; i++){
    if(g[i].getXL() == "E1"){
      fe1 += gdrGeneralizedLorentzian(1,7.0,&g[i],ldpa,0.0);
    }
  }

  /*** look for the first M1, should be one */
  int m = 0;
  for(int i=0 ; i<MAX_GDR ; i++){
    if(g[i].getXL() == "M1"){
      fm1 += gdrStandardLorentzian(1,7.0,&g[i]);
      m = i;
      break;
    }
  }
  g[m].setSigma(17.01*pow(a,-0.878)*fe1/fm1);
}


/**********************************************************/
/*      Set GDR Parameters in Scope Here                  */
/**********************************************************/
void gdrParameterSave(GDR *gdr0, const unsigned int a)
{
  ngdr = 0;
  for(int i=0 ; i<MAX_GDR ; i++){
    if(gdr0[i].getXL() != "  ") gdr[ngdr++] = &gdr0[i];
  }

  /*** save width predicted by droplet model */
  gdrDropletGammaSave(a);
}


/**********************************************************/
/*      Reset GDR Normalization Factor                    */
/**********************************************************/
void gdrParameterReset()
{
  gStrengthFunction = 1.0;
}


/**********************************************************/
/*      Save Gamma by Droplet Model                       */
/*      W.D. Myers et al. Phys. Rev. C 15, 2032 (1977)    */
/**********************************************************/
void gdrDropletGammaSave(const unsigned int a)
{
  gDropletGamma = 36.43 * pow((double)a,-1.0/3.0);
}


/**********************************************************/
/*      Gamma-ray Transmission                            */
/**********************************************************/
double gdrGammaTransmission(const GammaMultipolarity gm, const double eg, const double a, const double ex)
{
  double tg = 0.0;
  int    l  = 0;
  std::string xl = "  ";

  switch(gm){
  case E1: xl = "E1"; l = 1; break;
  case M1: xl = "M1"; l = 1; break;
  case E2: xl = "E2"; l = 2; break;
  case M2: xl = "M2"; l = 2; break;
  case E3: xl = "E3"; l = 3; break;
  case M3: xl = "M3"; l = 3; break;
  case E4: xl = "E4"; l = 4; break;
  case M4: xl = "M4"; l = 4; break;
  default:                   break;
  }
  if(l == 0) return(0.0);
 
  /*** first, if gSTF data table is given, interpolate data */
  if(opt.readphotoabsorption) tg = gdrGammaTransmissionInterpolate(gm,eg);

  /*** usual case, or when table interpolation retuned zero */
  if(tg == 0.0){
    tg = 0.0;
    for(int i=0 ; i<ngdr ; i++){
      if(gdr[i]->getXL() == xl){

        GammaProfile c = gdr[i]->getProfile();

        if(c == SL){
          tg += gdrStandardLorentzian(l,eg,gdr[i]);
        }
        else if(c == GL){
          tg += gdrGeneralizedLorentzian(l,eg,gdr[i],a,ex);
        }
        else if(c == ML){
          tg += gdrModifiedLorentzian(l,eg,gdr[i],a,ex);
        }
      }
    }
  }

  /*** transmission = 2 pi x E^{2L+1} gSTR */
  tg *= PI2*pow(eg,2*l+1.0)*gStrengthFunction;

  /*** if normalization factor is given, apply that factor */
  double tgs = parmGetFactor(parmTJ,gammaray);

  if(tgs > 0.0) tg *= tgs;

//  if(xl == "E1") cout << eg << " " << tg << endl;

  return(tg);
}


/**********************************************************/
/*      Brink-Axel Standard Lorentzian Function           */
/**********************************************************/
double gdrStandardLorentzian(const int l, const double eg, GDR *gam)
{
  double c,f;

  if((eg == 0.0) && (l > 0)) return(0.0);

  c = gStrNormFactor / (2*l+1.0);
  f = c * gam->getSigma() * gam->getWidth() * pow(eg,2.0-2.0*l)
        * gdrProfileLorentzian(eg,gam->getEnergy(),gam->getWidth());

   return(f);
}


/**********************************************************/
/*      Kopecky-Uhl Generalized Lorentzian Function       */
/**********************************************************/
double gdrGeneralizedLorentzian(const int l, const double eg, GDR *gam, const double a, const double ex)
{
  double c,g,f,t;

  c = gStrNormFactor / (2*l+1.0);
  t = (ex < 0.0) ? 0.0 : sqrt(ex/a);
  g = (l==1) ? gdrTempDependGamma(eg,gam->getEnergy(),gam->getWidth(),t) 
             : gam->getWidth();
  f = c * gam->getSigma() * gam->getWidth() * pow(eg,2.0-2.0*l)
      *( gdrProfileLorentzian(eg,gam->getEnergy(),g) 
      +( (l==1) ? 0.7*gdrTempDependGamma(0.0,gam->getEnergy(),gam->getWidth(),t)
                  /(gam->getEnergy() * gam->getEnergy() * gam->getEnergy())
                : 0.0 )  );
  return(f);
}


/**********************************************************/
/*      Plujko Modified Lorentzian Function               */
/**********************************************************/
double gdrModifiedLorentzian(const int l, const double eg, GDR *gam, const double a, const double ex)
{
  double c,g,f,t,p;

  c = gStrNormFactor / (2*l+1.0);
  t = (ex < 0.0) ? 0.0 : sqrt(ex/a);
  p = (t > 0.0) ? 1.0/(1.0 - exp(-eg/t)) : 1.0;

  g = (l==1) ? gdrTempDependGamma2(eg,gam->getEnergy(),gam->getWidth(),ex) 
             : gam->getWidth();

  f = c * p * gam->getSigma() * gam->getWidth() * pow(eg,2.0-2.0*l)
        * gdrProfileLorentzian(eg,gam->getEnergy(),g);

  return(f);
}


/**********************************************************/
/*      Renormalize E1 Strength Function for neutron      */
/**********************************************************/
double  gdrRenormE1StrengthFunction(const int pt, const int jt, const int a, const double d0, const double sn, Nucleus *n, double **tg)
{
  /*** call this only once */
  if(gStrengthFunction != 1.0) return(0.0);
  if(sn <= 0.0) return(0.0);

  /*** copy GDR parameters to local */
  gdrParameterSave(n->gdr,n->za.getA());

  /*** retrieve input <Gg>/D0 */
  double gdinput = parmGetValue(parmGSTR);

  /*** save bin parameters for re-binning the continuum */
  double ebin = n->de;
  double emax = n->max_energy;

  n->de = 0.02;               // 20 keV bin
  n->max_energy = sn + 0.001; // dummy 1-keV neutron incident energy

  statSetupEnergyBin(n);
  statSetupLevelDensity(n,&n->ldp);

  /*** prepare gamma-ray transmission coefficients */
  statStoreContinuumGammaTransmission(0,tg,n);

  /*** for s-wave neutron, J = I+1/2, s1 and s2 are dummies */
  double s1[1], s2[1], tsum = 0.0;
  if(jt == 0){
    tsum = specTransitionGamma(sumall,0,pt,1   ,tg,0.0,n,s1,s2);
  }else{
    tsum = specTransitionGamma(sumall,0,pt,jt+1,tg,0.0,n,s1,s2)
         + specTransitionGamma(sumall,0,pt,jt-1,tg,0.0,n,s1,s2);
  }

  double gammag = tsum * d0 /PI2; // calculated <Gamma_g>

  /*** no normalization
       if neutron separation energy is inside the discrete level region */
  if( (n->ncont > 0) && (tsum > 0.0) && gdinput != 0.0){
    /*** if negative Gg is given, renormalize by systematics */
    if(gdinput < 0.0){

      /*** Gamma_g systematics for Gg = -1 case */
      if(gdinput == -1.0){
        gammag = gdrGammagSystematics(a);
      }
      /*** direct Gamma_g input */
      else{
        gammag = fabs(gdinput);
      }
    }
    else if(gdinput > 0.0) gammag = gdinput * d0;

    /*** save Tg normalization constant */
    gStrengthFunction = PI2 * gammag / d0 / tsum;

    /*** ad hoc limit to the normalization to avoid unphysical value */
    if(gStrengthFunction > 5.0) gStrengthFunction = 5.0;
    else if(gStrengthFunction < 0.3) gStrengthFunction = 0.3;
  }

  n->de = ebin;
  n->max_energy =emax;
  statSetupEnergyBin(n);
  statSetupLevelDensity(n,&n->ldp);

  return(gammag);
}


/**********************************************************/
/*      Temperature dependent Gamma, by Kadmenskij        */
/**********************************************************/
double gdrTempDependGamma(const double eg, const double er, const double w, const double t)
{
  if(er == 0.0) return(0.0);
  return( w*(eg*eg + 4.0*PI*PI*t*t)/(er*er) );
}


/**********************************************************/
/*      Temperature dependent Gamma, by Plujko            */
/**********************************************************/
double gdrTempDependGamma2(const double eg, const double er, const double w, const double ex)
{
  if(er == 0.0) return(0.0);

  const double c = 5.42e-04;

  double gammac =  c * er * (eg + ex);

  double k0 = 0.3;
  double kr = (w - c*er*er)/gDropletGamma;
  double ks = (eg >= 2.0*er) ? k0 : kr + (k0 - kr)*fabs((eg - er)/er);

  double gammaf = ks * gDropletGamma;

  return(gammac + gammaf);
}


/**********************************************************/
/*      Lorentzian Profile                                */
/**********************************************************/
double gdrProfileLorentzian(const double eg, const double er, const double w)
{
  if((er == eg) && (eg*w == 0.0)) return(0.0);
  return (eg*w /( (er*er-eg*eg)*(er*er-eg*eg) + eg*eg*w*w )); 
}


/**********************************************************/
/*      Gamma_gamma Approximation                         */
/**********************************************************/
double gdrGammagSystematics(int a)
{
  double gg = 2.9666e-03 * pow((double)a,-2.0801);
  return(gg);
}


/******************************************************************************/
/*  gdr.cpp                                                                   */
/*        global GDR parameter library                                        */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "gdr.h"

static void gdrRatio (GDR *, GDR *, std::string, const double);

static const double SigmaRatio = 0.0008;

/**********************************************************/
/*      E1                                                */
/*      M. Herman et al. EMPIRE-3.2 Malta Manual          */
/*      Eq. (4.108) -- (4.116)                            */
/**********************************************************/
void gdrE1(const double mass, const double beta, GDR *g)
{
  double energy  = (49.336+7.34*beta)*pow(mass,-0.2409);
  double width   = 0.3 * energy;
  double sigma   = 10.6 * mass / width;
  std::string XL = "E1";
  g->setGDR(XL,energy,width,sigma,GL);
}

void gdrE1DoubleHump0(const double mass, const double beta, GDR *g)
{
  if(beta <= 0.064){
    gdrE1(mass,beta,g);
  }
  else{
    double energy  = 50.0*pow(mass,-0.232) * exp(-0.946*beta);
    double width   = (0.282-0.263*beta) * energy;
    double sigma   = 3.48  * mass / width;
    std::string XL = "E1";
    g->setGDR(XL,energy,width,sigma,GL);
  }
}

void gdrE1DoubleHump1(const double mass, const double beta, GDR *g)
{
  if(beta <= 0.064){
    g->clear();
  }
  else{
    double energy  = 50.0*pow(mass,-0.232);
    double width   = (0.35 -0.14 *beta) * energy;
    //double sigma  = 8.26  * mass / width;
    double sigma   = 1.464 * pow(mass,4.0/3.0) / width;
    std::string XL = "E1";
    g->setGDR(XL,energy,width,sigma,GL);
  }
}


/**********************************************************/
/*      M1                                                */
/**********************************************************/
void gdrM1(const double mass, GDR *g)
{
  double energy  = 41.0/pow(mass,1/3.0);
  double width   = 4.0;
  double sigma   = 1.0;
  std::string XL = "M1";
  g->setGDR(XL,energy,width,sigma,SL);
}


/**********************************************************/
/*      M1 Scissors                                       */
/**********************************************************/
void gdrM1scissors(const double mass, const double beta, GDR *g)
{
  double energy  = 80.0*fabs(beta)/pow(mass,1/3.0);
  double width   = 1.5;
  double sigma   = 42.4*beta*beta / width;
  std::string XL = "M1";
  g->setGDR(XL,energy,width,sigma,SL);
}


/**********************************************************/
/*      E2                                                */
/**********************************************************/
void gdrE2(const double z, const double mass, GDR *g)
{
  double a3 = 1./pow(mass,1/3.0);
  double energy  = 63.0*a3;
  double width   = 6.11 - 0.012*mass;
  double sigma   = 1.5e-4 * z * z * energy * energy*a3 / width;
  std::string XL = "E2";
  g->setGDR(XL,energy,width,sigma,SL);
}



/**********************************************************/
/*      M2                                                */
/**********************************************************/
void gdrM2(GDR *gm1, GDR *g)
{
  gdrRatio(gm1,g,"M2",SigmaRatio);
}


/**********************************************************/
/*      E3                                                */
/**********************************************************/
void gdrE3(GDR *ge2, GDR *g)
{
  gdrRatio(ge2,g,"E3",SigmaRatio);
}


/**********************************************************/
/*      M3                                                */
/**********************************************************/
void gdrM3(GDR *gm2, GDR *g)
{
  gdrRatio(gm2,g,"M3",SigmaRatio);
}


/**********************************************************/
/*      E4                                                */
/**********************************************************/
void gdrE4(GDR *ge3, GDR *g)
{
  gdrRatio(ge3,g,"E4",SigmaRatio);
}


/**********************************************************/
/*      Setting by Relative Strength                      */
/**********************************************************/
void gdrRatio(GDR *src, GDR *dst, std::string x, const double r)
{
  dst->setGDR(x, src->getEnergy(), src->getWidth(), src->getSigma()*r, src->getProfile());
}


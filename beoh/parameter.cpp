/******************************************************************************/
/*  parameter.cpp                                                             */
/*        parameter adjustment                                                */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "parameter.h"
#include "mt19937.h"

static double parmGaussianDistribution  (double, double);

/**********************************************************/
/*      Check range of given parameter                    */
/**********************************************************/
void    parmCheckFactor()
{
  /*** If width is given, Gaussian distribution
       with the averave of 1.0 and the width of w */
  for(int i=0 ; i < adj.nparm() ; i++){

    double w = adj.parm[i].getwidth();

    if(w > 0.0){
      double f = adj.parm[i].getfactor();
      /*** infinite loop to skip negative factor */
      for(;;){
        double g = parmGaussianDistribution(f, w);
        if(g > 0.0){
          adj.parm[i].setfactor(g); 
          adj.parm[i].setwidth(0.0);
          break;
        }
      }
    }
  }
}


/**********************************************************/
/*      Type: no argument                                 */
/**********************************************************/
double  parmGetFactor(const ParmType t)
{
  int    k = adj.getindex(t);
  double f = (k >= 0) ? adj.parm[k].getfactor() : 1.0;
  return(f);
}


/**********************************************************/
/*      Type: require particle channel                    */
/**********************************************************/
double  parmGetFactor(const ParmType t, const Particle p)
{
  int    k = adj.getindex(t,p);
  double f = (k >= 0) ? adj.parm[k].getfactor() : 1.0;
  return(f);
}


/**********************************************************/
/*      Type: require ZA number                           */
/**********************************************************/
double  parmGetFactor(const ParmType t, const ZAnumber za)
{
  int    k = adj.getindex(t,za);
  double f = (k >= 0) ? adj.parm[k].getfactor() : 1.0;
  return(f);
}


/**********************************************************/
/*      Type: require ZA number and paricle channel       */
/**********************************************************/
double  parmGetFactor(const ParmType t, const ZAnumber za, const Particle p)
{
  int    k = adj.getindex(t,za,p);
  double f = (k >= 0) ? adj.parm[k].getfactor() : 1.0;
  return(f);
}


/**********************************************************/
/*      Return Double: no argument                        */
/**********************************************************/
double  parmGetValue(const ParmType t)
{
  int    k = adj.getindex(t);
  double f = (k >= 0) ? adj.parm[k].getfactor() : 0.0;
  return(f);
}


/**********************************************************/
/*      Return Int: require ZA number                     */
/**********************************************************/
int     parmGetValue(const ParmType t, const ZAnumber za)
{
  int    k = adj.getindex(t,za);
  int    f = (k >= 0) ? (int)adj.parm[k].getfactor() : 0;
  return(f);
}


/**********************************************************/
/*      Generate Gaussian Distribution                    */
/**********************************************************/
double  parmGaussianDistribution(const double m, const double s)
{
  static double x1 = 0.0, x2 = 0.0;
  static bool   flag = false;
  double a,b,r;

  if(flag){
    flag = false;
    return(s*x2+m);
  }
  else{
    do{
      a = 2.0*genrand_real3() - 1.0;
      b = 2.0*genrand_real3() - 1.0;
      r  = a*a+b*b;
    }while (r>=1.0 || r==0.0);

    r = sqrt(-2.0*log(r)/r);
    x1 = a*r;
    x2 = b*r;
    flag = true;
  }

  return(s*x1+m);
}

/******************************************************************************/
/*  skyrmeparameter.cpp                                                       */
/*        Definition of Skyrme Parameters                                     */
/******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "HF.h"

void SkyrmeParameter(std::string skname, Skyrme *sp)
{
  double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
  double x0 = 0.0, x1 = 0.0, x2 = 0.0, x3 = 0.0;
  double w  = 0.0, a  = 0.0, v0 = 0.0;

  if(skname == "SIII"){ // Nucl. Phys. A264, 197
    t0 = -1128.75;
    t1 =   395.0;
    t2 =   -95.0;
    t3 = 14000.0;
    x0 =     0.45;
    x1 =     0.0;
    x2 =     0.0;
    x3 =     1.0;
    w  =   120.0;
    a  =     1.0;
    v0 =     0.0;
  }
  else if(skname == "SkM*"){
    t0 = -2645.0;
    t1 =   410.0;
    t2 =  -135.0;
    t3 = 15595.0;
    x0 =     0.09;
    x1 =     0.0;
    x2 =     0.0;
    x3 =     0.0;
    w  =   130.0;
    a  =     1.0/6.0;
    v0 =     0.0;
  }
  else if(skname == "SLy4"){ // Nucl. Phys. A635, 231
    t0 =  -2488.91;
    t1 =    486.82;
    t2 =   -546.39;
    t3 =  13777.0;
    x0 =      0.834;
    x1 =     -0.344;
    x2 =     -1.0;
    x3 =      1.354;
    w  =    123.0;
    a  =      1.0/6.0;
    v0 =      0.0;
  }
  else{
    std::cerr << "given Skyrme force " << skname << " not found" << std::endl;
    exit(-1);
  }


  /*** D. Vautherin, PRC 5, 626 (1972), Appendix A */
  sp->b1 =   t0 * (1.0 + x0/2.0) / 2.0;  // coeff. for rho^2
  sp->b2 =  -t0 * (x0 + 0.5) / 2.0;      // coeff. for rhon^2 + rhop^2
  sp->b3 =  (t1 * (1.0 + x1/2.0) + t2 * (1.0 + x2/2.0)) / 4.0;  // 1/4 (t1+t2) tau
  sp->b4 = -(t1 * (0.5 + x1) - t2 * (0.5 + x2)) / 4.0;          // 1/4 (t2-t1) tau_p
  sp->b5 = -(3.0 * t1 * (1.0 + x1/2.0) - t2 * (1.0 + x2/2.0)) / 16.0; // -1/8 (3t1-t2) div2 rho
  sp->b6 =  (3.0 * t1 * (0.5 + x1) + t2 * (0.5 + x2)) / 16.0;         // 1/16 (3t1+t2) div2 rho_p
  sp->b7 =  t3 * (1.0 + x3/2.0) / 12.0;  // for alpha
  sp->b8 = -t3 * (0.5 + x3) / 12.0;      // for alpha
  sp->b9 = -w / 2.0;  // W0/2 for divJ

  sp->alpha = a;
  sp->w = w;
  sp->v0 = v0;
/*
  cout << sp->b1 << endl;
  cout << sp->b2 << endl;
  cout << sp->b3 << endl;
  cout << sp->b4 << endl;
  cout << sp->b5 << endl;
  cout << sp->b6 << endl;
  cout << sp->b7 << endl;
  cout << sp->b8 << endl;
  cout << sp->b9 << endl;
*/
}

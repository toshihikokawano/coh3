/***********************************************************/
/*      modified Koning-Delaroche                          */
/*      by B-spline function                               */
/***********************************************************/

#include <iostream>
#include <fstream>
#include <cmath>

#include "spline.h"
#include "omplib.h"

#define SPLINE_FROM_FILE
#ifdef SPLINE_FROM_FILE
static const char *filename = "spline.dat";
#endif

unsigned int OMPspline(double e, int at, int zt, int ai, int zi, Optical *omp)
{
  unsigned int pf = OMPKoningDelaroche(e,at,zt,ai,zi,omp);
//unsigned int pf = OMPKunieda(e,at,zt,ai,zi,omp);

#ifdef SPLINE_FROM_FILE

  std::ifstream fs(filename);
  double emax  = 0.0;
  int    ndata = 0;
  double x,fv,fw,frv,frw,fav,faw;

  fs >> ndata;
  Spline spv(ndata),spw(ndata),sprv(ndata);
  Spline sprw(ndata),spav(ndata),spaw(ndata);

  for(int i=0 ; i<ndata ; i++){
    fs >> x;
    fs >> fv; fs >> frv; fs >> fav;
    fs >> fw; fs >> frw; fs >> faw;

    spv.Setdat(x,fv);
    spw.Setdat(x,fw);
    sprv.Setdat(x,frv);
    sprw.Setdat(x,frw);
    spav.Setdat(x,fav);
    spaw.Setdat(x,faw);
    if(emax < x) emax = x;
  }
  fs.close();
#else
  double emax  = 20.0;
  int    ndata = 5;
  Spline spv(ndata),spw(ndata),sprv(ndata),sprw(ndata);

  spv.Setdat( 0.0, 1.000); spw.Setdat( 0.0, 1.000); sprv.Setdat( 0.0, 1.0);  sprw.Setdat( 0.0, 1.0);
  spv.Setdat( 5.0, 1.050); spw.Setdat( 5.0, 1.050); sprv.Setdat( 5.0, 1.0);  sprw.Setdat( 5.0, 1.0);
  spv.Setdat(10.0, 1.100); spw.Setdat(10.0, 1.100); sprv.Setdat(10.0, 1.0);  sprw.Setdat(10.0, 1.0);
  spv.Setdat(15.0, 1.000); spw.Setdat(15.0, 1.000); sprv.Setdat(15.0, 1.0);  sprw.Setdat(15.0, 1.0);
  spv.Setdat(20.0, 0.950); spw.Setdat(20.0, 0.950); sprv.Setdat(20.0, 1.0);  sprw.Setdat(20.0, 1.0);
#endif

  if(e < emax){
    splineSupport(ndata,&spv);
    splineSupport(ndata,&spw);
    splineSupport(ndata,&sprv);
    splineSupport(ndata,&sprw);
    splineSupport(ndata,&spav);
    splineSupport(ndata,&spaw);
    omp->v1  *= splineInterpolation(ndata,e,&spv);
    omp->ws1 *= splineInterpolation(ndata,e,&spw);
    omp->r0  *= splineInterpolation(ndata,e,&sprv);
    omp->rs  *= splineInterpolation(ndata,e,&sprw);
    omp->a0  *= splineInterpolation(ndata,e,&spav);
    omp->as  *= splineInterpolation(ndata,e,&spaw);
  }

  return(pf);
}

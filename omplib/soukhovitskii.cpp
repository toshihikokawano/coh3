/*************************************************/
/*     E. Soukhovitskii, S. Chiba, J.Y. Lee,     */
/*     O. Iwamoto, T. Fukahori                   */
/*     J. Nucl.Phys. G 30, 905 (2004)            */
/*************************************************/

#include <cmath>

#include "omplib.h"


unsigned int OMPSoukhovitskii(double e,int at, int zt, int zi, Optical *omp)
{
  double v0,v1,v2,va,vd,wsd,wsa,wvd,vso,wsod,wis,wiv,wi0,wiso,
         lmdv,lmds,lmdso,cvso,cwso,phi;

  double a13 = pow((double)at,1.0/3.0);
  double eps = 1.0-2.0*(double)zt/(double)at;

  double ef = OMPFermiEnergy(at,zt,zi);
  double ex = e - ef;

  v0   =  -41.45;
  v1   =    0.03;
  v2   =    0.000205;
  vd   =   92.44;
  va   =   -0.06667;
  wsd  =   17.38;
  wsa  =    0.03833;
  wvd  =   14.74;
  vso  =    5.86;
  wsod =   -3.1;
  wis  =   11.79;
  wiv  =   81.63;
  wi0  =  100.0;
  wiso =  160.00;

  lmdv  =   0.0039075;
  lmds  =   0.01759;
  lmdso =   0.005;
  cvso  =  10.5;
  cwso  =  24.0;

/* common geometry part */
  omp->r0   = 1.245 * (1 - 0.05*(ex*ex/(ex*ex + wi0*wi0)) );
  omp->rs   = 1.2080;
  omp->rv   = 1.2476;
  omp->rvso = 1.1213;
  omp->rwso = omp->rvso;

  omp->rc   = (zi==0) ? 0.0: 1.2643;

  omp->a0   = 0.660 + 0.000253*e;
  omp->as   = 0.614;
  omp->av   = 0.594;
  omp->avso = 0.59;
  omp->awso = omp->avso;

  phi = (zi == 0) ? 0.0 :
                 (- v1 - 2*v2*ex + vd*lmdv*exp(-lmdv*ex))
                *(1.0 + ((zi==1)? 1: -1)*cvso*eps/(v0 + va*(at-232) + vd))
                * 0.90*(double)zt/a13;

  omp->v1   = (v0 + va*(at-232)+ v1*ex + v2*ex*ex + vd*exp(-lmdv*ex))
              *(1.0 + ((zi==1)? 1: -1)*cvso*eps/(v0 + va*(at-232) + vd)) + phi;
  omp->v2   =  0.0;
  omp->v3   =  0.0;

  omp->ws1  = (wsd + wsa*(at-232) + ((zi==1)? 1: -1)*cwso*eps)
              *exp(-lmds*ex) * ex*ex/(ex*ex+wis*wis);
  omp->ws2  = 0.0;
  omp->ws3  = 0.0;

  omp->wv1  = wvd * ex*ex/(ex*ex+wiv*wiv);
  omp->wv2  = 0.0;
  omp->wv3  = 0.0;

  omp->vso1 = vso*exp(-lmdso*ex);
  omp->vso2 = 0.0;

  omp->wso1 = wsod *ex*ex/(ex*ex+wiso*wiso);
  omp->wso2 = 0.0;

  return(0x0013);
}


/*************************************************/
/*     E. Soukhovitskii, R. Capote               */
/*     J. M. Quesada, S. Chiba                   */
/*     Phys. Rev. C 72, 024604 (2005)            */
/*************************************************/
unsigned int OMPSoukhovitskii2005(double e,int at, int zt, int zi, Optical *omp)
{
  double epst = 1.0-2.0*(double)zt/(double)at;
  double epsp = -1.0 + 2.0*(double)zi;

  double ef  = OMPFermiEnergy(at,zt,zi);
  double ex  = e - ef;
  double ex2 = ex*ex;

  double v0, lambdahf, cviso, av, bv, w0, bs, cs, cwiso, ea,// ccoul,
         vso, lambdaso, wso, bso;

  v0        =  49.97;
  lambdahf  =   0.01004;
  cviso     =  15.9;
  av        =  12.04;
  bv        =  81.36;
  ea        = 385.0;
  w0        =  17.20;
  bs        =  11.19;
  cs        =   0.01361;
  cwiso     =  23.5;
  vso       =   5.75;
  lambdaso  =   0.005;
  wso       =  -3.1;
  bso       = 160.0;
//ccoul     =   1.3;

  /*** Pu-239, Iwamoto */
  if( (zt == 94) && (at == 239) ){ 
    v0        =  50.054;
    w0        =  17.1463;
  }   


  omp->r0   = omp->rv   = 1.2568;
  omp->rs   = omp->r0s  = 1.1803;
  omp->rvso = omp->rwso = 1.1214;
  omp->rc   = 1.2452;
  
  omp->a0   = omp->av   = 0.633 ;
  omp->as   = omp->a0s  = 0.601 ;
  omp->avso = omp->awso = 0.590 ;
//omp->ac   = 0.545 ;
  
    
  double ahf = v0 + epst*epsp * cviso;
  double as  = w0 + epst*epsp * cwiso;

  omp->v1   = ahf * exp(-lambdahf*ex);
  omp->wv1  = av * ex2/(ex2 + bv*bv);
  omp->ws1  = as * ex2/(ex2 + bs*bs) * exp(-cs * abs(ex));

  omp->vso1 = vso * exp(-lambdaso * ex);
  omp->wso1 = wso * ex2/(ex2 + bso*bso);

  omp->vs1   = OMPdeltaSurface(ex,0.0,as,bs,cs);
  omp->v1   += OMPdeltaVolume(ex,0.0,av,bv) + OMPasymmetricVolume(e,ef,ea,av,bv);
  omp->vso1 += OMPdeltaVolume(ex,0.0,wso,bso);

  return(0x0013);
}


/*************************************************/
/*    modified Soukhovitskii                     */
/*        fitted to resonances                   */
/*************************************************/
unsigned int OMPModSoukhovitskii(double e,int at, int zt, int zi, Optical *omp)
{
  unsigned int id = OMPSoukhovitskii(e,at,zt,zi,omp);

/* for U238 */
// if(e<0.735) omp->ws1 = 2.37;
// if(e<1.130) omp->ws1 = 2.59;


/* for Np237 */
// if(e<1.130) omp->ws1 = 2.59;

/* for Am241 */
// omp->ws1 *= 0.709;

/* for Pu239 */
//   if(e<1.130) omp->ws1 = 2.15;
  omp->ws1 *= 0.7;

   return(id);
}

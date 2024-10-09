/******************************************************************************/
/**     Wahl's Systematics, LA-13928                                         **/
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "wahl.h"
#include "physicalconstant.h"

static double SigmaSystematics_Wahl (const int, const int, const double);
static double DeltaSystematics_Wahl (const int, const int, const double, const double);
static double FractSystematics_Wahl (const int, const int, const int, const double);

static inline double ParameterFunction3 (const int, double *, const double);
static inline double ParameterFunction4 (const int, double *, const double);

static void   ParameterFarWing (const double, const double, ZpModelData *, ZpModelParameter *);
static void   ParameterWing (const double, const double, ZpModelData *, ZpModelParameter *);
static void   ParameterPeak (const double, const double, ZpModelData *, ZpModelParameter *);
static void   ParameterSymm (const double, const double, ZpModelData *, ZpModelParameter *);

static int    FractionalYield(const int, const int, const int, double, ZpModelParameter *, double *);
static double EvenOddTerm(const int, const int, double, double);

static void   SetParameterPeak(const int, const int, const double, ZpModelData *);
static void   SetParameterSymm(const int, const int, const double, ZpModelData *);
static void   SetParameterWing(const int, const int, const double, ZpModelData *);
static void   SetRegionBoundary(const int, const int, ZpModelData *);
static double ApMax(const int, const int, const double, const double);

static inline double ParameterFunction(const int, const int, double *, double);
static inline double ParameterInterpolate(const int, const int, double *, double *, const double);

static double region_low  = 10.0;
static double region_high = 20.0;
static int    region_idx  =    0; // 0: below region_low, 1: between low and high, 2: above high

static const int nparm = 5;


/**********************************************************/
/*      Wahl's mass chain yield                           */
/*      -----------------------                           */
/*      zf, af: for fissining nucleus Z and A             */
/*      ex: CN excitation energy, en < 0 for SF           */
/**********************************************************/
void WahlChainYield(const int zf, const int af, const double ex, double *fs, double *fd, double *fr)
{
  double sigma[7], delta[7], fract[7];

  /* mass center, Eq. (1b) */
  double am = af * 0.5;

  /* Gaussian widths */
  for(int k=0 ; k<7 ; k++) sigma[k] = SigmaSystematics_Wahl(k,zf,ex);

  /* mass splits */
  for(int k=0 ; k<7 ; k++) delta[k] = DeltaSystematics_Wahl(k,zf,am,ex);

  /* Gaussian intensities */
  for(int k=0 ; k<7 ; k++) fract[k] = FractSystematics_Wahl(k,zf,af,ex);

  /* normalize the total fract */
  double sy = fract[1] + fract[2] + fract[3] + fract[5] + fract[6];
  fract[0] = fract[4] = (200.0 - sy) * 0.5;

  fs[0] = sigma[4];  fd[0] = delta[4];  fr[0] = fract[4] * 0.01;
  fs[1] = sigma[3];  fd[1] = delta[3];  fr[1] = fract[3] * 0.01;
  fs[2] = sigma[6];  fd[2] = delta[6];  fr[2] = fract[6] * 0.01;
  fs[3] = sigma[2];  fd[3] = delta[2];  fr[3] = fract[2] * 0.01;
}


/**********************************************************/
/*      Total Number of Neutrons                          */
/**********************************************************/
double NuSystematics_Wahl(const int zf, double ex)
{
  double p[6], nu = 0.0;
  for(int i=0 ; i<6 ; i++) p[i] = 0.0;

  p[0] =  1.563;
  p[1] = 16.66;
  p[2] = -0.00804;
  p[3] =  0.0918;
  
  nu = ParameterFunction3(zf,p,ex);

  return nu;
}


/**********************************************************/
/*      Gausian Widths                                    */
/**********************************************************/
double SigmaSystematics_Wahl(const int k, const int zf, const double ex)
{
  double p[6], sigma = 0.0;
  for(int i=0 ; i<6 ; i++) p[i] = 0.0;

  switch(k){
  case 0:
  case 2:
  case 4:
    p[0] =  2.808;
    p[1] =  8.685;
    p[2] = -0.0454;
    p[3] =  0.372;
    p[4] = -0.620;
    p[5] = -0.0122;
    sigma = ParameterFunction3(zf,p,ex);
    if((k == 2) && (sigma < 8.6)) sigma = 8.6;
    break;
  case 1:
  case 3:
    p[0] =  2.45;
    sigma = ParameterFunction4(zf,p,ex);
    break;

  case 5:
  case 6:
    p[0] =  3.17;
    p[3] =  0.303;
    sigma = ParameterFunction4(zf,p,ex);
    break;
  default:
    break;
  }

  return sigma;
}


/**********************************************************/
/*      mass splits                                       */
/**********************************************************/
double DeltaSystematics_Wahl(const int k, const int zf, const double am, const double ex)
{
  double p[6], delta = 0.0;
  for(int i=0 ; i<6 ; i++) p[i] = 0.0;

  switch(k){
  case 2:
    delta = 0.0;
    break;
  case 0:
  case 4:
    p[0] =  25.34;
    p[1] =  18.55;
    p[2] =  -0.0402;
    p[3] =  -1.220;
    p[4] =  -1.732;
    delta = ParameterFunction3(zf,p,ex);
    if(k == 0) delta *= -1.0;
    break;
  case 1:
  case 3:
    p[0] = 136.66;
    p[1] =  -0.177;
    p[3] =   0.060;
    p[4] =  -0.038;
    delta = ParameterFunction4(zf,p,ex) - am;
    if(k == 1) delta *= -1.0;
    break;
  case 5:
  case 6:
    p[0] =  30.31;
    delta = ParameterFunction4(zf,p,ex);
    if(k == 5) delta *= -1.0;
    break;
  default:
    break;
  }

  return delta;
}


/**********************************************************/
/*      Gaussian Peak Values                              */
/**********************************************************/
double FractSystematics_Wahl(const int k, const int zf, const int af, const double ex)
{
  double p[6], fract = 0.0;

  for(int i=0 ; i<6 ; i++) p[i] = 0.0;

  switch(k){
  case 2:
    if(ex < 11.96) fract = 4.060 * exp(0.470*(ex - 11.96));
    else{
      double t = -0.030+ 0.0050*(af - 236);
      fract = 4.060 + 86.02 * (1.0 - exp(t*(ex - 11.96)));
    }
    break;
  case 1:
  case 3:
    p[0] =  43.00;
    p[1] =  -1.91;
    p[3] =  -3.41;
    fract = ParameterFunction4(zf,p,ex);
    if(fract < 0.0) fract = 0.0;
    break;
  case 5:
  case 6:
    if((zf < 93) || (ex > 20.0)) fract = 0.0;
    else{
      p[0] =   6.80;
      if(ex > 8.0) fract = p[0] - (p[0]/12.0) * (ex - 8.0);
      else         fract = ParameterFunction4(zf,p,ex);
      if(zf == 93) fract *= 0.5;
    }
    break;
  default:
    break;
  }

  return fract;
}


/**********************************************************/
/*      Wahl's functional form in Table 1, Eq.(3)and (4)  */
/**********************************************************/
inline double ParameterFunction3(const int zf, double *p, const double pe)
{
  double dz = zf - 92.0;
  double p1 = p[0] + p[3] * dz;
  double p2 = p[1] + p[4] * dz;
  double p3 = p[2] + p[5] * dz;

  return( p1 + (p2 - p1)*(1.0 - exp(p3 * pe)) );
}

inline double ParameterFunction4(const int zf, double *p, const double pe)
{
  double dz = zf - 92.0;
  double p1 = p[0] + p[3] * dz;
  double p2 = p[1] + p[4] * dz;

  return( p1 + p2 * pe );
}


/**********************************************************/
/*      Wahl's Zp model                                   */
/*      -----------------------                           */
/*      zf, af: for fissining nucleus Z and A             */
/*      ex: CN excitation energy, en < 0 for SF           */
/**********************************************************/
void WahlZpModel(const int zf, const int af, const int mass_first, const int nmass, const int ncharge, const double ex, double **fz, int *z0, double *fzn)
{
  ZpModelData      zd;
  ZpModelParameter zp;

  if(     ex < region_low ) region_idx = 0;
  else if(ex > region_high) region_idx = 2;
  else                      region_idx = 1;

  region_idx = 0;

  /* store parameters */
  SetParameterPeak(zf, af, ex, &zd);
  SetParameterSymm(zf, af, ex, &zd);
  SetParameterWing(zf, af, ex, &zd);

  /* store A-boundary */
  SetRegionBoundary(zf, af, &zd);

  for(int i=0 ; i<nmass ; i++){
    int    a  = i + mass_first;
    double ap = a;
    double ah = af - ap;

    if(ap < zd.B[0])      ParameterFarWing(ap,ah,&zd,&zp);
    else if(ap < zd.B[1]) ParameterWing(ap,ah,&zd,&zp);
    else if(ap < zd.B[2]) ParameterPeak(ap,ah,&zd,&zp);
    else if(ap < zd.B[3]) ParameterSymm(ap,ah,&zd,&zp);
    else if(ap < zd.B[4]) ParameterPeak(ap,ah,&zd,&zp);
    else if(ap < zd.B[5]) ParameterWing(ap,ah,&zd,&zp);
    else                  ParameterFarWing(ap,ah,&zd,&zp);

    /* modify FZ and FN values by inputs */
    if(zp.FZ > 1.0 && fzn[0] != 0.0) zp.FZ = 1.0 + (zp.FZ - 1.0)*fzn[0];
    if(zp.FN > 1.0 && fzn[1] != 0.0) zp.FN = 1.0 + (zp.FN - 1.0)*fzn[1];

    /* modify Symm region, connect DeltaZ linearly between B2 and B3 */
    // if( (zd.B[2] <= ap) && (ap < zd.B[3]) ){
    //   ZpModelParameter zp0, zp1;
    //   ParameterPeak(zd.B[2],af - zd.B[2],&zd,&zp0);
    //   ParameterPeak(zd.B[3],af - zd.B[3],&zd,&zp1);

    //   zp.DeltaZ = (zp1.DeltaZ - zp0.DeltaZ)/(zd.B[3] - zd.B[2]) * (ap - zd.B[2]) + zp0.DeltaZ;
    // }

    z0[i] = FractionalYield(zf,af,ncharge,ap,&zp,fz[i]);
  }
}


/**********************************************************/
/*      Region Boundaries                                 */
/**********************************************************/
void SetRegionBoundary(const int zf, const int af, ZpModelData *d)
{
  double am = ApMax(zf, af, d->DeltaZmax, d->d50);

  d->B[0] = 70.0;
  d->B[1] = 77.0 + 0.036*(af - 236);
  d->B[3] = (d->DeltaZmax - d->DeltaZ140 + am * d->d50 + 140.0 * d->dDeltaZ)/(d->d50 + d->dDeltaZ);
  d->B[2] = af - d->B[3];
  d->B[4] = af - d->B[1];
  d->B[5] = af - d->B[0];
  d->Ba   = af - am;
  d->Bb   = am;
}


double ApMax(const int zf, const int af, const double dzmax, const double d50)
{
  double f1 = (250.0 - af)/(250.0 - 236.0);
  if(f1 < 0.0) f1 = 0.0;
  else if(f1 > 1.0) f1 = 1.0;

  double f2 = 1.0 - f1;
  double r  = (double)af / (double)zf;
  double ak1 = 50.0*r - dzmax/d50;
  double ak2 = (50.0 - dzmax)*r;

  return(f1 * ak1 + f2 * ak2);
}


/**********************************************************/
/*      Sigma, Delta, Fz(A), Fn(A) in Each Region         */
/**********************************************************/
void ParameterFarWing(const double ap, const double ah, ZpModelData *d, ZpModelParameter *p)
{
  p->SigmaZ = d->SigmaZ140 + d->dSigmaZ * (d->B[4] - 140.0);
  p->DeltaZ = d->DeltaZ140 + d->dDeltaZ * (d->B[4] - 140.0); if(ap < ah) p->DeltaZ *= -1.0;
  p->FZ     = d->FZ140;
  p->FN     = d->FN140;
}


void ParameterWing(const double ap, const double ah, ZpModelData *d, ZpModelParameter *p)
{
  double deltaB5 = d->DeltaZ140 + d->dDeltaZ * (d->B[4] - 140.0); 
  double sigmaB5 = d->SigmaZ140 + d->dSigmaZ * (d->B[4] - 140.0); 

  if(ap < ah){
    p->DeltaZ = - deltaB5 - d->dDeltaZW * (d->B[4] - ah);
    p->SigmaZ = sigmaB5   + d->dSigmaZW * (d->B[1] - ap);
    p->FZ     = d->FZ140 + d->dFZW      * (d->B[1] - ap);
    p->FN     = d->FN140 + d->dFNW      * (d->B[1] - ap);
  }
  else{
    p->DeltaZ = deltaB5   - d->dDeltaZW * (ap - d->B[4]);
    p->SigmaZ = sigmaB5   + d->dSigmaZW * (ap - d->B[4]);
    p->FZ     = d->FZ140 + d->dFZW      * (ap - d->B[4]);
    p->FN     = d->FN140 + d->dFNW      * (ap - d->B[4]);
  }
}


void ParameterPeak(const double ap, const double ah, ZpModelData *d, ZpModelParameter *p)
{
  double a = (ap < ah) ? ah : ap;

  p->DeltaZ = d->DeltaZ140 + d->dDeltaZ * (a - 140.0); if(ap < ah) p->DeltaZ *= -1.0;
  p->SigmaZ = d->SigmaZ140 + d->dSigmaZ * (a - 140.0);

  p->FZ = d->FZ140 + d->dFZ * (a - 140.0);
  p->FN = d->FN140;
}


void ParameterSymm(const double ap, const double ah, ZpModelData *d, ZpModelParameter *p)
{
  if(ap < d->Ba){
    p->DeltaZ = - d->DeltaZmax + d->d50*(ah - d->Bb);
    p->SigmaZ = d->SigmaZ50;
  }
  else if(ap < d->Bb){
    p->DeltaZ = - d->DeltaZmax +  2.0*d->DeltaZmax * (d->Bb - ah) / (d->Bb - d->Ba);
    p->SigmaZ = d->SigmaZ140 - d->dSigmaZ * (140.0 - d->Bb);
  }
  else{
    p->DeltaZ = d->DeltaZmax - d->d50*(ap - d->Bb);
    p->SigmaZ = d->SigmaZ50;
  }

  p->FZ = 1.0;
  p->FN = 1.0;
}


/**********************************************************/
/*      Fractional Independent Yields, FI                 */
/**********************************************************/
int FractionalYield(const int zf, const int af, const int ncharge, double ap, ZpModelParameter *p, double *f)
{
  double x0 = 1.0/(p->SigmaZ * sqrt(2.0));
  double zp = ap * (double)zf/(double)af + p->DeltaZ;

  /* calculate fraction from Zp-charge_range to Zp+charge_range */
  int charge_range = (ncharge - 1)/2;
  int z0 = (int)zp - charge_range;

  double sum = 0.0;
  for(int j=0 ; j<ncharge ; j++){
    double z = (double)(z0 + j);
    double v = (z - zp + 0.5) * x0;
    double w = (z - zp - 0.5) * x0;
    double fa = EvenOddTerm(z,ap,p->FZ,p->FN);

    f[j] = 0.5 * fa * (erf(v) - erf(w));
    sum += f[j];
  }

  /* renormalize fraction */
  if(sum > 0.0){
    for(int j=0 ; j< ncharge ; j++) f[j] /= sum;
  }

  return z0;
}


double EvenOddTerm(const int zf, const int af, double fz, double fn)
{
  double f = 1.0;

  int iz = zf%2;
  int in = (af - zf)%2;

  if     ((iz == 0) && (in == 0)) f = fz * fn;
  else if((iz == 0) && (in == 1)) f = fz / fn;
  else if((iz == 1) && (in == 0)) f = fn / fz;
  else                            f = 1.0/(fz * fn);

  return f;
}


/**********************************************************/
/*      Store Model Parameters in Each Region             */
/**********************************************************/
void SetParameterPeak(const int zf, const int af, const double ex, ZpModelData *d)
{
  static double val1[] = {
      0.566 ,  0.0   ,  0.0064,  0.0109,  0.0    ,
     -0.487 ,  0.0   ,  0.0180,  0.0   , -0.00203,
      1.207 ,  0.0   , -0.0420,  0.0   ,  0.0022 ,
      1.076 ,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
     -0.0038,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
     -0.0080,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
      0.0030,  0.0   ,  0.0   ,  0.0   ,  0.0    };

  static double val2[] = {
      0.542 ,  1.310 ,  0.033 ,  0.0   , -0.005  ,
     -0.428 ,  0.0   ,  0.0   ,  0.164 , -0.0116 ,
      1.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
      1.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
      0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
      0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
      0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    };

  int i = 0;
  if(region_idx != 1){
    double * val = (region_idx == 0) ? &val1[0] : &val2[0];
    d->SigmaZ140 = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->DeltaZ140 = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->FZ140     = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->FN140     = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->dSigmaZ   = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->dDeltaZ   = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->dFZ       = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
  }
  else{
    int i = 0;
    d->SigmaZ140 = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->DeltaZ140 = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->FZ140     = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->FN140     = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->dSigmaZ   = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->dDeltaZ   = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->dFZ       = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
  }
}


void SetParameterSymm(const int zf, const int af, const double ex, ZpModelData *d)
{
  static double val1[] = {
      0.191 ,  0.0   , -0.0076,  0.0   ,  0.0    ,
      0.356 ,  0.060 ,  0.0   ,  0.0   ,  0.0    ,
      0.699 ,  0.0   ,  0.0   ,  0.0   ,  0.0    };

  static double val2[] = {
      0.191 ,  0.0   , -0.0076,  0.0   ,  0.0    ,
      0.542 ,  1.310 ,  0.033 ,  0.0   , -0.005  ,
      0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    };


  int i = 0;
  if(region_idx != 1){
    double * val = (region_idx == 0) ? &val1[0] : &val2[0];
    d->d50       = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->SigmaZ50  = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->DeltaZmax = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
  }
  else{
    d->d50       = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->SigmaZ50  = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->DeltaZmax = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
  }

  if(d->d50 < 0.0) d->d50 = 0.0;
}


void SetParameterWing(const int zf, const int af, const double ex, ZpModelData *d)
{
  static double val1[] = {
     -0.045 ,  0.0094,  0.0   ,  0.0   ,  0.0    ,
      0.0   , -0.0045,  0.0   ,  0.0   ,  0.0    ,
      0.159 , -0.028 ,  0.0   ,  0.0   ,  0.0    ,
      0.039 ,  0.0   ,  0.0   ,  0.0   ,  0.0    };

  static double val2[] = {
      0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
      0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
      0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
      0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    };


  int i = 0;
  if(region_idx != 1){
    double * val = (region_idx == 0) ? &val1[0] : &val2[0];
    d->dSigmaZW  = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->dDeltaZW  = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->dFZW      = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
    d->dFNW      = ParameterFunction(zf,af,&val[nparm*i],ex); i++;
  }
  else{
    d->dSigmaZW  = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->dDeltaZW  = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->dFZW      = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
    d->dFNW      = ParameterInterpolate(zf,af,&val1[nparm*i],&val2[nparm*i],ex); i++;
  }
}


/**********************************************************/
/*      Parameter Functional Form                         */
/**********************************************************/
inline double ParameterFunction(const int zf, const int af, double *p, double pe)
{
  double dz = zf - 92.0;

  double x  = 0.0;
  if(region_idx == 0){ 
    double da = af - 236.0;
    double de = pe - 6.551;
    x = p[0] + p[1]*dz + p[2]*da + p[3]*de + p[4]*da*da;
  }
  else{
    double p1 = p[0] + p[2]*dz;
    double p2 = p[1] + p[3]*dz;
    x = p1 + (p2 - p1)*(1.0 - exp(-p[4]*pe));
  }

  return x;
}


/**********************************************************/
/*      Parameter Functional Form in Transitional Area    */
/**********************************************************/
inline double ParameterInterpolate(const int zf, const int af, double *p0, double *p1, const double pe)
{
  const double a = 1.0;

  region_idx = 0; double x0 = ParameterFunction(zf,af,p0,pe);
  region_idx = 2; double x1 = ParameterFunction(zf,af,p1,pe);
  region_idx = 1;

  double rm = 0.5 * (region_high + region_low);
  double f = 1/(1.0+exp(-(pe-rm)/a));
  double x  = x0 * f + x1 * (1.0 - f);

  return x;
}

#include <iostream>
#include <cmath>

#include "ffp.h"

static void   Gaussian_U235 (const double, double *, double *, double *);
static void   Gaussian_U238 (const double, double *, double *, double *);
static void   Gaussian_Pu239 (const double, double *, double *, double *);
static void   Gaussian_Cf252 (const double, double *, double *, double *);

static double TKEenergy_U235 (const double, const int, const double);
static double TKEenergy_U238 (const double, const int, const double);
static double TKEenergy_Pu239 (const double, const int, const double);
static double TKEenergy_Cf252 ();

static double TKEmass_U235 (int, double);
static double TKEmass_U238 (int, double);
static double TKEmass_Pu239 (int, double);
static double TKEmass_Cf252 (const int);

static double TKEwidth_U235 (const bool);
static double TKEwidth_U238 (const bool);
static double TKEwidth_Pu239 (const bool, const int);
static double TKEwidth_Cf252 (const bool);

static void   YAmodel_CGMF (const double, const double, double *, double *, double *, double *);
static double TKE_Systematics (const int, const int);


static double TKEsave = 0.0;

/**********************************************************/
/*      Gaussian Fitted Y(A)                              */
/**********************************************************/
void FFPGaussianParameters(const bool sf, const int zf, const int af, const double ex, const double sn, double *fs, double *fd, double *fr)
{
  for(int i=0 ; i<4 ; i++) fr[i] = fd[i] = fs[i] = 0.0;

  double en = ex - sn;
  if(en < 0.0) en = 0.0;

  /* spontaneous fission */
  if(sf){
    /* Cf252 */
    if((zf == 98) && (af == 252)) Gaussian_Cf252(0.0,fs,fd,fr);
  }

  /* neutron induced fission */
  else{
    /* U235 */
    if(zf == 92){
      switch(af){
      case 230:
      case 231:
      case 232:
      case 233:
      case 234:
      case 235:
      case 236: Gaussian_U235(en,fs,fd,fr); break;
      case 237:
      case 238:
      case 239: Gaussian_U238(en,fs,fd,fr); break;
      default: break;
      }
    }

    /* Pu239 */
    else if(zf == 94){
      switch(af){
      case 236:
      case 237:
      case 238:
      case 239:
      case 240:
      case 241:
      case 242:
      case 243: Gaussian_Pu239(en,fs,fd,fr); break;
      default: break;
      }
    }
  }
}


/**********************************************************/
/*      TKE as a function of A from input parameters      */
/**********************************************************/
double FFPConstructTKEA(const int af, const int a, double *tkea)
{
  double tke = 0.0;
  double am = 0.5 * af;
  double ah = a;

  double da = ah - am;
  if(da < 0.1) ah = af - a;

  tke = (tkea[0] - tkea[1]*ah) * (1.0 - tkea[2]*exp(-(da*da)/tkea[3]));

  return tke;
}


/**********************************************************/
/*     sigma_TKE as a function of A from input parameters */
/**********************************************************/
double FFPConstructSigmaTKEA(const int af, const int a, double *stkea)
{
  double stke = 0.0;
  double am = 0.5 * af;
  double ah = a;

  if (ah - am < 0.1) ah = af - a;

  //double A0 = stkea[0];

  double x = ah - am;

  stke = stkea[0] + stkea[1]*exp(-stkea[2]*x*x);

  //stke = stkea[1];
  //double x = ah - A0;
  //for (int i=2; i<=9; i++){
  //  stke += stkea[i] * x;
  //  x *= ah - A0;
  //}
  return stke;
}


/**********************************************************/
/*      TKE as a function of pre-A                        */
/**********************************************************/
double FFPSystematics_TKE_A(const bool sf, const int zf, const int af, const int a, const double t0)
{
  double tke = 0.0;
  double am  = af * 0.5;
  double ah  = ((double)a >= am) ? (double)a : (double)(af - a);

  if(sf){
    if((zf == 98) && (af == 252)){
      tke = TKEmass_Cf252(ah);
    }
  }
  else{
    if(zf == 92){
      if((230 <= af) && (af <= 236)) tke = TKEmass_U235(ah,am);
      else if((237 <= af) && (af <= 239)) tke = TKEmass_U238(ah,am);
    }
    else if(zf == 94){
      if((236 <= af) && (af <= 243)) tke = TKEmass_Pu239(ah,am);
    }
  }

  /* if not given in systematics, set the average <TKE> */
  if(tke == 0.0){
    tke = (t0 > 0.0) ? t0 : TKE_Systematics(zf,af);
  }

  TKEsave = tke;
  return tke;
}


/**********************************************************/
/*      TKE distribution width                            */
/**********************************************************/
double FFPSystematics_TKE_A_Width(const bool sf, const int zf, const int af, const int a)
{
  bool constwidth = false;
  double dtke = 0.0;
  double am  = af * 0.5;
  double ah  = ((double)a >= am) ? (double)a : (double)(af - a);

  if(sf){
    if((zf == 98) && (af == 252)){
      dtke = TKEwidth_Cf252(constwidth);
    }
  }
  else{
    if(zf == 92){
      if((230 <= af) && (af <= 236)) dtke = TKEwidth_U235(constwidth);
      else if((237 <= af) && (af <= 239)) dtke = TKEwidth_U238(constwidth);
    }
    else if(zf == 94){
      if((236 <= af) && (af <= 243)) dtke = TKEwidth_Pu239(constwidth,ah);
    }
  }

  return dtke;
}


/**********************************************************/
/*      TKE as a function of neutron incident energy      */
/**********************************************************/
double FFPSystematics_TKE_En(const bool sf, const int zf, const int af, const double ex, const double sn)
{
  double tke = 0.0;

  double en = ex - sn;
  if(en < 0.0) en = 0.0;

  double tke0 = TKE_Systematics(zf,af);

  if(sf){
    if((zf == 98) && (af == 252)) tke = TKEenergy_Cf252();
    else tke = tke0;
  }
  else{
    if(zf == 92){
      if((233 <= af) && (af <= 236)) tke = TKEenergy_U235(en,af,tke0);
      else if((237 <= af) && (af <= 239)) tke = TKEenergy_U238(en,af,tke0);
    }
    else if(zf == 94){
      if((236 <= af) && (af <= 243)) tke = TKEenergy_Pu239(en,af,tke0);
    }
  }

  return tke;
}


/**********************************************************/
/*       U235                                             */
/**********************************************************/
/* Shin Okumura, JNST 55, 1009 (2018)
   En < 10 keV, use the thermal yield
   fast range, fitted to D'yachenko data (1969) by P. Jaffke, updated by TK */
void Gaussian_U235(const double en, double *fs, double *fd, double *fr)
{
/* original */
  fr[0] = 0.78464   ; fd[0] = 22.997  ; fs[0] = 4.8278;
  fr[1] = 0.20323   ; fd[1] = 15.634  ; fs[1] = 2.7277;
  fr[2] = 0.0       ; fd[2] =  0.0    ; fs[2] = 0.0;
  fr[3] = 0.002916  ; fd[3] =  0.0    ; fs[3] = 8.6;

/* adjusted, 2020/5/10
  fr[0] = 8.1929E-01; fd[0] = 2.3116E+01; fs[0] = 5.0795E+00;
  fr[1] = 1.9853E-01; fd[1] = 1.5254E+01; fs[1] = 2.9359E+00;
  fr[2] = 0.00000   ; fd[2] =  0.000    ; fs[2] = 0.0000;
  fr[3] = 0.00291   ; fd[3] =  0.000    ; fs[3] = 8.6000;
*/
  if(en >= 0.01){
    fd[0] = fd[0] - 0.259156 * en;
    fd[1] = fd[1] - 0.1996   * en;
    fs[0] = fs[0] + 0.166617 * en;

    fr[0] = fr[0] + 0.0161376 * en;
    fr[1] = fr[1] - 0.0156474 * en;
    fr[3] = fr[3] + 0.004 * en;
    /* re-adjust F0 */
    fr[0] = 1.0 - fr[1] - 0.5 * fr[3];
  }
}


/* determined from experimental data of D. Duke (2014) */
double TKEenergy_U235(const double en, const int af, const double tke0)
{
  double tke = 0.0;
  if(af == 236){
//    double p1 = 171.239;
//    double p2 = 0.18;
//    double p3 = 0.0043;
//    double p4 = 0.32296;
//    tke = (p1 - p2 * en) * (1.0 - p3 * exp(-en / p4));

    if(en < 2.2778) tke = 1.7009e+02;
    else            tke = 1.7050e+02 -0.18 * en;

  }
  else{
    /*** rough estimate of TKE slope from 235U data */
    double tke1 = -0.18 - (236 - af) * 0.02;
    tke = tke0 + tke1 * en;
  }

  return tke;
}


/* U235 determined from Hambsch experimental data
   although A=236, we shift this functional form to the other A numbers */
double TKEmass_U235(int ah, double am)
{
  double p1 = 335.284;
  double p2 = 1.17404;
  double p3 = 0.187596;
  double p4 = 69.0834;

  double da = ah - am;
  if(da < 0.1) ah = am + 1.0;
  double tke = (p1 - p2 * ah) * (1.0 - p3 * exp(-da * da / p4));

  return tke;
}


double TKEwidth_U235(const bool cw)
{
  double dtke = 0.0;
  if(cw) dtke = 8.0;
  else   dtke = TKEsave * 0.04;

  return dtke;
}


/**********************************************************/
/*      U238                          from Okumura (2022) */
/**********************************************************/
/* data from Vives 2000 */
void Gaussian_U238(const double en, double *fs, double *fd, double *fr)
{
  static double p[] = {-2.167787, 5.0323, 135.16,  -0.09, 3.3868, 0.0142, -2.224051, -5.1629, 142.20,  -0.16, 5.5624, 0.1048, 10.0092, 0.0153};
  YAmodel_CGMF(239.0,en,p,fs,fd,fr);
}


double TKEenergy_U238(const double en, const int af, const double tke0)
{
  double tke = 0.0;
  if(af == 239){
    double p1 = 171.072;
    double p2 = -0.311464;
    tke = p1 + p2 * en;
  }
  else{
    double tke1 = -0.320;
    tke = tke0 + tke1 * en;
  }

  return tke;
}


double TKEmass_U238(int ah, double am)
{
  double p1 = 348.371;
  double p2 = 1.274;
  double p3 = 0.1800;
  double p4 = 59.199;

  double da = ah - am;
  if(da < 0.1) ah = am + 1.0;
  double tke = (p1 - p2 * ah) * (1.0 - p3 * exp(-da * da / p4));

  return tke;
}


double TKEwidth_U238(const bool cw)
{
  double dtke = 0.0;
  if(cw) dtke = 8.0;
  else   dtke = TKEsave * 0.04;

  return dtke;
}



/**********************************************************/
/*      Pu239                                             */
/**********************************************************/
/* data from Akimov 1971, Surin 1972, Wagemans 1984, Schillebeeckx 1992, 
   Nishio 1995, Tsuchiya 2000 */
void Gaussian_Pu239(const double en, double *fs, double *fd, double *fr)
{
  static double p[] ={-25.369127, 29.9818, 135.11,   0.13, 3.8465, 0.0689, -25.258746, -30.0000, 141.35,   0.20, 6.5176, 0.0324, 9.9823, 0.0580};
  YAmodel_CGMF(240.0,en,p,fs,fd,fr);
}


double TKEenergy_Pu239(const double en, const int af, const double tke0)
{
  double tke = 0.0;
  if(af == 240){
    double p1 = 178.21;
    double p2 = -0.3409;
    tke = p1 + p2 * en;
  }
  else{
    tke = tke0 - 0.3409 * en;
  }

  return tke;
}


/* parameter by A. Lovell */
double TKEmass_Pu239(int ah, double am)
{
  double p1 = 341.921;
  double p2 = 1.17193;
  double p3 = 0.16595;
  double p4 = 46.8726;
  
  double da = ah - am;
  if(da < 0.1) ah = am + 1.0;
  double tke = (p1 - p2 * ah) * (1.0 - p3 * exp(-da * da / p4));

  return tke;
}


double TKEwidth_Pu239(const bool cw, const int ah)
{
  static double a0 = 128.0;
  static double p[] = {  7.5837e+00,  1.0168e-01, -1.6588e-02,  3.9178e-04,  0.0000e+00,
                         0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00};
  double dtke = 0.0;

  if(cw) dtke = 8.0;
  else{
    dtke = p[0];
    double x = ah - a0;
    for (int i=1; i<=8; i++){
      dtke += p[i] * x;
      x *= ah - a0;
    }
  }

  return dtke;
}


/**********************************************************/
/*      Cf252                                             */
/**********************************************************/
/* Y(A) from Hambsch 1997 */
void Gaussian_Cf252(double en, double *fs, double *fd, double *fr)
{
  static double p[] = {-1.000000, -1.6706, 141.78,   0.00, 5.7679, 0.0000, -1.000000, 1.6706, 146.26,   0.00, 7.9062, 0.0000, 10.0496, 0.0000};

  YAmodel_CGMF(252.0,en,p,fs,fd,fr);
}


/* A.B. Smith Phys. Rev. 102, 813 (1956) */
double TKEenergy_Cf252()
{
  double tke = 186.0;
  return tke;
}


/* TKE determined from Gook 2014 data by P. Jaffke */
double TKEmass_Cf252(const int ah)
{
  static double a0 = 132.0;
  static double p[] = {  1.9084e+02, -1.6273e-01, -9.2138e-02,  1.2956e-02, -1.0775e-03,
                         4.6069e-05, -1.0307e-06,  1.1537e-08, -5.0992e-11};
  double tke = p[0];

  double x = ah - a0;
  for(int i=1 ; i<=8 ; i++){
    tke += p[i] * x;
    x *= ah - a0;
  }

  return tke;
}


double TKEwidth_Cf252(bool cw)
{
  double dtke = 0.0;
  if(cw) dtke = 8.0;
  else   dtke = TKEsave * 0.04;
  return dtke;
}


/**********************************************************/
/*      CGMF Parameterization of Y(A)                     */
/*      P. Jaffke's parameter from CGMF and               */
/*      updated by A. Lovell                              */
/**********************************************************/
/* G(A,E) = w(E)/sqrt(2*PI*sigma(E)^2) * (exp [-(A-mu(E))^2/(2*sigma(E)^2)])
   dependence on the incident neutron energy (E):
   w(E) = 1/(1+exp[(E-w_a)/w_b])
   mu(E) = mu_a + mu_b*E
   sigma(E) = sigma_a + sigma_b*E */
void YAmodel_CGMF(const double af, const double en, double *p, double *fs, double *fd, double *fr)
{
  fr[0] = 1.0/(1.0 + exp((en - p[0]) / p[1]));
  fr[1] = 1.0/(1.0 + exp((en - p[6]) / p[7]));
  fr[2] = 0.0;
  fr[3] = 2.0 - 2.0*(fr[0] + fr[1]);
  
  fd[0] = p[ 2] + p[ 3] * en - af/2.0;
  fd[1] = p[ 8] + p[ 9] * en - af/2.0;
  fd[2] = 0.0;
  fd[3] = 0.0;

  fs[0] = p[ 4] + p[ 5] * en;
  fs[1] = p[10] + p[11] * en;
  fs[2] = 0.0;
  fs[3] = p[12] + p[13] * en;
}


/**********************************************************/
/*      TKE systematics by Viola                          */
/**********************************************************/
double TKE_Systematics(const int zf, const int af)
{
  /*** V.E. Viola, K. Kwiatkowski, M. Walker,  Phys. Rev. C 31, 1550 (1985) */
  double a13 = pow((double)af,1/3.0);
  double tke = 0.1189 * zf * zf / a13 + 7.3;
  return tke;
}



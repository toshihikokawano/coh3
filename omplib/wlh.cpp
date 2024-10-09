/***********************************************************/
/*      T.R. Whitehead, Y. Lim, and J.W. Hold              */
/*      arXiv:2009.08436v3 [nucl-th] 27 Aug 2021           */
/*      OMPs fit from ab initio chiral interactions        */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPWLH1(double e, int at, int zt, int ai, int zi, Optical *omp)
{

   double eps=1.0-2.0*(double)zt/(double)at;
   double a3 = pow((double)at,1.0/3.0);
   double omp_parms[48] = { 0., 0., 5.71581e+01,  3.28055e-01,  1.93092e-04,  1.27824e-06,
         2.09515e+01,  2.67401e-01,  7.47553e-04,  1.25258e+00,
         3.21132e-01,  6.70967e-04,  4.24950e-06,  7.25804e-01,
         4.16904e-04,  3.86419e-06,  2.14543e-01,  7.58107e-01,
         3.25537e+00,  3.12599e-01,  6.25048e-04,  9.99957e+00,
         7.54228e-02,  5.59196e-01,  6.76906e+01,  6.89906e-01,
         7.64380e+01,  9.33670e-01,  2.57898e-06,  2.19943e-01,
         5.94204e-01,  6.05321e+00,  1.21647e-02,  4.84285e-04,
         2.18277e+00,  4.58225e-02,  3.02741e+00,  1.40682e-01,
         1.26379e+00,  1.48432e+00,  4.20908e-03,  8.20927e-01,
         1.00135e+01,  8.54929e-03,  1.22918e+00,  8.17408e-01,
         7.74917e-01,  3.17690e-04};

   if(ai != 1 || zi != 0) return(0);

   // real volume depth, radius, diffuseness
   omp->v1 = omp_parms[2] - omp_parms[6]*eps;
   omp->v2 = -omp_parms[3] + omp_parms[7]*eps;
   omp->v3 = omp_parms[4] - omp_parms[8]*eps;
   omp->v4 = omp_parms[5];

   omp->vs1 = 0.0;
   omp->vs2 = 0.0;
   omp->vs3 = 0.0;

   omp->r0 = omp_parms[9] - omp_parms[11]*e +omp_parms[12]*e*e - omp_parms[10]/a3;
   omp->r0s = 0.0;
   omp->a0 = omp_parms[13] + omp_parms[14]*e - omp_parms[15]*e*e - (omp_parms[16]-omp_parms[17]*eps)*eps;
   omp->a0s = 0.0;

   // imaginary volume depth, radius, diffuseness
   omp->wv1 = omp_parms[18] - omp_parms[21]*eps;
   omp->wv2 = omp_parms[19] - omp_parms[22]*eps;
   omp->wv3 = -omp_parms[20];

   omp->rv = omp_parms[23] + (omp_parms[24] + omp_parms[25]*(double)at)/(omp_parms[26] + (double)at + omp_parms[27]*e) + omp_parms[28]*e*e;
   omp->av = omp_parms[29] - omp_parms[30]*e/(-omp_parms[31] - e) + (omp_parms[32] - omp_parms[33]*e)*eps;

   // imaginary surface depth, radius, diffuseness
   omp->ws1 = omp_parms[34] - omp_parms[36]*eps;
   omp->ws2 = -omp_parms[35] + omp_parms[37]*eps;
   omp->ws3 = 0.0;

   if (e > 40.) {
      omp->ws1 = 0.0;
      omp->ws2 = 0.0;
   }

   omp->rs = omp_parms[38] - omp_parms[40]*e - omp_parms[39]/a3;
   omp->as = omp_parms[41];

   // spin orbit parameters
   omp->vso1 = omp_parms[42] - omp_parms[43]*(double)at;
   omp->vso2 = 0.0;
   omp->vso3 = 0.0;
 
   omp->wso1 = 0.0;
   omp->wso2 = 0.0;
   omp->wso3 = 0.0;

   omp->rvso = omp_parms[44] - omp_parms[45]/a3;
   omp->rwso = 0.0;
   omp->avso = omp_parms[46] - omp_parms[47]*(double)at;
   omp->awso = 0.0;

   return(0x0013);

}

unsigned int OMPWLH2(double e, int at, int zt, int ai, int zi, Optical *omp)
{

   double eps=1.0-2.0*(double)zt/(double)at;
   double a3 = pow((double)at,1.0/3.0);
   double omp_parms[48] = { 0., 0., 5.40220e+01,  2.99006e-01,  1.41452e-04,  1.66736e-06,
         2.08434e+01,  2.65492e-01,  6.62695e-04,  1.29798e+00,
         3.97229e-01,  5.40576e-04,  1.97630e-06,  7.21033e-01,
         3.94364e-04,  1.74962e-06,  1.70648e-01,  6.24857e-01,
         2.21505e+00,  2.61652e-01,  5.81812e-04,  7.60384e+00,
         4.76306e-02,  4.65605e-01,  8.85857e+01,  8.17832e-01,
         8.52889e+01,  8.32774e-01,  2.55027e-06,  2.04605e-01,
         6.30732e-01,  8.33116e+00,  7.10580e-02,  8.42787e-04,
         1.87864e+00,  4.04378e-02,  3.29653e+00,  1.53675e-01,
         1.30567e+00,  1.47385e+00,  3.91707e-03,  8.33677e-01,
         9.39118e+00,  8.02939e-03,  1.25967e+00,  8.26821e-01,
         7.80790e-01,  3.74290e-04};

   if(ai != 1 || zi != 0) return(0);

   // real volume depth, radius, diffuseness
   omp->v1 = omp_parms[2] - omp_parms[6]*eps;
   omp->v2 = -omp_parms[3] + omp_parms[7]*eps;
   omp->v3 = omp_parms[4] - omp_parms[8]*eps;
   omp->v4 = omp_parms[5];

   omp->vs1 = 0.0;
   omp->vs2 = 0.0;
   omp->vs3 = 0.0;

   omp->r0 = omp_parms[9] - omp_parms[11]*e +omp_parms[12]*e*e - omp_parms[10]/a3;
   omp->r0s = 0.0;
   omp->a0 = omp_parms[13] + omp_parms[14]*e - omp_parms[15]*e*e - (omp_parms[16]-omp_parms[17]*eps)*eps;
   omp->a0s = 0.0;

   // imaginary volume depth, radius, diffuseness
   omp->wv1 = omp_parms[18] - omp_parms[21]*eps;
   omp->wv2 = omp_parms[19] - omp_parms[22]*eps;
   omp->wv3 = -omp_parms[20];

   omp->rv = omp_parms[23] + (omp_parms[24] + omp_parms[25]*(double)at)/(omp_parms[26] + (double)at + omp_parms[27]*e) + omp_parms[28]*e*e;
   omp->av = omp_parms[29] - omp_parms[30]*e/(-omp_parms[31] - e) + (omp_parms[32] - omp_parms[33]*e)*eps;

   // imaginary surface depth, radius, diffuseness
   omp->ws1 = omp_parms[34] - omp_parms[36]*eps;
   omp->ws2 = -omp_parms[35] + omp_parms[37]*eps;
   omp->ws3 = 0.0;

   if (e > 40.) {
      omp->ws1 = 0.0;
      omp->ws2 = 0.0;
   }

   omp->rs = omp_parms[38] - omp_parms[40]*e - omp_parms[39]/a3;
   omp->as = omp_parms[41];

   // spin orbit parameters
   omp->vso1 = omp_parms[42] - omp_parms[43]*(double)at;
   omp->vso2 = 0.0;
   omp->vso3 = 0.0;
 
   omp->wso1 = 0.0;
   omp->wso2 = 0.0;
   omp->wso3 = 0.0;

   omp->rvso = omp_parms[44] - omp_parms[45]/a3;
   omp->rwso = 0.0;
   omp->avso = omp_parms[46] - omp_parms[47]*(double)at;
   omp->awso = 0.0;

   return(0x0013);

}

unsigned int OMPWLH3(double e, int at, int zt, int ai, int zi, Optical *omp)
{

   double eps=1.0-2.0*(double)zt/(double)at;
   double a3 = pow((double)at,1.0/3.0);
   double omp_parms[48] = { 0., 0., 5.22564e+01,  2.62373e-01,  2.05452e-04,  4.39954e-07,
         2.25296e+01,  2.38283e-01,  6.31830e-04,  1.24546e+00,
         2.33753e-01,  6.55036e-04,  5.23030e-06,  7.14386e-01,
         5.33923e-04,  5.24113e-06,  2.98017e-01,  8.17472e-01,
         1.82420e+00,  1.51790e-01,  3.81671e-05,  7.01318e+00,
        -1.01862e-03,  6.08065e-01,  6.73344e+01,  6.28708e-01,
         6.80065e+01,  9.60325e-01,  2.99417e-06,  1.89499e-01,
         6.10633e-01,  1.05308e+01,  8.35152e-02,  7.51665e-04,
         1.03690e+00,  5.08780e-03,  8.87440e-01,  3.52080e-02,
         1.26186e+00,  1.19301e+00,  3.58826e-03,  8.33927e-01,
         1.00957e+01,  8.49747e-03,  1.23265e+00,  7.92527e-01,
         7.69466e-01,  3.71410e-04};

   if(ai != 1 || zi != 0) return(0);

   // real volume depth, radius, diffuseness
   omp->v1 = omp_parms[2] - omp_parms[6]*eps;
   omp->v2 = -omp_parms[3] + omp_parms[7]*eps;
   omp->v3 = omp_parms[4] - omp_parms[8]*eps;
   omp->v4 = omp_parms[5];

   omp->vs1 = 0.0;
   omp->vs2 = 0.0;
   omp->vs3 = 0.0;

   omp->r0 = omp_parms[9] - omp_parms[11]*e +omp_parms[12]*e*e - omp_parms[10]/a3;
   omp->r0s = 0.0;
   omp->a0 = omp_parms[13] + omp_parms[14]*e - omp_parms[15]*e*e - (omp_parms[16]-omp_parms[17]*eps)*eps;
   omp->a0s = 0.0;

   // imaginary volume depth, radius, diffuseness
   omp->wv1 = omp_parms[18] - omp_parms[21]*eps;
   omp->wv2 = omp_parms[19] - omp_parms[22]*eps;
   omp->wv3 = -omp_parms[20];

   omp->rv = omp_parms[23] + (omp_parms[24] + omp_parms[25]*(double)at)/(omp_parms[26] + (double)at + omp_parms[27]*e) + omp_parms[28]*e*e;
   omp->av = omp_parms[29] - omp_parms[30]*e/(-omp_parms[31] - e) + (omp_parms[32] - omp_parms[33]*e)*eps;

   // imaginary surface depth, radius, diffuseness
   omp->ws1 = omp_parms[34] - omp_parms[36]*eps;
   omp->ws2 = -omp_parms[35] + omp_parms[37]*eps;
   omp->ws3 = 0.0;

   if (e > 40.) {
      omp->ws1 = 0.0;
      omp->ws2 = 0.0;
   }

   omp->rs = omp_parms[38] - omp_parms[40]*e - omp_parms[39]/a3;
   omp->as = omp_parms[41];

   // spin orbit parameters
   omp->vso1 = omp_parms[42] - omp_parms[43]*(double)at;
   omp->vso2 = 0.0;
   omp->vso3 = 0.0;
 
   omp->wso1 = 0.0;
   omp->wso2 = 0.0;
   omp->wso3 = 0.0;

   omp->rvso = omp_parms[44] - omp_parms[45]/a3;
   omp->rwso = 0.0;
   omp->avso = omp_parms[46] - omp_parms[47]*(double)at;
   omp->awso = 0.0;

   return(0x0013);

}

unsigned int OMPWLH4(double e, int at, int zt, int ai, int zi, Optical *omp)
{

   double eps=1.0-2.0*(double)zt/(double)at;
   double a3 = pow((double)at,1.0/3.0);
   double omp_parms[48] = { 0., 0., 4.89638e+01,  2.84317e-01, -8.84999e-04,  5.26989e-06,
         1.99247e+01,  3.23092e-01,  1.62662e-03,  1.35913e+00,
         3.07903e-01,  1.31969e-03,  1.57760e-05,  7.71034e-01,
         1.50701e-03,  2.03980e-05,  3.50893e-01,  1.09001e+00,
         3.27876e+00,  2.45623e-01,  7.36418e-04,  8.63975e+00,
        -2.67419e-04,  5.00823e-01,  1.10825e+02,  8.51822e-01,
         1.13970e+02,  7.69755e-01,  4.16265e-06,  2.65264e-01,
         5.66004e-01,  4.64734e+00,  4.74981e-02,  1.09281e-04,
         1.79793e+00,  6.41075e-02,  2.95102e+00,  2.26054e-01,
         1.34526e+00,  1.70187e+00,  2.02361e-03,  8.84760e-01,
         8.86358e+00,  8.46551e-03,  1.34759e+00,  9.75889e-01,
         8.54430e-01,  3.52910e-04};

   if(ai != 1 || zi != 0) return(0);

   // real volume depth, radius, diffuseness
   omp->v1 = omp_parms[2] - omp_parms[6]*eps;
   omp->v2 = -omp_parms[3] + omp_parms[7]*eps;
   omp->v3 = omp_parms[4] - omp_parms[8]*eps;
   omp->v4 = omp_parms[5];

   omp->vs1 = 0.0;
   omp->vs2 = 0.0;
   omp->vs3 = 0.0;

   omp->r0 = omp_parms[9] - omp_parms[11]*e +omp_parms[12]*e*e - omp_parms[10]/a3;
   omp->r0s = 0.0;
   omp->a0 = omp_parms[13] + omp_parms[14]*e - omp_parms[15]*e*e - (omp_parms[16]-omp_parms[17]*eps)*eps;
   omp->a0s = 0.0;

   // imaginary volume depth, radius, diffuseness
   omp->wv1 = omp_parms[18] - omp_parms[21]*eps;
   omp->wv2 = omp_parms[19] - omp_parms[22]*eps;
   omp->wv3 = -omp_parms[20];

   omp->rv = omp_parms[23] + (omp_parms[24] + omp_parms[25]*(double)at)/(omp_parms[26] + (double)at + omp_parms[27]*e) + omp_parms[28]*e*e;
   omp->av = omp_parms[29] - omp_parms[30]*e/(-omp_parms[31] - e) + (omp_parms[32] - omp_parms[33]*e)*eps;

   // imaginary surface depth, radius, diffuseness
   omp->ws1 = omp_parms[34] - omp_parms[36]*eps;
   omp->ws2 = -omp_parms[35] + omp_parms[37]*eps;
   omp->ws3 = 0.0;

   if (e > 40.) {
      omp->ws1 = 0.0;
      omp->ws2 = 0.0;
   }

   omp->rs = omp_parms[38] - omp_parms[40]*e - omp_parms[39]/a3;
   omp->as = omp_parms[41];

   // spin orbit parameters
   omp->vso1 = omp_parms[42] - omp_parms[43]*(double)at;
   omp->vso2 = 0.0;
   omp->vso3 = 0.0;
 
   omp->wso1 = 0.0;
   omp->wso2 = 0.0;
   omp->wso3 = 0.0;

   omp->rvso = omp_parms[44] - omp_parms[45]/a3;
   omp->rwso = 0.0;
   omp->avso = omp_parms[46] - omp_parms[47]*(double)at;
   omp->awso = 0.0;

   return(0x0013);

}

unsigned int OMPWLH5(double e, int at, int zt, int ai, int zi, Optical *omp)
{

   double eps=1.0-2.0*(double)zt/(double)at;
   double a3 = pow((double)at,1.0/3.0);
   double omp_parms[48] = { 0., 0., 5.12707e+01,  2.55931e-01, -9.42339e-04,  4.62002e-06,
         2.09509e+01,  3.30513e-01,  1.65245e-03,  1.34302e+00,
         3.36900e-01,  1.24823e-03,  1.52620e-05,  7.72197e-01,
         1.24094e-03,  1.83160e-05,  3.30972e-01,  1.03029e+00,
         3.70016e+00,  2.80893e-01,  5.25183e-04,  1.00883e+01,
         8.73254e-03,  6.46281e-01,  6.95997e+01,  6.85829e-01,
         9.13809e+01,  1.05856e+00,  5.74254e-06,  3.12911e-01,
         5.33470e-01,  4.41333e+00,  5.41500e-02,  7.73086e-04,
         1.42733e+00,  4.63125e-02,  2.04585e+00,  1.70510e-01,
         1.27478e+00,  1.56846e+00,  1.88844e-03,  8.78288e-01,
         9.59626e+00,  9.03525e-03,  1.33041e+00,  9.53738e-01,
         8.49260e-01,  3.44680e-04};

   if(ai != 1 || zi != 0) return(0);

   // real volume depth, radius, diffuseness
   omp->v1 = omp_parms[2] - omp_parms[6]*eps;
   omp->v2 = -omp_parms[3] + omp_parms[7]*eps;
   omp->v3 = omp_parms[4] - omp_parms[8]*eps;
   omp->v4 = omp_parms[5];

   omp->vs1 = 0.0;
   omp->vs2 = 0.0;
   omp->vs3 = 0.0;

   omp->r0 = omp_parms[9] - omp_parms[11]*e +omp_parms[12]*e*e - omp_parms[10]/a3;
   omp->r0s = 0.0;
   omp->a0 = omp_parms[13] + omp_parms[14]*e - omp_parms[15]*e*e - (omp_parms[16]-omp_parms[17]*eps)*eps;
   omp->a0s = 0.0;

   // imaginary volume depth, radius, diffuseness
   omp->wv1 = omp_parms[18] - omp_parms[21]*eps;
   omp->wv2 = omp_parms[19] - omp_parms[22]*eps;
   omp->wv3 = -omp_parms[20];

   omp->rv = omp_parms[23] + (omp_parms[24] + omp_parms[25]*(double)at)/(omp_parms[26] + (double)at + omp_parms[27]*e) + omp_parms[28]*e*e;
   omp->av = omp_parms[29] - omp_parms[30]*e/(-omp_parms[31] - e) + (omp_parms[32] - omp_parms[33]*e)*eps;

   // imaginary surface depth, radius, diffuseness
   omp->ws1 = omp_parms[34] - omp_parms[36]*eps;
   omp->ws2 = -omp_parms[35] + omp_parms[37]*eps;
   omp->ws3 = 0.0;

   if (e > 40.) {
      omp->ws1 = 0.0;
      omp->ws2 = 0.0;
   }

   omp->rs = omp_parms[38] - omp_parms[40]*e - omp_parms[39]/a3;
   omp->as = omp_parms[41];

   // spin orbit parameters
   omp->vso1 = omp_parms[42] - omp_parms[43]*(double)at;
   omp->vso2 = 0.0;
   omp->vso3 = 0.0;
 
   omp->wso1 = 0.0;
   omp->wso2 = 0.0;
   omp->wso3 = 0.0;

   omp->rvso = omp_parms[44] - omp_parms[45]/a3;
   omp->rwso = 0.0;
   omp->avso = omp_parms[46] - omp_parms[47]*(double)at;
   omp->awso = 0.0;

   return(0x0013);

}

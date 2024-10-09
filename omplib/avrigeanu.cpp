/***********************************************************/
/*      V. Avrigeanu, P. Hodgson, M. Avrigeanu             */
/*      Phys. Rev. C49, 2136 (1994)                        */
/***********************************************************/

#include <cmath>

#include "omplib.h"

unsigned int OMPAvrigeanu(double e, int at, int zt, int ai, int zi, Optical *omp)
{
   double a13;

   if(ai != 4 || zi != 2) return(0);

   a13 = pow((double)at,1.0/3.0);

   omp->v1  =101.1 + 6.051*zt/a13;    omp->v2  = -0.248;  omp->v3  =  0.0;
   if(e > 73.0){
     omp->wv1 = 26.82- 1.706*a13   ;  omp->wv2 =  0.006;  omp->wv3 =  0.0;
   }else{
     omp->wv1 = 12.64- 1.706*a13   ;  omp->wv2 =  0.200;  omp->wv3 =  0.0;
   }

   omp->ws1 =  0.0;    omp->ws2 =  0.0;    omp->ws3 =  0.0;
   omp->vso1=  0.0;    omp->vso2=  0.0;    omp->vso3=  0.0;

   omp->r0  =  1.245;  omp->a0  =  0.817 - 0.0085*a13;
   omp->rv  =  1.570;  omp->av  =  0.692 - 0.020 *a13;
   omp->rs  =  0.0;    omp->as  =  0.0;
   omp->rc  =  1.3;
/*
    This Coulomb radius was taken from L.W.Put, Nucl.Phys. A291, 93 (1977)
*/
   return(0x0011);
}


/***********************************************************/
/*      M. Avrigeanu, A.C. Obreja, F.L. Roman,             */
/*      V. Avrigeanu, W. von Oertzen                       */
/*      At. Data Nucl. Data Tables 95, 501 (2009)          */
/*                                                         */
/*      Regional OMP, 50 < A < 125, E < 50MeV              */
/***********************************************************/
unsigned int OMPAvrigeanu2009(double e, int at, int zt, int ai, int zi, Optical *omp)
{
  double a13,e1,e2,e3,rb;

  if(ai!=4 || zi!=2) return(0);

  a13=pow((double)at,1.0/3.0);

  rb = 2.66 + 1.36*a13;
  e2 = (2.59 + 10.4/(double)at)*zt/rb;
  e1 = -3.03 - 0.762*a13 + 1.24*e2;
  e3 = 23.6 + 0.181*zt/a13;


  if(e < e3){
    omp->v1  =168.0 + 0.733*zt/a13;  omp->v2  = -2.64;
  }else{
    omp->v1  =116.5 + 0.337*zt/a13;  omp->v2  = -0.453;
  }

  omp->wv1 = 2.73- 2.88*a13;  omp->wv2 =  1.11;

  if(e < e1){
   omp->ws1 =  4.0;                         omp->ws2 =  0.0;
  }else if(e < e2){
   omp->ws1 =  22.2 + 4.57*a13 - 7.446*e2;  omp->ws2 =  6.0;
  }else{
   omp->ws1 =  22.2 + 4.57*a13;             omp->ws2 = -1.446;
  }

  omp->r0  =  (e<25.0) ? 1.18 + 0.012*e : 1.48;
  omp->a0  =  (e<e2) ? 0.671 + 0.0012*at + (0.0094 - 0.0042*a13)*e2 :
                       fmax(0.671 + 0.0012*at + (0.0094 - 0.0042*a13)*e,0.55);

  omp->rv  =  1.34;
  omp->av  =  0.50;
  omp->rs  =  1.52;
  omp->as  =  0.729 - 0.074*a13;

  omp->rc  =  1.3;

   return(0x0013);
}


/***********************************************************/
/*      V. Avrigeanu, M. Avrigeanu, C. Manailescu          */
/*      Phys. Rev. C 90, 044612 (2014)                     */
/***********************************************************/
unsigned int OMPAvrigeanu2014(double e, int at, int zt, int ai, int zi, Optical *omp)
{
  double a13,e1,e2,e3,e4,rb;

  if(ai!=4 || zi!=2) return(0);

  a13=pow((double)at,1.0/3.0);

  rb = 2.66 + 1.36*a13;
  e2 = (2.59 + 10.4/(double)at)*zt/rb;
  e1 = -3.03 - 0.762*a13 + 1.24*e2;
  e3 = 22.2 + 0.181*zt/a13;
  e4 = 29.1 - 0.22*zt/a13;


  if(e < e3){
    omp->v1  =165.0 + 0.733*zt/a13;  omp->v2  = -2.64;
  }else{
    omp->v1  =116.5 + 0.337*zt/a13;  omp->v2  = -0.453;
  }

  omp->wv1 = 2.73- 2.88*a13;  omp->wv2 =  1.11;

  if(e <= e1){
   omp->ws1 =  4.0;                         omp->ws2 =  0.0;
  }else if(e <= e2){
   omp->ws1 =  22.2 + 4.57*a13 - 7.446*e2;  omp->ws2 =  6.0;
  }else{
   omp->ws1 =  22.2 + 4.57*a13;             omp->ws2 = -1.446;
  }

  omp->r0  =  (e <= 25.0) ? 1.18 + 0.012*e : 1.48;
  if(e <= e2){
    omp->a0  =  0.631 + (0.016 - 0.001*e2)*zt/a13;
  }else if(e <= e4){
    omp->a0  =  0.631 + 0.016*zt/a13 - (0.001*zt/a13)*e;
  }else{
    omp->a0  =  0.684 - 0.016*zt/a13 - (0.0026 - 0.00026*zt/a13)*e;
  }

  omp->rv  =  1.34;
  omp->av  =  0.50;
  
  omp->rs  =  1.52;
  if((152 < at) && (at < 190)) omp->rs = fmax(1.74 - 0.01*e, 1.52);

  omp->as  =  0.729 - 0.074*a13;

  omp->rc  =  1.3;

   return(0x0013);
}



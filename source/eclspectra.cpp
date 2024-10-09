/******************************************************************************/
/*  eclspectra.cpp                                                            */
/*        exclusive spectra calculation by integral method                    */
/******************************************************************************/

#include <iostream>
#include <cstdlib>

#include "structur.h"
#include "eclipse.h"

/*** skip loop if probability is less than this value */
static const double gProbabilityCut = 0.0;

/*** include/exclude 4th neutron emission */
static const bool include4n = true;

/*** calculation monitor */
static const bool monitor = false;

static void eclGammaSpec (const int, const int, const int, const double, double ***, double *, double *);
//static void eclGammaCapt (const int, const int, const int, const double, double ***, double *, double *);


/**********************************************************/
/*      Calculation of Particle Emission Spectra          */
/*      -----                                             */
/*           The most complicated subroutine in CoH       */
/*           You are not supposed to understand this.     */
/**********************************************************/
void eclSpectra(const bool ddx, const bool fns, const int cm, int *ktot, int **cidx, double ****ptbl, EXSpectra *dat)
{
  const int c0 = 0;
  int ncut = -1;

  /*** special case for the capture gamma-ray */
  if(ddx) ncut = dat[c0].getNc();
//eclGammaCapt(ncut,0,ktot[c0],1.0,ptbl[c0],dat[c0].spec[0],dat[c0].glin);
  eclGammaSpec(ncut,0,ktot[c0],1.0,ptbl[c0],dat[c0].spec[0],dat[c0].glin);

  /*** first gamma-ray before particle emission */
  for(int k0=0 ; k0<=ktot[c0] ; k0++){

    if(monitor) std::cerr << "monitor " << k0 << " / " << ktot[c0] << std::endl;

    if(ptbl[c0][k0][0][k0] < 0.0) continue;

    double p0 = 1.0;
    if(k0 != 0) p0 = ptbl[c0][0][0][k0];
    if(p0 <= gProbabilityCut) continue;
    if(!fns) p0 *= 1.0 - dat[c0].pfis[k0]; // fission probability at (c0,k0) subtracted

    /*** first particle emission */
    for(int i0=1 ; i0<cm ; i0++){ // i0 : decay channel from the first CN
      int c1 = cidx[c0][i0];      // c1 : CN index after particle emission
      if(c1 < 0) continue;

      /*** sum over residual excitation */
      for(int j0=k0 ; j0<=ktot[c1] ; j0++){
        double q0 = ptbl[c0][k0][i0][j0];
        if(q0 <= gProbabilityCut) continue;

        /*** see if nucleus after particle emission emits one more particle */
        double s1 = 0.0;
        if(fns) s1= dat[c1].pfis[j0];
        else{
          if(ptbl[c1][j0][0][j0] < 0.0)  s1 = 1.0;
          else{
            for(int k1=j0 ; k1<=ktot[c1] ; k1++) s1 += ptbl[c1][j0][0][k1];
          }
        }
        /*** probabilities :
               p0 : probability to arrive at (c0,k0)
               q0 : probability of i0-kind particle emission from (c0,k0) to (c1,j0)
               s1 : probability of gamma-decay from (c1,j0) state, or (c1,j0) fissions [fns case]
             assumption :
               if (c1,j0) emits gamma, the final production is c1. */
        s1 *= p0*q0;

        /*** exclude binary reactions for which the final states are
             in the discrete levels. we can know this by
               k0 = 0, which means that the transition is from the initial state, and
               j0 > Nc, the final state is below the continuum boundary. */
        if(ddx){
          ncut = dat[c1].getNc();
          if(dat[c1].getNl() == 1){
            dat[c1].spec[i0][j0-k0] += s1;
            eclGammaSpec(-1,j0,ktot[c1],s1,ptbl[c1],dat[c1].spec[0],dat[c1].glin);
          }
          else{
            eclGammaSpec(ncut,j0,ktot[c1],s1,ptbl[c1],dat[c1].spec[0],dat[c1].glin);
            /*** binary reaction, residual is in a discrete level */
            if( (k0 == 0) && (j0 >= ncut-1) ) dat[c1].sbin[j0-k0] += s1;
            /*** residual is in a continuum */
            else dat[c1].spec[i0][j0-k0] += s1;
          }
        }
        else{
          dat[c1].spec[i0][j0-k0] += s1;
          eclGammaSpec(-1,j0,ktot[c1],s1,ptbl[c1],dat[c1].spec[0],dat[c1].glin);
        }

        if(k0 != 0) dat[c1].spec[0][k0] += s1;

        /*** second gamma-ray before second particle emission */
        for(int k1=j0 ; k1<=ktot[c1] ; k1++){

          if(ptbl[c1][k1][0][k1] < 0.0) continue;

          double p1 = 1.0;
          if(k1 != j0){
            p1 = ptbl[c1][j0][0][k1];
            if(p1 <= gProbabilityCut) continue;
          }
          if(!fns) p1 *= 1.0 - dat[c1].pfis[k1];


          /*** second particle emission */
          for(int i1=1 ; i1<cm ; i1++){ // i1 : decay channel from the second CN
            int c2 = cidx[c1][i1];      // c2 : CN index after particle emission
            if(c2 < 0) continue;

            /*** sum over residual excitation */
            for(int j1=k1 ; j1<=ktot[c2] ; j1++){
              double q1 = ptbl[c1][k1][i1][j1];
              if(q1 <= gProbabilityCut) continue;

              double s2 = 0.0;
              if(fns) s2 = dat[c2].pfis[j1];
              else{
                if(ptbl[c2][j1][0][j1] < 0.0)  s2 = 1.0;
                else{
                  for(int k2=j1 ; k2<=ktot[c2] ; k2++) s2 += ptbl[c2][j1][0][k2];
                }
              }

              s2 *= p0*q0*p1*q1;
              dat[c2].spec[i0][j0-k0] += s2;
              dat[c2].spec[i1][j1-k1] += s2;

              if(k0 !=  0) dat[c2].spec[0][k0   ] += s2;
              if(k1 != j0) dat[c2].spec[0][k1-j0] += s2;

              if(ddx) ncut = dat[c2].getNc();
              eclGammaSpec(ncut,j1,ktot[c2],s2,ptbl[c2],dat[c2].spec[0],dat[c2].glin);


              /*** third gamma-ray before third particle emission */
              for(int k2=j1 ; k2<=ktot[c2] ; k2++){
                if(ptbl[c2][k2][0][k2] < 0.0) continue;

                double p2 = 1.0;
                if(k2 != j1){
                  p2 = ptbl[c2][j1][0][k2];
                  if(p2 <= gProbabilityCut) continue;
                }
                if(!fns) p2 *= 1.0 - dat[c2].pfis[k2];


                /*** third particle emission */
                for(int i2=1 ; i2<cm ; i2++){ // i2 : decay channel from the third CN
                  int c3 = cidx[c2][i2];      // c3 : CN index after particle emission
                  if(c3 < 0) continue;

                  /*** the third particle must be for (n,3n) only */
                  if( (i0 != 1) || (i1 != 1) || (i2 != 1) ) continue;

                  /*** sum over residual excitation */
                  for(int j2=k2 ; j2<=ktot[c3] ; j2++){
                    double q2 = ptbl[c2][k2][i2][j2];
                    if(q2 <= gProbabilityCut) continue;

                    double s3 = 0.0;
                    if(fns) s3 = dat[c3].pfis[j2];
                    else{
                      if(ptbl[c3][j2][0][j2] < 0.0)  s3 = 1.0;
                      else{
                        for(int k3=j2 ; k3<=ktot[c3] ; k3++) s3 += ptbl[c3][j2][0][k3];
                      }
                    }
                    s3 *= p0*q0*p1*q1*p2*q2;
                    dat[c3].spec[i0][j0-k0] += s3;
                    dat[c3].spec[i1][j1-k1] += s3;
                    dat[c3].spec[i2][j2-k2] += s3;

                    if(k0 !=  0) dat[c3].spec[0][k0   ] += s3;
                    if(k1 != j0) dat[c3].spec[0][k1-j0] += s3;
                    if(k2 != j1) dat[c3].spec[0][k2-j1] += s3;

                    if(ddx) ncut = dat[c3].getNc();
                    eclGammaSpec(ncut,j2,ktot[c3],s3,ptbl[c3],dat[c3].spec[0],dat[c3].glin);

                    if(!include4n) continue;
                    /*** 4th gamma-ray before 4th particle emission */
                    for(int k3=j2 ; k3<=ktot[c3] ; k3++){
                      if(ptbl[c3][k3][0][k3] < 0.0) continue;

                      double p3 = 1.0;
                      if(k3 != j2){
                        p3 = ptbl[c3][j2][0][k3];
                        if(p3 <= gProbabilityCut) continue;
                      }
                      if(!fns) p3 *= 1.0 - dat[c3].pfis[k3];

                      /*** 4th particle emission */
                      for(int i3=1 ; i3<cm ; i3++){ // i3 : decay channel from the 4th CN
                        int c4 = cidx[c3][i3];      // c4 : CN index after particle emission
                        if(c4 < 0) continue;

                        /*** the 4th particle must be for (n,4n) only */
                        if( (i0 != 1) || (i1 != 1) || (i2 != 1) || (i3 != 1) ) continue;

                        /*** sum over residual excitation */
                        for(int j3=k3 ; j3<=ktot[c4] ; j3++){
                          double q3 = ptbl[c3][k3][i3][j3];
                          if(q3 <= gProbabilityCut) continue;

                          double s4 = 0.0;
                          if(fns) s4 = dat[c4].pfis[j3];
                          else{
                            if(ptbl[c4][j3][0][j3] < 0.0)  s4 = 1.0;
                            else{
                              for(int k4=j3 ; k4<=ktot[c4] ; k4++) s4 += ptbl[c4][j3][0][k4];
                            }
                          }
                          s4 *= p0*q0*p1*q1*p2*q2*p3*q3;
                          dat[c4].spec[i0][j0-k0] += s4;
                          dat[c4].spec[i1][j1-k1] += s4;
                          dat[c4].spec[i2][j2-k2] += s4;
                          dat[c4].spec[i3][j3-k3] += s4;

                          if(k0 !=  0) dat[c4].spec[0][k0   ] += s4;
                          if(k1 != j0) dat[c4].spec[0][k1-j0] += s4;
                          if(k2 != j1) dat[c4].spec[0][k2-j1] += s4;
                          if(k3 != j2) dat[c4].spec[0][k3-j2] += s4;

                          if(ddx) ncut = dat[c4].getNc();
                          eclGammaSpec(ncut,j3,ktot[c4],s4,ptbl[c4],dat[c4].spec[0],dat[c4].glin);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


/**********************************************************/
/*      Calculation of Gamma-Ray Spectra                  */
/**********************************************************/
void eclGammaSpec(const int ncut, const int k0, const int ktot, const double f, double ***ptbl, double *gray, double *glin)
{
  if(f == 0.0) return;
  if(k0 == ktot) return;

  double *gpop;
  gpop = new double [ktot+1];

  for(int k1=k0 ; k1<=ktot ; k1++) gpop[k1] = 0.0;

  /*** first generation */
  for(int k1=k0 ; k1<=ktot ; k1++){
    double p0 = ptbl[k0][0][k1];
    if(p0 <= gProbabilityCut) continue;

    /*** g0 is a population at k1, fed from k0, multiplied by 
         a probability of which k1 decays by gamma emission */
    double g0 = 0.0;
    for(int k2=k1 ; k2<=ktot ; k2++) g0 += ptbl[k1][0][k2];
    g0 *= f*p0;

    if( (ncut >= 0) && (k0 >= ncut-1) ) glin[k1-k0] += g0;
    else                                gray[k1-k0] += g0;
    gpop[k1] = g0;
  }

  /*** for all second generation and above, no particle assumed */
  for(int k1=k0+1 ; k1<ktot ; k1++){

    double g0 = 0.0;
    for(int k2=k1 ; k2<=ktot ; k2++) g0 += ptbl[k1][0][k2];
    if(g0 <= gProbabilityCut) continue;

    g0 = gpop[k1] /g0;

    for(int k2=k1 ; k2<=ktot ; k2++){
      if(ptbl[k2][0][k2] < 0.0) continue;

      double g1 = g0 * ptbl[k1][0][k2];
      if(g1 <= gProbabilityCut) continue;

      if( (ncut >= 0) && (k1 >= ncut-1) ) glin[k2-k1] += g1;
      else                                gray[k2-k1] += g1;
      gpop[k2] += g1;
    }
  }

  delete [] gpop;
}


/*** multiplicity summed up to five */
/*
void eclGammaCapt(int ncut, int k0, int ktot, double f, double ***ptbl, double *gray, double *glin)
{
  if(f == 0.0) return;

  for(int k1=k0 ; k1<=ktot ; k1++){

    if(ptbl[k1][0][k1] < 0.0) continue;

    double p1 = ptbl[k0][0][k1];
    if(p1 == 0.0) continue;

    double x1 = f*p1;
    double g1 = 0.0;
    for(int m1=k1 ; m1<=ktot ; m1++) g1 += ptbl[k1][0][m1];
    gray[k1-k0] += x1*g1;

    for(int k2=k1+1 ; k2<=ktot ; k2++){

      if(ptbl[k2][0][k2] < 0.0) continue;

      double p2 = ptbl[k1][0][k2];
      if(p2 == 0.0) continue;

      double x2 = x1*p2;
      double g2 = 0.0;
      for(int m2=k2 ; m2<=ktot ; m2++) g2 += ptbl[k2][0][m2];

      if( (ncut >= 0) && (k1 > ncut) ) glin[k2-k1] += x2*g2;
      else                             gray[k2-k1] += x2*g2;

      for(int k3=k2+1 ; k3<=ktot ; k3++){

        if(ptbl[k3][0][k3] < 0.0) continue;

        double p3 = ptbl[k2][0][k3];
        if(p3 == 0.0) continue;

        double x3 = x2*p3;
        double g3 = 0.0;
        for(int m3=k3 ; m3<=ktot ; m3++) g3 += ptbl[k3][0][m3];

        if( (ncut >= 0) && (k2 > ncut) ) glin[k3-k2] += x3*g3;
        else                             gray[k3-k2] += x3*g3;

        for(int k4=k3+1 ; k4<=ktot ; k4++){

          if(ptbl[k4][0][k4] < 0.0) continue;

          double p4 = ptbl[k3][0][k4];
          if(p4 == 0.0) continue;

          double x4 = x3*p4;
          double g4 = 0.0;
          for(int m4=k4 ; m4<=ktot ; m4++) g4 += ptbl[k4][0][m4];

          if( (ncut >= 0) && (k3 > ncut) ) glin[k4-k3] += x4*g4;
          else                             gray[k4-k3] += x4*g4;

          for(int k5=k4 ; k5<=ktot ; k5++){

            if(ptbl[k5][0][k5] < 0.0) continue;

            double p5 = ptbl[k4][0][k5];
            if(p5 == 0.0) continue;

            double x5 = x4*p5;
            if( (ncut >= 0) && (k4 > ncut) ) glin[k5-k4] += x5;
            else                             gray[k5-k4] += x5;
          }
        }
      }
    }
  }
}
*/


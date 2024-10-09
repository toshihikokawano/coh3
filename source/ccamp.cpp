/******************************************************************************/
/*  ccamp.cpp                                                                 */
/*        scattering amplitudes for coupled-channels method                   */
/******************************************************************************/

#include <cmath>

#include "optical.h"
#include "coupling.h"
#include "etc.h"
#include "terminate.h"

static inline double abs_square(std::complex<double> x)
{ return( x.real() * x.real() + x.imag() * x.imag() ); }

extern double *fact;

static std::complex<double> ***amp;


/**********************************************************/
/*      Scattering Amplitude Memory Allocation            */
/*      --------                                          */
/*      index of amp: [level][angle][M (m0,M0,m1,M1)]     */
/**********************************************************/
void ccAllocateAmplitute(Collective *col)
{
  try{
    amp = new std::complex<double> ** [col->nlevel];

    int c0 = (int)(2*fabs(col->lev[0].spin));
    int m0 = 3*(c0 + 1);

    for(int i=0 ; i<col->nlevel ; i++){
      amp[i] = new std::complex<double> * [MAX_ANGDIST];

      int c1 = (int)(2*fabs(col->lev[i].spin));
      int m1 = 3*(c1 + 1);

      for(int j=0 ; j<MAX_ANGDIST ; j++){
        amp[i][j] = new std::complex<double> [m0*m1];
        for(int k=0 ; k<m0*m1 ; k++) amp[i][j][k] = std::complex<double>(0.0,0.0);
      }
    }
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("ccAllocateAmplitude");
  }
}


/**********************************************************/
/*      Release Scattering Amplitude Memory               */
/**********************************************************/
void ccFreeAmplitude(Collective *col)
{
  for(int i=0 ; i<col->nlevel ; i++){
    for(int j=0 ; j<MAX_ANGDIST ; j++){
      delete [] amp[i][j];
    }
    delete [] amp[i];
  }
  delete [] amp;
}


/**********************************************************/
/*      Scattering Amplitude for CC model                 */
/*      Satchler, Direct Reactions, Eq. (5.24a)           */
/**********************************************************/
void ccScatteringAmplitude(const int jj, const int nmax, double *th, Collective *col, CCdata *cdt, std::complex<double> *c)
{
  unsigned long idx;

  int is0 = (int)(2*col->pspin);
  int is1 = is0;        // inelastic scattering only, s0=s1

  /*** for target level */
  for(int n0=0 ; n0<nmax ; n0++){
    if(cdt[n0].level != col->target_index) continue;

    int i0 = cdt[n0].chn.st2;
    int l0 = (int)(2.0*cdt[n0].chn.l);
    int j0 = cdt[n0].chn.j2;
    int k0 = cdt[n0].level;

    std::complex<double> siga(0.0,0.0);
    if(col->lev[k0].coulomb.eta != 0.0)
      siga = omCoulombPhaseFactor(cdt[n0].chn.l,col->lev[k0].coulomb);

    for(int m0=-is0 ; m0<=is0 ; m0+=2){     int m0k  = (m0 +is0)/2; // projectile m
      for(int mm0=-i0 ; mm0<=i0 ; mm0+=2){  int mm0k = (mm0+i0 )/2; // target M
        int mp0 = m0+mm0;

        double cg0 = clebsh_gordan(l0,is0,0,m0,j0)*clebsh_gordan(j0,i0,m0,mm0,jj)*sqrt(l0+1.0); 
        if(cg0 == 0.0) continue;

        /*** for each discrete level */
        for(int n1=0 ; n1<nmax ; n1++){

          int i1 = cdt[n1].chn.st2;
          int l1 = (int)(2.0*cdt[n1].chn.l);
          int j1 = cdt[n1].chn.j2;
          int k1 = cdt[n1].level;

          int idc= n0*nmax+n1;

          std::complex<double> sigb(0.0,0.0);
          std::complex<double> a(0.0,0.0);

          if(col->lev[k1].coulomb.eta != 0.0){
            std::complex<double> sig(0.0,0.0);
            sigb = omCoulombPhaseFactor(cdt[n1].chn.l,col->lev[k1].coulomb);
            sig = siga * sigb;
            a   = sig  * c[idc];
          }
          else
            a = c[idc];

          for(int m1=-is1 ; m1<=is1 ; m1+=2){     int m1k  = (m1 +is1)/2;
            for(int mm1=-i1 ; mm1<=i1 ; mm1+=2){  int mm1k = (mm1+i1 )/2;
              int mp1 = m1+mm1;
              int ml1 = mp0-mp1;
              int m2  = std::abs(ml1/2);
              double x1  = ((ml1<0) ? ( (m2%2==0) ? 1 : -1 ) : 1)
                           *sqrt(exp(fact[cdt[n1].chn.l-m2]-fact[cdt[n1].chn.l+m2]));

              idx = (unsigned int)(  (i1+1)*((is1+1)*((i0+1)*m0k+mm0k)+m1k)+mm1k  );

              double cg1 = clebsh_gordan(l1,is1,ml1,m1,j1)*clebsh_gordan(j1,i1,ml1+m1,mm1,jj)*sqrt(l1+1.0);
              if(cg1 == 0.0) continue;

              double x2 = x1*cg0*cg1;

              for(int j=0 ; j<MAX_ANGDIST ; j++){
                double x3 = x2*assocLegendrePol(cdt[n1].chn.l,m2,th[j]);
                amp[k1][j][idx] += x3*a;
              }
            }
          }
        }
      }
    }
  }
}


/**********************************************************/
/*      Angular Ditribution for Coupled-Channels          */
/**********************************************************/
double ccAngularDistribution(const int k, const int j, double *th, Collective *col)
{
  unsigned long idx;
  double cx[3];

  double wave  = col->lev[col->target_index].wave_number;
  double sfact = 1.0/((2*col->pspin+1)*(2*fabs(col->lev[col->target_index].spin)+1.0));

  int is0 = (int)(2.0*col->pspin);
  int is1 = is0;
  int k0  = col->target_index;
  int i0  = (int)(2.0*fabs(col->lev[k0].spin));
  int i1  = (int)(2.0*fabs(col->lev[k].spin));

  cx[0] = cx[1] = cx[2] = 0.0; // cx[1] saved for polarization

  std::complex<double> f(0.0,0.0);
  if((k == col->target_index) && (col->lev[k0].coulomb.eta != 0.0)){
    f = omRutherfordAmplitude(th[j],col->lev[k0].wave_number,col->lev[k0].coulomb);
    cx[2] = abs_square(f);
  }

  for(int m0=-is0 ; m0<=is0 ; m0+=2){        int m0k  = (m0 +is0)/2;
    for(int mm0=-i0 ; mm0<=i0 ; mm0+=2){     int mm0k = (mm0+i0 )/2;
      for(int m1=-is1 ; m1<=is1 ; m1+=2){    int m1k  = (m1 +is1)/2;
        for(int mm1=-i1 ; mm1<=i1 ; mm1+=2){ int mm1k = (mm1+i1 )/2;
          idx = (unsigned int)(  (i1+1)*((is1+1)*((i0+1)*m0k+mm0k)+m1k)+mm1k  );
          std::complex<double> a = amp[k][j][idx] / wave;
          if((m0 == m1) && (mm0 == mm1)) a += f;

          cx[0] += abs_square(a);
        }
      }
    }
  }

  return(sfact * cx[0]);
}


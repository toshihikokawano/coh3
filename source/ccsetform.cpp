/******************************************************************************/
/*  ccsetform.cpp                                                             */
/*        calculate radial part potentials for coupled-channels               */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "optical.h"
#include "etc.h"
#include "coupling.h"

#define GAUSS_LEGENDRE20
#include "gauss.h"

static void   ccDeformedRadius (Optical *, Optical *, const double);
static double ccPotentialRadialCoulombRot (const int, const double, Optical *);
static double ccPotentialRadialCoulombVib (const int, const int, const double, Optical *);
static double ccChargedistRatio (const int, double *);

#undef  DEBUG

/**********************************************************/
/*     Deformed Optical Potential Radial Form             */
/*     Rotational Excitation by Legendre Exapnsion        */
/**********************************************************/
void ccPotentialFormRot(const int nclev, const int bmax, double *beta, const int zzprod, Optical *omp, double mu, Potential *pc0, Potential **pc1)
{
  double c0 = zzprod * PERMITTIV * COULOMBSQ;
  double c1 = 2.0*mu*AMUNIT/VLIGHTSQ/HBARSQ;
  double c2 = c1*CSPO;
  double c4 = sqrt(PI4);

  double cratio = (zzprod == 0) ? 0.0 : ccChargedistRatio(bmax,beta);
  for(int n=0 ; n<nclev ; n++){

    Optical dom = omp[n]; // Set the OMP

    for(int l=0 ; l<MAX_LAMBDA ; l++){
      int lambda = 2*l;

      double c3 = (zzprod == 0) ? 0.0 : c0 *3.0 * cratio / ((2*lambda+1.0)*pow(omp[n].Rc,3.0));
      double y  = sqrt( (2*lambda+1.0) ) * c1;

      /***  Tamura Eq.(16), except no sqrt(4pi) factor, to make the diagonal 
            the same as usual optical potential */
      for(int i=0 ; i<=pc0[n].n_match ; i++){
        double r  = (double)i * pc0[n].width;

        double sr = 0.0, si = 0.0;

        for(int j=0 ; j<MAX_GAUSSLEG ; j++){
          double cosx1 = 0.5*(1+gaussleg_x[j]);
          double cosx2 = 0.5*(1-gaussleg_x[j]);

          /*** Radius = R(1 + sum_L beta_L sqrt((2L+1)/4pi) P_L(cos theta) */
          double rx1 = 0.0;
          double rx2 = 0.0;
          for(int l=2 ; l<=bmax ; l+=2){
            /*** index of beta: 0 = beta2, 1 = beta4, 2 = beta6 */
            double x = beta[l/2-1]*sqrt( (2*l+1.0)/PI4 );
            rx1 += x*legendre(l,cosx1);
            rx2 += x*legendre(l,cosx2);
          }

          ccDeformedRadius(&dom,&omp[n],rx1);
          std::complex<double> p1 = omPotentialRadialMean(r,&dom);
          double  v1 = (zzprod !=0 ) ? c3 * ccPotentialRadialCoulombRot(lambda,r,&dom) : 0.0;

          ccDeformedRadius(&dom,&omp[n],rx2);
          std::complex<double> p2 = omPotentialRadialMean(r,&dom);
          double  v2 = (zzprod !=0 ) ? c3 * ccPotentialRadialCoulombRot(lambda,r,&dom) : 0.0;

          double s1 = y*legendre(lambda,cosx1);
          double s2 = y*legendre(lambda,cosx2);

          /*** deformed Coulomb potential subtracted here */
          sr += 0.5*(s1*(p1.real()- v1) + s2*(p2.real()- v2))*gaussleg_a[j];
          si += 0.5*(s1*(p1.imag()    ) + s2*(p2.imag()    ))*gaussleg_a[j];
        }

        pc1[n][l].mean_field[i] = c4 * std::complex<double>(sr,si);
        pc1[n][l].spin_orbit[i] = c2 * omPotentialRadialSpo(r,&dom);
        pc1[n][l].coulomb[i]    = (zzprod != 0) ? c0*c1* omPotentialRadialCoulomb(r,&dom) : 0.0;

        if(l == 0){
          pc0[n].mean_field[i] = std::complex<double>(sr,si);
          pc0[n].spin_orbit[i] = pc1[n][l].spin_orbit[i];
          pc0[n].coulomb[i]    = pc1[n][l].coulomb[i];
        }
      }
    }
  }
#ifdef DEBUG
  for(int n=0 ; n<nclev ; n++){
    for(int i=0 ; i<=pc0[n].n_match ; i++){
      std::cout << std::setw(12) << (double)i * pc0[n].width;
      std::cout << std::setw(12) << pc0[n].mean_field[i].real();
      std::cout << std::setw(12) << pc0[n].mean_field[i].imag();

      for(int l=1 ; l<MAX_LAMBDA ; l++){
        std::cout << std::setw(12) << pc1[n][l].mean_field[i].real();
        std::cout << std::setw(12) << pc1[n][l].mean_field[i].imag();
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
#endif
}


/**********************************************************/
/*     Deformed Optical Potential Radial Form             */
/*     Phonon Excitation by Power Series Expansion        */
/**********************************************************/
void ccPotentialFormVib(const int nclev, const int zzprod, Optical *omp, const double mu, Potential *pc0, Potential **pc1, Potential **pc2)
{
  double c1 = 2.0*mu*AMUNIT/VLIGHTSQ/HBARSQ;
  double c2 = c1*CSPO;
  double c3 = c1*zzprod * PERMITTIV * COULOMBSQ;

  for(int n=0 ; n<nclev ; n++){
    for(int i=0 ; i<=pc0[n].n_match ; i++){
      double r  = (double)i * pc0[n].width;

      /*** usual optical potential */
      pc0[n].mean_field[i] = c1 * omPotentialRadialMean(r,&omp[n]);
      pc0[n].spin_orbit[i] = c2 * omPotentialRadialSpo(r,&omp[n]);
      pc0[n].coulomb[i]    = (zzprod != 0) ? c3 * omPotentialRadialCoulomb(r,omp) : 0.0;

      pc0[n].mean_field[i] -= pc0[n].coulomb[i];

      std::complex<double> v1 = c1 * omPotentialRadialFirstDeriv(r,&omp[n]);
      std::complex<double> v2 = c1 * omPotentialRadialSecondDeriv(r,&omp[n]);

      for(int l=0 ; l<MAX_LAMBDA ; l++){

        /*** first and second order derivatives */
        double vc1 = (zzprod !=0 ) ? c3 * ccPotentialRadialCoulombVib(1,l,r,omp) : 0.0;
        double vc2 = (zzprod !=0 ) ? c3 * ccPotentialRadialCoulombVib(2,l,r,omp) : 0.0;

        pc1[n][l].mean_field[i] = v1 - vc1;
        pc2[n][l].mean_field[i] = v2 - vc2;
      }
    }
  }
#ifdef DEBUG
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(3);

  for(int n=0 ; n<nclev ; n++){
    for(int i=0 ; i<=pc0[n].n_match ; i++){
      std::cout << std::setw(11) << (double)i * pc0[n].width;
      std::cout << std::setw(11) << pc0[n].mean_field[i].real();
      std::cout << std::setw(11) << pc0[n].mean_field[i].imag();
      
      for(int l=0 ; l<MAX_LAMBDA ; l++){
        std::cout << std::setw(11) << pc1[n][l].mean_field[i].real();
        std::cout << std::setw(11) << pc1[n][l].mean_field[i].imag();

        std::cout << std::setw(11) << pc2[n][l].mean_field[i].real();
        std::cout << std::setw(11) << pc2[n][l].mean_field[i].imag();
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
#endif
}


/**********************************************************/
/*     Charge radius ratio for deformed Coulomb           */
/**********************************************************/
double ccChargedistRatio(const int bmax, double *beta)
{
  double x1,x2,x3,b1,b2,b3;

  double vol = PI4/3.0;
  double c   = SQRT_1_PI/6.0;

  for(int l1=2 ; l1<=bmax ; l1+=2){     int k1 = l1/2 - 1;
    vol += beta[k1]*beta[k1];
    b1   = beta[k1];
    x1   = sqrt(l1+1.0);

    for(int l2=2 ; l2<=bmax ; l2+=2){   int k2 = l2/2 - 1;
      b2 = b1 * beta[k2];
      x2 = x1 * sqrt(l2+1.0);

      for(int l3=2 ; l3<=bmax ; l3+=2){ int k3 = l3/2 - 1;
          b3 = b2 * beta[k3];
          double w = wigner_3j(l1,l2,l3,0,0,0);
          x3 = x2 * b3 * sqrt(l3+1.0) * w * w;
          vol += c*x3;
      }
    }
  }

  return(PI4/(3.0*vol));
}


/**********************************************************/
/*     change potential radius                            */
/**********************************************************/
void ccDeformedRadius(Optical *dom, Optical *omp, const double rx)
{
  dom->R0  = omp->R0* (1.0+rx);
  dom->R0s = omp->R0s*(1.0+rx);
  dom->Rv  = omp->Rv* (1.0+rx);
  dom->Rs  = omp->Rs* (1.0+rx);
  dom->Rc  = omp->Rc* (1.0+rx);
}


/**********************************************************/
/*     deformed Coulomb potential                         */
/*     R.H. Bassel, R.M. Drisko, G.R. Satchler,           */
/*     ORNL-3240 (1962)                                   */
/**********************************************************/
double ccPotentialRadialCoulombRot(const int l, const double r, Optical *omp)
{
  double vc = 0.0;

  if(r == 0.0){
    vc = (l == 0) ? omp->Rc*omp->Rc/2.0 : 0.0;
  }else if(r <= omp->Rc){
    vc = r*r/(l+3.0) + ((l==2) ? pow(r,(double)l)*log(omp->Rc/r) : 
                                 r*r/(l-2.0)*(1.0-pow(r/omp->Rc,l-2.0)));
  }else{
    vc = r*r*pow(omp->Rc/r,l+3.0)/(l+3.0);
  }
  return(vc);
}


double ccPotentialRadialCoulombVib(const int n, const int l, const double r, Optical *omp)
{
  double vc = 0.0;

  if(r == 0.0){
    vc = 0.0;
  }
  else{
    double c = 3.0/(2*l + 1.0);
    if(r <= omp->Rc){
      vc = c * pow(r,(double)l) * pow(omp->Rc,-(l+1.0));
      if(n == 2) vc *= 1.0-l;
    }
    else{
      vc = c * pow(omp->Rc,(double)l) * pow(r,-(l+1.0));
      if(n == 2) vc *= l + 2.0;
    }
  }
  return(vc);
}



/******************************************************************************/
/*  HObasis.cpp                                                               */
/*        Harmonic Oscillator basis for cylindrical coordinate                */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "HObasis.h"

static void   HOSphericalBasis (SphericalBasis *);
static double OverlapFunction  (const int, const int, const int, Basis *, SphericalBasis *, Operator *);
static double OverlapFunctionA (const int, const int, const double);
static double OverlapFunctionB (const int, const int, const int, const int, const int);

extern double *fact;

#undef DEBUG

/**********************************************************/
/*      Cylindrical Harmonic Oscillator Basis             */
/*      --------                                          */
/*      generate single particle states                   */
/*      based on the cylindrical HO                       */
/**********************************************************/
int HOCylindricalBasis(Basis *bs)
{
  /*** scaling constants */
  bs->bp  = bs->b * pow(bs->q, 1.0/6.0);
  bs->bz  = bs->b / pow(bs->q, 1.0/3.0);
  bs->bp2 = bs->bp * bs->bp;
  bs->bz2 = bs->bz * bs->bz;

  double q13 = pow(bs->q, 1.0/3.0);
  double q23 = q13*q13;
  double x   = (bs->getN0() + 2.0)/q13 - 1.0 - 0.5/bs->q;

  bs->nzmax = (int)(pow(bs->q, 2.0/3.0) * (bs->getN0()+2.0) - bs->q - 0.5);
  bs->nrmax = (int)(0.5 * x);
  bs->lambdamax = (int)x;

  double e0 = bs->getN0() + 2.0;
  int omega_max = 4*bs->nrmax+3;

  /*** count number of total vectors and blocks (same omega) */
  bs->nv = 0;
  bs->nb = 0;
  for(int m=1 ; m<=omega_max ; m+=2){
    int l = (m+1)/2;
    bool found = false;
    for(int is=-1 ; is<=1 ; is+=2){
      for(int ir=0 ; ir<=bs->nrmax ; ir++){
        for(int iz=0 ; iz<=bs->nzmax ; iz++){
          for(int il=l-1 ; il<=l ; il++){
            if((2*il+is) == m){
              double e1 = (2.0*ir + il + 1.0)*q13 + (iz + 0.5)/q23;
              if(e1 <= e0){
                bs->nv ++;
                found = true;
              }
            }
          }
        }
      }
    }
    if(found) bs->nb ++;
  }

  if( (bs->nv > 0) && (bs->nb > 0) ) bs->memalloc();
  else return(-1);

  int i  = 0;
  int nc = 0;
  bs->index[0] = 0;
  for(int m=1 ; m<=omega_max ; m+=2){
    int l = (m+1)/2;
    int k = 0;
    for(int is=-1 ; is<=1 ; is+=2){
      for(int ir=0 ; ir<=bs->nrmax ; ir++){
        for(int iz=0 ; iz<=bs->nzmax ; iz++){
          for(int il=l-1 ; il<=l ; il++){
            if((2*il+is) == m){
              double e1 = (2.0*ir + il + 1.0)*q13 + (iz + 0.5)/q23;
              if(e1 <= e0){
                int p = ((iz + 2*ir + il)%2 == 0) ? 1 : -1;
                bs->state[i].set(ir,iz,il,is,m,p);
                i++;
                k++;
              }
            }
          }
        }
      }
    }
    if(k > 0){
      bs->omega[nc] = m;                   // omega of the block
      bs->ndata[nc] = k;                   // number of states in the block
      bs->index[nc+1] = bs->index[nc] + k; // block starting number
      nc++;
    }
  }

#ifdef DEBUG
  std::cout << "# Basis Parameters"<< std::endl;
  std::cout << "#   N0 = " << std::setw(2) << bs->getN0() << std::endl;
  std::cout << "#   b  = " << std::setw(12) << bs->b << std::endl;
  std::cout << "#   q  = " << std::setw(12) << bs->q << std::endl;

  std::cout << "#  Nr max = " << std::setw(12) << bs->nrmax << std::endl;
  std::cout << "#  Nz max = " << std::setw(12) << bs->nzmax << std::endl;
  std::cout << "#  L  max = " << std::setw(12) << bs->lambdamax << std::endl;

  std::cout << "# Basis States       " << std::setw(4) << bs->nv << std::endl;
  std::cout << "# Diagonal Blocks    " << std::setw(4) << bs->nb << std::endl;
  std::cout << "# Size of Each Block ";
  for(int k=0 ; k<bs->nb ; k++) std::cout << std::setw(4) << bs->ndata[k];
  std::cout << std::endl;
  std::cout << "# Max Jz Value       " << std::setw(4) << bs->state[bs->nv-1].getJ() << "/2" << std::endl;


  i = 0;
  for(int k=0 ; k<bs->nb ; k++){
    std::cout << std::setw(4) << k;
    std::cout << "  N " << std::setw(4) << bs->ndata[k];
    std::cout << "  Omega " << std::setw(3) << bs->omega[k]<<"/2";
    std::cout << "  Index " << std::setw(4) << bs->index[k] << std::setw(4) << bs->index[k+1];
    std::cout << std::endl;

    for(int j=0 ; j<bs->ndata[k] ; j++){
      std::cout << " " << std::setw(4) << i;
      std::cout << " : nz" << std::setw(3) << bs->state[i].getNz();
      std::cout << " : nr" << std::setw(3) << bs->state[i].getNr();
      std::cout << " : ";
      std::cout << std::setw(3) << bs->state[i].getL();
      std::cout << std::setw(3) << bs->state[i].getS();
      std::cout << std::setw(3) << bs->state[i].getJ();
      std::cout << std::setw(3) << bs->state[i].getP() << std::endl;
      i++;
    }
  }
#endif
  return(0);
}



/**********************************************************/
/*      Conversion from Cylindrical to Spherical          */
/*      ----------                                        */
/*      L. Bonneau, et al. Phys. Rev. C 75, 054618 (2007) */
/**********************************************************/
void HOSphericalExpansion(Basis *dbase, SPEnergy *spe, Operator *H, HFInterface *hfsp)
{
  const int    sp_nmax =  5;
  const int    sp_lmax = 15;
  const double sp_emax = 10.0;
  const double eps     = 1.0e-03;
  const double epsv2   = 1.0e-05;

  /*** Spherical Harmonics Oscillator basis */
  SphericalBasis sbase;
  sbase.setLimit(sp_nmax, sp_lmax);
  HOSphericalBasis(&sbase);

  /*** for all single-particle state in Cylindrical HO */
  hfsp->n = 0;
  for(int k=0 ; k<dbase->nv ; k++){

    if(spe->energy[k] > sp_emax) break;
    if(k >= HF_MAX_SPSTATE) break;
    if(spe->v2[k] < epsv2) break;

    int k1 = spe->subindex[k];    // index inside block
    int kb = spe->block[k];       // block index

    int omega1 = dbase->omega[kb]; // 2*Omega in deformed HO

    double piave = 0.0;
    double pitot = 0.0;

    /*** states in Spherical HO */
    int nstate = 0, lmax = 0;
    for(int i=0 ; i < sbase.nv ; i++){

      if((sbase.getM(i) != omega1) || (sbase.getJ(i) < omega1)) continue;

      double sum = OverlapFunction(i,k1,kb,dbase,&sbase,H);

      double sumsq = sum*sum;
      piave += sumsq * ( (sbase.getL(i)%2 == 0) ? 1.0 : -1.0);
      pitot += sumsq;

      /*** save state information in the HFInterface object */
      if(fabs(sum) > eps){
        hfsp->state[k].set(sbase.getN(i),sbase.getL(i),sbase.getJ(i),sum);
        nstate ++;
        if(sbase.getL(i) > lmax) lmax = sbase.getL(i);
        if(nstate >= HF_MAX_EXPANSION) break;
      }
    }

    /*** average parity */
    spe->parity[k] = (pitot > 0.0) ? piave / pitot : 0.0;

    /*** save s.p. energy and v2 in the HFInterface object */
    hfsp->state[k].energy = spe->energy[k];
    hfsp->state[k].parity = (spe->parity[k] > 0.0) ? 1 : -1;
    hfsp->state[k].v2     = spe->v2[k];
    hfsp->state[k].k2     = omega1;
    hfsp->state[k].lmax   = lmax;
    hfsp->n++;
  }
  hfsp->b = 1.0/dbase->b;
}


/**********************************************************/
/*      Sum of Overlap Functions                          */
/**********************************************************/
double OverlapFunction(const int i, const int k1, const int kb, Basis *dbase, SphericalBasis *sbase, Operator *H)
{
  int omega0 = sbase->getM(i);
  int nperp0 = 2*sbase->getN(i) + sbase->getL(i);

  double sum = 0.0;
  for(int k2 = 0 ; k2<dbase->ndata[kb] ; k2++){
    int k0 = dbase->index[kb] + k2;

    int nz = dbase->state[k0].getNz();
    int nr = dbase->state[k0].getNr();
    int nl = dbase->state[k0].getL();  // Lambda
    int ns = dbase->state[k0].getS();  // 2 Sigma

    int nperp2 = 2*nr + nl;
    int dnperp = nperp0 - nperp2;

    if( (nz + nperp2 - nperp0)%2 !=0 ) continue;
    if(nperp2 > nperp0) continue;
    if(nl > sbase->getL(i)) continue;

    double x = sbase->getL(i) + 0.5 + sbase->getS(i) * ns * 0.5*omega0;
    double g = sqrt(x / (2.0*sbase->getL(i) + 1.0));
    if(sbase->getS(i) == -1) g *= -ns;

    double a = OverlapFunctionA(nz,dnperp,dbase->q);
    double b = OverlapFunctionB(dnperp,nperp2,nl,sbase->getN(i),sbase->getL(i));
    x = a * b * g * H[kb].vector[k1][k2];

#ifdef DEBUG
    if(x*x > 0.01){
       std::cout << "    " << std::setw(4) << k1;
       std::cout << std::setw(4) << nz << std::setw(4) << nperp2;
       std::cout << std::setw(4) << nl << std::setw(4) << ns;
       std::cout << std::setw(12) << H[kb].vector[k1][k2];
       std::cout << std::setw(4) << sbase->getN(i) << std::setw(4) << sbase->getL(i);
       std::cout << std::setw(4) << 2*(sbase->getJ(i) + sbase->getL(i))-3;
       std::cout << std::setw(12) << a << std::setw(12) << b;
       std::cout << std::setw(12) << g << std::setw(12) << x << std::endl;
    }
#endif

    sum += x;
  }
  return(sum);
}


/**********************************************************/
/*      Overlap Function Eq. (A2)                         */
/*      A_{n n'}(q)                                       */
/**********************************************************/
void HOSphericalBasis(SphericalBasis *sbase)
{
  /*** count total number of states */
  sbase->nv = 0;
  for(int l=0 ; l <=sbase->getLmax() ; l++){
    for(int s=-1 ; s<=1 ; s+=2){
      int j = 2*l + s;
      if(j < 0) continue;
      for(int m=1 ; m<=j ; m+=2){
        for(int n=0 ; n<=sbase->getNmax() ; n++){
          sbase->nv ++;
          continue;
        }
      }
    }
  }

  sbase->memalloc();

  int nc = 0;
  for(int l=0 ; l <=sbase->getLmax() ; l++){
    for(int s=-1 ; s<=1 ; s+=2){
      int j = 2*l + s;
      if(j < 0) continue;
      for(int m=1 ; m<=j ; m+=2){
        for(int n=0 ; n<=sbase->getNmax() ; n++){
          sbase->state[nc++].set(n,l,s,j,m);
        }
      }
    }
  }

#ifdef DEBUG
  for(int k2=0 ; k2<sbase->nv ; k2++){
    std::cout << std::setw(5) << k2;
    std::cout << std::setw(5) << sbase->state[k2].getN();
    std::cout << std::setw(5) << sbase->state[k2].getL();
    std::cout << std::setw(5) << sbase->state[k2].getS();
    std::cout << std::setw(5) << sbase->state[k2].getJ();
    std::cout << std::setw(5) << sbase->state[k2].getM() << std::endl;
  }
#endif
}


/**********************************************************/
/*      Overlap Function Eq. (A2)                         */
/*      A_{n n'}(q)                                       */
/**********************************************************/
double OverlapFunctionA(const int n1, const int n2, const double q)
{
  double a = 0.0;

  if((n1-n2)%2 != 0){ return(a); }

  if(fabs(q-1.0) <= 1.0e-9){
    a = 0.0;
    if(n1 == n2) a = 1.0;
    return(a);
  }

  int m_min = std::max(0,(n1 - n2) / 2);
  int m_max = (n1 - n1%2) / 2;

  double x1 = (q-1.0) / (4.0 * sqrt(q));

  double sum = 0.0;
  for(int m = m_min ; m<=m_max ; m++){
    double p  = (m%2 == 0) ? 1 : -1;
    double x2 = pow(x1,2*m);
    double x3 = exp( fact[m] + fact[m+(n2-n1)/2] + fact[n1-2*m] );
    sum += p*x2/x3;
  }

  double c1 = sqrt( pow(2.0, (double)(n1-n2)) * exp(fact[n1] + fact[n2]) );
  double c2 = pow( (q-1.0) / (q+1.0), (n2-n1)/2.0 );
  double c3 = pow( 2.0 * sqrt(q) / (q+1.0), n1+0.5 );

  a = c1 * c2 * c3 * sum;

  return(a);
}


/**********************************************************/
/*      Overlap Function Eq.(A5)                          */
/*      B_{nz,nperp,Lambda;n,l}                           */
/**********************************************************/
double OverlapFunctionB(const int nz, const int np, const int lambda, const int n, const int l)
{
  double b = 0.0;

  if( (np-lambda)%2 != 0 ) return(b);
  if( (nz+np) != (2*n+l) ) return(b);

  int alpha = (np + lambda)/2; // Eq. (A6)
  int beta  = (np - lambda)/2; // Eq. (A7)

  double sum = 0.0;
  for(int k=0 ; k<=l ; k++){
    int k2 = k + n - beta;
    int k1 = nz - 2*k2;

    if(k1<0 || k2<0) continue;

    double p  = (k%2 == 0) ? 1 : -1;
    double x1 = exp(fact[2*(l-k)] + fact[k+n]);
    double x2 = exp(fact[k] + fact[l-k] + fact[k1] + fact[k2]);
    sum += p * x1 / x2;
  }

  double c1 = ((alpha+n)%2 == 0) ? 1.0 : -1.0;
  double c2 = pow(2.0, (double)n);
  double c3 = exp(fact[alpha] - fact[beta]);
  double c4 = exp(fact[nz] + fact[l-lambda] + fact[n+l]) * (2*l + 1.0);
  double c5 = pow(2.0, (double)nz);
  double c6 = exp(fact[l+lambda] + fact[n] + fact[2*(n+l)+1]); 
  b = sum * c1 * c2 * sqrt(c3 * c4 / (c5 * c6));

  return(b);
}



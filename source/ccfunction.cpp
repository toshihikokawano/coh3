/******************************************************************************/
/*  ccfunc.cpp                                                                */
/*        special functions in coupled channels calculation                   */
/******************************************************************************/

#include <cmath>

#include "optical.h"
#include "coupling.h"

static inline int sign(double x){
  return( (((x)<0.0) ? -1: 1) );
}

static double ccGeometricalFactor (const int, const int, Chspin *, Chspin *);
static void   ccCouplingMatrix (const NuclearModel, const int, const int, double **, CCdata *, Collective *);

/**********************************************************/
/*      Maximum possible number of channels               */
/*      Tamura, Rev. Mod. Phys., 37, 679 (1965), Eq.(18)  */
/**********************************************************/
int ccNumberOfChannels(Collective *col)
{
  int n = 0;
  int i = (int)(2.0*col->pspin);
  for(int k=0 ; k<col->nlevel ; k++){
    int    p = ( (int)(2*col->lev[k].spin)%2 !=0 ) ? 1:0;
    double s = fabs(col->lev[k].spin);
    switch(i){
    case 0: if(p) n += (int)(s+0.5);
            else  n += (int)(s+1.0);
            break;
    case 1:       n += (int)(2.0*s+1.0);
            break;
    case 2: if(p) n += (int)(3.0*s+1.5);
            else  n += (int)(3.0*s+2.0);
            break;
    default:      n += (int)(2.0*s+1.0);
            break;
    }
  }
  return(n);
}


/**********************************************************/
/*      Setup Channels Those Coupled to the Same J Pi     */
/**********************************************************/
int ccSetChannel(const int jj, const int p, const int lmax, Collective *col, CCdata *cdt)
{
  int n = 0;
  int pi = sign(2*col->pspin);

  for(int k=0 ; k<col->nlevel ; k++){
    int pt = sign(2*col->lev[k].spin);
    int st = (int)fabs(2*col->lev[k].spin);

    for(int l=((p == pi*pt) ? 0 : 1) ; l<=lmax ; l+=2){
      for(int j=0 ; j<(int)(2*col->pspin + 1) ; j++){
        int j2=(int)(2.0*(l + col->pspin - j));
        if( j2 < 0 ) continue;

        int jjmin1 = std::abs(j2-st);
        int jjmax1 =          j2+st ;

        if( (jjmin1 <= jj) && (jj <= jjmax1)){
          cdt[n].level   = k;
          cdt[n].chn.si2 = 2*col->pspin;
          cdt[n].chn.st2 = st;
          cdt[n].chn.pi  = pi;
          cdt[n].chn.pt  = pt;
          cdt[n].chn.l   = l;
          cdt[n].chn.j2  = j2;
          cdt[n].lev     = &col->lev[k];
          cdt[n].xl_term = l*(l+1.0);
          cdt[n].so_term = j2/2.0*(j2/2.0+1.0) - cdt[n].xl_term - col->pspin*(col->pspin+1.0);

          cdt[n].open    = (col->lev[k].wave_number > 0.0) ? true : false;

          n++;
        }
      }
    }
   }
   return(n);
}


/**********************************************************/
/*      Nuclear Matrix Elements for Rotational Nucleus    */
/*      Tamura Eq.(40), (41)                              */
/**********************************************************/
void ccMatrixElementRot(const int k2, Collective *col)
{
  for(int i=0 ; i<col->nlevel ; i++){
    int si2 = (int)(2.0*fabs(col->lev[i].spin));
    double yi = sqrt(si2 + 1.0);

    for(int j=0 ; j<=i ; j++){
      int sj2 = (int)(2.0*fabs(col->lev[j].spin));
      double yj = sqrt(sj2 + 1.0);

      int ij = i*(i+1)/2+j;

      /*** lambda is limited to even numbers */
      for(int l=0 ; l<MAX_LAMBDA ; l++){
        int lambda = 2*l;
        col->qmatrix[l][ij] = clebsh_gordan(si2, 2*lambda, k2, 0, sj2) * yi / yj;
      }
    }
  }
}


/**********************************************************/
/*      Nuclear Matrix Elements for Vibrational Nucleus   */
/*      Tamura Eq.(36), (37)                              */
/**********************************************************/
void ccMatrixElementVib(Collective *col)
{
  for(int i=0 ; i<col->nlevel ; i++){
    int pi  = col->lev[i].phonon;
    int si2 = (int)(2.0*fabs(col->lev[i].spin));
    int pti = sign(col->lev[i].spin);

    for(int j=0 ; j<=i ; j++){
      int pj  = col->lev[j].phonon;
      int sj2 = (int)(2.0*fabs(col->lev[j].spin));
      double yj = sqrt(sj2 + 1.0);

      int ij = i*(i+1)/2+j;
      col->pindex[ij] = 0;

      for(int lambda=0 ; lambda<MAX_LAMBDA ; lambda++){
        int lambda2 = lambda * 2;
        double q = 0.0;
        /*** ground state : one-phonon */
        if((pj == 0) && (pi == 1)){
          if((lambda == 2) && (si2 ==  4)){
            q = col->beta[0];
            col->pindex[ij] = 1;
          }
          else if((lambda == 3) && (si2*pti == -6)){
            q = -col->beta[1];
            col->pindex[ij] = 1;
          }
        }

        /*** one-quadrupole-phonon : two-quadrupole-phonon */
        else if((pj == 1) && (pi == 2)){
          if((sj2 == 4) && (lambda == 2)){
            q = col->beta[0] * sqrt(2.0*(si2+1.0)/5.0);
            col->pindex[ij] = 1;
          }
        }

        /*** ground state : two-quadrupole-phonon */
        else if((pj == 0) && (pi == 2)){
          if(lambda2 == si2){
            q = col->beta[0] * col->beta[0] * clebsh_gordan(4,4,0,0,lambda2) / sqrt(PI2);
            col->pindex[ij] = 2;
          }
        }
        /*** one-quadrupole-phonon : one-octupole-two-phonon */
        else if((pj == 1) && (pi == 1)){
          if((sj2 == 4) && (si2*pti == -6)){
            q = col->beta[0] * col->beta[1] * clebsh_gordan(4,6,0,0,lambda2) / sqrt(PI4);
            col->pindex[ij] = 2;
          }
        }
        col->qmatrix[lambda][ij] = q / yj;
      }
    }
  }
}


/**********************************************************/
/*      Coupling Potential                                */
/**********************************************************/
void ccCouplingPotential(const NuclearModel model, const int m, const int jj, double **cv, std::complex<double> **Vpot, CCdata *cdt, Collective *col, Potential *pc0, Potential **pc1, Potential **pc2)
{
  std::complex<double> vi(0.0,0.0), vj(0.0,0.0);

  /*** A x Q, Tamura Eq. (27) */
  ccCouplingMatrix(model,m,jj,cv,cdt,col);

  for(int k=0 ; k<=pc0[0].n_match ; k++){
    for(int i=0 ; i<m ; i++){    int ni = cdt[i].level;
      for(int j=0 ; j<=i ; j++){ int nj = cdt[j].level;

        int id1 =  i*( i+1)/2 +  j;
        int id2 = ni*(ni+1)/2 + nj;

        Vpot[k][id1] = std::complex<double>(0.0,0.0);

        for(int l=1 ; l<MAX_LAMBDA ; l++){
          if(cv[l][id1] == 0.0) continue;

          /*** rotational model, Tamura Eq. (14),(15) */
          if(model == rotation){
            vi = pc1[ni][l].mean_field[k];
            vj = pc1[nj][l].mean_field[k];
          }
          /*** vibrational model, Tamura Eq. (13) */
          else{
            if(col->pindex[id2] == 2){       // second order
              vi = pc2[ni][l].mean_field[k];
              vj = pc2[nj][l].mean_field[k];
            }
            else{                            // first order
              vi = pc1[ni][l].mean_field[k];
              vj = pc1[nj][l].mean_field[k];
            }
          }
          Vpot[k][id1] -= 0.5*(vi+vj) * cv[l][id1];

        }
      }

      /*** diagonal elements */
      int id1 = i*(i+1)/2+i;
      Vpot[k][id1] += - pc0[ni].mean_field[k]
                     + cdt[i].xl_term * pc0[ni].r2inv[k] - cdt[i].lev->wavesq
                     - cdt[i].so_term * pc0[ni].spin_orbit[k];
    }
  }
}


/**********************************************************/
/*      Set Coupling Matrix                               */
/**********************************************************/
void ccCouplingMatrix(const NuclearModel model, const int m, const int jj, double **cv, CCdata *cdt, Collective *col)
{
  for(int i=0 ; i<m ; i++){    int ni = cdt[i].level;
    for(int j=0 ; j<=i ; j++){ int nj = cdt[j].level; 
      int id1 =  i*( i+1)/2 +  j;
      int id2 = ni*(ni+1)/2 + nj;

      for(int l=0 ; l<MAX_LAMBDA ; l++) cv[l][id1] = 0.0;
      for(int l=0 ; l<MAX_LAMBDA ; l++){
        int lambda = l;
        if(model == rotation) lambda *= 2;
        if(col->qmatrix[l][id2] == 0.0) continue;

        double a = ccGeometricalFactor(jj,lambda,&cdt[i].chn,&cdt[j].chn);
        cv[l][id1] = a * col->qmatrix[l][id2];
      }
    }
  }
}


/**********************************************************/
/*      Channel Coupling Geometrical Coefficients         */
/*      Satchler Eq.(5.59), Tamura Eq.(28)                */
/**********************************************************/
double ccGeometricalFactor(const int jj, const int lambda, Chspin *ci, Chspin *cj)
{
  /*** 
       (-)^{J-I0-s} i^{l0-l1-lambda}
       <l0,l1,0,0 | lambda,0>
       W(j0,j1,l0,l1;lambda,s) = W(l0,l1,j0,j1;lambda,s) = W(l0,j0,l1,j1;s,lambda)
       W(j0,j1,I0,I1;lambda,J) = W(j0,I0,j1,I1;J,lambda)
  */
  int lambda2 = 2 * lambda;
  int li2 = 2 * ci->l;
  int lj2 = 2 * cj->l;
  double xi = sqrt((li2 + 1.0) * (ci->j2 + 1.0));
  double xj = sqrt((lj2 + 1.0) * (cj->j2 + 1.0));
  double yj = sqrt(cj->st2 + 1.0);

  /*** parity(L) = i^L, parity(2J) = i^{2J} = (-)^J  */

  double a = parity(jj - ci->st2 - ci->si2) * parity(ci->l - cj->l - lambda)
           * clebsh_gordan(li2, lj2, 0, 0, lambda2)
           * racah(ci->j2, cj->j2, li2, lj2, lambda2, ci->si2)
           * racah(ci->j2, cj->j2, ci->st2, cj->st2, lambda2, jj)
           * xi * xj * yj
           * parity(lambda); // this phase needed to make Satchler = Tamura

  /*** Tamura Eq. (28), with a factor of sqrt{2I'+1} */
/*
  double a = parity(jj - cj->st2 - ci->si2 + ci->l + 3*cj->l)
           * clebsh_gordan(li2, lj2, 0, 0, lambda2)
           * racah(ci->j2, ci->st2, cj->j2, cj->st2, jj, lambda2)
           * racah(li2, ci->j2, lj2, cj->j2, ci->si2, lambda2)
           * xi * xj * yj;
*/
  
  return(a / sqrt(PI4));
}


/**********************************************************/
/*      Coulomb Functions for Coupled-Equations           */
/**********************************************************/
int ccExternalFunction(const double rm, Wavefunc *wfn, Collective *col)
{
   int    lmax0=0, j=0, jmax=0;

   for(int k=0 ; k<col->nlevel ; k++){
     j = (int)(2.0*(fabs(col->lev[k].spin)+col->pspin)); if(jmax < j) jmax=j;
   }
   lmax0 = ((j%2) == 0) ? j/2 : j/2+1;

   double rhom = rm*col->lev[0].wave_number;
   lmax0 = omExternalFunction(lmax0,rhom,col->lev[0].coulomb.eta,&wfn[0]);

   for(int k=1 ; k<col->nlevel ; k++){
     rhom = rm*col->lev[k].wave_number;

     if(rhom < 0.0){
       omExternalClosed(lmax0,-rhom,col->lev[k].coulomb.eta,&wfn[k]);
     }
     else{
       omExternalFunction(lmax0,rhom,col->lev[k].coulomb.eta,&wfn[k]);
     }
   }

   return(lmax0);
}


/**********************************************************/
/*      Check if Target State is in Rotatioanl Band       */
/**********************************************************/
int ccFindTargetLevel(const double ex, Collective *col)
{
  int  kt = -1;
  for(int k=0 ; k<col->nlevel ; k++){
    if( col->lev[k].excitation == ex ){
      col->target_index = kt = k;
      break;
    }
  }
  return(kt);
}


/**********************************************************/
/*      K Number in the Rotatioanl Band, Smallest Spin    */
/**********************************************************/
int ccBandheadSpin(Collective *col)
{
   double st0 = fabs(col->lev[0].spin);
   for(int k=1 ; k<col->nlevel ; k++){
     double st1 = fabs(col->lev[k].spin);
     if( st1 < st0 ){ st0=st1; }
   }
   return((int)(2*st0));
}



/******************************************************************************/
/**     POLYCALC, Solution of Least-Squares Equation                         **/
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "polysq.h"

const double EPS = 1.0e-72;

const int ERR_NEG = -1;  /* Matrix not positiv */
const int ERR_PIV = -2;  /* Pivot zero         */
const int ERR_INP = -3;  /* Data error         */

static int inverse         (double *, const int);
static int matrix_choleski (double *, const int);
static int matrix_inverse  (double *, const int);


/**********************************************************/
/*      Solve Least-Squares Equation                      */
/*                                                        */
/*   n: number of data points                             */
/*   m: number of parameters                              */
/*   p: prameter vector [m]                               */
/*   x: prameter covariance [m*(m+1)/2]                   */
/*   y: data vector [n]                                   */
/*   v: data covariance [n*(n+1)/2]                       */
/*   f: sensitivity matrix [n*m]                          */
/**********************************************************/
double LSQCalc(const int n, const int m,  double *y, double *p, double *v, double *x, double *f)
{
  double *work1,*work2;

  // V = V^{-1}
  if( (inverse(v,n)) != 0 ) return(-1.0);

  work1 = new double [m*n];
  work2 = new double [n*m];

  // W = F' V^{-1}
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<n ; j++){
      work1[i*n+j] = 0.0;
      for(int k=0 ; k<n ; k++){
        int jk = (j<=k) ? k*(k+1)/2+j : j*(j+1)/2+k;
        work1[i*n+j] +=  f[k*m+i] * v[jk];
      }
    }
  }

  // X = W F = F' V^{-1} F
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<=i ; j++){
      int ij = (i+1)*i/2+j;
      x[ij] = 0.0;
      for(int k=0 ; k<n ; k++) x[ij] += work1[i*n+k] * f[k*m+j];
    }
  }

  // X = X^{-1}
  if( (inverse(x,m)) != 0 ){
    delete [] work1;
    delete [] work2;
    return(-1.0);
  }

  // 
  for(int i=0 ; i<m ; i++){
    work2[i] = 0.0;
    for(int j=0 ; j<n ; j++) work2[i] += work1[i*n+j]*y[j];
  }

  for(int i=0 ; i<m ; i++){
    p[i] = 0.0;
    for(int j=0 ; j<m ; j++){
      int ij = (j<=i) ? i*(i+1)/2+j : j*(j+1)/2+i;
      p[i] += x[ij]*work2[j];
    }
  }

  delete [] work1;
  delete [] work2;

  return(0.0);
}


/**********************************************************/
/*      Solve Generalized Least-Squares Equation          */
/*                                                        */
/*  p1 = p0 + X0 F'(F X0 F' + V)^(-1) {y-f(p0)}           */
/*  X1 = X0 - X0 F'(F X0 F' + V) F X0                     */
/*                                                        */
/*   n: number of data points                             */
/*   m: number of parameters                              */
/*  p0: prior prameter vector [m]                         */
/*  p1: posterior prameter vector [m]                     */
/*   x: prameter covariance [m*(m+1)/2]                   */
/*   y: data vector [n]                                   */
/*   v: data covariance [n*(n+1)/2]                       */
/*   f: sensitivity matrix [n*m]                          */
/*      (0,0)   (0,1)   ... (0,m-1)                       */
/*      (1,0)   (1,1)   ... (1,m-1)                       */
/*      (n-1,0) (n-1,1) ... (n-1,m-1)                     */
/**********************************************************/
int LSQGeneral(const int n, const int m, double *p0, double *p1, double *x, double *y, double *v, double *f)
{
  double *w, xx;
  w = new double [m*n];

  /***  W = X F' */
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<n ; j++){
      xx = 0.0;
      for(int k=0 ; k<m ; k++){
        xx += x[((k<=i) ? i*(i+1)/2+k : k*(k+1)/2+i)]*f[j*m+k];
      }
      w[i*n+j] = xx;
    }
  }
  /*** V <- F X F'+ V =F W + V */
  for(int i=0 ;  i<n ; i++){
    for(int j=0 ; j<=i ; j++){
      xx = 0.0;
      for(int k=0 ; k<m ; k++) xx += f[i*m+k] * w[k*n+j];
      v[i*(i+1)/2+j]+=xx;
    }
  }
  /*** V <- V^(-1) */
  if(inverse(v,n) != 0) return -1;

  /*** F = W V^(-1) */
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<n ; j++){
      xx = 0.0;
      for(int k=0 ; k<n ; k++){
        xx += w[i*n+k]*v[((k<=j) ? j*(j+1)/2+k : k*(k+1)/2+j)];
      }
      f[i*n+j] = xx;
    }
  }

  /*** p1 = p0 + F y */
  for(int i=0 ; i<m ; i++){
    xx = 0.0;
    for(int j=0 ; j<n ; j++) xx += f[i*n+j] * y[j];
    p1[i] = p0[i] + xx;
  }

  /*** X0 <- X0 -  F W' */
  for(int i=0 ; i<m ; i++){
    for(int j=0 ; j<=i ; j++){
      xx = 0.0;
      for(int k=0 ; k<n ; k++) xx += f[i*n+k] * w[j*n+k];
      x[i*(i+1)/2+j] -= xx;
    }
  }

   delete [] w;
   return 0;
}


/**********************************************************/
/*      Levenberg-Marquardt Method                        */
/*                                                        */
/*   d: damping facotr                                    */
/*  p0: initial prameter vector [m]                       */
/*  p1: output vector [m]                                 */
/*   y: data vector [n]                                   */
/*   w: error vector [n]                                  */
/*   z: Jacobian matrix [n*m]                             */
/*      (0,0)   (0,1)   ... (0,m-1)                       */
/*      (1,0)   (1,1)   ... (1,m-1)                       */
/*      (n-1,0) (n-1,1) ... (n-1,m-1)                     */
/**********************************************************/
int LSQMarquardt(const int n, const int m, const double d, double *p0, double *p1, double *y, double *w, double *z)
{
  double *v;
  v = new double [m*(m+1)/2];

  /*** J' W J */
  for(int j0=0 ; j0<m ; j0++){
    for(int j1=0 ; j1<=j0 ; j1++){
      double x = 0.0;
      for(int i=0 ; i<n ; i++) x += z[j0+i*m] * z[j1+i*m] * w[i];
      v[j0*(j0+1)/2 + j1] = x;
    }
  }

  /*** add lambda diag(J' J), and invert */
  for(int j0=0 ; j0<m ; j0++) v[j0*(j0+1)/2 + j0] *= 1.0 + d;

  if(inverse(v,m) != 0) return -1;

  /*** [J' J + lambda x diag (J' J)]^{-1} J' */
  for(int j0=0 ; j0<m ; j0++){

    p1[j0] = p0[j0];
    for(int i=0 ; i<n ; i++){
      double x = 0.0;
      for(int j1=0 ; j1<m ; j1++){
        int k = (j1 <= j0) ? j0*(j0+1)/2+j1 : j1*(j1+1)/2+j0;
        x += v[k]*z[j1+i*m];
      }
      p1[j0] += x*y[i]*w[i];
    }
  }

  delete [] v;
  return 0;
}


/**********************************************************/
/*      Matrix Inversion                                  */
/**********************************************************/
int inverse(double *a, const int n)
{
   if(n < 1) return(ERR_INP);

   int c = matrix_choleski(a,n);

   if(c != 0) return(c);
   else c = matrix_inverse(a,n);

   return(c);
}

int matrix_choleski(double *a, const int n)
{
   double  x1,x2,x4,x3,sum;
   int c=0,n2,l,i,jk,ik,i1,jj,ij,lj,j;

   if(a[0] == 0.0) return(ERR_PIV);
   if(a[0] < 0.0) c=ERR_NEG;
   a[0]=1.0/a[0];
   if(n == 1) return(c);

   x1=a[0]*a[1];
   x2=a[1]*x1;
   a[1]=x1;
   x4=a[2];
   x3=x4-x2;
   if( fabs(x3) < fabs(x4*EPS) ) return(ERR_PIV);
   if (x3<0.0) c=ERR_NEG;
   a[2]=1.0/x3;
   if (n==2) return(c);

   n2=n-2;
   l=4;
   for(i=1;i<=n2;i++){
       jk=2;
       ik=l;
       for(j=1;j<=i;j++){
           sum=0.0;
           lj=l+j-1;
           for(ij=l;ij<=lj;ij++){
               sum+=a[ij-1]*a[jk-1];
               jk++;
           }
           jk++;
           ik++;
           a[ik-1]-=sum;
       }
       i1=i+2;
       jj=1;
       ij=l;
       sum=0.0;
       for(j=2;j<=i1;j++){
           x1=a[ij-1]*a[jj-1];
           sum+=a[ij-1]*x1;
           a[ij-1]=x1;
           jj+=j;
           ij++;
       }
       x4=a[jj-1];
       x3=x4-sum;
       if ( fabs(x3)<fabs(x4*EPS) ) return(ERR_PIV);
       if(x3 < 0.0) c=ERR_NEG;
       a[jj-1]=1.0/x3;
       l=jj+1;
   }
   return(c);
}


int matrix_inverse(double *a, const int n)
{
   double sum,x1;
   int i,ij,ije,ijsp,kj,kjs,j,kjsp,ik,ki,k,ii,i1;

   int n1 = 0;
   int c  = 0;

   if(a[0] < 0) c=ERR_NEG;
   if(n == 1) return(c);

   a[1]=-a[1];
   if (n!=2){
       n1=n-1;
       ije= 2;
       for(i=2;i<=n1;i++){
           ij=ije+2;
           ije+=i+1;
           kjs=0;
           for(j=2;j<=i;j++){
               sum=a[ij-1];
               ij++;
               kjsp=j;
               kjs+=kjsp;
               kj= kjs;
               for(ik=ij;ik<=ije;ik++){
                   sum+=a[ik-1]*a[kj-1];
                   kj+=kjsp;
                   kjsp++;
               }
               a[kj-1]=-sum;
           }
           a[ije-1]=-a[ije-1];
       }
   }

   ii= 1;
   ki= 2;
   sum=a[0];
   for(k=2;k<=n;k++){
       ii+=k;
       if(a[ii-1] < 0) c=ERR_NEG;
       x1=a[ki-1]*a[ii-1];
       sum+=a[ki-1]*x1;
       a[ki-1]=x1;
       ki+=k;
   }
   a[0]=sum;
   if(n == 2) return(c);
   ij=1;
   for(i=2;i<=n1;i++){
       i1=i+1;
       for(j=2;j<=i;j++){
           ij++;
           kj=ij;
           sum=a[kj-1];
           ijsp=i1-j;
           for(k=i;k<=n1;k++){
               kj+=k;
               ki=kj+ijsp;
               sum+=a[kj-1]*a[ki-1];
           }
           a[ij-1]=sum;
       }
       ij++;
       ki=ij+i;
       ii=ij;
       sum=a[ii-1];
       for(k=i1;k<=n;k++){
           ii+=k;
           x1=a[ki-1]*a[ii-1];
           sum+=a[ki-1]*x1;
           a[ki-1]=x1;
           ki+=k;
       }
       a[ij-1]=sum;
   }
   return(c);
}


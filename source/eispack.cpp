/******************************************************************************/
/**     Subroutines from EISPACK                                             **/
/**                                                                          **/
/**     Original EISPACK FORTRAN routines developed at ANL                   **/
/**     Four sub-routines to diagonalize real diagonal                       **/
/**     and complex Hermite matrices were converted into C++                 **/
/**     They include tred2, tqli, tql2, and htribk                           **/
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "eispack.h"


static int EISPACK_tred2 (const int, double **, double *, double *);
static int EISPACK_tqli (const int, double *, double *, double **);

static int EISPACK_htridi (const int, double **, double **, double *, double *, double *, double **);
static int EISPACK_tql2 (const int, double *, double *, double **, int);
static int EISPACK_htribk (const int, double **, double **, double **, double **, double **);


int EISPACKUnitaryDiag(const int n, double *a,  double *w, double **z)
{
  double **v, *d, *e;
  try{
    e  = new double [n+1];
    v  = new double * [n+1];
    for(int i=0 ; i<n+1 ; i++){
      v[i] = new double [n+1];
    }
  }
  catch(std::bad_alloc &e){ return(-1); }

  d = &w[-1];

  for(int i=1 ; i<=n ; i++){
    for(int j=1 ; j<=n ; j++){
      int k = (j <= i) ? i*(i-1)/2 + j-1 : j*(j-1)/2 + i-1;
      v[i][j] =  a[k];
    }
    d[i] = 0.0;
  }

  EISPACK_tred2(n,v,d,e);
  EISPACK_tqli(n,d,e,v);

  for(int i=0 ; i<n ; i++){
    for(int j=0 ; j<n ; j++)  z[j][i] = v[i+1][j+1];
  }

  for(int i=0 ; i<n+1 ; i++){
    delete [] v[i];
  }
  delete [] v;
  delete [] e;

  return(0);
}


int EISPACKHermiteDiag(const int n, std::complex<double> *a, double *w, std::complex<double> *z)
{
  double **ar, **ai, *f1, *f2, *f3[2], **zr, **zi, *x;

  try{
    ar = new double * [n+1];
    ai = new double * [n+1];
    zr = new double * [n+1];
    zi = new double * [n+1];
    f1 = new double [n+1];
    f2 = new double [n+1];
    x  = new double [n+1];
    for(int i=0 ; i<n+1 ; i++){
      ar[i] = new double [n+1];
      ai[i] = new double [n+1];
      zr[i] = new double [n+1];
      zi[i] = new double [n+1];
    }
    f3[0] = new double [n+1];
    f3[1] = new double [n+1];
  }
  catch(std::bad_alloc &e){ return(-1); }

  for(int i=0 ; i<n ; i++){
    int p = i+1;
    for(int j=0 ; j<n ; j++){
      int ij = i*n+j;
      int q = j+1;
      ar[p][q] = a[ij].real();
      ai[p][q] = a[ij].imag();
      zr[p][q] = (i == j) ? 1.0 : 0.0;
      zi[p][q] = 0.0;
    }
  }


  EISPACK_htridi(n,ar,ai,x,f1,f2,f3);
  if(EISPACK_tql2(n,x,f1,zr,1) != 0) return(-1);
  EISPACK_htribk(n,ar,ai,f3,zr,zi);

  for(int i=1 ; i<=n ; i++){
    w[i-1] = x[i];
    for(int j=1 ; j<=n ; j++){
      int ij = (i-1)*n+j-1;
      z[ij] = std::complex<double>(zr[i][j],-zi[i][j]);
    }
  }

  for(int i=0 ; i<n+1 ; i++){
    delete [] ar[i];
    delete [] ai[i];
    delete [] zr[i];
    delete [] zi[i];
  }
  delete [] ar;
  delete [] ai;
  delete [] zr;
  delete [] zi;
  delete [] f1;
  delete [] f2;
  delete [] f3[0];
  delete [] f3[1];
  delete [] x;

  return(0);
}


int EISPACK_tred2(const int n, double **a, double *d, double *e)
{
  int l,k,j,i;
  double scale,hh,h,g,f;

  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)  scale += fabs(a[i][k]);
      if (scale == 0.0)   e[i]=a[i][l];
      else {
        for (k=1;k<=l;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g = f>0 ? -sqrt(h) : sqrt(h);
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=1;j<=l;j++) {
          /* Next statement can be omitted if eigenvectors not wanted */
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=1;k<=j;k++)    g += a[j][k]*a[i][k];
          for (k=j+1;k<=l;k++)  g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=1;j<=l;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=1;k<=j;k++)  a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  /* Next statement can be omitted if eigenvectors not wanted */
  d[1]=0.0;
  e[1]=0.0;
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i]=a[i][i]; */
  for (i=1;i<=n;i++) {
    l=i-1;
    if (d[i]) {
      for (j=1;j<=l;j++) {
        g=0.0;
        for (k=1;k<=l;k++)   g += a[i][k]*a[k][j];
        for (k=1;k<=l;k++)   a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
  }
  return(0);
}


int EISPACK_tqli(const int n, double *d, double *e, double **z)
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;

  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
        dd = fabs(d[m]) + fabs(d[m+1]);
        if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
        if (iter++ == 30){
          std::cerr << "Too many iterations in TQLI";
          return(-1);
        }
        g=(d[l+1]-d[l])/(2.0*e[l]);
        r=sqrt((g*g)+1.0);
        g=d[m]-d[l]+e[l]/(g+ r*((g>=0.0) ? 1.0 : -1.0));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=sqrt(f*f+g*g));
          if (r==0.0) {
            d[i+1] -= p;
            e[m]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;
          /* Next loop can be omitted if eigenvectors not wanted */
          for (k=1;k<=n;k++) {
            f=z[k][i+1];
            z[k][i+1]=s*z[k][i]+c*f;
            z[k][i]=c*z[k][i]-s*f;
          }
        }
        d[l] -= p;
        e[l]=g;
        e[m]=0.0;
      }
    } while (m != l);
  }
  return(0);
}


int EISPACK_htridi(const int n, double **ar, double **ai, double *d, double *e, double *e2, double **tau)
{
  double scale,f,fi,g,gi,h,hh,si;

  tau[0][n] = 1.0;
  tau[1][n] = 0.0;
  for(int i=1 ; i<=n ; i++) d[i] = ar[i][i];

  for(int ii=1 ; ii<=n ; ii++){

    int i = n + 1 - ii;
    int l = i - 1;
    h = 0.0;
    scale = 0.0;

    if(l >= 1){
      for(int k = 1 ;  k<=l ; k++){
        scale += fabs(ar[i][k]) + fabs(ai[i][k]);
      }

      if(scale != 0.0) goto L140;
      tau[0][l] = 1.0;
      tau[1][l] = 0.0;
    }

    e[i] = 0.0;
    e2[i] = 0.0;
    goto L290;
      
L140:
    for(int k = 1 ; k<=l ; k++){
      ar[i][k] = ar[i][k] / scale;
      ai[i][k] = ai[i][k] / scale;
      h += ar[i][k] * ar[i][k] + ai[i][k] * ai[i][k];
    }
    e2[i] = scale * scale * h;
    g = sqrt(h);
    e[i] = scale * g;
    f = sqrt(ar[i][l]*ar[i][l] + ai[i][l]*ai[i][l]);

    if(f == 0.0){
      tau[0][l] = -tau[0][i];
      si = tau[1][i];
      ar[i][l] = g;
    }
    else{
      tau[0][l] = (ai[i][l] * tau[1][i] - ar[i][l] * tau[0][i]) / f;
      si = (ar[i][l] * tau[1][i] + ai[i][l] * tau[0][i]) / f;
      h += f * g;
      g = 1.0 + g / f;
      ar[i][l] = g * ar[i][l];
      ai[i][l] = g * ai[i][l];

      if(l == 1) goto L270;

      f = 0.0;
    }

    for(int j=1 ; j<=l ; j++){
      g  = 0.0;
      gi = 0.0;
      for(int k=1 ; k<=j ; k++){
        g  +=   ar[j][k] * ar[i][k] + ai[j][k] * ai[i][k];
        gi += - ar[j][k] * ai[i][k] + ai[j][k] * ar[i][k];
      }

      int jp1 = j + 1;
      if(l >= jp1){
        for(int k=jp1 ; k<=l ; k++){
          g  +=   ar[k][j] * ar[i][k] - ai[k][j] * ai[i][k];
          gi += - ar[k][j] * ai[i][k] - ai[k][j] * ar[i][k];
        }
      }

      e[j] = g / h;
      tau[1][j] = gi / h;
      f += e[j] * ar[i][j] - tau[1][j] * ai[i][j];
    }

    hh = f / (h + h);
    for(int j=1 ; j<=l ; j++){
      f = ar[i][j];
      g = e[j] - hh * f;
      e[j] = g;
      fi = -ai[i][j];
      gi = tau[1][j] - hh * fi;
      tau[1][j] = -gi;

      for(int k=1 ; k<=j ; k++){
        ar[j][k] = ar[j][k] - f*e[k] - g*ar[i][k] + fi*tau[1][k] + gi*ai[i][k];
        ai[j][k] = ai[j][k] - f*tau[1][k] - g*ai[i][k] - fi*e[k] - gi*ar[i][k];
      }
    }
L270:
    for(int k=1 ; k<=l ; k++){
      ar[i][k] = scale * ar[i][k];
      ai[i][k] = scale * ai[i][k];
    }
    tau[1][l] = -si;
L290:
    hh = d[i];
    d[i] = ar[i][i];
    ar[i][i] = hh;
    ai[i][i] = scale * sqrt(h);
  }

  return(0);
}


int EISPACK_tql2(const int n, double *d, double *e, double **z, int job)
{
  int    m=0, mml=0, l1=0, l2=0;
  double c, f, g, h,  p, r, s, c2, c3, s2,  dl1, el1,  tst1, tst2;
  s2 = c3 = 0.0; // to make compiler quiet

  int ierr = 0;
  if(n == 1) return(ierr);

  for (int i=2 ; i<=n ; i++)  e[i-1] = e[i];

  f = 0.0;
  tst1 = tst2 = 0.0;
  e[n] = 0.0;

  for(int l=1 ; l<=n ; l++){

    int j = 0;
    h = fabs(d[l]) + fabs(e[l]);
    if(tst1 < h) tst1 = h;

    for(m=l ; m<=n ; m++){
      tst2 = tst1 + fabs(e[m]);
      if (tst2 == tst1) break;
    }

    if(m == l){
      d[l] += f;
      continue;
    }

    do{
      if(j == 30){
        ierr = l;
        return(ierr);
      }
      j++;

      l1 = l + 1;
      l2 = l1 + 1;
      g = d[l];
      p = (d[l1] - g) / (2.0 * e[l]);
      r = sqrt(p*p + 1.0);
      d[l]  = e[l] / (p + r*((p>=0.0) ? 1.0 : -1.0));
      d[l1] = e[l] * (p + r*((p>=0.0) ? 1.0 : -1.0));
      dl1 = d[l1];
      h = g - d[l];

      if(l2 <= n){
        for(int i=l2 ; i<=n ; i++)  d[i] = d[i] - h;
      }

      f += h;
      p = d[m];
      c = 1.0;
      c2 = c;
      el1 = e[l1];
      s = 0.0;
      mml = m - l;

      for(int ii=1 ; ii<=mml ; ii++){
        c3 = c2;
        c2 = c;
        s2 = s;
        int i = m - ii;
        g = c * e[i];
        h = c * p;
        r = sqrt(p*p + e[i]*e[i]);
        e[i+1] = s * r;
        s = e[i] / r;
        c = p / r;
        p = c * d[i] - s * g;
        d[i+1] = h + s * (c * g + s * d[i]);

        if(job == 1){
          for(int k=1 ; k<=n ; k++){
            h = z[k][i+1];
            z[k][i+1] = s * z[k][i] + c * h;
            z[k][i]   = c * z[k][i] - s * h;
          }
        }
      }

      p = -s * s2 * c3 * el1 * e[l] / dl1;
      e[l] = s * p;
      d[l] = c * p;
      tst2 = tst1 + fabs(e[l]);
    }while(tst2 > tst1);

    d[l] += f;
  }


  for(int ii=2 ; ii<=n ; ii++){
    int i = ii - 1;
    int k = i;
    p = d[i];

    for(int j=ii ; j<=n ; j++){
      if(d[j] >= p) continue;
      k = j;
      p = d[j];
    }

    if(k == i) continue;
    d[k] = d[i];
    d[i] = p;

    for(int j=1 ; j<=n ; j++){
      p = z[j][i];
      z[j][i] = z[j][k];
      z[j][k] = p;
    }
  }

  return(ierr);
}


int EISPACK_htribk(const int n, double **ar, double **ai, double **tau, double **zr, double **zi)
{
  if (n == 0) return(0);

  for (int k=1; k<=n; k++){
    for(int j=1 ; j<=n ; j++){
      zi[k][j] = -zr[k][j] * tau[1][k];
      zr[k][j] *= tau[0][k];
    }
  }

  for(int i=2 ; i<=n ; i++){
    int l = i - 1;
    double h = ai[i][i];
    if(h == 0.0) continue;

    for(int j=1 ; j<=n ; j++){
      double sr = 0.0;
      double si = 0.0;

      for(int k=1 ; k<=l ; k++){
        sr += ar[i][k] * zr[k][j] - ai[i][k] * zi[k][j];
        si += ar[i][k] * zi[k][j] + ai[i][k] * zr[k][j];
      }

      sr = sr / h / h;
      si = si / h / h;

      for(int k=1 ; k<=l ; k++){
        zr[k][j] = zr[k][j] - sr * ar[i][k] - si * ai[i][k];
        zi[k][j] = zi[k][j] - si * ar[i][k] + sr * ai[i][k];
      }
    }
  }

  return (0);
}

// class and functions for B-spline function

/**************************************/
/*      B-Spline Interpolation        */
/**************************************/
class Spline{
 private:
    int       nmax      ;     /* maximal data points   */
    int       nd        ;     /* number of data points */
    double    *xdata    ;     /* x knots               */
    double    *ydata    ;     /* y knots               */
    double    *spline   ;     /* B-spline at (x,y)     */

    bool inrange(int n, int m){
      if(0<=n && n<m) return true;
      else            return false;
    }
 public:
    Spline(int n){
      nmax   = n;
      nd     = 0;
      xdata  = new double [n];
      ydata  = new double [n];
      spline = new double [n];
    }
    ~Spline(){
      delete [] xdata;
      delete [] ydata;
      delete [] spline;
    }
    void Setdat(double x, double y){
      if(inrange(nd,nmax)){
        xdata[nd] = x;
        ydata[nd] = y;
        nd++;
      }
    }

    double GetSpline(int i){
      if(inrange(i,nd)) return spline[i];
      else return 0.0;
    }
    void SetSpline(int i, double s){
      if(inrange(i,nd)){ spline[i] = s; }
    }
    double GetXdat(int i){
      if(inrange(i,nd)) return xdata[i];
      else return 0.0;
    }
    double GetYdat(int i){
      if(inrange(i,nd)) return ydata[i];
      else return 0.0;
    }
};

/**************************************/
/*      spline.cpp                    */
/**************************************/
double  splineInterpolation (const int, const double, Spline *);
void    splineSupport (const int, Spline *);


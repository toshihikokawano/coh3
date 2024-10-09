const int DefaultXsizeLimit = 1000;
const int DefaultYsizeLimit =  100;

/**************************************/
/*      Two-Dim Data Array            */
/**************************************/
class DataTable{
  private:
    int    xsizelimit;        // X-data size limit
    int    ysizelimit;        // Y-data size limit
    int    xsize;             // actual X-data size
    int    ysize;             // actual Y-data size
    bool   arrayalloc;        // flag for data allocation
  public:
    double *xdata;            // data buffer for X[xsize]
    double **ydata;           // data buffer for Y[xsize x ysize]
    double p;                 // temporal x data to be copied into actual X-data
    double *q;                // temporal y data

    DataTable(){
      arrayalloc = false;
      xsize = 0;
      ysize = 0;
      setDefaultLimit();
    }

    ~DataTable(){
      if(arrayalloc){
        for(int i=0 ; i<xsize ; i++) delete [] ydata[i];
        delete [] ydata;
        delete [] xdata;
        delete [] q;
        arrayalloc = false;
      }
    }

    void setLimit(const int n, const int m){
      xsizelimit = n;
      ysizelimit = m;
    }

    void setDefaultLimit(){
      xsizelimit = DefaultXsizeLimit;
      ysizelimit = DefaultYsizeLimit;
    }

    void memalloc(const int n, const int m){
      if(!arrayalloc){
        /*** when size limit is given */
        if((xsizelimit > 0) && (xsizelimit < n)) xsize = xsizelimit;
        else xsize = n;

        if((ysizelimit > 0) && (ysizelimit < m)) ysize = ysizelimit;
        else ysize = m;

        /*** allocate 2-dim data */
        xdata = new double [xsize];
        ydata = new double * [xsize];
        for(int i=0 ; i<xsize ; i++) ydata[i] = new double [ysize];

        q = new double [ysize];

        clearTable();
        arrayalloc = true;
      }
    }

    void clearTable(){
      for(int i=0 ; i<xsize ; i++){
        xdata[i] = 0.0;
        for(int j=0 ; j<ysize ; j++) ydata[i][j] = 0.0;
      }
    }

    void clearTemp(){
      p = 0.0;
      for(int j=0 ; j<ysize ; j++) q[j] = 0.0;
    }


    int getXsizeLimit(){ return xsizelimit; }
    int getYsizeLimit(){ return ysizelimit; }
  
    int getXsize(){ return xsize; }
    int getYsize(){ return ysize; }

    bool setData(const int n){
      bool flag = true;
      if((0 <= n) && (n < xsize)){
        xdata[n] = p;
        for(int i=0 ; i<ysize ; i++) ydata[n][i] = q[i];
      }
      else flag = false;
      return flag;
    }
  
    bool CheckRange(double x){
      bool flag = true;
      if( (x < xdata[0]) || (xdata[xsize-1] < x) ) flag = false;
      return flag;
    }
};


/**************************************/
/*      datatable.cpp                 */
/**************************************/
int dataTableRead (const std::string, DataTable *);
bool dataTableInterpolate (const double, DataTable *);
void dataTableWrite (DataTable *);


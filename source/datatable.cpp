/******************************************************************************/
/*  datatable.cpp                                                             */
/*        reading tabulated data                                              */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>

#include "datatable.h"


/**********************************************************/
/*      Read Tablulated Data in File                      */
/*      ----                                              */
/*      External data are given in a tabulated format     */
/*      # comment                                         */
/*      X(1)  Y(1,1) Y(1,2) Y(1,3) ...                    */
/*      X(2)  Y(2,1) Y(2,2) Y(2,3) ...                    */
/*      ...                                               */
/**********************************************************/
int dataTableRead(const std::string file, DataTable *tab)
{
  std::ostringstream os;
  std::ifstream      fp;
  std::string        line, item;

  /*** try current directry first */
  fp.open(&file[0]);
  if(!fp) return(-1);

  /*** set data size limit when not given */
  if(tab->getXsizeLimit() == 0 || tab->getYsizeLimit()  == 0) tab->setDefaultLimit();

  int nx = 0, ny = tab->getYsizeLimit();
  bool firstcall = true;
  const char *delim = " ";

  /*** pre-scan the data file to determine the data size */
  while(getline(fp,line)){
    if(line[0] == '#') continue;
    if(line.length() == 0) continue;

    if(firstcall){
      char *tok = strtok(&line[0],delim);
      if(tok != nullptr){
        for(int j=0 ; j<tab->getYsizeLimit() ; j++){
          tok = strtok(nullptr,delim);
          if(tok == nullptr){ ny = j; break; }
        }
        firstcall = false;
      }
    }

    std::istringstream ss(line);
    ss >> tab->p;
    if(tab->p == 0.0) break;
    nx ++;
  }
  fp.close();
  if(nx == 0 || ny == 0) return(-2);

  /*** allocate data buffer */
  tab->memalloc(nx,ny);
  
  /*** store all data in the buffer */
  fp.open(&file[0]);
  int k = 0;
  while(getline(fp,line)){
    if(line[0] == '#') continue;
    if(line.length() == 0) continue;

    std::istringstream ss(line);
    ss >> tab->p;
    for(int j=0 ; j<ny ; j++){
      ss >> item;
      tab->q[j] = atof(item.c_str());
    }

    tab->setData(k); // (p,q) will be copied to (xdata,ydata) in the k-th row
    k++;
  }
  fp.close();

  return(0);
}


/**********************************************************/
/*      Interpolate Data                                  */
/**********************************************************/
bool dataTableInterpolate(const double x0, DataTable *tab)
{
  bool binsearch = true;
  int k0 = 0;
  int k1 = tab->getXsize() - 1;

  tab->clearTemp(); // reset (p,q) output

  bool checkrange = tab->CheckRange(x0);

  /*** check if x0 is inside the range, then use the lowest/highest point */
  if(!checkrange){
    int k = (x0 < tab->xdata[k0]) ? k0 : k1;
    tab->p = x0;
    for(int j=0 ; j<tab->getYsize() ; j++) tab->q[j] = tab->ydata[k][j];
  }

  /*** general case */
  else{
    int k = k0;

    /*** find the interval xdata[k] <= x0 < xdata[k+1] */
    if(binsearch){ // binary search
      int c0 = k0;
      int c1 = k1;

      while(1){
        int c2 = (c0 + c1)/2;
        if(x0 < tab->xdata[c2]) c1 = c2;
        else c0 = c2;

        if((c1 - c0) == 1){ k = c0; break; }
      }
    }
    else{ // regular search
      for(k=k0 ; k<k1 ; k++){
        if( (tab->xdata[k] <= x0) && (x0 < tab->xdata[k+1]) ) break;
      }
    }

    /*** linear interpolation */
    tab->p = x0;
    double dx = tab->xdata[k+1] - tab->xdata[k];
    double x1 = x0 - tab->xdata[k];

    for(int j=0 ; j<tab->getYsize() ; j++){
      double dy = tab->ydata[k+1][j] - tab->ydata[k][j];
      tab->q[j] = dy / dx * x1 + tab->ydata[k][j];
    }
  }

  return(checkrange);
}


/**********************************************************/
/*      Print Data in a Table Format                      */
/**********************************************************/
void dataTableWrite(DataTable *tab)
{
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(5);

  for(int i=0 ; i<tab->getXsize() ; i++){
    std::cout << std::setw(13) << tab->xdata[i];
    for(int j=0 ; j<tab->getYsize() ; j++){
      std::cout << std::setw(13) << tab->ydata[i][j];
    }
    std::cout << std::endl;
  }
}

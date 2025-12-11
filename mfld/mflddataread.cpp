/******************************************************************************/
/*  mflddataread.cpp                                                          */
/*        functions to read input parameters                                  */
/******************************************************************************/

#include <iostream>
#include <fstream>

#include "mfld.h"
#include "elements.h"
#include "terminate.h"
#include "global.h"

static const unsigned int W_HEAD = 10;  // width of keyword

static char *   readExtractData (const int);
static int      readFgetOneline (char *);
static inline int  readElementToZ (char *);

static std::string head;
static char inputdata[256];

/**********************************************************/
/*      Read System Parameter                             */
/**********************************************************/
int readSystem(char *s, SystemData *dat, PrintData *prn, double *eps)
{
  while(1){
    if( readFgetOneline(s) < 0 ) break;

    if(head == "nucleus  :" || head == "target   :"){
      int z = readElementToZ(readExtractData(0));
      int a = atoi(          readExtractData(1));
      dat->target.setZA(z,a);

    }else if(head == "emax     :"){
      dat->emax = atof(readExtractData(0));

    }else if(head == "bin      :"){
      dat->de = atof(readExtractData(0));

    }else if(head == "eaverage :"){
      dat->eave = atof(readExtractData(0));

    }else if(head == "spmax    :"){
      dat->spmax = atof(readExtractData(0));

    }else if(head == "phmax    :"){
      dat->phmax = atoi(readExtractData(0));
      if(dat->phmax > dat->psize) dat->phmax = dat->psize;

    }else if(head == "epsilon  :"){
      eps[0]    = atof(readExtractData(0));
      eps[1]    = atof(readExtractData(1));
      eps[2]    = atof(readExtractData(2));
      eps[3]    = atof(readExtractData(3));
      eps[4]    = atof(readExtractData(4));
      eps[5]    = atof(readExtractData(5));

    }else if(head == "sharpcut :"){
      dat->sharpcutoff = true;

    }else if(head == "printLD  :"){
      std::string item = readExtractData(0);
      if     (item == "TotalDensity"      ){ prn->TotalDensity = true; }
      else if(item == "PartialDensity"    ){ prn->PartialDensity = true; }
      else if(item == "JDependentDensity" ){ prn->JDependentDensity = true; }
      else if(item == "BothParity"        ){ prn->BothParity = true; }
      else if(item == "Plotting"          ){ prn->Plotting = true; }
      else if(item == "SpinDistribution"  ){ prn->SpinDistribution = true; }
      else if(item == "SpinCutOff"        ){ prn->SpinCutOff = true; }
      else if(item == "AverageConfig"     ){ prn->AverageConfig = true; }
      else{
        message << "no such print option [" << item << "]";
        cohTerminateCode("readSystem");
      }

    }else if(head == "printFRDM:"){
      std::string item = readExtractData(0);
      if(item == "spstate"){
        pmf.spstate = true;
      }
      else if(item == "shiftedstate"){
        pmf.shiftedstate = true;
      }
      else if(item == "bcs"){
        pmf.bcs = true;
      }
      else if(item == "energy"){
        pmf.energy = true;
      }
      else if(item == "expansion"){
        pmf.expansion = true;
      }
      else if(item == "deformation"){
        pmf.deformation = true;
      }
      else{
        message << "no such print option [" << item << "]";
        cohTerminateCode("readSystem");
      }

    }else{
      message << "unknown key-word [" << head << "]";
      cohTerminateCode("readSystem");
    }
  }

  return(0);
}


/**********************************************************/
/*      Skip comment-line and read 1-line                 */
/**********************************************************/
static std::string line;

int readFgetOneline(char *s)
{
  /*** read one line, skip comment or blank line */
  while(1){
    std::cin.getline(s,255);
    if(std::cin.eof() != 0) return(-1);
    if(s[0] == '#') continue;
    if(strlen(s) <= 1) continue;
    break;
  }

  /*** copy data line to str, extract line head */
  line = (std::string)s;
  head = line.substr(0, W_HEAD);
  if(head.size() < W_HEAD){
    std::ostringstream os;
    message << "key-word [" << head << "] too short";
    cohTerminateCode("readFgetOneLine");
  }

  return(1);
}


/**********************************************************/
/*      Extract Data from Line                            */
/**********************************************************/
char * readExtractData(const int n)
{
  inputdata[0] = '\0';

  if( strlen(line.c_str()) > W_HEAD ){
    std::string work = line;
    char *str = &work[W_HEAD];
    char *tok = NULL;
    const char *delim = " ";

    tok = strtok(str,delim);

    if(tok != NULL){
      for(int i=1 ; i<=n ; i++){
        tok = strtok(NULL,delim);
        if(tok == NULL) break;
      }
    }
    if(tok != NULL) strncpy(&inputdata[0],tok,sizeof(inputdata)-1);
  }

  return(&inputdata[0]);
}


/**********************************************************/
/*      Convert Element Name into Z Number                */
/**********************************************************/
int readElementToZ(char *d)
{
  int z = element_getZ(d);

  if(z <=0 ){
    message << "element " << d << " not found";
    cohTerminateCode("readElementToZ");
  }

  return(z);
}



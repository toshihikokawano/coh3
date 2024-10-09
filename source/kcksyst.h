// define level density parameter file location
// function to read level density parameters

#ifndef __DIR_H__
#define __DIR_H_
#include "dir.h"
#endif

/**************************************/
/*      kcksyst.cpp                   */
/**************************************/
int     kckDataRead (ZAnumber *, LevelDensity *);
double  kckAsymptoticLevelDensity (const double);
double  kckSpinCutoff (const double);
double  kckTemperature (const double, const double);
double  kckE0 (const double, const double, const double);

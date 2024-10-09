// define discrete level data file location
// prototype of function to read disclete level data

#ifndef __DIR_H__
#define __DIR_H_
#include "dir.h"
#endif


/**************************************/
/*      Selection for Nmax            */
/**************************************/
typedef enum {normal=0, extended=1, reassign=2, all=3} MaxLevelCtl;


/**************************************/
/*      ripl2levels.cpp               */
/**************************************/
int     riplReadDiscreteLevels (ZAnumber *, Level *, const MaxLevelCtl);

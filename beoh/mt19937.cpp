#include <random>
#include "mt19937.h"

static std::mt19937 MersenneTwister;
static const double pow32 = 4294967296.0;


void genrand_init(unsigned long seed)
{
  std::mt19937 genrand(seed);
}

/* [0,1] interval */
double genrand_real1(void)
{
  double x = MersenneTwister();
  return(x / (pow32 - 1.0) ); 
}

/* [0,1) interval */
double genrand_real2(void)
{
  double x = MersenneTwister();
  return(x / pow32 ); 
}

/* (0,1) interval */
double genrand_real3(void)
{
  double x = MersenneTwister() + 0.5;
  return(x / pow32 ); 
}


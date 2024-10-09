// incident energy grid for MACS calculation

const int N_REACTION_RATE = 12;
const int N_MACS_EGRID = 30;
const int N_MACS_GGRID = 30;


#ifndef COH_TOPLEVEL
static double macs_egrid[N_MACS_EGRID]={
 0.00001, 0.0001,  0.001,  0.005,  0.010,  0.020,  0.040,  0.060,  0.080,  0.100,
   0.120,  0.150,  0.200,  0.300,  0.400,  0.500,  0.600,  0.800,  1.000,  1.200,
   1.500,  2.000,  2.500,  3.000,  4.000,  5.000,  6.000,  7.000,  8.000, 10.000};

static double macs_ggrid[N_MACS_GGRID]={
   1.000,  1.500,  2.000,  2.500,  3.000,  3.500,  4.000,  4.500,  5.000,  5.500,
   6.000,  6.500,  7.000,  7.500,  8.000,  8.500,  9.000,  9.500, 10.000, 10.500,
  11.000, 12.000, 13.000, 14.000, 15.000, 16.000, 17.000, 18.000, 19.000, 20.000};
#endif


/**************************************/
/*      macs.cpp                      */
/**************************************/
void    macsStoreCrossSection (const int, const int);
void    macsReactionRate (const int, const Particle, double *);

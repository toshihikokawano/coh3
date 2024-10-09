// functions for data output in main calculations

#ifndef __OUTPREP_H__
#define __OUTPREP_H__
#include "outprep.h"
#endif

/**********************************************************/
/*      Some Print Out Options                            */
/**********************************************************/

/*** print a banner at the top */
#undef PRINT_BANNER
//#define PRINT_BANNER

/*** print discrete levels on fission barriers */
#undef PRINT_BANDLEVELS


/**************************************/
/*      output.cpp                    */
/**************************************/
void    outSetSigmaReaction (const int);
void    outZA (ZAnumber *);
void    outTitle (char *);
void    outBanner (void);
void    outSystem  (System *, bool);
void    outTargetState (const int, const int, double);
void    outCompound (const int, Pdata *);
void    outLevelDensity (const int, double);
void    outGDR (const bool, const int, const double);
void    outFissionBarrier (const int);
void    outCrossSection (const int, double, double, double);
void    outReaction (const int, const int, double);
void    outTotalResidual (CumulativeResidualProduct *);
void    outIsomericRatio (CumulativeResidualProduct *);
void    outFission (const int);
void    outParticleProduction (const int, Channel *, double **);
void    outAngularDistribution (const int, int, int, int, ZAnumber *);
void    outLegendreCoefficient (const int, int, int, int, ZAnumber *, double ***);
void    outSpectrum (const int, double **, Nucleus *);
void    outSpectrumFineGamma (double *, GammaProduction *, Nucleus *);
void    outSpectrumSum (const int, double **, Nucleus *);
void    outPrimaryGammaSpectrum (const int, double **, Nucleus *);
void    outDiscreteGamma (GammaProduction *);
void    outDiscreteLevelPopulation (double, Nucleus *);
void    outPopulation (Nucleus *);
void    outGammaCascade (const int, Nucleus *);
void    outIsomerProduction (Nucleus *);
void    outDSD (double, double, double, Dcapt *);
void    outExcitonParameter (double **, double **, double **);
void    outChanceFispecParameter (const int,FChance *, Nucleus *);
void    outChanceFispec (const int, double *, double **, FChance *, Nucleus *);
void    outChanceFispecEnergy (double *);
void    outFissionNeutronEnergy (const int, double, double, double *, FNSpec *);
void    outFissionNeutronEnergyTable (const int, double *, FNSpec *);
void    outFissionNeutronSpectrum (const int, double *, double *, double);
void    outParameter (void);


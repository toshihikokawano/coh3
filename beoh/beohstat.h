/**************************************/
/*      beohspectra.cpp               */
/**************************************/
void    beohspectra (System *, Pdata *, Transmission **, Transmission **, double **,  Spectra *, GammaProduction *);
void    beohspectraLAB (System *, Pdata *, Transmission **, Transmission **, double **,  Spectra *, double, double);

/**************************************/
/*      beohspecmc.cpp                */
/**************************************/
void    beohspectraMC (Pdata *, Transmission **, Transmission **, double **,  Spectra *, const unsigned long);

/**************************************/
/*      beohprofile.cpp               */
/**************************************/
void    beohStrengthProfile (Beta *, Beta *, BetaProfile *);

/**************************************/
/*      beohpopinit.cpp               */
/**************************************/
void    beohInitialPopulation (const double, const int, BetaProfile *, Nucleus *);
void    beohInitialPopulationStat (const double, const double, const int, const double, const double, const double, Nucleus *);
int     beohStoreLevelExcite (const double, Nucleus *);
int     beohPopinitJmax (const double, Nucleus *);
void    beohClearInitialPopulation (Nucleus *);

/**************************************/
/*      beohneutrino.cpp              */
/**************************************/
void    beohElectronSpectrum (const double, double *, double *, Nucleus *, double *, double *);

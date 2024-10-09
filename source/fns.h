// functions for fission neutron spectrum model calculations

const double MAX_SPECTRUM_ENERGY = 30.0;

/**************************************/
/*      fns.cpp                       */
/**************************************/
void    fnsFissionNeutronModel (System *, Pdata *, double *, FNSpec *);


/**************************************/
/*      fnscalc.cpp                   */
/**************************************/
void    fnsMadlandNixModel (const int, const int, double *, double **, FChance *);
void    fnsInverseCrossSection (Pdata *, ZAnumber *, ZAnumber *);
double  fnsAverageSpectrumEnergy (const int, double *, double *);
void    fnsNormalizeSpectrum (const int, double *, double *);
void    fnsNormalizeExcitationDist (const int, double *, double *);


/**************************************/
/*      fnsprefis.cpp                 */
/**************************************/
double  fnsPreFissionSpectrum (const int, double *, double *, double *, Nucleus *);

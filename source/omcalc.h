// functions for optical model calculations

/**************************************/
/*      omcalc.cpp                    */
/**************************************/
int     omCalc (const double, Pdata *, ZAnumber *, const double, double *, CrossSection *);

/**************************************/
/*      cccalc.cpp                    */
/**************************************/
int     ccCalc (const double, const double, Pdata *, ZAnumber *, const double, Direct *, double **, CrossSection *);
int     ccPreset (const double, const double, Pdata *, ZAnumber *, const double, Direct *);
void    ccMain (const double, const double, Particle, ZAnumber *, double **, CrossSection *);
void    ccCleanUp (void);
int     ccHermiteMatrix (const int, const int);


/**************************************/
/*      dwcalc.cpp                    */
/**************************************/
int     dwbaCalc (const double, const double, Pdata *, ZAnumber *, const double, Direct *, CrossSection *);


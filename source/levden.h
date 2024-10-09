// functions for level density formulas

/**************************************/
/*      levden.cpp                    */
/**************************************/

double  ldLevelDensity (const double, const double, LevelDensity *);
double  ldDensityParameter (const double, const double, LevelDensity *);
double  ldSpinDistribution (const double, const double, const double, const double);
double  ldParityDistribution (const int);
double  ldParityDistributionEdepend (const int, const double, const double);
double  ldLevelSpinCutoff (const int, Level *);
void    ldTGconnect (const double, const double, const double, LevelDensity *);
void    ldTextrapolate (const double, const double, LevelDensity *);
void    ldGextrapolate (const double, LevelDensity *);
double  ldFermiGas (const double, const double, const double);
double  ldConstantTemperature (const double, const double, const double);
double  ldSpinCutoff (double, const double, const double, const double);
double  ldShellCorrection (const double, const double, const double, const double);


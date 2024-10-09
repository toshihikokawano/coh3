// functions for the least-squares method

/**************************************/
/*     POLYSQ.CPP                     */
/**************************************/
int     LSQPolynomial (const int, const int, double *, double *, double *);
int     LSQLegendre (const bool, const int, const int, double *, double *, double *);

/**************************************/
/*     POLYCALC.CPP                   */
/**************************************/
double  LSQCalc (const int, const int, double *, double *, double *, double *, double *);

int     LSQGeneral (const int, const int, double *, double *, double *, double *, double *, double *);

int     LSQMarquardt (const int, const int, const double, double *, double *, double *, double *, double *);

// functions for data output in optical model calculations
// these functions are independent from output.cpp

/**************************************/
/*      omoutput.cpp                  */
/**************************************/
void    outOMP (const int, Optical *);
void    outOMPtable (const int, Optical *);
void    outCoupledState (const NuclearModel, const int, LevelData *);
void    outDeformation (const NuclearModel, const int, double *);
void    outLevelExcite (const int, const int, const double, const double, LevelData *, double *);
void    outRmatrix (const int, std::complex<double> *, std::complex<double> *);
void    outSmatrix (const int, int, std::complex<double> *);
void    outTransmission (const int, const int, const double, double *);
void    outCCSmatrix (const int, const int, const int, Collective *, CCdata *, std::complex<double> *);

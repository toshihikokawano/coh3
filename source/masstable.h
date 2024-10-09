// define mass excess data file location
// function to read mass table

class MassExcess{
 public:
  unsigned int za;    // Z*1000 + A
  float        mass;  // mass excess
};


/**************************************/
/*      masstable.cpp                 */
/**************************************/
void    massReadFile (std::string);
double  mass_excess (const int, const int);
double  mass_excess (const int, const int, bool *);


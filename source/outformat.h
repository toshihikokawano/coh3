// formatting output values

/*** lower limit of output values, enforce zero */
static const double output_eps = 1.0e-99;

static std::string particle_name[9]={
   "      gamma","    neutron","     proton","      alpha",
   "   deuteron","     triton","     helion","    fission",
   "    unknown"};

static std::string p_name = "gnpadthf ";

static std::string cline="#          ";
static std::string blank="           ";
static std::string dashl=" ----------";
static const int DisplayWidth = 99;


// eliminate super small numbers
static inline double lowfilter(double x)
{ if(fabs(x) < output_eps) return 0.0; else return(x); }

// fixed width (11) with 4 digits below decimal point
static inline void outVal(double x)
{ std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4) << std::setw(11) << x; }

// variable width (w)
static inline void outVal(int w, double x)
{ std::cout.setf(std::ios::scientific, std::ios::floatfield);
  int p = w - 8;
  if(p >= 1){ std::cout << std::setprecision(w-8) << std::setw(w) << x; }
  else      { std::cout << std::setw(w) << x; }}

static inline void outVal(int w, int x)
{ std::cout << std::dec << std::setw(w) << x; }

static inline void outVal(int w, unsigned int x)
{ std::cout << std::dec << std::setw(w) << x; }

// variable width (w) and precistion (p)
static inline void outVal(int w, int p, double x)
{ std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout << std::setw(w) << std::setprecision(p) << x; }

// new line
static inline void nl()
{  std::cout << std::endl; }


/**************************************/
/*      output.cpp                    */
/**************************************/
// they need to be separated from outupt.h because of omoutput.cpp
void outSectionHead (const char *);
double outRetrieveLabE (void);

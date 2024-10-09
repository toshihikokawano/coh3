// formatting of output values

/*** lower limit of output values, enforce zero */
static const double output_eps = 1.0e-99;


static std::string cline="#          ";
static std::string blank="           ";
static std::string dashl=" ----------";
static const int DisplayWidth = 99;

void outSectionHead (const char *);
double outRetrieveLabE (void);

static inline double lowfilter(double x)
{ if(fabs(x) < output_eps) return 0.0; else return(x); }

static inline void outVal(double x)
{ std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << std::setprecision(4) << std::setw(11) << x; }

static inline void outVal(int w, double x)
{ std::cout.setf(std::ios::scientific, std::ios::floatfield);
  int p = w - 8;
  if(p >= 1){ std::cout << std::setprecision(w-8) << std::setw(w) << x; }
  else      { std::cout << std::setw(w) << x; }}

static inline void outVal(int w, int x)
{ std::cout << std::dec << std::setw(w) << x; }

static inline void outVal(int w, unsigned int x)
{ std::cout << std::dec << std::setw(w) << x; }

static inline void outVal(int w, int p, double x)
{ std::cout.setf(std::ios::fixed, std::ios::floatfield);
  std::cout << std::setw(w) << std::setprecision(p) << x; }

static inline void nl()
{  std::cout << std::endl; }


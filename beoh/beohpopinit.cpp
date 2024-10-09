/******************************************************************************/
/*  popinit.cpp                                                               */
/*        store initial polulation into energy bin or level                   */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>

#include "beoh.h"
#include "statmodel.h"
#include "beohstat.h"
#include "levden.h"

static void beohPopinitSingle (const double, const int, const int, const double, Nucleus *);
static void beohPopinitLevden (const int, const double, const double, Nucleus *);
static void beohPopinitFile (const std::string, Nucleus *);
static void beohPopinitContDist (const int, const double, const double, const double, const double, Nucleus *, const double, const int);

#undef POPDEBUG
#ifdef POPDEBUG
static void beohPopCheck (Nucleus *);
static void beohPopProfile (Nucleus *);
#endif

static const std::string fname = "";

/**********************************************************/
/*      Store Initial Population for Beta-Decay           */
/**********************************************************/
void  beohInitialPopulation(const double j0, const int p0, BetaProfile *p, Nucleus *n)
{
  n->jmax = MAX_J-1;

  /*** divide strength by spin fraction */
  double m = (j0 <= 0.5) ? 2.0 : 3.0;
  for(int k=0 ; k<n->ncont ; k++){
    if(j0 <= 0.5){
      beohPopinitSingle(j0    ,p0,k,1.0,n);
      beohPopinitSingle(j0+1.0,p0,k,1.0,n);
    }
    else{
      beohPopinitSingle(j0-1.0,p0,k,1.0,n);
      beohPopinitSingle(j0    ,p0,k,1.0,n);
      beohPopinitSingle(j0+1.0,p0,k,1.0,n);
    }
    for(int j=0 ; j<=n->jmax ; j++){
      n->pop[k][j].even /= m;
      n->pop[k][j].odd  /= m;
    }
  }

  /*** multiply k0-th bin population by the strength profile */
  for(int k=0 ; k<n->ncont ; k++){
    for(int j=0;j<=n->jmax;j++){
      n->pop[k][j].even *= p->rcont[k];
      n->pop[k][j].odd  *= p->rcont[k];
    }
  }

  /*** discrete population */
  for(int i=0 ; i<n->ndisc ; i++) n->lpop[i] = p->rdisc[i];

#ifdef POPDEBUG
  beohPopCheck(n);
#endif
}


/**********************************************************/
/*      Store Initial Population for Statistical Decay    */
/**********************************************************/
void  beohInitialPopulationStat(const double q, const double j0, const int p0, const double sf, const double e0, const double ew, Nucleus *n)
{
  int jmax = MAX_J - 1;

  /* the excitation energy is distributed in the continuum */
  if(ew > 0.0){
    jmax = beohPopinitJmax(e0,n);
    beohPopinitContDist(jmax,q,sf,e0,ew,n,j0,p0);
  }
  /* all population at the top energy bin */
  else{
    if(fname != ""){
      /*** read data file, now optional */
      beohPopinitFile(fname,n);
    }
    else{
      /*** when initial spin and parity are not given, use level density spin distribution */
      if(p0 == 0){
        jmax = beohPopinitJmax(0.0,n);
        beohPopinitLevden(jmax,q,sf,n);
      }
      /*** otherwise, start with the given spin and parity */
      else{
        beohPopinitSingle(j0,p0,0,q,n);
      }
    }
   }

  n->jmax = jmax;

#ifdef POPDEBUG
  beohPopCheck(n);
//beohPopProfile(n);
#endif
}


/**********************************************************/
/*      Popolation Localized on a Single JP State         */
/**********************************************************/
void beohPopinitSingle(const double j0, const int p0, const int k0, double q, Nucleus *n)
{
  for(int j=0 ; j<MAX_J ; j++){
    if(j == (int)j0){
      /*** J is given but parity is not */
      if(p0 == -2){
        n->pop[k0][j].even += 0.5 * q;
        n->pop[k0][j].odd  += 0.5 * q;
      }
      /*** J and parity are both given */
      else if(p0 == -1){
        n->pop[k0][j].even  = 0.0;
        n->pop[k0][j].odd  += q;
      }
      else{
        n->pop[k0][j].even += q;
        n->pop[k0][j].odd   = 0.0;
      }
    }
  }
}


/**********************************************************/
/*      Population Distributed Over the Spin States       */
/**********************************************************/
void beohPopinitLevden(const int jmax, const double q, const double sf, Nucleus *n)
{
  double sd;
  double se = 0.0;
  double so = 0.0;

  for(int j=0 ; j<=jmax ; j++){
    sd = ldSpinDistribution(j+halfint(n->lev[0].spin),n->ldp.spin_cutoff*sf,n->ldp.sigma0,1.0);
    se += (n->pop[0][j].even = sd * ldParityDistribution( 1));
    so += (n->pop[0][j].odd  = sd * ldParityDistribution(-1));
  }

  sd = se + so;
  for(int j=0 ; j<=jmax ; j++){
    n->pop[0][j].even = q * n->pop[0][j].even / sd;
    n->pop[0][j].odd  = q * n->pop[0][j].odd  / sd;
  }
}


/**********************************************************/
/*      Population Distributed Over the Spin States       */
/**********************************************************/
void beohPopinitFile(const std::string fname, Nucleus *n)
{
  double sd;
  double se = 0.0;
  double so = 0.0;
  std::ifstream fp;

  fp.open(&fname[0]);

  int j = 0;
  std::string line;
  while(getline(fp,line)){
    if(j >= MAX_J) break;
    std::stringstream ss(&line[1]);  // skip comment #
    ss >> n->pop[0][j].even >> n->pop[0][j].odd;
    se += n->pop[0][j].even;
    so += n->pop[0][j].odd;
    j++;
  }
  fp.close();

  int jmax = j-1;

  sd = se + so;
  se = so = 0.0;
  for(int j=0 ; j<=jmax ; j++){
    n->pop[0][j].even /= sd;
    n->pop[0][j].odd  /= sd;

    se += n->pop[0][j].even;
    so += n->pop[0][j].odd;
  }
}


/**********************************************************/
/*      Population Distributed in the Continuum           */
/**********************************************************/
void beohPopinitContDist(const int jmax, const double q, const double sf, double e0, double ew, Nucleus *n, const double j0, const int p0)
{
  int k0 = 0;
  double sd, *g;
  Parity **p;

  g = new double [MAX_ENERGY_BIN];
  p = new Parity * [MAX_ENERGY_BIN];
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    p[k] = new Parity [jmax + 1];

    g[k] = 0.0;
    for(int j=0; j<=jmax ; j++) p[k][j].even = p[k][j].odd = 0.0;
  }
                                        

  /*** find bin-index for the lowest bin */
  int k1 = 0;
  double el = e0 - EXRANGE_FACTOR * ew;

  if(el < n->excitation[n->ncont-1]){
    k1 = n->ncont-1;
  }
  else{
    for(int k=k0 ; k<n->ncont-1 ; k++){
      if((n->excitation[k+1] < el) && (el <= n->excitation[k])){ k1 = k; break; }
    }
  }

  /*** Gaussian distribution */
  double sp = 0.0;
  for(int k=k0 ; k<=k1 ; k++){
    double de = n->excitation[k] - e0;
    g[k] = exp(-de*de/(2.0*ew*ew));
    sp += g[k];
  }
  if(sp > 0.0) for(int k=k0 ; k<=k1 ; k++) g[k] =  g[k] / sp;


  /*** population = Gaussian x spin-distribution */
  if(p0 == 0){
    sp = 0.0;
    for(int k=k0 ; k<=k1 ; k++){
      double se = 0.0;
      double so = 0.0;
      for(int j=0 ; j<=jmax ; j++){
        sd = ldSpinDistribution(j+halfint(n->lev[0].spin),n->ldp.spin_cutoff*sf,n->ldp.sigma0,1.0);
        se += (p[k][j].even = g[k] * sd * ldParityDistribution( 1));
        so += (p[k][j].odd  = g[k] * sd * ldParityDistribution(-1));
      }
      sp += se + so;
    }
    for(int k=k0 ; k<=k1 ; k++){
      for(int j=0 ; j<=jmax ; j++){
        n->pop[k][j].even += q * p[k][j].even / sp;
        n->pop[k][j].odd  += q * p[k][j].odd  / sp;
      }
    }
  }
  else{
    for(int k=k0 ; k<=k1 ; k++) beohPopinitSingle(j0,p0,k,g[k],n);
  }
  
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    delete [] p[k];
  }
  delete [] p;
  delete [] g;
}


/***********************************************************/
/*     Find Nearest Discrete Level, and  Store Population  */
/***********************************************************/
int  beohStoreLevelExcite(const double q, Nucleus *n)
{
  double d = 100.0, x;
  int k = 0;

  if( (n->ndisc == 1) && (n->ncont == 0) ) return(k);

  for(int i=0 ; i<n->ndisc-1 ; i++){
    x = fabs(n->max_energy - n->lev[i].energy);
    if(d>x){
      k = i;
      d = x;
    }
  }
  if(k>0) n->lpop[k] = q;

#ifdef POPDEBUG
  beohPopCheck(n);
#endif

  return(k+1);
}


/**********************************************************/
/*      Determine Jmax and save it                        */
/**********************************************************/
int beohPopinitJmax(const double e0, Nucleus *n)
{
  int jmax = MAX_J - 1;
  int k0 = 0, km = 0;

  /*** find bin-index for the mean energy */
  if(e0 > 0.0){
    for(int k=k0 ; k<n->ncont-1 ; k++){
      if((n->excitation[k+1] < e0) && (e0 <= n->excitation[k])){ km = k; break; }
    }
  }

  ldLevelDensity(n->excitation[km],(double)n->za.getA(),&n->ldp);

  /*** first determine the highest J without scaling */
  for(int j=0 ; j<MAX_J ; j++){
    double sd = ldSpinDistribution(j+halfint(n->lev[0].spin),n->ldp.spin_cutoff,n->ldp.sigma0,1.0);
    if(sd == 0.0){
      jmax = j-1;
      break;
    }
  }
  return jmax;
}


/**********************************************************/
/*      Fill Zeros in Population Array                    */
/**********************************************************/
void beohClearInitialPopulation(Nucleus *n)
{
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    for(int j=0 ; j<MAX_J ; j++) {
      n->pop[k][j].even = n->pop[k][j].odd = 0.0;
    }
  }
  for(int k=0 ; k<MAX_LEVELS ; k++) n->lpop[k] = 0.0;
}


#ifdef POPDEBUG
/**********************************************************/
/*      For Debugging                                     */
/**********************************************************/
void beohPopCheck(Nucleus *n)
{
  double s1 = 0.0, s2 = 0.0;
  const int jm = 10;

  std::cout << "#  ";
  std::cout << std::setw(3) << setfill('0') << n->za.getZ() << '-';
  std::cout << std::setw(3) << setfill('0') << n->za.getA() << std::endl;
  std::cout << setfill(' ');

  for(int k=0 ; k<n->ncont ; k++){
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout << std::setw(5) << k;
    std::cout << std::setw(7) << std::setprecision(3) << n->excitation[k];
    std::cout.setf(std::ios::scientific, std::ios::floatfield);

    double s3 = 0.0;
    for(int j=0 ; j<=n->jmax ; j++) {
      s1 += n->pop[k][j].even;
      s3 += n->pop[k][j].even;
      if(j < jm) std::cout << std::setprecision(2) << std::setw(9) << n->pop[k][j].even;
    }
    std::cout << std::endl;

    std::cout << "            ";
    for(int j=0 ; j<=n->jmax ; j++) {
      s1 += n->pop[k][j].odd;
      s3 += n->pop[k][j].odd;
      if(j < jm) std::cout << std::setprecision(2) << std::setw(9) << n->pop[k][j].odd;
    }
    std::cout << std::setprecision(4) << std::setw(11) << s3 << std::endl;
  }

  s2 = 0.0;
  for(int k=0 ; k<n->ndisc ; k++){
    s2 += n->lpop[k];
    std::cout << std::setw(5) << k
         << std::setprecision(4) << std::setw(13) << n->lev[k].energy
         << std::setw(13) << n->lpop[k] << std::endl;
  }

  std::cout << "# SUMcont " << s1 << std::endl;
  std::cout << "# SUMdisc " << s2 << std::endl;
  std::cout << "# SUMall  " << s1+s2 << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
/*
  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  for(int k=0 ; k<n->ncont ; k++){
    for(int j=0 ; j<=n->jmax ; j++) {
      std::cout << std::setw(5) << j << std::setprecision(4) << std::setw(13) << n->excitation[k] << std::setw(13)<< n->pop[k][j].even + n->pop[k][j].odd << std::endl;
    }
    std::cout << std::endl;
  }
*/
}

void beohPopProfile(Nucleus *n)
{
  double s1 = 0.0;
  for(int k=0 ; k<n->ncont ; k++){
    double s2 = 0.0;
    for(int j=0 ; j<MAX_J-1 ; j++) s2 += n->pop[k][j].even + n->pop[k][j].odd;

    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout << std::setw(7) << std::setprecision(3) << n->excitation[k];
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout << std::setprecision(4) << std::setw(11) << s2 << std::endl;

    s1 += s2;
  }
  std::cout << "# " << s1 << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

#endif


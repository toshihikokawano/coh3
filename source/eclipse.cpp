/******************************************************************************/
/*  eclipse.cpp                                                               */
/*      ExClusive LIght Particle SPEctra                                      */
/*      eclipse incorporated into main CoH source, Sept. 2012                 */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "structur.h"
#include "nucleus.h"
#include "global.h"
#include "eclipse.h"
#include "terminate.h"

static void   eclClearArray (const int, const int, const int);
static void   eclCountParticle (const int, const int);
static void   eclGammaMultiplicity (const int, const int, double, EXSpectra *);
static void   eclPEFraction (const int, double **, Nucleus *);
static void   eclAddBinaryReaction (const int, const int, EXSpectra *);

static double   ****ptbl;  // decay probability table
static int        **cidx;  // decay chain indices
static int         *ktot;  // total bin number
static double     **mtpl;  // number of emitted particls (multiplicity)
static EXSpectra    *dat;
static int          max_energybin = 0;

static bool memallocated = false;

#undef DEBUG

/**********************************************************/
/*      ECLIPSE Main Section                              */
/**********************************************************/
void eclCalc(System *sys, double **pe, const unsigned long nsim)
{
  const int nm = sys->max_compound;
  const int cm = sys->max_channel;
  double    de = sys->energy_bin;

  if(pex.ddx) ctl.ddx = true;

  //--------------------------------------------------------
  // Prepare decay probability matrix, and indices

  for(int c0=0 ; c0<nm ; c0++){
    ktot[c0] = ncl[c0].ntotal - 1;

    int nl = ncl[c0].ndisc;
    if(nl > NLEV) nl = NLEV;
    dat[c0].set(ncl[c0].ncont,nl,ncl[c0].max_energy,ncl[c0].lev[ ncl[c0].ndisc-1 ].energy);

    if(ctl.ddx){
      /*** check if binary reaction */
      dat[c0].setGcon(false); // means discrete population is by a particle only
      for(int c1=1 ; c1<cm ; c1++){
        /*** look for the case where the residual from the initial CN is C0 */
        if(c0 == ncl[0].cdt[c1].next){
          dat[c0].setGcon(true); // means discrete transition not in spectrum
          if(!opt.chargediscrete && (c1 != neutron)) dat[c0].setGcon(false);
          break;
        }
      }
    }

    for(int c1=0 ; c1<cm ; c1++) cidx[c0][c1] = ncl[c0].cdt[c1].next;
  }

  if(pex.ptable) eclOutPtable(nm,cm,ktot,cidx,ptbl);

  
  //--------------------------------------------------------
  // Fractions of preequilibrium emission

  eclPEFraction(cm,pe,&ncl[0]);


  //--------------------------------------------------------
  // Reconstruct energy spectra, perform MC simulation

  /*** count number of particles emitted */
  eclCountParticle(nm,cm);

  /*** exclusive pre-fission neutron spectrum calculation */
  if(ctl.fns){
    eclSpectra(false,ctl.fns,cm,ktot,cidx,ptbl,dat);
//  eclOutSpectra(nm,max_energybin,cm,de,mtpl,dat);
  }

  /*** Monte Carlo Simulation */
  else if(ctl.montecarlo){
    eclMC(nm,max_energybin,cm,ktot,cidx,ptbl,mtpl,nsim,de,dat);
    eclOutSpectra(nm,max_energybin,cm,de,mtpl,dat);
  }

  /*** deterministic method */
  else{
    eclSpectra(ctl.ddx,ctl.fns,cm,ktot,cidx,ptbl,dat);
    eclGammaMultiplicity(nm,cm,de,dat);
    if(ctl.ddx){
      /*** calculate DDX and gamma-ray multiplicities */
      eclAddBinaryReaction(nm,cm,dat);
//      eclOutSpectra(nm,max_energybin,cm,de,mtpl,dat);
      eclOutHead(nm,cm,sys->lab_energy);
      eclDDX(nm,max_energybin,cm,sys,mtpl,pe,dat);
    }
    else{
      eclOutSpectra(nm,max_energybin,cm,de,mtpl,dat);
    }
  }
}


/**********************************************************/
/*      Allocate Memory                                   */
/**********************************************************/
void eclAllocateMemory(const int nm, const int cm)
{
  if(memallocated) return;

  /*** determine maximum bin number for all compound nuclei */
  max_energybin = 0;
  for(int c=0 ; c<nm ; c++){
    for(int j=0 ; j<cm ; j++){
      if(!ncl[c].cdt[j].status) continue;
      int p = ncl[c].cdt[j].next;
      if(ncl[p].ntotal > max_energybin) max_energybin = ncl[p].ntotal;
    }
  }
  if(max_energybin == 0) return;

  max_energybin += 2;

  //--------------------------------------------------------
  // Memory Allocation
  //

  try{
    /*** calculated spectrum array */
    dat = new EXSpectra [nm];
    for(int n=0 ; n<nm ; n++){
      dat[n].memalloc(cm,max_energybin);
    }

    /*** decay probability from CoH calculation : */
    /*** ptbl[Parent Index][Parent Energy][Decay Channel][Daughter Energy] */
    ptbl = new double *** [nm];
    for(int n=0 ; n<nm ; n++){
      ptbl[n] = new double ** [max_energybin];

      for(int k=0 ; k<max_energybin ; k++){
        ptbl[n][k] = new double * [cm];
        for(int c=0 ; c<cm ; c++){
          ptbl[n][k][c] = new double [max_energybin];
        }
      }
    }
    /*** cidx[Parent Index][Decay Channel] = Daughter Index */
    cidx = new int * [nm];
    for(int i=0 ; i<nm ; i++) cidx[i] = new int [cm];

    /*** mtpl[Parent Index][Particle ID] = Number of Emitted Particles */
    mtpl = new double * [nm];
    for(int i=0 ; i<nm ; i++) mtpl[i] = new double [cm];

    /*** ktot[Parent Index] = Total Number of Bins */
    ktot = new int [nm];
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error" << e.what();
    cohTerminateCode("eclAllocateMemory");
  }

  eclClearArray(nm,max_energybin,cm);

  memallocated = true;
}


/**********************************************************/
/*      Free Allocated Memory                             */
/**********************************************************/
void eclDeleteAllocated(const int nm, const int cm)
{
  if(!memallocated) return;

  for(int n=0 ; n<nm ; n++){
    for(int k=0 ; k<max_energybin ; k++){
      for(int c=0 ; c<cm ; c++) delete [] ptbl[n][k][c];
      delete [] ptbl[n][k];
    }
    delete [] ptbl[n];
    delete [] cidx[n];
    delete [] mtpl[n];
  }
  delete [] ptbl;
  delete [] cidx;
  delete [] mtpl;
  delete [] ktot;
  delete [] dat;

  memallocated = false;
}


/**********************************************************/
/*      Initialize All Data Arrays                        */
/**********************************************************/
void eclClearArray(const int nm, const int km, const int cm)
{
  for(int c0=0 ; c0<nm ; c0++){
    for(int k0=0 ; k0<km ; k0++){
      for(int c1=0 ; c1<cm ; c1++){
        dat[c0].spec[c1][k0] = 0.0;
        for(int k1=0 ; k1<km ; k1++) ptbl[c0][k0][c1][k1] = 0.0;
      }

      // flag to check if table exists
      // this will be overridden when data are read in
      ptbl[c0][k0][0][k0] = -1.0;

      dat[c0].glin[k0] =
      dat[c0].sbin[k0] =
      dat[c0].pfis[k0] = 0.0;
    }
    for(int c1=0 ; c1<cm ; c1++){
      cidx[c0][c1] = -1;
      mtpl[c0][c1] = 0.0;
    }
  }
}


/**********************************************************/
/*      Generate Decay Probability Table                  */
/**********************************************************/
void  eclGenerateDecayTable(const int c0, const int k0, const int cm, double **dlp)
{
  int km = ncl[c0].ntotal;
  if(km == 0) return;

  /*** fission probablity of bin(c0,k0) */
  dat[c0].pfis[k0] = (dlp[0][0] >= 1.0) ? 1.0 : dlp[0][0];
  dlp[0][0] = 0.0;

  /*** calculate sum of population increment */
  double sum = 0.0;
  for(int k=k0 ; k<km ; k++){
    for(int j=0 ; j<cm ; j++){
      if( (k == 0) && (j == 0) ) continue;
      sum += dlp[j][k];
    }
  }

  /*** make dlp probabilities */
  if(sum > 0.0){
    for(int k=k0 ; k<km ; k++){
      for(int j=0 ; j<cm ; j++) dlp[j][k] /= sum;
    }
  }
  else{
    dlp[0][km-1] = 1.0;
  }

  /*** generate decay probability table */
  for(int k1=k0 ; k1<km ; k1++){
    for(int c1=0 ; c1<cm ; c1++) {
      ptbl[c0][k0][c1][k1] = dlp[c1][k1];
    }
  }

#ifdef DEBUG
  double s = 0.0;
  for(int k1=k0 ; k1<km ; k1++){
    for(int c1=0 ; c1<cm ; c1++){
      std::cout << std::setw(5) << c0 << std::setw(5) << k0;
      std::cout << std::setw(5) << c1 << std::setw(5) << k1;
      std::cout << std::setw(13) << ptbl[c0][k0][c1][k1];
      s += ptbl[c0][k0][c1][k1];
    }
    std::cout << std::endl;
  }
  std::cout <<"SUM " << s << std::endl;
#endif
}


/**********************************************************/
/*      Store Level Population                            */
/**********************************************************/
void  eclLevelPopulation(const int nm, double **lp0, double **lp1)
{
  for(int c=0 ; c<nm ; c++){
    for(int i=0 ; i<ncl[c].ndisc ; i++){
      if(i < NLEV){
        dat[c].lpop[i].particle = lp0[c][i];
        dat[c].lpop[i].gamma    = lp1[c][i];
      }
    }
  }
}


/**********************************************************/
/*      Count Number of Emitted Particles                 */
/**********************************************************/
void  eclCountParticle(const int nm, const int cm)
{
  for(int n=0 ; n<nm ; n++){
    for(int c=0 ; c<cm ; c++){
      if(cidx[n][c] < 0) continue;
      mtpl[cidx[n][c]][c]  = mtpl[n][c];
      mtpl[cidx[n][c]][c] += 1.0;
    }
  }
} 


/**********************************************************/
/*      Gamma-Ray Multiplicity                            */
/**********************************************************/
void eclGammaMultiplicity(const int nm, const int cm, double de, EXSpectra *dat)
{
  double *eav,*tot,*sum;
  eav = new double [cm];
  tot = new double [cm];
  sum = new double [cm];

  for(int n=0 ; n<nm ; n++){

    /*** total and average energies */
    for(int c=0 ; c<cm ; c++){
      sum[c] = 0.0;
      tot[c] = 0.0;

      /*** an implicit first point is zero, 
           e0 is shifted to dE/4 */
      double s0 = 0.0, s1 = 0.0;
      if(c == 0){
        s0 = dat[n].spec[0][0] + dat[n].glin[0];
        s1 = dat[n].spec[0][1] + dat[n].glin[1];
      }
      else{
        s0 = dat[n].spec[c][0];
        s1 = dat[n].spec[c][1];
        /*** for binary reactions, add sbin */
        if(n == c){
          s0 += dat[n].sbin[0];
          s1 += dat[n].sbin[1];
        }
      }

      double e0 = 1.0/4.0 *de;
      double e1 = de;

      /*** the first triangle in [0,E0], and the second trapezoid [E0,E1] */
      sum[c] = s0*e0/2.0 + (s0+s1)*(e1-e0)/2.0;

      /*** integral for SxE in this energy range */
      tot[c] = s0*e0*e0/3.0 + (s0*(e1+2.0*e0) + s1*(e0+2.0*e1))*(e1-e0)/6.0;

      for(int k=2 ; k<=ktot[n]+1 ; k++){
        if(c == 0){
          s0 = dat[n].spec[0][k-1] + dat[n].glin[k-1];
          s1 = dat[n].spec[0][k  ] + dat[n].glin[k  ];
        }
        else{
          s0 = dat[n].spec[c][k-1];
          s1 = dat[n].spec[c][k  ];
          if(n == c){
            s0 += dat[n].sbin[k-1];
            s1 += dat[n].sbin[k  ];
          }
        }

        e0 = (k-1)*de;
        e1 =  k   *de;
        if(k == ktot[n]+1) e1 = e0 + 0.5*de;

        sum[c] += (s0+s1)*(e1-e0)/2.0;
        tot[c] += (s0*(e1+2.0*e0) + s1*(e0+2.0*e1))*(e1-e0)/6.0;
      }
      eav[c] = (sum[c] > 0.0) ? tot[c]/sum[c] : 0.0;
    }

    /*** sum of average particle energies */
    double ep = 0.0;
    for(int c=1 ; c<cm ; c++){
      if(mtpl[n][c] != 0.0) ep += mtpl[n][c]*eav[c];
    }
    /*** multiplicity = average excitation / average gamma energy */
    mtpl[n][0] = (eav[0] > 0.0) ? (dat[n].getEm() - ep)/eav[0] : 0.0;
    if(mtpl[n][0] < 0.0) mtpl[n][0] = 0.0;

#ifdef DEBUG
    std::cout << std::setprecision(6) << std::setiosflags(std::ios::scientific);
    std::cout << "# " << std::setw(3) << n << std::setw(13) << mtpl[n][0] << std::setw(13) << dat[n].getEm() << std::setw(13)<< ep << std::endl;
    std::cout << "# Sum        ";
    for(int c=0 ; c<cm ; c++) std::cout << std::setw(13) << sum[c];
    std::cout << std::endl;

    std::cout << "# Etot       ";
    for(int c=0 ; c<cm ; c++) std::cout << std::setw(13) << tot[c];
    std::cout << std::endl;

    std::cout << "# Eave       ";
    for(int c=0 ; c<cm ; c++) std::cout << std::setw(13) << eav[c];
    std::cout << std::endl;
#endif
  }

  delete [] eav;
  delete [] tot;
  delete [] sum;
}


/**********************************************************/
/*      Preequilibrium Emission Fraction                  */
/**********************************************************/
void eclPEFraction(const int cm, double **pe, Nucleus *n)
{
  int kmax = (int)(n->max_energy/n->de) +1;

  for(int j=0 ; j<cm ; j++){
    if (!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    int k = (int)(ncl[p].max_energy/ncl[p].de);
    if(k > kmax) kmax = k;
  }

  /*** ratio to the total emissions */
  for(int j=1 ; j<cm ; j++){
    for(int k=0 ; k<kmax ; k++){
      pe[j][k] = (crx.spectra[j][k] > 0.0) ? pe[j][k]/crx.spectra[j][k] : 0.0;
    }
  }
}


/**********************************************************/
/*      Copy PreFission Neutron Spectra to Array          */
/**********************************************************/
void eclRetrievePrefissionSpectrum(const int c0, double *pf)
{
  for(int k=0 ; k<max_energybin; k++){
    pf[k] = dat[c0].spec[neutron][k];
  }
}


/**********************************************************/
/*      Copy Fission Probabilities to Array               */
/**********************************************************/
void eclRetrieveFissionProbability(const int c0, double *pfex)
{
  for(int k=0 ; k<max_energybin; k++){
    pfex[k] = (k <= ktot[c0]+1) ? dat[c0].pfis[k] : 0.0;
  }
}


/**********************************************************/
/*      Add Binary Reaction Component if not Separated    */
/**********************************************************/
void eclAddBinaryReaction(const int nm, const int cm, EXSpectra *dat)
{
  for(int n=0 ; n<nm ; n++){
    if(!dat[n].getGcon()){
      for(int k=0 ; k<=ktot[n]+1 ; k++){
        for(int c=1 ; c<cm ; c++){
          if(n == c) dat[n].spec[c][k] += dat[n].sbin[k];
        }
      }
    }
  }
}


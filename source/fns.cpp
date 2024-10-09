/******************************************************************************/
/*  fns.cpp                                                                   */
/*        prompt fission neutron spectrum calculation                         */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "nucleus.h"
#include "masstable.h"
#include "fns.h"
#include "fnssystematics.h"
#include "output.h"
#include "global.h"
#include "kcksyst.h"
#include "levden.h"
#include "eclipse.h"
#include "terminate.h"


static void   fnsCheckParameters (const int, FChance *);
static void   fnsEachChanceFissionSpectrum (Pdata *, FChance *, Nucleus *);
static void   fnsFissionFragmentZA (const int, const int, FChance *);
static void   fnsChanceFissionZA (const double, FChance *, FChance *);
static void   fnsFissionEnergy (const double, const double, const unsigned int, Nucleus *, FChance *);
static void   fnsLevelDensityParameter (Nucleus *, FChance *);
static void   fnsEnergySharing (FChance *);
static void   fnsMaximumTemperature (FChance *);
static int    fnsOutgoingEnergyGrid (const int, double *);

static inline double fnsAverageSn (const int, const int, const double);
static inline void   fnsCleanSpecArray (void);

static const int MAX_FNSPEC = 5;
static       int ne = 0;
static double *eout, *spec[MAX_FNSPEC];
/*** spec array (fc->cms has the same but for these averages)
           spec[0] : fission spectrum from light fragment
           spec[1] : fission spectrum from heavy fragment
           spec[2] : light and heavy average
           spec[3] : pre-fission (Hauser-Feshbach) spectrum
           spec[4] : total neutron spectrum */

#undef DEBUG_EXFISS
#ifdef DEBUG_EXFISS
static void fnsExfissCheck (const int, const double, double *, FNSpec *);
#endif


/**********************************************************/
/*      Main Madland-Nix Model Calculation                */
/**********************************************************/
void fnsFissionNeutronModel(System *sys, Pdata *pdt, double *chi, FNSpec *fns)
{
  double *pf = nullptr;

  fnsCheckParameters(fns->maxchance(), fns->fc);

  try{
    eout  = new double [MAX_ENERGY_BIN];
    pf    = new double [fns->maxchance()];
    for(int j=0 ; j<MAX_FNSPEC ; j++){
      spec[j] = new double [MAX_ENERGY_BIN];
    }
  }
  catch(std::bad_alloc &e){
    message << "memory allocation error " << e.what();
    cohTerminateCode("fnsFissionNeutronModel");
  }

  for(int j=0 ; j<MAX_ENERGY_BIN ; j++)  chi[j] = 0.0;

  /*** copy neutron data */
  Pdata pn = *pdt;
  if(fns->omindex == 0) fns->omindex = 6; // default OMP, Koning-Delaroche
  pn.omp = fns->omindex;                  // replace OM index

  /*** determine the outgoing energy grid */
  ne = fnsOutgoingEnergyGrid(1,eout);

  /*** check fragment ZA for the first chance fission */
  fnsFissionFragmentZA(sys->compound.getZ(),sys->compound.getA(),&fns->fc[0]);

  /*** total fission cross section */
  double ftot = 0.0;
  for(int i=0 ; i<sys->max_compound ; i++) ftot += crx.prod[i].fiss;

  /*** follow neutron emission chain */
  int cid = 0; // CN index
  int mch = 0; // chance fission
  double nutot = 0.0;
  do{
    /*** fission probability at each fission chance
         if ftot = 0, only the first chance fission is calculated */
    if(ftot == 0.0) pf[mch] = (mch == 0) ? 1.0 : 0.0;
    else            pf[mch] = crx.prod[cid].fiss / ftot;
    if(pf[mch] <= 0.0) break;

    fnsCleanSpecArray();

    /*** representative fission fragments */
    if(mch > 0) fnsChanceFissionZA(pn.mass_excess,&fns->fc[mch-1],&fns->fc[mch]);

    /*** normalized pre-fission evaporation neutrons */
    double ep = 0.0;
    if(mch > 0){
      eclRetrievePrefissionSpectrum(cid,spec[3]);
      ep = fnsPreFissionSpectrum(ne,eout,spec[3],spec[4],&ncl[cid]);
      ep *= mch;
    }
    fns->fc[mch].ecms[3] = ep;

    /*** ncl.popfis is fission probabilities at each extation energy, 
         here we calculate the average energy of that. 
         note: normaliation not really required but do it for plotting */
    fnsNormalizeExcitationDist(ncl[cid].ntotal,ncl[cid].excitation,ncl[cid].popfis);
    fns->fc[mch].exfiss = fnsAverageSpectrumEnergy(ncl[cid].ntotal,ncl[cid].excitation,ncl[cid].popfis);

    /*** TKE, total fission energy, Egamma, and TXE  */
    fnsFissionEnergy(pn.mass_excess,ep,sys->incident.za.getA(),&ncl[cid],&fns->fc[mch]);
    if(fns->fc[mch].txe < 0.0) break;

    /*** main Madland-Nix model calculation */
    fnsEachChanceFissionSpectrum(&pn,&fns->fc[mch],&ncl[cid]);

    /*** add CN and fission spectrum, weighted by the number of neutrons */
    double nt = fns->fc[mch].nubar + mch;
    for(int i=0 ; i<ne ; i++){
      spec[4][i] = (fns->fc[mch].nubar*spec[2][i] + mch*spec[3][i])/nt;
      chi[i]    += spec[4][i] * pf[mch] * nt;
    }
    nutot += nt * pf[mch];

    /*** spectrum LAB energy */
    fns->fc[mch].elab    = fnsAverageSpectrumEnergy(ne,eout,spec[4]);

    /*** average prefission neutron energy */
    fns->fc[mch].ecms[4] = (fns->fc[mch].nubar*fns->fc[mch].ecms[2] + mch*fns->fc[mch].ecms[3])/nt;

    if(prn.system)   outChanceFispecParameter(mch,&fns->fc[mch],&ncl[cid]);
    if(prn.xsection) outChanceFispecEnergy(fns->fc[mch].ecms);
    if(prn.spectra)  outChanceFispec(ne,eout,spec,&fns->fc[mch],&ncl[cid]);

    /*** move to the next compound after neutron emission */
    mch ++;
    if(mch >= fns->maxchance()) break;

    int cnext = ncl[cid].cdt[neutron].next;
    if(cnext < 0) break;
    cid = cnext;

  }while(cid < sys->max_compound);

  /*** renormalize Chi, and calculate average energy */
  fnsNormalizeSpectrum(ne,eout,chi);
  double eavrg = fnsAverageSpectrumEnergy(ne,eout,chi);

  if(prn.system)    outFissionNeutronEnergy(mch-1,eavrg,nutot,pf,fns);
  if(prn.spectra)   outFissionNeutronSpectrum(ne,eout,chi,fns->maxwell);
  if(prn.fisenergy) outFissionNeutronEnergyTable(mch-1,pf,fns);

#ifdef DEBUG_EXFISS
  fnsExfissCheck(mch, sys->lab_energy, pf, fns);
#endif


  delete [] eout;
  delete [] pf;
  for(int j=0 ; j<MAX_FNSPEC ; j++){
    delete [] spec[j];
  }
}


/**********************************************************/
/*      Check Input Parameters                            */
/**********************************************************/
void fnsCheckParameters(const int nf, FChance *fc)
{
  for(int i=0 ; i<nf ; i++){
    if(fc[i].rt < 0.0){
      message << i << "-th Rt parameter " << fc[i].rt << " negative";
      cohTerminateCode("fnsCheckParameters");
    }
    if(fc[i].nuratio < 0.0){
      message << i << "-th Nu-ratio parameter " << fc[i].nuratio << " negative";
      cohTerminateCode("fnsCheckParameters");
    }
    if(fc[i].tps < 1.0 || fc[i].tps >= 2.0){
      message << i << "-th T-distribution parameter " << fc[i].tps << " out of range";
      cohTerminateCode("fnsCheckParameters");
    }
    if(fc[i].anisotropy < 0.0){
      message << i << "-th anisotropy parameter " << fc[i].anisotropy << " negative";
      cohTerminateCode("fnsCheckParameters");
    }
  }
}


/**********************************************************/
/*      For Each Chance-Fission                           */
/**********************************************************/
void fnsEachChanceFissionSpectrum(Pdata *pn, FChance *fc, Nucleus *n)
{
  /*** level density parameters */
  fnsLevelDensityParameter(n, fc);

  /*** maximum temperature in each fragment */
  fnsMaximumTemperature(fc);

  /*** calculate inverse reaction cross sections */
  fnsInverseCrossSection(pn,&fc->lf.za,&fc->hf.za);

  /*** main MN model calculation, and normalize spec to SigF */
  fnsMadlandNixModel(n->za.getA(),ne,eout,spec,fc);

  /*** calculate nu-bar */
  fc->nubar = fnsNubar(fc->txe,fc->ecms[2],fc->meansep,fc->egamma);
}


/**********************************************************/
/*      Representative Fission Fragments                  */
/**********************************************************/
void fnsFissionFragmentZA(const int zt, const int at, FChance *fc)
{
  /*** initial case
       if not given, assume symmetric fission */
  if( (fc->lf.za.getZ() == 0) || (fc->hf.za.getZ() == 0) ||
      (fc->lf.za.getA() == 0) || (fc->hf.za.getA() == 0) ){

    int zl = (int)(zt/2.0);
    int zh = zl;
    if( (zl + zh) < zt ) zh++;

    int al = (int)(at/2.0);
    int ah = al;
    if( (al + ah) < at ) ah++;

    fc->lf.za.setZA(zl,al);
    fc->hf.za.setZA(zh,ah);
  }
}


/**********************************************************/
/*      Fission Fragments at Each Chance Fission          */
/**********************************************************/
void fnsChanceFissionZA(const double nmass, FChance *fc0, FChance *fc1)
{
  /*** if not given,
       each time, subtract a neutron from a fragment whose Sn is smaller */
  if( (fc1->lf.za.getZ() == 0) || (fc1->hf.za.getZ() == 0) ||
      (fc1->lf.za.getA() == 0) || (fc1->hf.za.getA() == 0) ){
    double sl =  mass_excess(fc0->lf.za.getZ(),fc0->lf.za.getA()-1) + nmass
               - mass_excess(fc0->lf.za.getZ(),fc0->lf.za.getA()  );
    double sh =  mass_excess(fc0->hf.za.getZ(),fc0->hf.za.getA()-1) + nmass
               - mass_excess(fc0->hf.za.getZ(),fc0->hf.za.getA()  );

    if(sl < sh){
      fc1->lf.za.setZA(fc0->lf.za.getZ(), fc0->lf.za.getA()-1);
      fc1->hf.za.setZA(fc0->hf.za.getZ(), fc0->hf.za.getA()  );
    }
    else{
      fc1->lf.za.setZA(fc0->lf.za.getZ(), fc0->lf.za.getA()  );
      fc1->hf.za.setZA(fc0->hf.za.getZ(), fc0->hf.za.getA()-1);
    }
  }
}


/**********************************************************/
/*      Fission Energy Balance and TXE                    */
/**********************************************************/
void fnsFissionEnergy(const double nmass, const double eavepf, const unsigned int ap, Nucleus *n, FChance *fc)
{
  /*** separation energy for calculating equivalent excitation energy */
  double sn = mass_excess(n->za.getZ(),n->za.getA()-1) + nmass - mass_excess(n->za.getZ(),n->za.getA());

  /*** for the photon-induced case and no multi-chance case, set Sn = 0 */
  if(ap == 0) sn = 0.0;

  /*** check TKE */
  if(fc->tke0 == 0.0){
    fc->tke0 = fnsTKESystematics(n->za.getZ(),n->za.getA());
    if(fc->tke1 == 0.0) fc->tke1 = fnsTKEEdepSystematics(n->za.getA());
    fc->tke = fnsTKEEnergyDependence(fc->tke0,fc->tke1,sn,n->ntotal,n->excitation,n->popfis);
  }
  else if(fc->tke0 == -1.0){
    fc->tke = fnsTKEEnergyDependenceCGMF(n->za.getZ(),n->za.getA(),sn,n->ntotal,n->excitation,n->popfis);
  }
  else{
    if(fc->tke1 != 0.0) fc->tke = fnsTKEEnergyDependence(fc->tke0,fc->tke1,sn,n->ntotal,n->excitation,n->popfis);
    else fc->tke = fc->tke0;
  }

  /*** check total energy release */
  if(fc->etotal == 0.0){
    double mcn = mass_excess(n->za.getZ()    , n->za.getA()    );
    double mlf = mass_excess(fc->lf.za.getZ(), fc->lf.za.getA());
    double mhf = mass_excess(fc->hf.za.getZ(), fc->hf.za.getA());
    fc->etotal = mcn - mlf - mhf;
  }

  /*** check Egamma */
  if(fc->egamma == 0.0) fc->egamma = fnsEgammaSystematics(n->za.getA());

  /*** average separation energy */
  if(fc->meansep == 0.0){
    double s1 = fnsAverageSn(fc->lf.za.getZ(),fc->lf.za.getA(),nmass);
    double s2 = fnsAverageSn(fc->hf.za.getZ(),fc->hf.za.getA(),nmass);
    fc->meansep = (s1+s2)*0.5;
  }

  /*** total excitation energy, TXE */
  fc->txe = fc->etotal + n->max_energy - fc->tke - eavepf;
}


/**********************************************************/
/*      Level Densities of FF at Each Chance Fission      */
/**********************************************************/
void fnsLevelDensityParameter(Nucleus *n, FChance *fc)
{
  /*** level density parameters for fission fragments */
  if(fc->lf.a == 0.0 || fc->hf.a == 0.0){
    kckDataRead(&fc->lf.za, &fc->lf.ldp);
    fc->lf.ldp.a      = kckAsymptoticLevelDensity(fc->lf.za.getA());
    fc->lf.ldp.sigma0 = sqrt(kckSpinCutoff(fc->lf.za.getA()));

    kckDataRead(&fc->hf.za, &fc->hf.ldp);
    fc->hf.ldp.a      = kckAsymptoticLevelDensity(fc->hf.za.getA());
    fc->hf.ldp.sigma0 = sqrt(kckSpinCutoff(fc->hf.za.getA()));

    fnsEnergySharing(fc);

    fc->ac = ldDensityParameter(n->max_energy,(double)n->za.getA(),&n->ldp);
  }
}


/**********************************************************/
/*      Calculate a-parameters for Given TXE and RT       */
/**********************************************************/
void fnsEnergySharing(FChance *fc)
{
  const double eps = 1e-6;

  double el = fc->txe * 0.5;
  double eh = fc->txe * 0.5;
  double e0 = el/eh;

  double al = ldDensityParameter(el, (double)fc->lf.za.getA(), &fc->lf.ldp);
  double ah = ldDensityParameter(eh, (double)fc->hf.za.getA(), &fc->hf.ldp);
  double p  = al/ah * fc->rt*fc->rt;

  eh = (fc->txe - (fc->lf.ldp.pairing_energy + fc->hf.ldp.pairing_energy))/(1.0+p) + fc->hf.ldp.pairing_energy;
  el = fc->txe - eh;

  double e1 = el/eh;
  for(int itr = 0; itr < 10 ; itr++){
    al = ldDensityParameter(el, (double)fc->lf.za.getA(), &fc->lf.ldp);
    ah = ldDensityParameter(eh, (double)fc->hf.za.getA(), &fc->hf.ldp);
    p  = al/ah * fc->rt*fc->rt;
    eh = (fc->txe - (fc->lf.ldp.pairing_energy + fc->hf.ldp.pairing_energy))/(1.0+p) + fc->hf.ldp.pairing_energy;
    el = fc->txe - eh;
    e1 = el/eh;

    if( fabs(e0-e1) < eps ){ break; }
    e0 = e1;
  }

  fc->lf.a = al;
  fc->hf.a = ah;
}


/**********************************************************/
/*      Maximum Temperature                               */
/**********************************************************/
void fnsMaximumTemperature(FChance *fc)
{
  fc->tmax    = sqrt(fc->txe/fc->ac);
  fc->hf.tmax = sqrt(fc->ac / (fc->lf.a*fc->rt*fc->rt + fc->hf.a)) * fc->tmax;
  fc->lf.tmax = fc->hf.tmax * fc->rt;
}


/**********************************************************/
/*      Secondary Neutron Energy Grid                     */
/**********************************************************/
int fnsOutgoingEnergyGrid(const int m, double *energy)
{
  double de = 1e-6, emax = MAX_SPECTRUM_ENERGY;
  int k = 0;

  if(m == 1){
    /*** generate finer energy grid from 1e-5 to 30 MeV */
    bool f = false;
    for(int order=-5 ; order<=2 ; order++){
      double base = pow(10.0,order);
      for(int i=0 ; i<90 ; i++){
        double e = base + de * i;
        energy[k++] = e;
        if( (e >= emax) || (k >= MAX_ENERGY_BIN) ){
          f = true;
          break;
        }
      }
      if(f) break;
      de *= 10.0;
    }
  }
  else if(m == 2){
    /*** outgoing energy grid used in Madland-Nix model */
    bool f = false;
    int  step = 30;
    double e1 = 10.0; // above E1, use const. grid
    for(int order=-5 ; order<=2 ; order++){
      double base = pow(10.0,order);
      de = base*9.0/step;
      for(int i=0 ; i<step ; i++){
        double e = base + de * i;
        energy[k++] = e;
        if( (e >= e1) || (k >= MAX_ENERGY_BIN) ){
          f = true;
          break;
        }
      }
      if(f) break;
    }
    de = 0.2;
    for(int i=1 ; ; i++){
      double e = e1 + i*de;
      energy[k++] = e;
      if( (e >= emax) || (k >= MAX_ENERGY_BIN) ) break;
    }
  }
  else{
    /*** up to 100keV,  step of 1,2, and 5 */
    for(int order=-5 ; order<-1 ; order++){
      double base = pow(10.0,order);
      energy[k++] = base;
      energy[k++] = base*2;
      energy[k++] = base*5;
    }
    double e = 0.1;
    de = 0.1;
    for(int i=0 ; i<300 ; i++){
      energy[k++] = e;
      if( (e >= emax) || (k >= MAX_ENERGY_BIN) ) break;
      e += de;
    }
  }

  return(k);
}


/**********************************************************/
/*      Average Neutron Separation Energy                 */
/**********************************************************/
inline double fnsAverageSn(const int z, const int a, const double mneutron)
{
  /*** average Sn, S2n, S3n, and S4n */
  double m0 = mass_excess(z,a  );
  double m1 = mass_excess(z,a-1);  double s1 = m1 + mneutron - m0;
  double m2 = mass_excess(z,a-2);  double s2 = m2 + mneutron - m1;
  double m3 = mass_excess(z,a-3);  double s3 = m3 + mneutron - m2;
  double m4 = mass_excess(z,a-4);  double s4 = m4 + mneutron - m3;

  return( (s1+s2+s3+s4)/4. );
}


/**********************************************************/
/*      Set Zero in Spectrum Array                        */
/**********************************************************/
inline void fnsCleanSpecArray(void)
{
  for(int i=0 ; i<MAX_FNSPEC ; i++){
    for(int j=0 ; j<MAX_ENERGY_BIN ; j++) spec[i][j] = 0.0;
  }
}


/**********************************************************/
/*      Debug Routine                                     */
/**********************************************************/
#ifdef DEBUG_EXFISS
#include <iomanip>
void fnsExfissCheck (const int m, const double elab, double *pf, FNSpec *fns)
{
  int n = ncl[0].ntotal;
  double *x[m], *y[m];

  for(int i=0 ; i<m ; i++){
    x[i] = new double [n];
    y[i] = new double [n];
    for(int j=0 ; j<n ; j++) x[i][j] = y[i][j] = 0.0;
  }

  int c = 0;
  for(int i=0 ; i<m ; i++){
    if(pf[i] <= 0.0) break;
    for(int j=0 ; j<ncl[c].ntotal ; j++){
      x[i][j] = ncl[c].excitation[j];
      y[i][j] = ncl[c].popfis[j];
    }
    c = ncl[c].cdt[neutron].next;
    if(c < 0) break;
  }

  std::cout.setf(std::ios::scientific, std::ios::floatfield);
  std::cout << setprecision(5);

  std::cout << "# " << elab << std::endl;
  std::cout << "#";
  for(int i=0 ; i<m ; i++){
    std::cout << std::setw(13) << pf[i];
    std::cout << std::setw(13) << fns->fc[i].exfiss;
  }
  std::cout << std::endl;

  for(int j=0 ; j<n ; j++){
    for(int i=0 ; i<m ; i++){
      std::cout << std::setw(13) << x[i][j];
      std::cout << std::setw(13) << y[i][j];
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;


  for(int i=0 ; i<m ; i++){  
    delete [] x[i];
    delete [] y[i];
  }
}
#endif

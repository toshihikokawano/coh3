/******************************************************************************/
/*  output.cpp                                                                */
/*        main data output                                                    */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "physicalconstant.h"
#include "structur.h"
#include "output.h"
#include "parameter.h"
#include "nucleus.h"
#include "elements.h"
#include "global.h"


/**********************************************************/
/*      Global Parameters                                 */
/**********************************************************/

extern std::string version;

static std::string particle_name[9]={
   "      gamma","    neutron","     proton","      alpha",
   "   deuteron","     triton","     helion","    fission",
   "    unknown"};

static std::string p_name = "gnpadth";

static double labE = 0.0;   // laboratory energy
static double sigR = 0.0;   // total reaction cross section
static ZAnumber targZA(0,0);

void outSetSigmaReaction(const int n)
{
  sigR = 0.0;
  for(int i=0 ; i<n ; i++) if(crx.prod[i].xsec > 0.0) sigR += crx.prod[i].xsec;
  if(sigR == 0.0) sigR = crx.reaction;
}

double outRetrieveLabE(void){ return labE; }

#include "outbanner.h"
#include "outformat.h"


/**********************************************************/
/*      Header                                            */
/**********************************************************/
void outTitle(char *str)
{
#ifdef PRINT_BANNER
  outBanner1();
#endif
  outSectionHead(&version[0]);
  std::cout <<cline;
  std::cout << &str[10] << std::endl;
}


/**********************************************************/
/*      Print Main Banner                                 */
/**********************************************************/
void outBanner()
{
  static int val[] = {0,0,2,0,2,2,4,5,6,7,8,9,10};
  time_t c;

  /*** get today */
  time(&c);
  struct tm *p = localtime(&c);
  int year  = p->tm_year + 1900;
  int month = p->tm_mon  + 1;
  int day   = p->tm_mday;

  /*** calculate moon age */
  int a = ( (year - 11)%19 ) * 11;
  int b = val[month];
  int x = (a + b + day)%30;

  /*** print banner */
  outSectionHead(&version[0]);
  if     (x <  8)  std::cout <<banner1 << std::endl;
  else if(x < 15)  std::cout <<banner2 << std::endl;
  else if(x < 23)  std::cout <<banner3 << std::endl;
  else             std::cout <<banner4 << std::endl;

  std::string bar = "#";
  for(int i=0 ; i<DisplayWidth-12 ; i++) bar += "#";
  std::cout << bar << " Moonage" << std::setw(3) << x << std::endl;
}


/**********************************************************/
/*      Write Z and A                                     */
/**********************************************************/
void outZA(ZAnumber *za)
{
  char element[3];

  if(za->getZ() >= N_ELEMENTS){
    element[0] = element[1] = ' ';
  }
  else{
    strncpy(element,element_name[za->getZ()].c_str(),2);
    if(strlen(element_name[za->getZ()].c_str()) == 1) element[1] = ' ';
  }
  element[2] = '\0';

  std::cout << "  ";
  std::cout << std::setw(3) << std::setfill('0') << za->getZ() << '-';
  std::cout << std::setw(3) << std::setfill('0') << za->getA();
  std::cout << std::setw(2) << element;
  std::cout << std::setfill(' ');
}


/**********************************************************/
/*      Write Header of Each Section                      */
/**********************************************************/
void outSectionHead(const char *sec)
{
  int n = strlen(sec);

  std::string bar1 = "#";

  if(labE == 0.0){
    for(int i=0 ; i<DisplayWidth-1 ; i++) bar1 += "#";
    std::cout << bar1 << std::endl;
  }else{
    for(int i=0 ; i<DisplayWidth-27 ; i++) bar1 += ".";
    std::cout << bar1 << "/";  outZA(&targZA);  std::cout << " // ";  outVal(10,6,labE); nl();
  }

  std::string bar2 = "#";
  for(int i=0 ; i<DisplayWidth-n-1 ; i++) bar2 += " ";
  std::cout << bar2 << sec << std::endl;
}


/**********************************************************/
/*      System Parameters                                 */
/**********************************************************/
void outSystem(System *sys, bool flag)
{
  /*** save Z, A, E for the first time */
  targZA.setZA(sys->target.getZ(),sys->target.getA());
  labE  = sys->lab_energy;
  if(!flag) return;

  outSectionHead("SYSTEM PARAMETERS");
  std::cout << cline << "  AtomicNum    MassNum" << std::endl;
  std::cout << "#    Target"; outVal(11,sys->target.getZ()); outVal(11,sys->target.getA()); nl();
  std::cout << "#  Incident"; outVal(11,sys->incident.za.getZ()); outVal(11,sys->incident.za.getA()); nl();
  std::cout << cline << "  Ecms[MeV]  Elab[MeV] ReduceMass WaveNumber" << std::endl;
  std::cout << "     Energy"; outVal(sys->cms_energy);   outVal(sys->lab_energy);
                         outVal(sys->reduced_mass); outVal(sys->wave_number); nl();
}


/**********************************************************/
/*      Target State                                      */
/**********************************************************/
void outTargetState(const int targid, const int targlev, double ex)
{
  outSectionHead("TARGET STATE");
  std::cout << cline << " Excitation  Spin" << std::endl;

  outZA(&ncl[targid].za);
  outVal(ex);
  outVal(5,1,ncl[targid].lev[targlev].spin);
  char p = (ncl[targid].lev[targlev].parity < 0) ? '-' : '+';
  std::cout << p << std::endl;
}


/**********************************************************/
/*      Compound Nucleus, Z and A, Binding Energies       */
/**********************************************************/
void outCompound(const int n, Pdata *pdt)
{
  outSectionHead("COMPOUND NUCLEUS DATA");
  std::cout << cline << "    Z     A  MassExess    Ex(max) Separation     N-cont    N-level" << std::endl;

  for(int i=0 ; i<n ; i++){
    outZA(&ncl[i].za);
    std::cout << blank;
    outVal(11,4,ncl[i].mass_excess); outVal(11,4,ncl[i].max_energy);
    std::cout << blank;
    outVal(11,ncl[i].ncont); outVal(11,ncl[i].ndisc); nl();
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(ncl[i].cdt[j].status){
        std::cout << particle_name[j];
        outZA(&ncl[ ncl[i].cdt[j].next ].za);
        outVal(11,4,pdt[j].mass_excess);
        std::cout << blank;
        outVal(11,4,ncl[i].cdt[j].binding_energy); nl();
      }
    }
  }
}


/**********************************************************/
/*      Level Density Parameters of Each Compound         */
/**********************************************************/
void outLevelDensity(const int n, double d0)
{
  outSectionHead("LEVEL DENSITY PARAMETERS");
  if(d0 > 0.0)
    std::cout << cline <<"       a       Pairing     T          E0         Em         Sigma2     Eshell         D0" << std::endl;
  else
    std::cout << cline <<"       a       Pairing     T          E0         Em         Sigma2     Eshell" << std::endl;
  for(int i=0 ; i<n ; i++){
    if(ncl[i].ncont > 0){
      /*** check redundancy */
      bool found = false;
      for(int j=0 ; j<i ; j++){
        if(ncl[i].za == ncl[j].za){
          found = true;
          break;
        }
      }
      if(found) continue;

      outZA(&ncl[i].za);
      outVal(11,4,ncl[i].ldp.a);
      outVal(11,4,ncl[i].ldp.pairing_energy);
      outVal(11,4,ncl[i].ldp.temperature);
      outVal(11,4,ncl[i].ldp.E0);
      outVal(11,4,ncl[i].ldp.match_energy);
      outVal(11,4,ncl[i].ldp.sigma0*ncl[i].ldp.sigma0);
      outVal(11,6,ncl[i].ldp.shell_correct);
      if( (i == 0) && (d0 > 0.0) ) outVal(d0);
      nl();
    }
  }
}


/**********************************************************/
/*      GDR Parameter                                     */
/**********************************************************/
void outGDR(const bool printall, const int n, const double gg)
{
  outSectionHead("GIANT DIPOLE RESONANCE DATA");
  if(gg > 0.0){
    std::cout << cline << "Energy[MeV] Width[MeV] Sigma0[mb] <Gam>[MeV] Profile" << std::endl;
  }
  else{
    std::cout << cline << "Energy[MeV] Width[MeV] Sigma0[mb]            Profile" << std::endl;
  }

  for(int i=0 ; i<n ; i++){
    if(printall){ outZA(&ncl[i].za); nl(); }

    for(int k=0 ; k<MAX_GDR ; k++){
      if(ncl[i].gdr[k].getEnergy() > 0.0){
        std::cout << "         " << ncl[i].gdr[k].getEM() << std::setw(1) << ncl[i].gdr[k].getL();
        std::string prof;
        if(ncl[i].gdr[k].getProfile() == EX){
          prof = "External";
          std::cout << blank << blank << blank << blank << prof;
        }
        else{
          outVal(11,4,ncl[i].gdr[k].getEnergy());
          outVal(11,4,ncl[i].gdr[k].getWidth());
          outVal(11,6,ncl[i].gdr[k].getSigma());

          switch(ncl[i].gdr[k].getProfile()){
          case SL: prof = "SL"; break;
          case GL: prof = "GL"; break;
          case ML: prof = "ML"; break;
          default: prof = "  "; break;
          }
          std::cout << blank << "      " << prof;
        }
        nl();
      }
    }
    if((gg > 0.0) && (i == 0)){
      std::cout << "   Average Gamma Width" << blank << blank;
      outVal(gg); nl();
    }

    if(!printall) break;
  }
}


/**********************************************************/
/*      Fission Barriers                                  */
/**********************************************************/
void outFissionBarrier(const int n)
{
  outSectionHead("FISSION BARRIER PARAMETERS");

  std::cout << cline << blank << "Height[MeV] Width[MeV] KBand[MeV]       Elmax[MeV]" << std::endl;
  for(int i=0 ; i<n ; i++){
    for(int m=0 ; m<MAX_HUMP ; m++){
      if(ncl[i].fissile){
        if(ncl[i].fission->barrier[m].height > 0.0){
          outZA(&ncl[i].za);
          if     (m == 0) std::cout << " 1stBarrier";
          else if(m == 1) std::cout << " 2ndBarrier";
          else if(m == 2) std::cout << " 3rdBarrier";
          else            std::cout << blank;
          outVal(11,4,ncl[i].fission->barrier[m].height);
          outVal(11,4,ncl[i].fission->barrier[m].curvature);
          std::cout << blank << "      ";
          outVal(ncl[i].fission->barrier[m].elmax); nl();

          for(int k=0 ; k<ncl[i].fission->barrier[m].nband ; k++){
            std::cout << blank << blank << blank << blank;
            outVal(ncl[i].fission->barrier[m].kband[k].excitation);
            outVal(5,1,ncl[i].fission->barrier[m].kband[k].k2/2.0);
            char p = (ncl[i].fission->barrier[m].kband[k].parity < 0) ? '-' : '+';
            std::cout << p << std::endl;

            /*** to print band levels */
#ifdef PRINT_BANDLEVELS
            int l0 = ((ncl[i].fission->barrier[m].kband[k].k2 == 0)
                   && (ncl[i].fission->barrier[m].kband[k].parity < 0)) ? 1 : 0;
            int ls = (ncl[i].fission->barrier[m].kband[k].k2 == 0) ? 2 : 1;
            double xk = ncl[i].fission->barrier[m].kband[k].k2/2.0;
            double hb = ncl[i].fission->barrier[m].inertia;
            double e0 = ncl[i].fission->barrier[m].kband[k].excitation;
            for(int l=l0 ; ; l+=ls){
              double xj = xk + l;
              double e = (xj*(xj+1.0)-xk*(xk+1.0))* hb + e0;
              if(e > ncl[i].fission->barrier[m].elmax) break;
              std::cout << blank << blank << blank << blank << blank;
              outVal(5,1,xj);
              std::cout << p;
              outVal(e); nl();
            }
#endif
          }
        }
      }
    }
  }

  if(ncl[0].fissile && ncl[0].fission->potwell[0].height > 0.0){
    nl();
    std::cout << cline << blank << "Height[MeV] Width[MeV] Aborb[MeV]" << std::endl;
    for(int i=0 ; i<n ; i++){
      for(int m=0 ; m<MAX_HUMP-1 ; m++){
        if(ncl[i].fissile){
          if(ncl[i].fission->potwell[m].height > 0.0){
            outZA(&ncl[i].za);
            if     (m == 0) std::cout << " Class-II  ";
            else if(m == 1) std::cout << " Class-III ";
            else            std::cout << blank;
            outVal(11,4,ncl[i].fission->potwell[m].height);
            outVal(11,4,ncl[i].fission->potwell[m].curvature);
            outVal(11,4,ncl[i].fission->potwell[m].absorb); nl();
          }
        }
      }
    }
  }

  bool f1 = false, f2 = false;
  for(int i=0 ; i<n ; i++){
    if(ncl[i].fissile){
      if(ncl[i].fission->fisenhance.energy != 0.0) f1 = true;
      if( (ncl[i].fission->fisenhance.compfact1p != 0.0) || (ncl[i].fission->fisenhance.compfact2p != 0.0) ) f2 = true;
    }
  }

  if(f1){
    nl();
    std::cout << cline << blank << "Energy[MeV] Width[MeV]     Factor" << std::endl;
    for(int i=0 ; i<n ; i++){
      if(ncl[i].fissile){
        outZA(&ncl[i].za); std::cout << " EnhanceFct";
        outVal(11,4,ncl[i].fission->fisenhance.energy);
        outVal(11,4,ncl[i].fission->fisenhance.width);
        outVal(11,4,ncl[i].fission->fisenhance.peak); nl();
      }
    }
  }

  if(f2){
    nl();
    std::cout << cline << blank << "     Factor Damp[/MeV]" << std::endl;
    for(int i=0 ; i<n ; i++){
      if(ncl[i].fissile){

        if(ncl[i].fission->fisenhance.compfact1n > 0.0){
          outZA(&ncl[i].za); std::cout << " CmprsF(+)";
          outVal(11,4,ncl[i].fission->fisenhance.compfact1p);
          outVal(11,4,ncl[i].fission->fisenhance.compfact2p); nl();
          outZA(&ncl[i].za); std::cout << " CmprsF(-)";
          outVal(11,4,ncl[i].fission->fisenhance.compfact1n);
          outVal(11,4,ncl[i].fission->fisenhance.compfact2n); nl();
        }
        else{
          outZA(&ncl[i].za); std::cout << " CmprsFact";
          outVal(11,4,ncl[i].fission->fisenhance.compfact1p);
          outVal(11,4,ncl[i].fission->fisenhance.compfact2p); nl();
        }
      }
    }
  }

}


/**********************************************************/
/*      Cross Section                                     */
/**********************************************************/
void outCrossSection(const int ctl, double s1, double s2, double s3)
{
  outSectionHead("TOTAL REACTION CROSS SECTION");
  if(ctl == 0){
    std::cout << cline << "      Total    Elastic   Reaction" << std::endl;
    std::cout << "  Sigma[mb]";
    outVal(s1); outVal(s2); outVal(s3); nl();
  }
  else if(ctl == 1){
    std::cout << "  Sigma[mb]"; outVal(s3); nl();
  }
  else if(ctl == 2){
    std::cout << cline << "        GDR    Quasi-d   Reaction" << std::endl;
    std::cout << "  Sigma[mb]";
    outVal(s1); outVal(s2); outVal(s3); nl();
  }
}


/**********************************************************/
/*      Reaction Cross Section                            */
/**********************************************************/
void outReaction(const int n, const int targid, double ex)
{
  outSectionHead("INDIVIDUAL GROUND STATE PRODUCTION");
  std::cout << cline << "   Residual  Prod.[mb]  Qval[MeV]   Particles" << std::endl;
  std::cout << cline << blank << blank << blank << "  ";
  for(int j=1 ; j<(int)p_name.length() ; j++){
    std::cout << std::setw(2)<< p_name.substr(j,1);
  }
  nl();

  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].xsec > 0){
      std::cout << " Production";
      outZA(&ncl[i].za);
      outVal(lowfilter(crx.prod[i].xsec));
      outVal(11,5, lowfilter(ex + ncl[i].max_energy - ncl[targid].max_energy));
      std::cout << "  ";

      for(int j=1 ; j<MAX_CHANNEL ; j++) std::cout << std::setw(2) << (int)crx.prod[i].par[j];
      nl();
    }
  }
  std::cout << blank << "        sum";
  outVal(lowfilter(sigR)); nl();
}


/**********************************************************/
/*      Total Residual Nucleus Production Cross Section   */
/**********************************************************/
void outTotalResidual(CumulativeResidualProduct *res)
{
  outSectionHead("STABLE AND LONG-LIVED STATE PRODUCTIONS");
  if(res->getNcurrent() == 0) return;

  std::cout << "#       Residual    Ex[MeV]     Jpi     T(1/2)[s]  Prod.[mb] Meta" << std::endl;

  int c = 1;
  double sum = 0.0;
  for(unsigned int z=res->zmin ; z<=res->zmax ; z++){
    for(unsigned int a=res->amin ; a<=res->amax ; a++){
      ZAnumber za(z,a);
      for(int i=0 ; i<res->getNcurrent() ; i++){
        if(res->rp[i].za != za) continue;

        outVal(5,c++);
        outZA(&res->rp[i].za);
        std::cout << "  ";
        outVal(9,5,res->rp[i].energy);
        outVal(7,1,res->rp[i].spin);
        char p0 = (res->rp[i].parity < 0) ? '-' : '+';
        std::cout << p0 << "   ";
        if(res->rp[i].halflife < 0.0) std::cout << dashl;
        else outVal(res->rp[i].halflife);
        outVal(lowfilter(res->rp[i].production));
        if(res->rp[i].metaflag > 0) std::cout << std::setw(5) << res->rp[i].metaflag;
        nl();

        if(res->rp[i].metaflag == 0) sum += res->rp[i].production;
      }
    }
  }
  std::cout << cline << blank << blank << "         g.s.sum";
  outVal(lowfilter(sum)); nl();
  nl(); nl();
}


/**********************************************************/
/*      Isomeric Ratios for Long-Lived Nuclides           */
/**********************************************************/
void outIsomericRatio(CumulativeResidualProduct *res)
{
  outSectionHead("ISOMERIC RATIOS");
  if(res->getNcurrent() == 0) return;

  std::cout << "#       Residual    Ex[MeV]     Jpi     T(1/2)[s] Ratio      m/(g-m)" << std::endl;

  int c = 1;
  for(unsigned int z=res->zmin ; z<=res->zmax ; z++){
    for(unsigned int a=res->amin ; a<=res->amax ; a++){
      ZAnumber za(z,a);

      /*** look at isomer */
      for(int i1=0 ; i1<res->getNcurrent() ; i1++){
        if(res->rp[i1].za != za) continue;
        if(res->rp[i1].metaflag == 0) continue;

      /*** look for ground state of the same ZA */
        for(int i0=0 ; i0<res->getNcurrent() ; i0++){
          if((res->rp[i0].za == za) && (res->rp[i0].metaflag == 0)){

            /*** production of gs includes meta-production */
            double d  = res->rp[i0].production - res->rp[i1].production;
            double r0 = (res->rp[i0].production == 0.0) ? 0.0 : res->rp[i1].production / res->rp[i0].production;
            double r1 = (d == 0.0) ? 0.0 : res->rp[i1].production / d;

            outVal(5,c++);
            outZA(&za);
            std::cout << "  ";
            outVal(9,5,res->rp[i1].energy);
            outVal(7,1,res->rp[i1].spin);
            char p0 = (res->rp[i1].parity < 0) ? '-' : '+';
            std::cout << p0 << "   ";
            outVal(res->rp[i1].halflife);
            outVal(lowfilter(r0));
            outVal(lowfilter(r1));
            nl();
          }
        }
      }
    }
  }
  nl(); nl();
}


/**********************************************************/
/*      Fission Cross Section                             */
/**********************************************************/
void outFission(const int n)
{
  outSectionHead("FISSON CROSS SECTION");
  std::cout << cline << "   Nucleus Fission[mb] Fiss.Prob." << std::endl;

  double sum = 0.0;
  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].fiss > 0) sum += crx.prod[i].fiss;
  }

  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].fiss > 0){
      int c = crx.prod[i].par[neutron] + 1;
      std::cout << " FisChance" << std::setw(1) << c;
      outZA(&ncl[i].za);
      outVal(lowfilter(crx.prod[i].fiss));
      outVal(11,3,crx.prod[i].fiss/sum*100.0);
      nl();
    }
  }
  std::cout << blank << " TotFission";
  outVal(lowfilter(sum)); nl();
  std::cout << blank << "  Fiss+Prod";
  outVal(lowfilter(sum+sigR)); nl();
}


/**********************************************************/
/*      Particle Production Cross Section                 */
/**********************************************************/
void outParticleProduction(const int n, Channel *cdt, double **spc)
{
  double sum[MAX_CHANNEL], ave[MAX_CHANNEL];

  for(int j=0 ; j<MAX_CHANNEL ; j++) sum[j] = ave[j] = 0.0;

  /*** for particles, use the ground-state production */
  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].xsec > 0){
      for(int j=1 ; j<MAX_CHANNEL ; j++){
        sum[j] += (int)crx.prod[i].par[j] * crx.prod[i].xsec;
      }
    }
  }

  /*** for gamma-ray, use the emission spectrum */
  int kmax = (int)(ncl[0].max_energy/ncl[0].de) +1;
  for(int k=0 ; k<kmax ; k++){
    sum[gammaray] += spc[gammaray][k] * ncl[0].de;
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(cdt[j].status) ave[j] += spc[j][k] * ncl[0].de * k * ncl[0].de;
    }
  }
  
  outSectionHead("PARTICLE PRODUCTION");
  std::cout << cline << "   Particle  Prod.[mb] Multiplic. Emean[MeV]" << std::endl;

  std::cout << " ParticlPrd" << particle_name[0];
  outVal(lowfilter(sum[0]));
  outVal(lowfilter(sum[0] / sigR));
  outVal(lowfilter(ave[0] / sum[0])); nl();

  for(int j=1 ; j<MAX_CHANNEL ; j++){
    if (!cdt[j].status) continue;

    /*** particle multiplicity and average energy */
    std::cout << " ParticlPrd" << particle_name[j];

    double mp = sum[j] / sigR;
    double em = (sum[j] == 0.0) ? 0.0 : ave[j] / sum[j];
    outVal(lowfilter(sum[j])); outVal(lowfilter(mp)); outVal(lowfilter(em)); nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Scattering Angular Distribution                   */
/**********************************************************/
void outAngularDistribution(const int ctl, int n0, int np, int step, ZAnumber *za)
{
  if(crx.theta[0] == 0.0) return;

  const int column = 8;
  int page = 1, m = 0;

  switch(ctl){
  case 0:
    outSectionHead("ELASTIC SCATTERING ANGULAR DISTRIBUTION");
    break;
  case 1:
    outSectionHead("CC ELASTIC AND INELASTIC SCATTERING ANGULAR DISTRIBUTION");
    break;
  case 2:
    outSectionHead("DWBA INELASTIC SCATTERING ANGULAR DISTRIBUTION");
    break;
  case 3:
    outSectionHead("COMPOUND PLUS DIRECT SCATTERING ANGULAR DISTRIBUTION");
    if(np > MAX_ANGDISTLEVELS) np=MAX_ANGDISTLEVELS;
    for(int k=np-1 ; k>=0 ; k--){
      if(crx.angdist[k][0] > 0.0){
        np = k+1;
        break;
      }
    }
    break;
  default:
    break;
  }
  
  if(np > 0) page = (np-n0-1)/column+1;

  for(int p=0 ; p<page ; p++){
    m = (np-n0 >= 1) ? std::min(np-n0-1,(p+1)*column-1) : 0;

    std::cout << cline; outZA(za);  std::cout << std::endl;
    std::cout << "#  Ang[deg]";
    for(int j=p*column+n0 ; j<=m+n0 ; j++) std::cout << std::setw(3) << j << " [mb/sr]";
    nl();

    for(int i=0 ; i<MAX_ANGDIST ; i++){
      int k = (int)crx.theta[i];
      if(k%step == 0){
        outVal(crx.theta[i]);
        for(int j=p*column+n0 ; j<=m+n0 ; j++) outVal(crx.angdist[j][i]);
        nl();
      }
    }
    nl();
    nl();
  }
}


/**********************************************************/
/*      Legendre Coefficients                             */
/**********************************************************/
void outLegendreCoefficient(const int ctl, int n0, int np, int id, ZAnumber *za, double ***cl)
{
  if(crx.theta[0] == 0.0) return;

  const int column = 8;
  int page = 1, m = 0;

  switch(ctl){
  case 0:
    outSectionHead("ELASTIC SCATTERING LEGENDRE COEFFICIENTS");
    break;
  case 1:
    outSectionHead("CC ELASTIC AND INELASTIC LEGENDRE COEFFICIENTS");
    break;
  case 2:
    outSectionHead("DWBA INELASTIC LEGENDRE COEFFICIENTS");
    break;
  case 3:
    outSectionHead("COMPOUND PLUS DIRECT REACTION LEGENDRE COEFFICIENTS");
    if(np > MAX_ANGDISTLEVELS) np=MAX_ANGDISTLEVELS;
    for(int k=np-1 ; k>=0 ; k--){
      if(crx.angdist[k][0] > 0.0){
        np = k+1;
        break;
      }
    }
    break;
  case 4:
    outSectionHead("COMPOUND REACTION LEGENDRE COEFFICIENTS IN CONTINUUM");
  default:
    break;
  }
  
  if(np > 0) page = (np-n0-1)/column+1;

  for(int p=0 ; p<page ; p++){
    if(np-n0 >= 1) m = std::min(np-n0-1,(p+1)*column-1);

    std::cout << cline; outZA(za);  std::cout << std::endl;
    std::cout << "#        L ";
    for(int j=p*column+n0 ; j<=m+n0 ; j++){
      if(ctl == 4) outVal(11,3,cl[id][j][MAX_J]);
      else         std::cout << std::setw(3) << j << "th level";
    }
    nl();

    int jmax = 0;
    for(int j=MAX_J-1 ; j>=0 ; j--){
      bool flag = false;
      for(int k=p*column+n0 ; k<=m+n0 ; k++){
        if(cl[id][k][j] != 0.0){ flag = true; break; }
      }
      if(flag){ jmax = j+1; break; }
    }

    std::cout.setf(std::ios::scientific, std::ios::floatfield);

    for(int j=0 ; j<=jmax ; j++){
      outVal(10,j);
      std::cout << " ";
      for(int k=p*column+n0 ; k<=m+n0 ; k++){
        double x0 = cl[id][k][0];
        double x1 = (x0 > 0.0) ? cl[id][k][j]/x0/(2.0*j+1.0) : 0.0;
        std::cout << std::setprecision(3) << std::setw(11) << lowfilter(x1);
      }
      nl();
    }
    nl();
    nl();
  }
}


/**********************************************************/
/*      Particle Emission Spectra                         */
/**********************************************************/
void outSpectrum(const int ctl, double **spc, Nucleus *n)
{
  int kmax = (int)(n->max_energy/n->de) +1;

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if (!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    int k = (int)(ncl[p].max_energy/ncl[p].de);
    if(k > kmax) kmax = k;
  }

  double gw = parmGetValue(parmBROD);
  if(gw > 0.0) kmax += (int)(3.0*gw/ncl[0].de);

  if     (ctl == 0) outSectionHead("SPECTRA FROM COMPOUND NUCLEUS");
  else if(ctl == 1) outSectionHead("PRECOMPOUND PARTICLE SPECTRA");
  else if(ctl == 2) outSectionHead("TOTAL PARTICLE EMISSION SPECTRA");

  std::cout << cline;
  if(ctl == 0) outZA(&n->za);
  nl();

  std::cout << "# Emin[MeV]  Emax[MeV]";
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    std::cout << std::setw(11) << particle_name[j];
  }
  std::cout << std::endl;

  for(int k=0 ; k<=kmax ; k++){
    if(k == 0){
      outVal(11,4,0.0);
      outVal(11,4,0.5*n->de);
    }
    else{
      outVal(11,4,((double)k-0.5)*n->de);
      outVal(11,4,((double)k+0.5)*n->de);
    }
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(!n->cdt[j].status) continue;
      outVal(lowfilter(spc[j][k]));
    }
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Finer Gamma-Ray Spectra                           */
/**********************************************************/
void outSpectrumFineGamma(double *gc, GammaProduction *gp, Nucleus *n)
{
  double dec = n->de;
  double del = FINE_GRID;

  /*** find the highest bin of the continuum spectrum */
  int mc;
  for(mc=MAX_ENERGY_BIN - 1 ; mc>=0.0 ; mc--) if(gc[mc] != 0.0) break;
  mc ++;

  int ml = (mc * dec) / del; // number of finer bins

  double el0, el1, ec0, ec1, gc0 = 0.0, gc1 = 0.0;
  double *gxc = new double [ml + 1]; // re-binned continuum spec
  double *gxl = new double [ml + 1]; // spectrum for discrete transitions

  /*** adjust continuum spectrum by subtracting discrete gammas */
  for(int kl=0 ; kl<=ml ; kl++) gxc[kl] = 0.0;

  int mc0 = 0;
  for(int kc=0 ; kc<mc-1 ; kc++){
    ec0 = (kc == 0) ? 0.0 : (kc-0.5)*dec;
    ec1 = (kc+0.5)*dec;
    /*** first copy the continuum spectrum x bin width = area */
    gc0 = gc[kc] * dec;

    /*** subtract lines */
    for(int i=0 ; i<gp->getN() ; i++){
      if( (ec0 <= gp->line[i].energy) && (gp->line[i].energy < ec1) ){
        gc0 -= gp->line[i].production;
        mc0 = kc; // remember the highest bin number that was corrected
      }
    }

    if(gc0 < 0.0) gc0 = 0.0; // avoid round-off error

    /*** re-bin the continuum same as the finer grid */
    for(int kl=0 ; kl<=ml ; kl++){
      el0 = (kl == 0) ? 0.0 : (kl-0.5)*del;
      el1 = (kl+0.5)*del;

      if( (ec0 <= el0) && (el1 <= ec1) ){ // when finer bin is inside wider bin
        gxc[kl] += gc0 * del / dec;
      }
      else if( (el0 <= ec0) && (ec0 < el1) ){ // when boundary overlaps
        if(kc != 0) gxc[kl] += gc1 * (ec0 - el0) / dec;
        gxc[kl] += gc0 * (el1 - ec0) / dec;
      }
    }
    gc1 = gc0; // previous data (kc-1)
  }
  mc0 ++;

  /*** bining discrete lines */
  for(int kl=0 ; kl<=ml ; kl++){
    el0 = (kl == 0) ? 0.0 : (kl-0.5)*del;
    el1 = (kl+0.5)*del;
    gxl[kl] = 0.0;
    for(int i=0 ; i<gp->getN() ; i++){
      if( (el0 <= gp->line[i].energy) && (gp->line[i].energy < el1) ) gxl[kl] += gp->line[i].production;
    }
    /*** convert them into /MeV unit */
    gxl[kl] /= del;
    gxc[kl] /= del;
  }

  int mx0 = ml; // highest bin where discrete line data exist
  for(mx0=ml ; mx0>=0 ; mx0--) if(gxl[mx0] != 0.0) break;
  mx0 ++;


  outSectionHead("GAMMARAY SPECTRA FOR DISCRETE TRANSITIONS");

  std::cout << cline;
  outZA(&n->za);
  nl();

  std::cout << "# Emin[MeV]  Emax[MeV] Disc[mb/MeV] Cont[mb/MeV]  Sum[mb/MeV]" << std::endl;

  for(int kl=0 ; kl<=mx0 ; kl++){
    el0 = (kl == 0) ? 0.0 : (kl-0.5)*del;
    el1 = (kl+0.5)*del;
    outVal(11,5,el0);
    outVal(11,5,el1);
    outVal(13,lowfilter(gxl[kl]));
    outVal(13,lowfilter(gxc[kl]));
    outVal(13,lowfilter(gxl[kl] + gxc[kl]));
    nl();
  }

  for(int kc=mc0 ; kc<mc ; kc++){
    ec0 = (kc == mc0) ? el1 : (kc-0.5)*dec;
    ec1 = (kc+0.5)*dec;
    outVal(11,5,ec0);
    outVal(11,5,ec1);
    outVal(13,0.0);
    outVal(13,lowfilter(gc[kc]));
    outVal(13,lowfilter(gc[kc]));
    nl();
  }
  nl();
  nl();

  delete [] gxc;
  delete [] gxl;
}


/**********************************************************/
/*      Sum All Spectrum                                  */
/**********************************************************/
void outSpectrumSum(const int ctl, double **spc, Nucleus *n)
{
  double sum[MAX_CHANNEL],ave[MAX_CHANNEL];
  int    kmax = (int)(n->max_energy/n->de) +1;

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    sum[j] = ave[j] = 0.0;
    if(!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    int k = (int)(ncl[p].max_energy/ncl[p].de);
    if(k > kmax) kmax = k;
  }

  for(int k=0 ; k<kmax ; k++){
    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(n->cdt[j].status){
        sum[j] += spc[j][k] * n->de;
        ave[j] += spc[j][k] * n->de * k*n->de;
      }
    }
  }

  if(ctl == 0)       outSectionHead("SUM SPECTRA FOR COMPOUND");
  else if(ctl == 1)  outSectionHead("SUM OF PRECOMPOUND EMISSION");
  else if(ctl == 2)  outSectionHead("TOTAL PARTICLE EMISSION SPECTRA");

  std::cout << cline;
  if(ctl == 2) std::cout << blank;
  else         outZA(&n->za);

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    std::cout << std::setw(11) << particle_name[j];
  }
  nl();

  std::cout << " SpcSum[mb]" << blank;
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    outVal(lowfilter(sum[j]));
  }
  nl();

  std::cout << " Emean[MeV]" << blank;
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    if(sum[j] > 0.0) ave[j] = ave[j]/sum[j];
    else             ave[j] = 0.0;
    outVal(lowfilter(ave[j]));
  }
  nl();

  if(ctl == 1){
    std::cout << cline << std::endl;
    std::cout << cline << " TotalPreEq   Compound   Fraction" << std::endl;
    std::cout << "  Sigma[mb]";
    outVal(crx.preeq); outVal(crx.reaction); outVal(crx.preeq/crx.reaction); nl();
  }
}


/**********************************************************/
/*      Primary Gamma-Ray Spectrum                        */
/**********************************************************/
void outPrimaryGammaSpectrum(const int ng, double **spc, Nucleus *n)
{
  outSectionHead("PRIMARY GAMMA-RAY SPECTRUM FOR CAPTURE");

  std::cout << cline; outZA(&n->za);  std::cout << std::endl;
  std::cout << "#            Egam[MeV]  Sigma[mb]" << std::endl;

  for(int i=ng-1 ; i>=0 ; i--){
    std::cout << blank;
    outVal(          spc[0][i] );
    outVal(lowfilter(spc[1][i]));
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      All Discrete Gamma Lines                          */
/**********************************************************/
void outDiscreteGamma(GammaProduction *gp)
{
  outSectionHead("DISCRETE GAMMA PRODUCTION");

  std::cout << "#       Z   A  Energy[MeV]   Production" << std::endl;
  for(int i=0 ; i<gp->getN() ; i++){
    outVal(5,i+1);
    outVal(4,gp->line[i].za.getZ());
    outVal(4,gp->line[i].za.getA());
    outVal(13,7,gp->line[i].energy);
    outVal(13,gp->line[i].production);
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Discrete Level Population for Binary Reaction     */
/**********************************************************/
void outDiscreteLevelPopulation(double extot, Nucleus *n)
{
  for(int j=1 ; j<MAX_CHANNEL ; j++){
    if (!n->cdt[j].status) continue;
    int p = n->cdt[j].next;
    if(ncl[p].ndisc == 0) continue;

    outSectionHead("DISCRETE LEVEL POPULATION");
    std::cout << cline; outZA(&ncl[p].za); std::cout << std::endl;
    std::cout << "#   Ex[MeV]     JP      Ecms[MeV]  Prod.[mb]" << std::endl;

    double sum = 0.0;
    for(int i0=0 ; i0<ncl[p].ndisc ; i0++){
      double e = extot - ncl[p].lev[i0].energy - n->cdt[j].binding_energy;
      if(e < 0.0) continue;

      outVal(11,4,ncl[p].lev[i0].energy);
      outVal(7, 1,ncl[p].lev[i0].spin);
      char p0 = (ncl[p].lev[i0].parity < 0) ? '-' : '+';
      std::cout << p0 << "   ";
      outVal(e);
      outVal(lowfilter(ncl[p].lpop[i0])); nl();

      sum += ncl[p].lpop[i0];
    }
    std::cout << blank << blank << "        sum";
    outVal(lowfilter(sum)); nl();
  }
}


/**********************************************************/
/*      Population of Continuum in Compound Nucleus       */
/**********************************************************/
void outPopulation(Nucleus *n)
{
  const int jmax = 6;

  if(n->ncont == 0) return;

  outSectionHead("CONTINUUM POPULATION");
  std::cout << cline; outZA(&n->za);  std::cout << std::endl;
  std::cout << "#   Ex[MeV]";
  for(int j=0 ; j<jmax ; j++) std::cout << "     "<< std::setw(3) << j << "   ";
  std::cout << "   Total   " << std::endl;

  for(int k=0 ; k<n->ncont ; k++){
    outVal(n->excitation[k]);

    double sum = 0.0;
    for(int j=0 ; j<=n->jmax ; j++) sum += n->pop[k][j].even+n->pop[k][j].odd;
    for(int j=0 ; j<jmax ; j++){
      double x = n->pop[k][j].even+n->pop[k][j].odd;
      outVal(lowfilter(x));
    }
    outVal(lowfilter(sum)); nl();
  }
}


/**********************************************************/
/*      Gamma Cascade                                     */
/**********************************************************/
void outGammaCascade(const int c0, Nucleus *n)
{
  outSectionHead("GAMMAS FROM DISCRETE TRANSITION");
  std::cout << cline;
  for(int j=1 ; j<MAX_CHANNEL ; j++) std::cout << std::setw(2) << (int)crx.prod[c0].par[j];
  std::cout << blank;
  outZA(&n->za); nl();
  std::cout << "#   Ex[MeV]     Jpi    FinalState";
  std::cout << "  Branching Production    Egamma Production   T(1/2)[s]" << std::endl;

  double tisom = parmGetValue(parmISOM);
  bool   gcut  = (tisom > 0.0) ? true : false;

  for(int i0=n->ndisc-1 ; i0>=0 ; i0--){
    if(n->lpop[i0] < output_eps && i0 > 0) continue;

    std::cout << std::setw(3) << i0;
    outVal(8,4,n->lev[i0].energy);
    outVal(7,1,n->lev[i0].spin);
    char p0 = (n->lev[i0].parity < 0) ? '-' : '+';
    std::cout << p0 << "   " << blank << blank;
    outVal(lowfilter(n->lpop[i0]));
    if((n->lev[i0].halflife > 0.0) && (lowfilter(n->lpop[i0]) > 0.0)){
      std::cout << blank <<  blank;
      outVal(n->lev[i0].halflife);
    }
    nl();

    if(i0 == 0) continue;
    if( gcut && (n->lev[i0].halflife > tisom) ) continue;
    for(int j=0 ; j<n->lev[i0].ngamma ; j++){
      int    i1 = n->lev[i0].fstate[j];
      double eg = n->lev[i0].energy-n->lev[i1].energy;
      double sg = n->lev[i0].branch[j]*n->lpop[i0];
      if(opt.internalconversion) sg *= n->lev[i0].gratio[j];

      std::cout << blank << blank << std::setw(3) << i1;
      outVal(8,4,n->lev[i1].energy);
      outVal(11,4,n->lev[i0].branch[j]);
      std::cout << "          ";
      outVal(11,4,eg);
      outVal(sg); nl();
    }
  }
}


/**********************************************************/
/*      Isomer Production                                 */
/**********************************************************/
void outIsomerProduction(Nucleus *n)
{
  bool   found = false;
  for(int i0=1 ; i0<n->ndisc ; i0++){
    if(n->lev[i0].halflife > thalfmin && lowfilter(n->lpop[i0]) > 0.0) found = true;
  }
  if(!found) return;

  outSectionHead("ISOMERIC STATE PRODUCTION");
  std::cout << cline; outZA(&n->za); std::cout << std::endl;
  std::cout << "#   Ex[MeV]     Jpi     T(1/2)[s]  Prod.[mb]" << std::endl;

  for(int i0=1 ; i0<n->ndisc ; i0++){
    if(n->lev[i0].halflife <= thalfmin) continue;

    std::cout << std::setw(3) << i0;
    outVal(8,4,n->lev[i0].energy);
    outVal(7,1,n->lev[i0].spin);
    char p0 = (n->lev[i0].parity < 0) ? '-' : '+';
    std::cout << p0 << "   ";
    outVal(n->lev[i0].halflife);
    outVal(lowfilter(n->lpop[i0])); nl();
  }
}


/**********************************************************/
/*      DSD Capture Cross Section                         */
/**********************************************************/
void outDSD(double sigt, double sigd, double sigc, Dcapt *dsd)
{
  outSectionHead("DSD CROSS SECTION");
  std::cout << cline << "       Real       Imag   NonLocal" << std::endl;
  std::cout << "    V1[MeV]";
  outVal(11,3, dsd->v1); outVal(11,3, dsd->w1); outVal(11,5, dsd->nonloc); nl();
  std::cout << cline << std::endl;
  std::cout << cline << "      Total     Direct Semidirect" << std::endl;
  std::cout << "  Sigma[mb]";
  outVal(sigt); outVal(sigd); outVal(sigc); nl();
}


/**********************************************************/
/*      Exciton Model Parameters                          */
/**********************************************************/
void outExcitonParameter(double **spd, double **m2, double **m2r)
{
  const int npmax = 5;

  outSectionHead("PRE-COMPOUND PARAMETERS");
  std::cout << cline << "   particle   gz[/MeV]   gn[/MeV] delta[MeV]  Vdep[MeV]" << std::endl;

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(ncl[0].cdt[j].status){
      std::cout << " PreEq     " << particle_name[j];
      outVal(spd[j][0]); outVal(spd[j][1]); outVal(spd[j][2]); outVal(spd[j][3]); nl();
    }
  }

  std::cout << cline << blank;
  for(int i=0 ; i<npmax ; i++) std::cout << std::setw(11) << i;
  std::cout << std::endl;
  std::cout <<  " M2 [/hbar]"; outVal(0); outVal(m2[0][0]); nl();

  std::cout <<  " M2R[/hbar]" << blank;
  outVal(m2r[0][0]); nl();

  for(int i=1 ; i<npmax ; i++){
    std::cout << blank << std::setw(11) << i;
    for(int j=0 ; j<=i ; j++) outVal(m2[i][j]);
    nl();

    std::cout << blank << blank;
    for(int j=0 ; j<=i ; j++) outVal(m2r[i][j]);
    nl();
  }
}


/**********************************************************/
/*      Fission Neutron Spectrum Model Parameters         */
/**********************************************************/
void outChanceFispecParameter(const int m, FChance *fc, Nucleus *n)
{
  outSectionHead("MULTICHANCE FISSION NEUTRON MODEL PARAMETERS");

  std::cout << "#";
  switch(m){
  case  0: std::cout << " 1st"; break;
  case  1: std::cout << " 2nd"; break;
  case  2: std::cout << " 3rd"; break;
  default: std::cout << " "<<std::setw(1) << m << "th"; break;
  }

  std::cout <<"          LightFF    HeavyFF   Compound" << std::endl;
  std::cout << blank;
  outZA(&fc->lf.za);
  outZA(&fc->hf.za);
  outZA(&n->za);
  nl();

  std::cout << "  LDP[/MeV]";
  outVal(11,4,fc->lf.a);
  outVal(11,4,fc->hf.a);
  outVal(11,4,fc->ac);
  nl();

  std::cout << "  Tmax[MeV]";
  outVal(11,4,fc->lf.tmax);
  outVal(11,4,fc->hf.tmax);
  outVal(11,4,fc->tmax);
  nl();

  std::cout << cline << std::endl;
  std::cout << cline << " Efiss[MeV]   TKE[MeV]   TXE[MeV]  Egam[MeV] Exfis[MeV]" << std::endl;
  std::cout << "     Energy";
  outVal(fc->etotal); outVal(fc->tke); outVal(fc->txe); outVal(fc->egamma); outVal(fc->exfiss);
  nl();

  std::cout << cline << std::endl;
  std::cout << cline << "         Rt   Nu-Ratio     tpdf_s anisotropy" << std::endl;
  std::cout << "  Parameter";
  outVal(11,4,fc->rt); outVal(11,4,fc->nuratio);
  outVal(11,4,fc->tps); outVal(11,4,fc->anisotropy);
  nl();
}


/**********************************************************/
/*      Multi-Chance Fission Neutron Spectrum             */
/**********************************************************/
void outChanceFispec(const int m, double *e, double **s, FChance *fc, Nucleus *n)
{
  outSectionHead("MULTICHANCE FISSION NEUTRON SPECTRUM");
  std::cout << cline; outZA(&fc->lf.za); outZA(&fc->hf.za);
  std::cout << blank; outZA(&n->za); nl();
  std::cout << "#    E[MeV]    LightFF    HeavyFF    Average   Compound      Total" << std::endl;

  for(int k=0 ; k<m ; k++){
    outVal(e[k]);
    for(int i=0 ; i<5 ; i++) outVal(lowfilter(s[i][k]));
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Multi-Chance Fission Neutron Average Energy       */
/**********************************************************/
void outChanceFispecEnergy(double *ecms)
{
  outSectionHead("MULTICHANCE FISSION NEUTRON ENERGIES");

  std::cout << cline << "    LightFF    HeavyFF    Average   Compound      Total" << std::endl;
  std::cout << "  Ecms[MeV]";
  for(int i=0 ; i<5 ; i++) outVal(11,4,ecms[i]);
  nl();
}


/**********************************************************/
/*      Total Fission Neutron Energy (<E> and Nu-bar)     */
/**********************************************************/
void outFissionNeutronEnergy(const int mc, double e, double nu, double *pf, FNSpec *fns)
{
  outSectionHead("AVERAGE FISSION NEUTRON ENERGY AND MULTIPLICITY");

  std::cout << cline << "  FissProb. Epref[MeV]     Nu-bar Efiss[MeV]   TKE[MeV]   TXE[MeV]  Egam[MeV] Exfis[MeV]" << std::endl;

  const int neave = 5;
  double eave[neave];
  for(int j=0 ; j<neave ; j++) eave[j] = 0.0;

  for(int i=0 ; i<=mc ; i++){
    std::cout << " FisChance" << std::setw(1) <<i+1;
    outVal(pf[i]);
    outVal(11,4,fns->fc[i].ecms[3]);
    outVal(11,4,fns->fc[i].nubar + i);
    outVal(11,4,fns->fc[i].etotal);
    outVal(11,4,fns->fc[i].tke);
    outVal(11,4,fns->fc[i].txe);
    outVal(11,4,fns->fc[i].egamma);
    outVal(11,4,fns->fc[i].exfiss);
    nl();

    int j = 0;
    eave[j++] += pf[i] * fns->fc[i].etotal;
    eave[j++] += pf[i] * fns->fc[i].tke;
    eave[j++] += pf[i] * fns->fc[i].txe;
    eave[j++] += pf[i] * fns->fc[i].egamma;
    eave[j]   += pf[i] * fns->fc[i].exfiss;
  }

  std::cout << "    Average" << blank;
  outVal(11,4,e);
  outVal(11,4,nu);
  for(int j=0 ; j<neave ; j++) outVal(11,4,eave[j]);
  nl();
}


void outFissionNeutronEnergyTable(const int mc, double *pf, FNSpec *fns)
{
  outVal(labE);
  for(int i=0 ; i<4 ; i ++){
    if(i <= mc) outVal(pf[i]); else{ outVal(11,4,0.0); }
  }
  for(int i=0 ; i<4 ; i ++){
    if(i <= mc) outVal(fns->fc[i].tke); else{ outVal(11,4,0.0); }
  }
  for(int i=0 ; i<4 ; i ++){
    if(i <= mc) outVal(fns->fc[i].exfiss); else{ outVal(11,4,0.0); }
  }
  for(int i=0 ; i<4 ; i ++){
    if(i <= mc) outVal(fns->fc[i].ecms[3]); else{ outVal(11,4,0.0); }
  }
  nl();
}


/**********************************************************/
/*      Total Fission Neutron Spectrum (Chi)              */
/**********************************************************/
void outFissionNeutronSpectrum(const int n, double *e, double *s, double t)
{
  outSectionHead("FISSION NEUTRON SPECTRUM");

  if(t > 0.0){
    std::cout << "#    E[MeV]  Chi[/MeV]   /Maxwell" << std::endl;
    double c = 2.0/sqrt(PI * t*t*t);
    for(int k=0 ; k<n ; k++){
      outVal(e[k]);
      outVal(lowfilter(s[k]));
      double y = c * exp(-e[k]/t) * sqrt(e[k]);
      outVal(lowfilter(s[k] / y));
      nl();
    }
  }
  else{
    std::cout << "#    E[MeV]  Chi[/MeV]" << std::endl;
    for(int k=0 ; k<n ; k++){
      outVal(e[k]);
      outVal(lowfilter(s[k]));
      nl();
    }
  }

  double x0 = 0.0, x1 = 0.0, de;
  for(int k=0 ; k<n-1 ; k++){
    de = e[k+1] - e[k];
    x0 += s[k] * de;
    x1 += s[k] * e[k] * de;
  }
  x1 = (x0 > 0.0) ? x1 / x0 : 0.0;
  nl();

  std::cout << "#   Average";
  outVal(x1);
  nl();
  nl();
}


/**********************************************************/
/*      Adjustable Model Parameters                       */
/**********************************************************/
void outParameter()
{
  if(adj.nparm() <= 0) return;

  int k = 0;
  for(int i=0 ; i<adj.nparm() ; i++){
    if(adj.parm[i].type == parmLCUT) continue;
    if(adj.parm[i].type == parmGSTR) continue;
    if(adj.parm[i].type == parmISOM) continue;
    if(adj.parm[i].type == parmD0SR) continue;
    if(adj.parm[i].type == parmESET) continue;
    k++;
  }
  if(k == 0) return;

  outSectionHead("ADJUSTED PARAMETERS");
  std::cout << "# Parameter    Nucleus   Particle     Factor" << std::endl;

  for(int i=0 ; i<adj.nparm() ; i++){
    if(adj.parm[i].type == parmLCUT) continue;
    if(adj.parm[i].type == parmGSTR) continue;
    if(adj.parm[i].type == parmISOM) continue;
    if(adj.parm[i].type == parmD0SR) continue;
    if(adj.parm[i].type == parmESET) continue;

    int      t  = adj.parm[i].type;
    int      p  = adj.parm[i].particle; if(p > 8) p = 8;
    ZAnumber za = adj.parm[i].za;
    double   f  = adj.parm[i].getfactor();

    std::cout << "    " << pname[t] << "   ";
    switch(t){
    case parmM2  :
    case parmM2R :
    case parmKO  :
    case parmDSDV:
    case parmGDRM1:
    case parmCOLL:
                    std::cout << blank << blank; 
                    break;
    case parmOV  :
    case parmOW  :
    case parmORV :
    case parmORW :
    case parmOAV :
    case parmOAW :
    case parmOC  :
    case parmSD  :  std::cout << blank << particle_name[p];
                    break;
    case parmLD  :
    case parmSPIN:
    case parmPAIR:
    case parmPDST:
    case parmFL1 :
    case parmFL2 :
    case parmFL3 :
                    outZA(&za);   std::cout << blank;
                    break;
    case parmTJ  :  outZA(&za);   std::cout << particle_name[p];
                    break;
    default      :  break;
    }
    outVal(f); nl();
  }
}


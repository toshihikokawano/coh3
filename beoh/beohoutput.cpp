/******************************************************************************/
/*  output.cpp                                                                */
/*        main data output                                                    */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "beoh.h"
#include "beohoutput.h"
#include "nucleus.h"
#include "global.h"
#include "elements.h"
#include "outformat.h"


/**********************************************************/
/*      Global Parameters                                 */
/**********************************************************/

extern std::string version;

static std::string particle_name[7]={
   "      gamma","    neutron","     proton","      alpha",
   "   deuteron","     triton","     helion"};

static std::string p_name = "gnpadth";
static double Qbeta = 0.0;
static ZAnumber targZA(0,0);


/**********************************************************/
/*      Print Banner                                      */
/**********************************************************/
void outBanner()
{
  outSectionHead(&version[0]);
  std::cout
    <<"#    oooooooooo.                     ooooo   ooooo\n"
    <<"#     888     Y8b                     888     888      In spring, it's the dawn.\n"
    <<"#     888     d8P  .oooo.   .ooooo.   888     888      The mountain rim is slowly getting light,\n"
    <<"#     888ooooooB  d88  88b d88   88b  888ooooo888      and thin purple cloud wisps there.\n"
    <<"#     888     Y8b 888oo888 888   888  888     888\n"
    <<"#     888     d8P 888      888   888  888     888                                The Pillow Book\n"
    <<"#    o8888bood8P   Y8bodP'  Y8boodP' o888o   o888o                                  Sei Shonagon\n"
    <<"#\n";
}


/**********************************************************/
/*      Header                                            */
/**********************************************************/
void outTitle(char *str)
{
  outBanner();
  outSectionHead(&version[0]);
  std::cout <<cline;
  std::cout << &str[10] << std::endl;
}


/**********************************************************/
/*      Write  Z and A                                    */
/**********************************************************/
void outZA(ZAnumber *za)
{
  char element[3];

  if(za->getZ() >= N_ELEMENTS){
    element[0] = element[1] = ' ';
  }
  else{
    strncpy(element,element_name[za->getZ()].c_str(),2);
    if(strlen(element_name[za->getZ()].c_str())==1) element[1] = ' ';
  }
  element[2] = '\0';

  std::cout << "  ";
  std::cout << std::setw(3) << std::setfill('0') << za->getZ() << '-';
  std::cout << std::setw(3) << std::setfill('0') << za->getA();
  std::cout << std::setw(2) << element;
  std::cout << std::setfill(' ');
}


/**********************************************************/
/*      Write  Header of Each Section                     */
/**********************************************************/
void outSectionHead(const char *sec)
{
  int n = strlen(sec);
  const int DisplayWidth = 99;

  std::string bar1="#";

  if(Qbeta == 0.0){
    for(int i=0 ; i<DisplayWidth-1 ; i++) bar1 += "#";
    std::cout << bar1 << std::endl;
  }else{
    for(int i=0 ; i<DisplayWidth-27 ; i++) bar1 += ".";
    std::cout << bar1 << "/";  outZA(&targZA);  std::cout << "  //";  outVal(10,6,Qbeta); nl();
  }

  std::string bar2 = "#";
  for(int i=0 ; i<DisplayWidth-n-1 ; i++) bar2 += " ";
  std::cout << bar2 << sec << std::endl;
}


/**********************************************************/
/*      System Parameters                                 */
/**********************************************************/
void outSystem(CalcMode mode, System *sys)
{
  /*** save Z, A, E for the first time */
  targZA.setZA(sys->target.getZ(),sys->target.getA());
  Qbeta  = sys->ex_total;

  outSectionHead("SYSTEM PARAMETERS");
  if(mode == betadecay){
    std::cout << cline << "  AtomicNum    MassNum     Q-beta" << std::endl;
    std::cout << "# Precursor"; outVal(11,sys->target.getZ()); outVal(11,sys->target.getA()); outVal(11,4,sys->ex_total); nl();
    std::cout << "#  Compound"; outVal(11,sys->compound.getZ()); outVal(11,sys->compound.getA()); nl();
  }
  else if(mode == statdecay){
    std::cout << cline << "  AtomicNum    MassNum Excitation";
    if(sys->beta2 != 0.0) std::cout << "      beta2";
    nl();

    std::cout << "#  Compound"; outVal(11,sys->compound.getZ()); outVal(11,sys->compound.getA()); outVal(11,4,sys->ex_total);
    if(sys->beta2 != 0.0) outVal(11,4,sys->beta2);
    nl();
  }
  else if(mode == fissiondecay || mode == fissionspec || mode == cumulativeyield){
    std::cout << cline << "  AtomicNum    MassNum Excitation" << std::endl;
    std::cout << "#  Target  "; outVal(11,sys->target.getZ()); outVal(11,sys->target.getA()); nl();
    std::cout << "#  Compound"; outVal(11,sys->compound.getZ()); outVal(11,sys->compound.getA()); outVal(11,4,sys->ex_total); nl();
  }
}


/**********************************************************/
/*      Parameters fof Fission Fragment Decay Mode        */
/**********************************************************/
void outFissionFragment(const int af, FFragData *fdt)
{
  outSectionHead("FISSION FRAGMENT DECAY PARAMETER");

  std::cout << cline << "Yield Cutoff          "; outVal(11,fdt->ycutoff); nl();
  std::cout << cline; nl();

  std::cout << "# FisChance   Fraction      Rtemp SpinFactor   Zpfac(Z)   Zpfac(N)        TKE    Ex_fiss   Eprefiss" << std::endl;
  for(int n=0 ; n<fdt->getFissionChance() ; n++){
    std::cout << "# " << std::setw(2) << n+1 << std::setw(7) << af - n;
    outVal(11,fdt->mc[n].fraction);
    outVal(11,5,fdt->mc[n].rt);
    outVal(11,5,fdt->mc[n].spinfactor);
    outVal(11,5,fdt->mc[n].ZpFactor[0]);
    outVal(11,5,fdt->mc[n].ZpFactor[1]);
    outVal(11,fdt->mc[n].tke);
    outVal(11,fdt->mc[n].exfis);
    outVal(11,fdt->mc[n].eprefis);
    nl();
  }
  nl();

  std::cout << "# Gaussian Parameters    fraction      width     center"; nl();
  for(int n=0 ; n<fdt->getFissionChance() ; n++){
    for(int i=0 ; i<4 ; i++){
      if(i == 0)  std::cout << "# " << std::setw(9) << n+1;
      else std::cout << cline;
      std::cout << std::setw(11) << i+1;
      outVal(11,5,fdt->mc[n].GaussFract[i]);
      outVal(11,5,fdt->mc[n].GaussSigma[i]);
      outVal(11,5,fdt->mc[n].GaussDelta[i]); nl();
    }
  }
  nl();
}


/**********************************************************/
/*      Target State                                      */
/**********************************************************/
void outTargetState(NuclearStructure *nst, const int isostate)
{
  outSectionHead("PRECURSOR STATE");
  std::cout << cline << " Excitation  Spin" << std::endl;

  outZA(&nst->za);
  outVal(nst->lev[isostate].energy);
  outVal(5,1,nst->lev[isostate].spin);
  char p = (nst->lev[isostate].parity < 0) ? '-' : '+';
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
    if(ncl[i].ncont>0){
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
void outGDR(const bool printall, const int n)
{
  outSectionHead("GIANT DIPOLE RESONANCE DATA");

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

    if(!printall) break;
  }
}


/**********************************************************/
/*      Fission Barriers                                  */
/**********************************************************/
void outFissionBarrier(const int n)
{
  outSectionHead("FISSION BARRIER PARAMETERS");
  std::cout << cline << blank << "Height[MeV] Width[MeV] KBand[MeV]" << std::endl;
  for(int i=0 ; i<n ; i++){
    for(int m=0 ; m<MAX_HUMP ; m++){
      if(ncl[i].fissile){
        if(ncl[i].fission->barrier[m].height > 0.0){
          outZA(&ncl[i].za);
          if     (m==0) std::cout << " 1stBarrier";
          else if(m==1) std::cout << " 2ndBarrier";
          else if(m==2) std::cout << " 3rdBarrier";
          else          std::cout << blank;
          outVal(11,4,ncl[i].fission->barrier[m].height);
          outVal(11,4,ncl[i].fission->barrier[m].curvature); nl();

          for(int k=0 ; k<ncl[i].fission->barrier[m].nband ; k++){
            std::cout << blank << blank << blank << blank;
            outVal(ncl[i].fission->barrier[m].kband[k].excitation);
            outVal(5,1,ncl[i].fission->barrier[m].kband[k].k2/2.0);
            char p = (ncl[i].fission->barrier[m].kband[k].parity < 0) ? '-' : '+';
            std::cout << p << std::endl;
          }
        }
      }
    }
  }
}


/**********************************************************/
/*      Beta Decay Strength                               */
/**********************************************************/
void outBeta(ZAnumber *res, Beta *gts, Beta *ens)
{
  outSectionHead("BETA DECAY FINAL STATE");
  std::cout << "# Daughter     GTstate ENSDFlevel" << std::endl;
  outZA(res);
  outVal(gts->nstate);
  outVal(ens->nstate);
  std::cout << std::endl;
}


/**********************************************************/
/*      Beta Decay Strength Distribution                  */
/**********************************************************/
void outBetaProfile(BetaProfile *bpf)
{
  outSectionHead("BETA STRENGTH DISTRIBUTION");

  std::cout << "# DataFrom " << bpf->source << std::endl;
  std::cout << "# Fraction    Discrete  Continuum" << std::endl;
  std::cout << blank;
  outVal(11,4,bpf->tdisc); outVal(11,4,bpf->tcont);
  std::cout << std::endl;

/*
  if(bpf->ndisc > 0){
    std::cout << "# Discrete" << std::endl;
    std::cout << "#     LeveL     Spin   Population" << std::endl;
    for(int i=0 ; i<bpf->ndisc ; i++){
      outVal(11,4,bpf->lev[i].energy);
      outVal(8,1,bpf->lev[i].spin);
      std::cout << ((bpf->lev[i].parity < 0) ? " - " : " + " );
      outVal(bpf->rdisc[i]);
      std::cout << std::endl;
    }
  }
  if(bpf->ncont > 0){
    std::cout << "# Continuum" << std::endl;
    std::cout << "#      Emin       Emax Population" << std::endl;
    for(int i=bpf->ncont-1 ; i>=0 ; i--){
      outVal(11,3,bpf->excitation[i+1]);
      outVal(11,3,bpf->excitation[i]);
      outVal(bpf->rcont[i]);
      std::cout << std::endl;
    }
  }
*/
}


/**********************************************************/
/*      Ground State Production Rate                      */
/**********************************************************/
void outGSProduction(const int n)
{
  outSectionHead("GROUND STATE PRODUCTION RATE");
  std::cout << cline << "   Residual  Prod.Rate   Particles" << std::endl;
  std::cout << cline << blank << blank << "  ";
  for(int j=1 ; j<(int)p_name.length() ; j++){
    std::cout << std::setw(2)<< p_name.substr(j,1);
  }
  nl();

  double sum=0.0;
  for(int i=0 ; i<n ; i++){
    std::cout << " Production";
    outZA(&ncl[i].za);
    outVal(lowfilter(crx.prod[i].xsec));
    std::cout << "  ";
    sum += crx.prod[i].xsec;

    for(int j=1 ; j<MAX_CHANNEL ; j++) std::cout << std::setw(2) << (int)crx.prod[i].par[j];
    nl();
  }
  std::cout << cline << "        sum";
  outVal(lowfilter(sum)); nl();
}


/**********************************************************/
/*      Fission Rate                                      */
/**********************************************************/
void outFission(int n)
{
  outSectionHead("FISSON RATE");
  std::cout << cline << "   Nucleus FissionRate Fiss.Prob." << std::endl;

  double sum = 0.0;
  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].fiss>0) sum += crx.prod[i].fiss;
  }

  for(int i=0 ; i<n ; i++){
    if(crx.prod[i].fiss>0){
      std::cout << " FisChance" << std::setw(1) <<i+1;
      outZA(&ncl[i].za);
      outVal(lowfilter(crx.prod[i].fiss));
      outVal(11,3,crx.prod[i].fiss/sum*100.0);
      nl();
    }
  }
  std::cout << blank << " TotFission";
  outVal(lowfilter(sum)); nl();
}


/**********************************************************/
/*      Particle Emission Spectra                         */
/**********************************************************/
void outSpectrum(const bool betacalc, const double de, double **spc, Nucleus *n)
{
  int k0 = beohZeroCut(spc);

  outSectionHead("TOTAL PARTICLE EMISSION SPECTRA");
  std::cout << "#      Emin       Emax"
       << "  Spectra / Decay / MeV" << std::endl;
  std::cout << "#     [MeV]      [MeV]";

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    std::cout << std::setw(11) << particle_name[j];
  }
  if(betacalc) std::cout << "   electron   neutrino";
  nl();

  for(int k=0 ; k<=k0 ; k++){
    if(k==0){
      outVal(11,4,0.0);
      outVal(11,4,0.5*de);
    }
    else{
      outVal(11,4,((double)k-0.5)*de);
      outVal(11,4,((double)k+0.5)*de);
    }

    for(int j=0 ; j<MAX_CHANNEL ; j++){
      if(!n->cdt[j].status) continue;
      outVal(lowfilter(spc[j][k]));
    }
    if(betacalc){
      outVal(lowfilter(spc[MAX_CHANNEL  ][k]));
      outVal(lowfilter(spc[MAX_CHANNEL+1][k]));
    }
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Sum of Spectra and Average Energies               */
/**********************************************************/
void outSpectrumSum(const bool betacalc, const double de, double **spc, Nucleus * n)
{
  double et[MAX_CHANNEL+2], ea[MAX_CHANNEL+2], em[MAX_CHANNEL+2];

  int k0 = beohZeroCut(spc);
  double e, e0, e1;

  int cm = MAX_CHANNEL;
  if(betacalc) cm += 2;

  for(int c=0 ; c<cm ; c++){
    em[c] = et[c] = 0.0;
    for(int k=0 ; k<=k0 ; k++){

      e0 = (k > 0)  ? (k - 0.5) * de : 0;
      e1 = (k + 0.5) * de;
      e  = (e0 + e1) * 0.5;

      em[c] += spc[c][k] * de;
      et[c] += spc[c][k] * de * e;
    }
    ea[c] = (em[c]>0.0) ? et[c]/em[c] : 0.0;
  }


  outSectionHead("SUM SPECTRA AND AVERAGE ENERGIES");

  std::cout << cline << blank;
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    std::cout << std::setw(11) << particle_name[j];
  }
  if(betacalc) std::cout << "   Electron   Neutrino";
  nl();

  std::cout << "      TotalEnergy[MeV]";
  double etot = 0.0;
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    outVal(lowfilter(et[j]));
    etot += et[j];
  }
  if(betacalc){
    outVal(lowfilter(et[MAX_CHANNEL  ]));  etot += et[MAX_CHANNEL  ];
    outVal(lowfilter(et[MAX_CHANNEL+1]));  etot += et[MAX_CHANNEL+1];
  }
  nl();

  std::cout << "    AverageEnergy[MeV]";
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    outVal(lowfilter(ea[j]));
  }
  if(betacalc){
    outVal(lowfilter(ea[MAX_CHANNEL  ]));
    outVal(lowfilter(ea[MAX_CHANNEL+1]));
  }
  nl();

  std::cout << "          Multiplicity";
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    if(!n->cdt[j].status) continue;
    outVal(lowfilter(em[j]));
  }
  if(betacalc){
    outVal(lowfilter(em[MAX_CHANNEL  ]));
    outVal(lowfilter(em[MAX_CHANNEL+1]));
  }
  nl();
  nl();

  std::cout << "    EnergyRelease[MeV]";
  outVal(etot); nl();
  nl();
  nl();
}


/**********************************************************/
/*      Neutron Emission Spectra in Lab Frame             */
/**********************************************************/
void outSpectrumLab(const int k0, const double de, double *spl)
{
  outSectionHead("NEUTRON EMISSION SPECTRA IN LAB FRAME");
  std::cout << "#      Emin       Emax"
       << "  Spectra / Decay / MeV" << std::endl;
  std::cout << "#     [MeV]      [MeV]";

  std::cout << std::setw(11) << particle_name[1];
  nl();

  for(int k=0 ; k<=k0 ; k++){
    if(k==0){
      outVal(11,4,0.0);
      outVal(11,4,0.5*de);
    }
    else{
      outVal(11,4,((double)k-0.5)*de);
      outVal(11,4,((double)k+0.5)*de);
    }

    outVal(lowfilter(spl[k]));
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Gamma Cascade                                     */
/**********************************************************/
void outGammaCascade(const double f, Nucleus *n)
{
  outSectionHead("GAMMAS FROM DISCRETE TRANSITION");
  std::cout << cline; outZA(&n->za); std::cout << std::endl;
  std::cout << "#   Ex[MeV]     Jpi    FinalState";
  std::cout << "  Branching Production    Egamma Production   T(1/2)[s]" << std::endl;

  for(int i0=n->ndisc-1 ; i0>=0 ; i0--){
    if(n->lpop[i0] < output_eps && i0 > 0) continue;

    std::cout << std::setw(3) << i0;
    outVal(8,4,n->lev[i0].energy);
    outVal(7,1,n->lev[i0].spin);
    char p0 = (n->lev[i0].parity < 0) ? '-' : '+';
    std::cout << p0 << "   " << blank << blank;
    outVal(lowfilter(n->lpop[i0] * f));
    if((n->lev[i0].halflife > 0.0) && (lowfilter(n->lpop[i0]) > 0.0)){
      std::cout << blank <<  blank;
      outVal(n->lev[i0].halflife);
    }
    nl();

    if(i0 == 0) continue;
    for(int j=0 ; j<n->lev[i0].ngamma ; j++){
      int    i1 = n->lev[i0].fstate[j];
      double eg = n->lev[i0].energy-n->lev[i1].energy;
      double sg = n->lev[i0].branch[j] * n->lpop[i0] * f;

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
/*      All Discrete Gamma Lines                          */
/**********************************************************/
void outDiscreteGamma(const double f, GammaProduction *gp)
{
  outSectionHead("DISCRETE GAMMA PRODUCTION");

  std::cout << "#       Z   A  Energy[MeV]   Production" << std::endl;
  for(int i=0 ; i<gp->getN() ; i++){
    outVal(5,i+1);
    outVal(4,gp->line[i].za.getZ());
    outVal(4,gp->line[i].za.getA());
    outVal(13,7,gp->line[i].energy);
    outVal(13,gp->line[i].production * f);
    nl();
  }
  nl();
  nl();
}


/**********************************************************/
/*      Finer Gamma-Ray Spectra                           */
/**********************************************************/
void outSpectrumFineGamma(double *gc, GammaProduction *gp, const double dec, const double del)
{
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

  std::cout << "# Emin[MeV]  Emax[MeV]   Disc[/MeV]   Cont[/MeV]    Sum[/MeV]" << std::endl;

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

  for(int kc=mc0 ; kc<=mc ; kc++){
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
/*      Total Residual Nucleus Production Cross Section   */
/**********************************************************/
void outTotalResidual(CumulativeResidualProduct *res)
{
  outSectionHead("STABLE AND LONG-LIVED STATE PRODUCTIONS");
  if(res->getNcurrent() == 0) return;

  std::cout << "#       Residual    Ex[MeV]     Jpi     T(1/2)[s] Prod.Ratio Meta" << std::endl;

  int c = 1;
  double sum0 = 0.0;
  double sum1 = 0.0;
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

        if(res->rp[i].metaflag == 0) sum0 += res->rp[i].production;
        sum1 += res->rp[i].production;
      }
    }
  }
  std::cout << cline << blank << blank << "         g.s.sum";
  outVal(lowfilter(sum0)); nl();
  std::cout << cline << blank << blank << "             sum";
  outVal(lowfilter(sum1)); nl();
  nl(); nl();
}


/**********************************************************/
/*      Isomeric Ratios for Long-Lived Nuclides           */
/**********************************************************/
void outIsomericRatio(CumulativeResidualProduct *res)
{
  outSectionHead("ISOMERIC RATIOS");
  if(res->getNcurrent() == 0) return;

  std::cout << "#       Residual    Ex[MeV]     Jpi     T(1/2)[s] m/(g+m)    m/g" << std::endl;

  int c = 1;
  for(unsigned int z=res->zmin ; z<=res->zmax ; z++){
    for(unsigned int a=res->amin ; a<=res->amax ; a++){
      ZAnumber za(z,a);

      /*** look for ground state of the same ZA */
      for(int i0=0 ; i0<res->getNcurrent() ; i0++){
        if(res->rp[i0].za != za) continue;
        if(res->rp[i0].metaflag != 0) continue;

        double pg = res->rp[i0].production;

        /*** sum g.s. and isomers */
        double pt = res->rp[i0].production;
        for(int i1=0 ; i1<res->getNcurrent() ; i1++){
          if((res->rp[i1].za == za) && (res->rp[i1].metaflag > 0)){
            pt += res->rp[i1].production;
          }
        }

        /*** look at the isomer with the same ZA */
        for(int i1=0 ; i1<res->getNcurrent() ; i1++){
          if((res->rp[i1].za != za) || (res->rp[i1].metaflag == 0)) continue;

          double r0 = (pt == 0.0) ? 0.0 : res->rp[i1].production / pt;
          double r1 = (pg == 0.0) ? 0.0 : res->rp[i1].production / pg;

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
  nl(); nl();
}


/**********************************************************/
/*      Print Fission Product Yield Data                  */
/**********************************************************/
void outFissionProductYield(const int nuktotal, Isotope *nuk)
{
  outSectionHead("INDEPENDENT AND CUMULATIVE FISSION YIELDS");
  if(nuktotal == 0) return;

  /*** find min/max Z and A */
  int z0 = 100, z1 = 0, a0 = 1000, a1 = 0;
  for(int k=0 ; k<nuktotal ; k++){
    int z = nuk[k].getZ();
    int a = nuk[k].getA();
    if(z < z0) z0 = z;
    if(z > z1) z1 = z;
    if(a < a0) a0 = a;
    if(a > a1) a1 = a;
  }

  std::cout << "#        Isotope M    I P      Indep.FPY     Cumul.FPY  T(1/2)[s]" << std::endl;

  int c = 1;
  double s0 = 0.0, s1 = 0.0;
  for(int z = z0 ; z <= z1 ; z++){
    for(int a = a0 ; a <= a1 ; a++){

      ZAnumber za(z,a);

      for(int k=0 ; k<nuktotal ; k++){
        if((unsigned int)z != nuk[k].getZ() || (unsigned int)a != nuk[k].getA()) continue;
        if((nuk[k].initial == 0.0) && (nuk[k].yield == 0.0)) continue;

        double t2 = (nuk[k].getLambda() == 0.0) ? 0.0 : log(2.0)/nuk[k].getLambda();

        outVal(5,c++);
        outZA(&za);
        std::cout << std::setw(2) << nuk[k].getM();
        outVal(5,1,nuk[k].getJ());
        if(nuk[k].getP() > 0) std::cout << " + ";
        else std::cout << " - ";
        outVal(14,nuk[k].initial);
        outVal(14,nuk[k].yield);
        outVal(11,t2);
        nl();
        s0 += nuk[k].initial;
        s1 += nuk[k].yield;
      }
    }
  }
  nl();

  std::cout << "#                    Total";
  outVal(14,s0);
  outVal(14,s1);
  nl(); nl();
}


/**********************************************************/
/*      Calculate and Print Mass and Chain Yield          */
/**********************************************************/
void outFissionProductChainYield(const int nuktotal, Isotope *nuk)
{
  outSectionHead("MASS AND CHAIN FISSION YIELDS");
  if(nuktotal == 0) return;

  /*** find min/max A */
  int a0 = 1000, a1 = 0;
  for(int k=0 ; k<nuktotal ; k++){
    int a = nuk[k].getA();
    if(a < a0) a0 = a;
    if(a > a1) a1 = a;
  }

  std::cout <<"#        A    IndepYield     MassYield    ChainYield" << std::endl;

  double s0, s1, s2, c0, c1, c2;
  /*** nuclide which has the longest half-life in a mass chain */
  s0 = s1 = s2 = 0.0;
  for(int a=a0 ; a<=a1 ; a++){

    c0 = c1 = c2 = 0.0;
    for(int k=0 ; k<nuktotal ; k++){
      if(nuk[k].getA() != (unsigned int) a) continue;

      /*** initial yields */
      c0 += nuk[k].initial;

      /*** just add all the yields of the same mass */
      c1 += nuk[k].yield;

      /*** find a stable nuclide in a decay chain */
      if(nuk[k].isstable()) c2 += nuk[k].yield;
    }

    std::cout << std::setw(10) << a;
    outVal(14,c0);
    outVal(14,c1);
    outVal(14,c2); nl();

    s0 += c0;
    s1 += c1;
    s2 += c2;
  }
  nl();

  std::cout << "#    Total";
  outVal(14,s0);
  outVal(14,s1);
  outVal(14,s2);
  nl(); nl();
}


/**********************************************************/
/*      Delayed Neutron Yield                             */
/**********************************************************/
void outDelayedNeutronYield(const int nuktotal, Isotope *nuk)
{
  outSectionHead("DELAYED NEUTRON YIELD");
  if(nuktotal == 0) return;

  /*** find min/max Z and A */
  int z0 = 100, z1 = 0, a0 = 1000, a1 = 0;
  for(int k=0 ; k<nuktotal ; k++){
    int z = nuk[k].getZ();
    int a = nuk[k].getA();
    if(z < z0) z0 = z;
    if(z > z1) z1 = z;
    if(a < a0) a0 = a;
    if(a > a1) a1 = a;
  }

  std::cout << "#     DN-Emitter M  T(1/2)[s]  DN-Branch     Cumul.FPY      DN-Yield" << std::endl;

  int c = 1;
  double dn = 0.0, dgroup6[6], lambda6[6];
  for(int i=0 ; i<6 ; i++) dgroup6[i] = lambda6[i] = 0.0;

  for(int z = z0 ; z <= z1 ; z++){
    for(int a = a0 ; a <= a1 ; a++){

      ZAnumber za(z,a);

      for(int k=0 ; k<nuktotal ; k++){
        if(nuk[k].isstable()) continue;
        if((unsigned int)z != nuk[k].getZ() || (unsigned int)a != nuk[k].getA()) continue;

        double t2 = (nuk[k].getLambda() == 0.0) ? 0.0 : log(2.0)/nuk[k].getLambda();

        /*** search for delayed neutron emitters */
        for(int m=0 ; m<nuk[k].mode.getNDM() ; m++){
          int j = nuk[k].mode.getNext(m);
          /*** beta-decay and delayed neutron emission, Zi+1 = Zj and Ai != Aj  */
          if( (nuk[k].getZ()+1 == nuk[j].getZ()) && (nuk[k].getA() != nuk[j].getA()) ){

            /*** (cumulative yield) x (number of neutrons) x (branching ratio) */
            double x = nuk[k].yield * (nuk[k].getA() - nuk[j].getA()) * nuk[k].mode.getBranch(m);

            outVal(5,c++);
            outZA(&za);
            std::cout << std::setw(2) << nuk[k].getM();
            outVal(11,t2);
            outVal(11,nuk[k].mode.getBranch(m));
            outVal(14,nuk[k].yield);
            outVal(14,x);
            nl();

            dn += x;

            /*** six group constant */
            int gid = 0;
            if(     t2 > 40.0) gid = 0;
            else if(t2 >  8.0) gid = 1;
            else if(t2 >  3.0) gid = 2;
            else if(t2 >  1.0) gid = 3;
            else if(t2 >  0.3) gid = 4;
            else               gid = 5;

            dgroup6[gid] += x;
            lambda6[gid] += x * nuk[k].getLambda();
          }
        }
      }
    }
  }

  /*** decay constant weighted average */
  for(int i=0 ; i<6 ; i++){
    if(dgroup6[i] > 0.0) lambda6[i] = lambda6[i] / dgroup6[i];
    else lambda6[i] = 0.0;
  }

  nl();
  std::cout << "#                          Total Delayed Neutron Yield";
  outVal(14,dn);
  nl();
  nl();

  std::cout << "#         DN-Group DecayConstant   Group-Yield" << std::endl;
  for(int i=0 ; i<6 ; i++){
    std::cout << std::setw(18) << i+1;
    outVal(14,lambda6[i]);
    outVal(14,dgroup6[i]);
    nl();
  }
  nl(); nl();
}


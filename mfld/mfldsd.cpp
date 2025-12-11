/******************************************************************************/
/*  mfldsd.cpp                                                                */
/*        Unperturbed Single Particle State                                   */
/******************************************************************************/

#include <iostream>
#include <iomanip>

#include "mfld.h"
#include "terminate.h"

static void SDPrepareSingleParticle (SystemData *, HFInterface *, int *, SingleParticle *, State **);
static void SDM2Average (const unsigned int, const unsigned int, SingleParticle *);
static void SDResetEnergy (SingleParticle *, State *, State *);
static void SDPHDensity (const int, MCount *, MCount *, MCount *);
static void SDConvolution (const int, SystemData *, MCount *, MCount *, Density *);
static void SDStateDensity (const int, Density *);
static void SDSingleShell (const int, const int, MCount *, double **);
static void SDBothShells (const int, const int, MCount *, MCount *, double **);
static void SDFix (const int, Density *);

const bool FIX_NEGATIVE = true;
const bool M2AVERAGE = false;
#undef DEBUG_MSTATE

/**********************************************************/
/*      Combinatorial Microscopic Level Density           */
/*      without Interaction                               */
/**********************************************************/
void SDUnperturbed(SystemData *sys, HFInterface *hfsp, Density *sd)
{
  const int pmax = sys->phmax;

  /*** prepare single particle states */
  SingleParticle sp[2];
  int            nl[4];
  State         *st[4];

  SDPrepareSingleParticle(sys,hfsp,nl,sp,st);

  if(M2AVERAGE){
    /*** calculate M2 average */
    SDM2Average(sys->target.getZ(),sys->target.getA(),sp);
    return;
  }

  /*** arrays to count the number of M-states */
  /*** mcp,h will have the number of 1p, 2p, 3p, ... configurations in each Z or N shell.
       mcp,h arrays will be re-used for different particle numbers. */
  MCount mcp, mch;
  mch.memalloc(sys->msize*pmax,sys->nsize);
  mcp.memalloc(sys->msize*pmax,sys->nsize);

  /*** mcz and mcn are the number of Np-Nh configurations for given N, upto pmax */
  MCount *mcz, *mcn;
  mcz = new MCount [pmax];
  mcn = new MCount [pmax];
  for(int p=0 ; p<pmax ; p++){
    mcz[p].memalloc(sys->msize*2*(p+1),sys->nsize);  mcz[p].clear();
    mcn[p].memalloc(sys->msize*2*(p+1),sys->nsize);  mcn[p].clear();
  }

  /*** for each particle and hole state, count the total number of configurations */
  for(int p=0 ; p<pmax ; p++){

    /*** neutron shell */
    mch.clear();
    mcp.clear();
    mcn[p].clear();

    SDCombination(p+1,sys,nl[0],st[0],&mch);       // number of (p+1)-hole states
    SDCombination(p+1,sys,nl[1],st[1],&mcp);       // number of (p+1)-particle states
    SDPHDensity(sys->nsize,&mch,&mcp,&mcn[p]);     // calculate particle-hole numbers

#ifdef DEBUG_MSTATE
    printMState("Nh",(p+1)%2,p,&mch);
    printMState("Np",(p+1)%2,p,&mcp);
#endif

    /*** proton shell */
    mch.clear();
    mcp.clear();
    mcz[p].clear();

    SDCombination(p+1,sys,nl[2],st[2],&mch);
    SDCombination(p+1,sys,nl[3],st[3],&mcp);
    SDPHDensity(sys->nsize,&mch,&mcp,&mcz[p]);

#ifdef DEBUG_MSTATE
    printMState("Zh",(p+1)%2,p,&mch);
    printMState("Zp",(p+1)%2,p,&mcp);
#endif

    message << "counting particle-hole states for p-h:" << p;
    cohNotice("SDUnperturbed");
  }

#ifdef DEBUG_MSTATE
  for(int p=0 ; p<pmax ; p++){
    printMState("N ",0,p,&mcn[p]);
    printMState("Z ",0,p,&mcz[p]);
  }
#endif

  /*** convolute Z and N shells */
  /*** the actual p-h configurations include different p-h states in Z or N shell.
       for example, 2p-2h can be 2p-2h in the same shell, or 1p-1h in both Z and N */
  SDConvolution(pmax,sys,mcz,mcn,sd);

  /*** convert M-scheme into J-scheme by rho(J) = rho(J=M) - rho(J=M+1) */
  SDStateDensity(pmax,sd);

  /*** remove wiggle */
  if(FIX_NEGATIVE) SDFix(pmax,sd);

  /*** free allocated memory */
  sp[proton ].memfree();
  sp[neutron].memfree();

  mch.memfree();
  mcp.memfree();
  for(int p=0 ; p<pmax ; p++){
    mcz[p].memfree();
    mcn[p].memfree();
  }
  delete [] mcz;
  delete [] mcn;

  for(int k=0 ; k<4 ; k++) delete [] st[k];
}


/**********************************************************/
/*      Prepare Single Particle State Arrays              */
/**********************************************************/
void SDPrepareSingleParticle(SystemData *sys, HFInterface *hfsp, int *nl, SingleParticle *sp, State **st)
{
  sp[neutron].memalloc(MAX_SPSTATE);
  sp[proton ].memalloc(MAX_SPSTATE);

  /*** prepare particle and hole states in each Z and N shell */
  for(int n=0 ; n<2 ; n++){
    SDStoreState((Particle)n,&hfsp[n],&sp[n]);
    PHCount(sys->sharpcutoff,&sp[n]);
  }

  /*** number of particle or hole levels in Z and N shells */
  /*** sp.h0, h1, p0, and p1 are the level index of lowest/highest hole/particle state.
       nl([nz]x[hp]) are the number of levels for N/Z shell Hole/Particle state.  */
  int k = 0;
  nl[k] = sp[neutron].h1 - sp[neutron].h0 + 1; st[k] = new State [nl[k]]; k++;
  nl[k] = sp[neutron].p1 - sp[neutron].p0 + 1; st[k] = new State [nl[k]]; k++;
  nl[k] = sp[proton ].h1 - sp[proton ].h0 + 1; st[k] = new State [nl[k]]; k++;
  nl[k] = sp[proton ].p1 - sp[proton ].p0 + 1; st[k] = new State [nl[k]]; k++;

  SDResetEnergy(&sp[neutron],st[0],st[1]);
  SDResetEnergy(&sp[proton ],st[2],st[3]);
}


/**********************************************************/
/*      Store Single Particle States and V2               */
/**********************************************************/
void SDStoreState(const Particle p, HFInterface *hfsp, SingleParticle *sp)
{
  const bool boundonly = false;

  const double ecut = 10.0;
  const double vcut = 1.0e-10;

  sp->lambda = hfsp->lambda;

  for(int k=0 ; k<hfsp->n ; k++){

    int    j2 = hfsp->state[k].k2;
    int    pt = (hfsp->state[k].parity > 0.0) ? 1 : -1;
    double e  = hfsp->state[k].energy;
    double v  = hfsp->state[k].v2;

    int t = 0;
    if(p == neutron) t = (v > 0.5) ? -1 :  1;
    else             t = (v > 0.5) ?  1 : -1;

    if(boundonly){
      if(e < 0.0){
        sp->set(-j2,pt,t,e,v);
        sp->set( j2,pt,t,e,v);
      }
      if(e > ecut) break;
      if(v < vcut) break;
    }
    else{
      sp->set(-j2,pt,t,e,v);
      sp->set( j2,pt,t,e,v);
    }
/*
    int c = sp->getNmax();
    std::cout << std::setw(4) << c-2;
    std::cout << std::setw(6) << sp->getM(c-2);
    std::cout << std::setw(3) << ((pt == 1) ? " + " : " - ");
    std::cout << std::setw(6) << sp->getT(c-2);
    std::cout << std::setw(12) << sp->getE(c-2);
    std::cout << std::setw(12) << sp->getV2(c-2) << std::endl;

    std::cout << std::setw(4) << c-1;
    std::cout << std::setw(6) << sp->getM(c-1);
    std::cout << std::setw(3) << ((pt == 1) ? " + " : " - ");
    std::cout << std::setw(6) << sp->getT(c-1);
    std::cout << std::setw(12) << sp->getE(c-1);
    std::cout << std::setw(12) << sp->getV2(c-1) << std::endl;
*/
  }
}


/**********************************************************/
/*      Store Single Particle States and V2               */
/**********************************************************/
void SDM2Average(const unsigned int z, const unsigned int a, SingleParticle *sp)
{
  double sn[2];
  double sm[2];

  for(int n=0 ; n<2 ; n++){

    sn[n] = 0.0;
    sm[n] = 0.0;

    double m2 = 0.0, xm2 = 0.0, xn = 0.0;

    for(int i=0 ; i<sp[n].getNmax() ; i++){

      if( fabs(sp[n].getE(i) - sp[n].lambda) < 5.0){

        m2   = sp[n].getM(i) * sp[n].getM(i) / 4.0;
        xm2 += sp[n].getV2(i) * m2;
        xn  += sp[n].getV2(i);

        // std::cout << std::setw(12) << sp[n].lambda;
        // std::cout << std::setw(12) << sp[n].getE(i);
        // std::cout << std::setw(12) << sp[n].getV2(i);
        // std::cout << std::setw(4) << sp[n].getM(i);
        // std::cout << std::setw(12) << xn;
        // std::cout << std::setw(12) << m2 << std::endl;
      }
    }
    sn[n] += xn;
    sm[n] += xm2;
  }

  std::cout << std::setw(5) << z;
  std::cout << std::setw(5) << a;
  std::cout << std::setw(12) << sm[0] / sn[0];
  std::cout << std::setw(12) << sm[1] / sn[1];
  std::cout << std::setw(12) << (sm[0] + sm[1]) / (sn[0] + sn[1]) << std::endl;
}


/**********************************************************/
/*      Shift Energies Relative to E-Fermi                */
/**********************************************************/
void SDResetEnergy(SingleParticle *p, State *sth, State *stp)
{
  /*** For hole states, Eh = E_Fermi - Es.p., and
       for particle states, Ep = E.s.p - E_Fermi
       to make the energies (almost) positive.
       However, the results are indepent of E_Fermi. */
  int c = 0;
  for(int i=p->h1 ; i>=p->h0 ; i--){
    double e = p->lambda - p->getE(i);
    sth[c++].set(p->getM(i),p->getP(i),p->getT(i),e,p->getV2(i));
  }

  c = 0;
  for(int i=p->p0 ; i<= p->p1 ; i++){
    double e = p->getE(i) - p->lambda;
    stp[c++].set(p->getM(i),p->getP(i),p->getT(i),e,1.0-p->getV2(i)); 
  }
}


/**********************************************************/
/*      Count Number of M-States for Np-Nh Configuration  */
/**********************************************************/
void SDCombination(const int n, SystemData *sys, const int nlev, State *lev, MCount *mc)
{
  MConfig cf;
  cf.memalloc(n);
  cf.init(sys->emax, sys->de);

  /*** many-fold loop by recursive call */
  int depth = 1;
  SDRecursive(depth,n+1,nlev,lev,&cf,mc);

  cf.memfree();
}


/**********************************************************/
/*      Recursive Call                                    */
/**********************************************************/
void SDRecursive(const int depth, const int n, const int nlev, State *lev, MConfig *cf, MCount *mc)
{
  /*** loop terminate when the depth is the same as N */
  if(depth == n){
    /*** energy bin index for the total N-particle (N-hole) states */
    int ke = (int)(cf->energy[depth-1]/cf->de);

    /*** for quasi-particles, the occupation probability is rounded to 0 - 9 */
    unsigned long g = round(cf->fraction[depth-1] * 10.0);

    int m = cf->mvalue[depth-1];

    if( (g > 0) && (ke >= 0) && (m >= 0) ){
      /*** because only positive M-states are recorded in MCount,
           M=0 should be doubled for the sum(M) = 0 case */
      unsigned long c = (m == 0) ? 2L : 1L;
      if(cf->parity[n-1] > 0) mc->m0[ke][m] += g * c;
      else                    mc->m1[ke][m] += g * c;
    }
  }
  else{
    int id0 = cf->id[depth-1] + 1;
    int id1 = nlev;
    int id2 = cf->step[depth];
    for(cf->id[depth] = id0 ; cf->id[depth] < id1 ; cf->id[depth] += id2){

      /*** cumulative excitaion energy, fraction, M, and parity */
      cf->energy[depth]   = cf->energy[depth-1]   + lev[cf->id[depth]].getE();
      if(cf->energy[depth] >= cf->emax) break;

      cf->fraction[depth] = cf->fraction[depth-1] * lev[cf->id[depth]].getV2();
      cf->mvalue[depth]   = cf->mvalue[depth-1]   + lev[cf->id[depth]].getM();
      cf->parity[depth]   = cf->parity[depth-1]   * lev[cf->id[depth]].getP();

      /*** go into one deeper level */
      SDRecursive(depth+1,n,nlev,lev,cf,mc);
    }
  }
}


/**********************************************************/
/*      Nparticle - Nhole Configuration State Density     */
/**********************************************************/
void SDPHDensity(const int nsize, MCount *mch, MCount *mcp, MCount *mc)
{
  int mmax = mc->getMsize();
  int mth  = mch->getMsize();
  int mtp  = mcp->getMsize();

  /*** instead of counting the number of Np-Nh states exactly,
       we bin the hole and particle states separately,
       then the number of Np-Nh is given by their convolution. */
  for(int nh=0 ; nh<nsize ; nh++){
    for(int np=0 ; np<nsize - nh ; np++){

      int ke = np + nh;

      for(int mh=-mth ; mh<=mth ; mh++){
        unsigned long gh0 = mch->m0[nh][abs(mh)];
        unsigned long gh1 = mch->m1[nh][abs(mh)];
        if( (gh0 == 0L) && (gh1 == 0L) ) continue;

        for(int mp=-mtp ; mp<=mtp ; mp++){
          int ms = mh + mp;
          if((ms < 0) || (ms > mmax)) continue;

          unsigned long gp0 = mcp->m0[np][abs(mp)];
          unsigned long gp1 = mcp->m1[np][abs(mp)];
          if( (gp0 == 0L) && (gp1 == 0L) ) continue;

          unsigned long g0 = gh0 * gp0 + gh1 * gp1;
          unsigned long g1 = gh0 * gp1 + gh1 * gp0;
          mc->m0[ke][ms] += g0;
          mc->m1[ke][ms] += g1;
        }
      }
    }
  }
}


/**********************************************************/
/*      Convolution of States to Calculate Level Density  */
/**********************************************************/
void SDConvolution(const int pmax, SystemData *sys, MCount *mcz, MCount *mcn, Density *sd)
{
  double *fbuf[2], *f[2];
  int fmax  = 2 * pmax * sys->msize;
  int fsize = 2 * pmax * (2*sys->msize+1);

  fbuf[even] = new double [fsize];  f[even] = &fbuf[even][fmax];
  fbuf[odd ] = new double [fsize];  f[odd ] = &fbuf[odd ][fmax];

  sd->memclear();

  /*** loop for energy bins */
  for(int n=0 ; n<sys->nsize ; n++){

    /*** loop for Np-Nh */
    for(int p=0 ; p<pmax ; p++){
      int pt = p+1; // actual total particle number = index + 1

      /*** f is the total number of intrinsic states at the n-th excitation energy */
      for(int m=-fmax ; m<=fmax ; m++) f[even][m] = f[odd][m] = 0.0;

      /*** loop for number of particles in Z-shell */
      for(int pz = 0 ; pz<=pt ; pz++){
        int pn = pt - pz; // particles in N-shell = total - Zp
        
        /*** N shell only */
        if(pz == 0){
          SDSingleShell(n,fmax,&mcn[pn-1],f);
        }
        /*** Z shell only */
        else if(pn == 0){
          SDSingleShell(n,fmax,&mcz[pz-1],f);
        }
        /*** convolution of N and Z configurations */
        else{
          SDBothShells(n,fmax,&mcz[pz-1],&mcn[pn-1],f);
        }
      }

/*
      std::cout << std::setw(3) << n << std::setw(3) << p << std::setw(3) << std::endl;
      for(int m=-10 ; m<=10 ; m+=2) std::cout << std::setw(10) << f[even][m];
      std::cout << std::endl;
      for(int m=-10 ; m<=10 ; m+=2) std::cout << std::setw(10) << f[odd][m];
      std::cout << std::endl;
*/
      int jmax = sys->jsize-1;
      if(jmax > 2*fmax) jmax = 2*fmax;

      /*** temporary ucopy number of states in n-th bin in sd */
      for(int j=0 ; j<=jmax ; j++){
        sd->r0[p][n][j] = f[even][j*2];
        sd->r1[p][n][j] = f[odd ][j*2];
      }
    }

    if(n%10 == 0){
      message << "comvoluting p-h states in the energy bin " << n << "/" << sys->nsize;
      cohNotice("SDConvolution");
    }
  }

  delete [] fbuf[even];
  delete [] fbuf[odd ];
}


/**********************************************************/
/*      Convert J-Scheme State Density from M-Scheme      */
/**********************************************************/
void SDStateDensity(const int pmax, Density *sd)
{
  for(int n=0 ; n<sd->getNsize() ; n++){

    for(int p=0 ; p<pmax ; p++){
      /*** state density = {rho(M=J) - rho(M=J+1)} / dE */
      for(int j=0 ; j<sd->getJsize()-2 ; j++){
        sd->r0[p][n][j] = (sd->r0[p][n][j] - sd->r0[p][n][j+1]) / sd->binwidth();
        sd->r1[p][n][j] = (sd->r1[p][n][j] - sd->r1[p][n][j+1]) / sd->binwidth();
      }
    }
  }
}


/**********************************************************/
/*      Np-Nh in Same Shell                               */
/**********************************************************/
void SDSingleShell(const int n, const int fmax, MCount *mc, double **f)
{
  int mt = mc->getMsize();
  for(int m=-mt ; m<=mt ; m++){
    if(abs(m) <= fmax){
      f[even][m] += 0.01 * mc->m0[n][abs(m)];
      f[odd ][m] += 0.01 * mc->m1[n][abs(m)];
    }
  }
}


/**********************************************************/
/*      (N-x)p-(N-x)h and xp-xh Configurations            */
/**********************************************************/
void SDBothShells(const int n, const int fmax, MCount *mcz, MCount *mcn, double **f)
{
  int mzt = mcz->getMsize();
  int mnt = mcn->getMsize();

  for(int kz=0 ; kz<=n ; kz++){
    /*** energy bins whose sum give the N-th bin eneregy */
    int kn = n - kz;

    for(int mn=-mnt ; mn<=mnt ; mn++){
      double fn0 = 0.01 * mcn->m0[kn][abs(mn)];
      double fn1 = 0.01 * mcn->m1[kn][abs(mn)];

      if((fn0 == 0.0) && (fn1 == 0.0)) continue;

      for(int mz=-mzt ; mz<=mzt ; mz++){

        double fz0 = 0.01 * mcz->m0[kz][abs(mz)];
        double fz1 = 0.01 * mcz->m1[kz][abs(mz)];

        if((fz0 == 0.0) && (fz1 == 0.0)) continue;

        int m = mn + mz;

        if(abs(m) <= fmax){
          f[even][m] += fz0 * fn0 + fz1 * fn1;
          f[odd ][m] += fz0 * fn1 + fz1 * fn0;
        }
      }
    }
  }
}


/**********************************************************/
/*      Fix Numerical Noise in Calculated Level Density   */
/**********************************************************/
void SDFix(const int pmax, Density *sd)
{
  const double eps = 1.0e-10;

  for(int p=0 ; p<pmax ; p++){
    for(int n=0 ; n<sd->getNsize() ; n++){
      for(int j=0 ; j<sd->getJsize() ; j++){

        if(abs(sd->r0[p][n][j]) < eps) sd->r0[p][n][j] = 0.0;
        if(abs(sd->r1[p][n][j]) < eps) sd->r1[p][n][j] = 0.0;

        if(sd->r0[p][n][j] < 0.0){
          double x1 = 0.0, x2 = 0.0;
          if(n > 0)                x1 = sd->r0[p][n-1][j] + sd->r0[p][n][j];
          if(n < sd->getNsize()-2) x2 = sd->r0[p][n+1][j] + sd->r0[p][n][j];

          sd->r0[p][n][j] = 0.0;
          if((n > 0) && (x1 >= 0.0))                     sd->r0[p][n-1][j] = x1;
          else if((n < sd->getNsize()-2) && (x2 >= 0.0)) sd->r0[p][n+1][j] = x2;
        }

        if(sd->r1[p][n][j] < 0.0){
          double x1 = 0.0, x2 = 0.0;
          if(n > 0)                x1 = sd->r1[p][n-1][j] + sd->r1[p][n][j];
          if(n < sd->getNsize()-2) x2 = sd->r1[p][n+1][j] + sd->r1[p][n][j];

          sd->r1[p][n][j] = 0.0;
          if((n > 0) && (x1 >= 0.0))                     sd->r1[p][n-1][j] = x1;
          else if((n < sd->getNsize()-2) && (x2 >= 0.0)) sd->r1[p][n+1][j] = x2;
        }
      }
    }
  }
}


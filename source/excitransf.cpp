/******************************************************************************/
/*  excitransf.cpp                                                            */
/*        pick-up and  stripping, in terms of exciton model                   */
/******************************************************************************/

#include <iostream>
#include <cmath>

#include "structur.h"
#include "nucleus.h"
#include "statmodel.h"
#include "exciton.h"

static double  preqStateDensityNt (const double, const double, SPdens *, Exconf *);


/**********************************************************/
/*      Transfer Reaction in Pre-Equilibrium Process      */
/**********************************************************/
void preqExcitonTransfer(System *sys, Pdata *pdt, Transmission **tc, Transmission **td, Spectra *spc)
{
  Preeq    prq;
  Exconf   exc;
  SPdens   spd[MAX_CHANNEL];

/************************************************/
/*      Parameter Setting                       */
/************************************************/

  /*** Total excitation energy */
  prq.ex_total   = sys->ex_total;

  /*** Parameters for state density */
  for(int id=0 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;
    int i = ncl[0].cdt[id].next;

    /*** Default single particle state density */
    spd[id].gz         = preqSingleStateDensity(id, ncl[i].za.getZ());
    spd[id].gn         = preqSingleStateDensity(id, ncl[i].za.getN());

    /*** Pairing energy */
    spd[id].pairing    = preqPairingDelta(ncl[i].za.getZ(),ncl[i].za.getA());

    /*** Potential well depth */
    spd[id].well_depth = preqPotentialDepth(sys->lab_energy,
                                            sys->incident.za.getZ(),ncl[i].za.getA());
  }
  prq.spd = spd;

  /*** Clear spc array */
  spc->memclear("pe");


/************************************************/
/*      Phenomenological Inputs                 */
/*      from PRECO-2006 and TALYS manual        */
/************************************************/

  double va  = 12.5 * pdt[sys->inc_id].za.getA();
  double c   = (sys->incident.pid == neutron) ? 5500.0 : 3800.0;
  double v1;
  if(     sys->incident.pid == neutron)        v1 =  7.0;
  else if(sys->incident.pid == proton   ||
          sys->incident.pid == deuteron ||
          sys->incident.pid == triton)         v1 = 17.0;
  else                                         v1 = 25.0;
  double vw  = v1+ (38.0-v1) * pow(sys->cms_energy,4.0)/
                              (pow(sys->cms_energy,4.0)+ pow(150.0,4.0));
  double xnt = 7.0*sqrt(sys->cms_energy/(double)sys->incident.za.getA())
              /(vw * (double)sys->target.getA() * (double)sys->target.getA());


/************************************************/
/*      For Each Ejectile                       */
/************************************************/

  double c1,c2,c3,d1,d2,d3,d4;

  c1 = (pdt[sys->inc_id].spin2+1.0)*(double)(sys->incident.za.getA()*sys->incident.za.getA());
  c2 = sys->incident.za.getA()/(sys->cms_energy + va);
  c3 = 2.0*(double)sys->target.getZ()/(double)sys->target.getA();

  for(int id=1 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;
    if(id == sys->inc_id) continue;
    if(sys->incident.za.getA()==1 && pdt[id].za.getA()==1) continue;

    int j = ncl[0].cdt[id].next;

    /*** Exciton numbers */
    int dz = sys->incident.za.getZ() - pdt[id].za.getZ();
    int dn = sys->incident.za.getN() - pdt[id].za.getN();
    int  n = sys->incident.za.getA() - pdt[id].za.getA();

    if(dz > 0){ exc.zp = dz;  exc.zh =   0; } // Za>Zb, proton-particle
    else      { exc.zp =  0;  exc.zh = -dz; } // Za<Zb, proton-hole

    if(dn > 0){ exc.np = dn;  exc.nh =   0; } // Na>Nb, neutron-particle
    else      { exc.np =  0;  exc.nh = -dn; } // Na<Nb, neutron-hole

    /*** Pick up, stripping, exchange */
    if(n < 0){     d1 = 1.0/(  80.0*     sys->cms_energy) ; n *= -1; }
    else if(n > 0) d1 = 1.0/( 580.0*sqrt(sys->cms_energy));
    else           d1 = 1.0/(1160.0*sqrt(sys->cms_energy));


    d2 = 1.0;
    if( (sys->incident.pid == alpha) && (pdt[id].za.getA() == 1) ) d2 = 12.0;
    else if( (sys->incident.za.getA() == 1) && (pdt[id].pid == alpha) ){
      if(sys->cms_energy < 20.0) d2 = 12.0;
      else d2 = 12.0 - 11.0*(sys->cms_energy-20.0)/sys->cms_energy;
    }
    d3 = (pdt[id].spin2+1.0)*pdt[id].za.getA();
    d4 = d1*d2*d3/c1
       * pow(c2,2.0*n)
       * pow(c/(double)ncl[j].za.getA(),(double)n)
       * pow(c3,2.0*(sys->incident.za.getZ()+2)*exc.zh+2.0*exc.np);
/*
    printf("%4d%4d%4d%4d  ",exc.zp,exc.zh,exc.np,exc.nh);
    printf(" % 11.4e % 11.4e % 11.4e\n",d2,d3,d4);
 */
    /*** Calculation for continuum region only */
    for(int k=1 ; k<ncl[j].ncont ; k++){
      double omega = preqStateDensityNt(xnt,ncl[j].excitation[k],&prq.spd[id],&exc);
      spc->pe[id][k] = d4 * tc[id][k].ecms * tc[id][k].sigr * omega;
/*
      printf("%3d %7.3f %10.3e %11.4e %11.4e\n",
              k,tc[id][k].ecms,tc[id][k].sigr,omega,spc->pe[id][k]);
*/
    }

    preqExcitonPopulation(id,tc[id],td[id],spc->pe[id]);
  }
}


/**********************************************************/
/*      State Density for Nucleon Transfer Reaction       */
/**********************************************************/
double preqStateDensityNt(const double c, const double ex, SPdens *spd, Exconf *exc)
{
  Exconf x;

  double w1 = 0.0;
  for(int i=0 ; i<=3 ; i++){
    for(int j=0 ; j<=(3-i) ; j++){
      double xnt = c*(  exc->np*exc->np +     exc->zp*exc->zp
                      + exc->nh*exc->nh + 1.5*exc->zh*exc->zh );

      x.set(exc->zp+i, exc->zh+i, exc->np+j, exc->nh+j);

      w1 += pow(xnt,(double)(i+j)) * preqStateDensity(ex,spd,&x);
    }
  }

  double w2 = 0.0;
  for(int i=0 ; i<=exc->zp ; i++){
    for(int j=0 ; j<=exc->zh ; j++){
      for(int k=0 ; k<=exc->np ; k++){
        for(int l=0 ; l<=exc->nh ; l++){
          if( (double)(i+j+k+l) > 0.5 ){
            x.set(exc->zp-i, exc->zh-j, exc->np-k, exc->nh-l);
            w2 += preqStateDensity(ex,spd,&x);
          }
        }
      }
    }
  }

  return(w1+w2);
}

/******************************************************************************/
/*  FRDMfileout.cpp                                                           */
/*        Print FRDM results on a file                                        */
/******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>

#include "FRDM.h"
#include "global.h"


/**********************************************************/
/*      Print FRDM Nuclear Shape                          */
/**********************************************************/
void FRDMWrite3DMesh(std::string filename, PolarCoordinate *p)
{
  bool polygon = false;
  std::ofstream fp;
  fp.open(&filename[0]);

  fp.setf(std::ios::scientific, std::ios::floatfield);
  fp << std::setprecision(4);

  if(polygon){
    for(int i = 0 ; i<p->nt-1 ; i++){
      for(int j = 0 ; j<p->np-1 ; j++){
        double x11 = p->getR(i,j) * p->getX(i,j);
        double y11 = p->getR(i,j) * p->getY(i,j);
        double z11 = p->getR(i,j) * p->getZ(i);

        double x12 = p->getR(i+1,j) * p->getX(i+1,j);
        double y12 = p->getR(i+1,j) * p->getY(i+1,j);
        double z12 = p->getR(i+1,j) * p->getZ(i+1);

        double x21 = p->getR(i,j+1) * p->getX(i,j+1);
        double y21 = p->getR(i,j+1) * p->getY(i,j+1);
        double z21 = p->getR(i,j+1) * p->getZ(i);

        double x22 = p->getR(i+1,j+1) * p->getX(i+1,j+1);
        double y22 = p->getR(i+1,j+1) * p->getY(i+1,j+1);
        double z22 = p->getR(i+1,j+1) * p->getZ(i+1);

        fp << std::setw(12) << x11 << std::setw(12) << y11 << std::setw(12) << z11;
        fp << std::setw(12) << x12 << std::setw(12) << y12 << std::setw(12) << z12;
        fp << std::setw(12) << x22 << std::setw(12) << y22 << std::setw(12) << z22;
        fp << std::setw(12) << x21 << std::setw(12) << y21 << std::setw(12) << z21 << std::endl;
      }
    }
  }
  else{
    for(int i = 0 ; i<p->nt ; i++){
      for(int j = 0 ; j<p->np ; j++){

        double x = p->getR(i,j) * p->getX(i,j);
        double y = p->getR(i,j) * p->getY(i,j);
        double z = p->getR(i,j) * p->getZ(i);

        fp << std::setw(12) << x;
        fp << std::setw(12) << y;
        fp << std::setw(12) << z << std::endl;
      }
      fp << std::endl;
    }
  }

  fp.close();
}


/**********************************************************/
/*      Print FRDM Microscopic Model Parameters           */
/**********************************************************/
void FRDMWriteParameters(const int z, const int a, LNParameter *bcs, SCParameter *scp)
{
  std::string fname = fileFRDMParameter;

  std::ofstream fp;
  fp.open(fname.c_str(),std::ios::out | std::ios::app);
  if(!fp) return;

  fp.setf(std::ios::scientific, std::ios::floatfield);
  fp << std::setprecision(4);

  fp << std::setw(4) << z << std::setw(4) << a;
  for(int p=0 ; p<=1 ; p++){
    fp << std::setw(12) << bcs[p].pairing_gap;
    fp << std::setw(12) << bcs[p].lambda2;
    fp << std::setw(12) << bcs[p].chemical_potential;
    fp << std::setw(12) << scp[p].chemical_potential;
    fp << std::setw(12) << scp[p].rho * 2.0;
  }

  double ep[2],es[2];
  for(int p=0 ; p<=1 ; p++){
    ep[p] = bcs[p].energy_micro - bcs[p].energy_average;
    es[p] = scp[p].energy_micro - scp[p].energy_average;
  }

  fp << std::setw(12) << ep[0] << std::setw(12) << ep[1] << std::setw(12) << ep[0]+ep[1];
  fp << std::setw(12) << es[0] << std::setw(12) << es[1] << std::setw(12) << es[0]+es[1];

  fp << std::endl;
}


/**********************************************************/
/*      Print FRDM Energies                               */
/**********************************************************/
void FRDMWriteEnergy(std::string fname, MFTSystem *sys, FRDMEnergy *e)
{
  std::ofstream fp;
  fp.open(fname.c_str(),std::ios::out | std::ios::app);
  if(!fp) return;

  fp.setf(std::ios::scientific, std::ios::floatfield);
  fp << std::setprecision(4);

  fp << std::setw(4) << sys->getZ()
     << std::setw(4) << sys->getA();

  fp << std::setw(12) << sys->getEps(2)
     << std::setw(12) << sys->getEps(3)
     << std::setw(12) << sys->getEps(4);

  fp << std::setw(12) << e->spherical
     << std::setw(12) << e->macro
     << std::setw(12) << e->micro
     << std::setw(12) << e->total
     << std::endl;
}

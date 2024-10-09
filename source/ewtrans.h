// transfer data between coupled-channels and Hauser-Feshbach

#ifdef CC_TOPLEVEL
Collective       gCol;
CCdata          *gCdt;
double          *gHermiteEigenvalue, *gChannelSigma, *gSPhase;
std::complex<double> *gHermiteEigenvector;
#endif

#ifndef CC_TOPLEVEL
extern Collective       gCol;
extern CCdata          *gCdt;
extern double          *gHermiteEigenvalue, *gChannelSigma, *gSPhase;
extern std::complex<double> *gHermiteEigenvector;
#endif

double NuSystematics_Wahl     (const int, double);

void   WahlZpModel            (const int, const int, const int, const int, const int, const double, double **, int *, double *);
void   WahlChainYield         (const int, const int, const double, double *, double *, double *);


class ZpModelParameter{
 public:
  /* parameters */
  double SigmaZ;
  double DeltaZ;
  double FZ;
  double FN;

  ZpModelParameter(){
    SigmaZ    = 0.0;
    DeltaZ    = 0.0;
    FZ        = 0.0;
    FN        = 0.0;
  }
};


class ZpModelData{
 public:
  /* A-boundary */
  double B[6];
  double Ba;
  double Bb;

  /* peak reagion*/
  double SigmaZ140;     // sigma_Z(140)
  double DeltaZ140;     // Delta_Z(140)
  double FZ140;         // F_Z(140)
  double FN140;         // F_N(140)
  double dSigmaZ;       //  dsigma_Z / dA
  double dDeltaZ;       //  dDelta_Z / dA
  double dFZ;           //  dF_Z / dA

  /* symmetry region */
  double d50;           //  d
  double SigmaZ50;      //  sigma_Z50
  double DeltaZmax;     //  Delta_Zmax

  /* wing region */
  double dSigmaZW;      //  dsigma_Z / dA in Wing
  double dDeltaZW;      //  dDelta_Z / dA in Wing
  double dFZW;          //  dF_Z / dA in Wing
  double dFNW;          //  dF_N / dA in Wing

  ZpModelData(){
    for(int i=0 ; i<6 ; i++) B[i] = 0.0;
    Ba = 0.0;
    Bb = 0.0;

    SigmaZ140 = 0.0;
    DeltaZ140 = 0.0;
    FZ140     = 0.0;
    FN140     = 0.0;
    dSigmaZ   = 0.0;
    dDeltaZ   = 0.0;
    dFZ       = 0.0;

    d50       = 0.0;
    SigmaZ50  = 0.0;
    DeltaZmax = 0.0;

    dSigmaZW  = 0.0;
    dDeltaZW  = 0.0;
    dFZW      = 0.0;
    dFNW      = 0.0;
  }
};

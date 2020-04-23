#include "createControlSolid.H"

int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);
nCorr=nCorr;

scalar convergenceTolerance(readScalar(stressControl.lookup("D")));
convergenceTolerance=convergenceTolerance;

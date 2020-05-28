//Solving Solid Continuity Equation

#include "updateRxnCoeff.H"

Info<< "\nCalculating epss" << endl;
fvScalarMatrix epssEqn
(       
        fvm::ddt(epss)
      + fvm::div(phis,epss)
      -  Krxn //units of 1/(seconds * control Vol)
);

epssEqn.solve();
epss.correctBoundaryConditions();
epss.max(0.0);
epss.min(0.999);

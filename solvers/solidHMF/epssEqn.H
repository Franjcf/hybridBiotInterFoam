//Solving Solid Continuity Equation

Info<< "\nCalculating epss" << endl;
fvScalarMatrix epssEqn
(
	fvm::ddt(epss)
      + fvm::div(phis,epss)
);

epssEqn.solve();
epss.correctBoundaryConditions();

epss.min(0.999);
epss.max(0.0);

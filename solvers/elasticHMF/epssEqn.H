//Update Solid Speed from Displacement Field
Us = fvc::ddt(D);
Us.correctBoundaryConditions();
phis = fvc::interpolate(Us) & mesh.Sf();  
Dsup = D*epss;

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

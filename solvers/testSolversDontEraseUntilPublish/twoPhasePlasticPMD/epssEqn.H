//Solving Solid Continuity Equation

Info<< "\nCalculating epss" << endl;

fvScalarMatrix epssEqn
(
	fvm::ddt(epss)
      + fvm::div(phip,epss)
	//add reaction here
);

epssEqn.solve();
epss.correctBoundaryConditions();

//epss.min(0.9999);

//updating the solid indicator function
forAll(epss,celli) 
{
	if(epss[celli]>0.01) //0.001 (good) //0.0001 (also Good)
	{		
	B[celli]=1; //in clay (use porous Pc models)
	}
	else
	{
	B[celli]=0; //in free pore (use interface curvature dependent Pc model)
	}
}

B.correctBoundaryConditions();

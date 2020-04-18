//Solving Solid Continuity Equation

//Defining reaction to just be at the interface (for biofilm experiments)

if(krxn.value()>0) //if biofilm is growing
{

forAll(epss,celli)
{
	if(epss[celli] >= epssEq.value()) 
        {
        Solid[celli]=1; 
        }
        else
        {
        Solid[celli]=0; 
        }
}

Solid.correctBoundaryConditions();

SolidBorder = fvc::grad(Solid);
//Bborder = fvc::grad(mag(fvc::grad(B)));

forAll (Krxn, celli)
    {
        if (mag(SolidBorder[celli]) != 0 && Solid[celli] ==0  )
           {Krxn[celli] = krxn.value();}
	else
	   {Krxn[celli] = 0;}
    }
}

if(krxn.value()<0) //if it is a decomposition reaction 
{

forAll(epss,celli)
{
        if(epss[celli] > 1e-20)
        {
        Solid[celli]=1;
        }
        else
        {
        Solid[celli]=0;
        }
}

Solid.correctBoundaryConditions();

SolidBorder = fvc::grad(Solid);
//Bborder = fvc::grad(mag(fvc::grad(B)));

forAll (Krxn, celli)
    {
            if (mag(SolidBorder[celli]) != 0 && Solid[celli] ==1  )
            {Krxn[celli] = krxn.value();}
            else
            {Krxn[celli] = 0;}
     }
}

Info<< "\nCalculating epss" << endl;
fvScalarMatrix epssEqn
(       
        fvm::ddt(epss)
      + fvm::div(phis,epss)
      -  Krxn //units of 1/(seconds * control Vol)
);

epssEqn.solve();
epss.correctBoundaryConditions();

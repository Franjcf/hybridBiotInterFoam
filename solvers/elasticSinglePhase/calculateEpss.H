//Updating speed and surface fluxes

Us = fvc::ddt(D);
Us.correctBoundaryConditions();
phis = fvc::interpolate(Us) & mesh.Sf();
phic=-fvc::interpolate(epss)*phis;

//Solving solid conservation equation 

Info<< "\nCalculating epss" << endl;
fvScalarMatrix epssEqn
(
	fvm::ddt(epss)
      + fvm::div(phis,epss)
      - (Krxn)
);

epssEqn.solve();
epss.correctBoundaryConditions();

//Calculating the superficial displacement and velocity averages. 
Dsup=D*epss;
UsSup=Us*epss;

//Defining reactions at the interfaces (for biofilm growth)

if(krxn.value()>0) //if growing reaction 
{

forAll(epss,celli)
{
        if(epss[celli] >= epssEq.value())
        {
        B[celli]=1;
        }
        else
        {
        B[celli]=0;
        }
}

B.correctBoundaryConditions();

Bborder = fvc::grad(B);

forAll (Krxn, celli)
    {
        if (mag(Bborder[celli]) != 0 && B[celli] ==0  )
           {Krxn[celli] = krxn.value();}
        else
           {Krxn[celli] = 0;}
    }
}

if(krxn.value()<0) // if decaying reaction
{

forAll(epss,celli)
{
        if(epss[celli] > 1e-20)
        {
        B[celli]=1;
        }
        else
        {
        B[celli]=0;
        }
}

B.correctBoundaryConditions();

Bborder = fvc::grad(B);

forAll (Krxn, celli)
    {
            if (mag(Bborder[celli]) != 0 && B[celli] ==1  )
            {Krxn[celli] = krxn.value();}
            else
            {Krxn[celli] = 0;}
     }
}


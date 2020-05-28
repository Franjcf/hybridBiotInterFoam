// Defining Permeability Models

Info<< "Reading intrinsic permeability field field K\n" << endl;
volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField k0("k0", K);

scalar useKozenyCarman(transportProperties.lookupOrDefault<scalar>("useKozenyCarman",1));

if (useKozenyCarman == 1)
{k0 = K*(eps*eps*eps)/(max((1-eps)*(1-eps),SMALL)) ;}

//Drag Coefficient Calculation
volScalarField Drag ("Drag", nuf/k0);

// A. Ferrari 2015 2.5 D drag coefficient
dimensionedScalar invFk0("Fk0", transportProperties.lookup("invFk0",0.));
volScalarField DragF("DragF", (1-Solid)*invFk0*nuf);

// Defining Reaction Variabels at Interface
Info<< "Reading reaction rate constant\n" << endl;
dimensionedScalar krxn
(
    "krxn",
    transportProperties.lookup("krxn")
);

Info<< "Reading equilibrium porosity\n" << endl;
dimensionedScalar epssEq
(
    "epssEq",
    transportProperties.lookup("epssEq")
);

Info<< "Reaction Rate Constant\n" << endl;
volScalarField Krxn
(   
    IOobject
    (   
        "Krxn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    krxn
);

// Labeling the places that have porous media
forAll(epss,celli) 
{
	if(epss[celli]>0.01) 
	{	
	Solid[celli]=1; 
	}
	else
	{
	Solid[celli]=0; 
	}
}

volVectorField SolidBorder(fvc::grad(Solid));

forAll(mesh.boundary(), patchi) 
{
	forAll(Solid.boundaryField()[patchi],facei)
	{
		if(epss.boundaryField()[patchi][facei]>0.01) 
		{Solid.boundaryFieldRef()[patchi][facei]=1;}
		else
		{Solid.boundaryFieldRef()[patchi][facei]=0;}
	}
}

surfaceScalarField Solidf(fvc::interpolate(Solid));

// Eliminating intermediate values of Solid Indicator Function
forAll(Solidf,facei)
{
	if(Solidf[facei]<1) // if <1 ==0.  
	{	
	Solidf[facei]=0; 
	}
}


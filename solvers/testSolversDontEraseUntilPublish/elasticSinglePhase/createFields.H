Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//Reading constants that parametrize the model

Info<< "Reading the fluid's kinematic viscocity 'nuf'\n" << endl;
dimensionedScalar nuf
(
    "nuf",
    transportProperties.lookup("nuf")
);

Info<< "Reading Kozeny-Carman Constant 'k0inv'\n" << endl;
dimensionedScalar k0inv
(
    "k0inv",
    transportProperties.lookup("k0inv")
);

Info<< "Reading gravity vector\n" << endl;
dimensionedVector g
(
    "g",
    transportProperties.lookup("g")
);

Info<< "Reading diffusivity of ions in fluid 'Df'\n" << endl;
dimensionedScalar Df
(
    "Df",
    transportProperties.lookup("Df")
);

Info<< "Reading diffusivity tourtuosity coefficient n\n" << endl;
dimensionedScalar n
(
    "n",
    transportProperties.lookup("n")
);

Info<< "Reading artificial clay diffusion (for Correction) 'Dc'\n" << endl;
dimensionedScalar Dc
(
    "Dc",
    transportProperties.lookup("Dc")
);

Info<< "Reading clay sheet thickness 'clayThickness'\n" << endl;
dimensionedScalar clayThickness
(
    "clayThickness",
    transportProperties.lookup("clayThickness")
);

Info<< "Reading the swelling force Coeff. 'sigma'\n" << endl;
dimensionedScalar sigma
(
    "sigma",
    transportProperties.lookup("sigma")
);

Info<< "Reading density of clay 'rhoc'\n" << endl;
dimensionedScalar rhoc
(
    "rhoc",
    transportProperties.lookup("rhoc")
);

Info<< "Reading density of the fluid 'rhof'\n" << endl;
dimensionedScalar rhof
(
    "rhof",
    transportProperties.lookup("rhof")
);

Info<< "Reading bulk viscosity of clay 'zeta'\n" << endl;
dimensionedScalar zeta
(
    "zeta",
    transportProperties.lookup("zeta")
);

Info<< "Reading equilibrium porosity\n" << endl;
dimensionedScalar epssEq
(
    "epssEq",
    transportProperties.lookup("epssEq")
);

Info<< "Reading reaction rate constant\n" << endl;
dimensionedScalar krxn
(
    "krxn",
    transportProperties.lookup("krxn")
);

Info<< "Define if we want to solve Convection Eqn.\n" << endl;
dimensionedScalar SolveConvectionEqn
(
    "SolveConvectionEqn",
    transportProperties.lookup("SolveConvectionEqn")
);

//Initializing Some Dimensionalized Scalar and Vectors////

dimensionedScalar ZeroValue("ZeroValue", dimensionSet(0,0,0,0,0,0,0), 0); 

dimensionedScalar ZeroLength("ZeroLength", dimensionSet(0,1,0,0,0,0,0), 0); 

dimensionedScalar ZeroP("ZeroP", dimensionSet(1,-1,-2,0,0,0,0), 0);

dimensionedScalar Linv("Linv", dimensionSet(0,-1,0,0,0,0,0), 0);

dimensionedScalar correction("correction", dimensionSet(0,0,0,0,0,0,0), 1e-20); 

/////////Creating fields in Mesh//////////////////////

Info<< "Reading field fluid pressure field 'p'\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading fluid velocity field 'U'\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading solid fraction field 'epss'\n" << endl;
volScalarField epss
(
    IOobject
    (
        "epss",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading initial porous field 'epssInitial'\n" << endl;
volScalarField epssInitial
(
    IOobject
    (
        "epssInitial",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT, 
        IOobject::NO_WRITE 
    ),
    epss
);

Info<< "Reading places where Solid Fraction has changed 'epssChange' \n" << endl;
volScalarField epssChange
(
    IOobject
    (
        "epssChange",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar ("epssChange", dimensionSet(0,0,0,0,0,0,0), 0)
);

Info<< "Reading fluid fraction field 'epsf'\n" << endl;
volScalarField epsf
(
    IOobject
    (
        "epsf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE 
    ),
    (1-epss)
);

Info<< "Reading total solid pressure 'Ps'\n" << endl;
volScalarField Ps
(
    IOobject
    (
        "Ps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    ZeroP
);

Info<< "Reading indicator variable B\n" << endl;
volScalarField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("B", dimensionSet(0,0,0,0,0,0,0), 0)
);

Info<< "Reading inverse permeability field 'kinv'\n" << endl;
volScalarField kinv
(
    IOobject
    (
        "kinv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE 
    ),
    k0inv*(epsf+epss) //just a placeholder for constant perm
    //(k0inv*(epss*epss)/((1-epss)*(1-epss)*(1-epss))) //for variable perm
);

Info<< "Reading Corrected Ion Diffusivity 'd'\n" << endl;
volScalarField d
(
    IOobject
    (
        "d",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE 
    ),
    Df*pow((1-epss),(n-1))
);

Info<< "Reading Ionic Concentration 'C'\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
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

#include "createPhi.H"
#include "createSurfaceFluxes.H"

Info<< "Reading Transport Momentum Source Term\n" << endl;
volVectorField SourceTransport
(
    IOobject
    (
        "SourceTransport",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	epsf*nuf*kinv*rhof/rhoc*(U)

);

Info<< "Reading Swelling Momentum Source Term\n" << endl;
volVectorField SourceSwelling
(
    IOobject
    (
        "SourceSwelling",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    epss*(sigma/rhoc)*fvc::grad(Ps)
);

volVectorField Bborder
(
    IOobject
    (
        "Bborder",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
      fvc::grad(B)
);

forAll(epss,celli) //P-pore only in clay
{
        if(epss[celli]>1e-20)
        {
        B[celli]=1; //in clay (use porous Pc models)
        }
        else
        {
         B[celli]=0; //in free pore (use large scale Pc model)
        }
}

///////////////////////

//Creating Mass Conservation parameters (can delete)

//Initial Values of Solid Fraction

//Epss of Places that gave changed
dimensionedScalar Epss0("Epss0", dimensionSet(0,0,0,0,0,0,0), 0); 
//Initial Epss of places that have changed
dimensionedScalar EpssIn("EpssIn", dimensionSet(0,0,0,0,0,0,0), 0);
dimensionedScalar localEpssChange("localEpssChange", dimensionSet(0,0,0,0,0,0,0), 0);

//for mass conservation indicators
dimensionedScalar initialEpss("initialEpss", dimensionSet(0,0,0,0,0,0,0), 0);
dimensionedScalar initialVol("initialEpss", dimensionSet(0,0,0,0,0,0,0), 0);
dimensionedScalar finalEpss("finalEpss", dimensionSet(0,0,0,0,0,0,0), 0); 
dimensionedScalar finalVol("finalEpss", dimensionSet(0,0,0,0,0,0,0), 0);

dimensionedScalar PercentEpssChange("PercentEpssChange", 0);
dimensionedScalar PercentVolChange("PercentVolChange", 0);

dimensionedScalar redCount("redCount", 0);

/////////////
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

//including the swelling pressure model
#include "createSwellingPressureModelFields.H"

//Setting permeabilities if we have constant permeability in the solid
forAll (epss,celli)
{
	if (epss[celli]<1e-1)
	{
		kinv[celli]=1e-20;
	}
	else 
	{kinv[celli]=k0inv.value();}
}

#include "createMRF.H"

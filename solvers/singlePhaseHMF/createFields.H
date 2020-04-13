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

Info<< "Reading 2D extra permeability term\n" << endl;
dimensionedScalar twoDkinv
(
    "twoDkinv",
    transportProperties.lookup("twoDkinv")
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

Info<< "Reading density of porous media 'rhos'\n" << endl;
dimensionedScalar rhos
(
    "rhos",
    transportProperties.lookup("rhos")
);

Info<< "Reading density of the fluid 'rhof'\n" << endl;
dimensionedScalar rhof
(
    "rhof",
    transportProperties.lookup("rhof")
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

Info<< "Reading porous media velocity field 'Us'\n" << endl;
volVectorField Us
(
    IOobject
    (
        "Us",
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

Info<< "Reading fluid fraction field 'eps'\n" << endl;
volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE   //change to no write for faster
    ),
    (1-epss)
);

Info<< "Reading Solid indicator variable\n" << endl;
volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Solid", dimensionSet(0,0,0,0,0,0,0), 0)
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
        IOobject::NO_WRITE 
    ),
    (k0inv*(epss*epss)/(max((1-epss)*(1-epss)*(1-epss),SMALL)))
    
);

Info<< "Reading Corrected Ion Diffusivity 'D'\n" << endl;
volScalarField D
(
    IOobject
    (
        "D",
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

singlePhaseTransportModel solid(Us, phis);

volVectorField SolidBorder
(
    IOobject
    (
        "SolidBorder",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
      fvc::grad(Solid)
);

//Setting solid phase indicator variable
forAll(epss,celli) 
{
	if(epss[celli]>1e-20) //0.01
        {
        Solid[celli]=1; 
        }
        else
        {
        Solid[celli]=0; 
        }
}

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

//Reading Max Courant Number for Porous Media
scalar maxCoUs = runTime.controlDict().lookupOrDefault<scalar>("maxCoUs", 1.0);
maxCoUs=maxCoUs;//this prevents the "unused variable" warning


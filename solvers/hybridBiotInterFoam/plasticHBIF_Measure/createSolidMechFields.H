//Creating Fields for Solid Mechanics 

Info<< "Reading Solid Velocity Field\n" << endl;
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

#include "createSurfaceFluxesSolid.H"

// Defining Solid Transport Model
singlePhaseTransportModel solid(Us, phis);     

//Reading Solid Density
dimensionedScalar rho_s ("rho_s",transportProperties.lookup("rho_s"));

// Reading solid Mechanics activator
scalar activateSolidMechanics (transportProperties.lookupOrDefault<scalar>("activateSolidMechanics",0.));
activateSolidMechanics = activateSolidMechanics; //prevents unused variable warning

// Reading Maximum Courant Number for Solid 
scalar maxCoUs = runTime.controlDict().lookupOrDefault<scalar>("maxCoUs", 1.0);
maxCoUs=maxCoUs; //prevents unused variable warning

Info<< "Reading Solid Velocity Field\n" << endl;
volVectorField CapTerm
(   
    IOobject
    (   
        "CapTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (Solid/rho_s)*
                         (
                           (eps*PcCoeff*dPcdS - eps*pc)*fvc::grad(alpha1)
                         +
                            eps*PcCoeff*(rho1-rho2)*g
                         )

);

Info<< "Reading Solid Velocity Field\n" << endl;
volVectorField SolidDrag
(
    IOobject
    (
        "SolidDrag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    epss*eps*Drag/rho_s*Us - eps*Drag/rho_s*U
);

Info<< "Reading Solid Velocity Field\n" << endl;
volVectorField Terzaghi
(
    IOobject
    (
        "Terzaghi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -(epss/rho_s)*fvc::grad(p)
);

Info<< "Reading Solid Velocity Field\n" << endl;
volScalarField StrainRateU
(
    IOobject
    (
        "StrainRateU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.414*mag(symm(fvc::grad(U)))
);

Info<< "Reading Solid Velocity Field\n" << endl;
volScalarField StrainRateUs
(
    IOobject
    (
        "StrainRateUs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    1.414*mag(symm(fvc::grad(Us)))
);


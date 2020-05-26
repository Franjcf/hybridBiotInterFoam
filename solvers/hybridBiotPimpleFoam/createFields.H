Info<< "Reading field p\n" << endl;
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

Info<< "Reading solid fraction field\n" << endl;     
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

//Defining Fluid Fraction Field (porosity)
volScalarField eps("eps",(1-epss));

// Defenition of Solid Indicator
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

#include "createPhi.H"

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

dimensionedScalar nuf
(
    "nuf",
    transportProperties.lookup("nuf")
);

dimensionedScalar rho
(
    "rho",
    transportProperties.lookup("rho")
);

// Defining Mass Flux
Info<< "Calculating superficial mass flux\n" << endl;
surfaceScalarField phiByEps
(
	IOobject
	(
        "phiByEps",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    phi*linearInterpolate(1/(eps+SMALL))
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p.name());

#include "createPorousMediaFields.H"
#include "createSolidMechFields.H"

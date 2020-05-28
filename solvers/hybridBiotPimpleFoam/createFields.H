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

// Reading Model Parameters //
Info<< "Reading the fluid's kinematic viscocity 'nuf'\n" << endl;
dimensionedScalar nuf
(
    "nuf",
    transportProperties.lookup("nuf")
);

Info<< "Reading density of fluid 'rhof'\n" << endl;
dimensionedScalar rhof
(
    "rhof",
    transportProperties.lookup("rhof")
);

Info<< "Reading density of solid 'rhos'\n" << endl;
dimensionedScalar rhos
(
    "rhos",
    transportProperties.lookup("rhos")
);

Info<< "Reading inverse of intrinsic permeability 'k0inv'\n" << endl;
dimensionedScalar k0inv
(
    "k0inv",
    transportProperties.lookup("k0inv")
);

Info<< "Reading Ferrari Drag term\n" << endl;
dimensionedScalar Fkinv
(
    "Fkinv",
    transportProperties.lookup("Fkinv")
);

Info<< "Reading Terzaghi Coeff. 'sigma'\n" << endl;
dimensionedScalar sigma
(
    "sigma",
    transportProperties.lookup("sigma")
);

// Creating Fields //

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

Info<< "Readin solid equilibrium indicator variable Solid\n" << endl;
volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar ("Solid", dimensionSet(0,0,0,0,0,0,0), 0)
);

volVectorField SolidBorder(fvc::grad(Solid));

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

#include "createPhi.H"
#include "createSolidMechFields.H"
#include "createTransportRxnFields.H"

//Initializing Terzaghi Stress 
volVectorField Terzaghi= -(sigma*rhof/rhos)*epss*fvc::grad(p);

//Pressure reference cells
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

//Reading Max Courant Number for Porous Media
scalar maxCoUs = runTime.controlDict().lookupOrDefault<scalar>("maxCoUs", 1.0);
maxCoUs=maxCoUs;




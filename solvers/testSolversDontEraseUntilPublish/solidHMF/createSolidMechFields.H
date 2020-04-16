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

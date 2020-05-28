// Solid Mechanics Fields //

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

#include "createSurfaceFluxes.H"

//Choosing solid viscosity model
singlePhaseTransportModel solid(Us, phis);

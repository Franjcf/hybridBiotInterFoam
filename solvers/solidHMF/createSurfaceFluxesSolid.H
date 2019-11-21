// Defining fluxes based on porous media velocity

surfaceScalarField phis
(
	IOobject
	(
        	"phis",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::NO_WRITE
    ),
    linearInterpolate(Us) & mesh.Sf()
);

surfaceScalarField phisEpss
(
	IOobject
	(
        	"phisEpss",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::NO_WRITE
    ),
    phis*linearInterpolate(epss)
);


//Fluxes based on porous media velocity

Info<< "Reading/calculating intrinsic solid face flux field 'phis'\n" << endl;
surfaceScalarField phis
(
	IOobject
	(
        	"phis",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::AUTO_WRITE
    ),
    linearInterpolate(Us) & mesh.Sf()
);

Info<< "Reading/calculating superficial solid face flux field 'phisEpss'\n" << endl;
surfaceScalarField phisEpss
(
	IOobject
	(
        	"phisEpss",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::AUTO_WRITE
    ),
    phis*linearInterpolate(epss)
);

//Modified flux based on fluid velocity and volume fraction
surfaceScalarField phiByEps
(
	IOobject
	(
        	"phiByEps",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::AUTO_WRITE
    ),
    phi*linearInterpolate(1/max(eps,SMALL))
);


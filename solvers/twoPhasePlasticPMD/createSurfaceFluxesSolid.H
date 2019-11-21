//Fluxes based on porous media velocity

Info<< "Reading/calculating intrinsic solid face flux field 'phip'\n" << endl;
surfaceScalarField phip
(
	IOobject
	(
        	"phip",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::NO_WRITE
    ),
    linearInterpolate(Up) & mesh.Sf()
);

Info<< "Reading/calculating superficial solid face flux field 'phipEpss'\n" << endl;
surfaceScalarField phipEpss
(
	IOobject
	(
        	"phipEpss",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::NO_WRITE
    ),
    phip*linearInterpolate(epss)
);

Info<< "Reading scalar transport correction field 'phic'\n" << endl;
surfaceScalarField phic
(
	IOobject
	(
        	"phic",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::NO_WRITE
    ),
      -phipEpss
);

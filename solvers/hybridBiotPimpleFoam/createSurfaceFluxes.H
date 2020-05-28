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

Info<< "Reading scalar transport correction field 'phic'\n" << endl;
surfaceScalarField phic
(
	IOobject
	(
        	"phic",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::AUTO_WRITE
    ),
      -phisEpss
);

//Modified flux based on fluid velocity and volume fraction

Info<< "Reading/calculating solid face flux field phiByEpsf\n" << endl;
surfaceScalarField phiByEpsf
(
	IOobject
	(
        	"phiByEpsf",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::AUTO_WRITE
    ),
    phi*linearInterpolate(1/max(epsf,SMALL))
);


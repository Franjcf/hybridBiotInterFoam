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
    phi*linearInterpolate(1/epsf)
);

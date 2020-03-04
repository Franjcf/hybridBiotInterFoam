//Modified flux based on fluid velocity and volume fraction

Info<< "Reading/calculating solid face flux field rhoPhiByEpsf\n" << endl;
surfaceScalarField rhoPhiByEpsf
(
	IOobject
	(
        	"rhoPhiByEpsf",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
        	IOobject::NO_WRITE
    ),
    rhoPhi*linearInterpolate(1/epsf)
);


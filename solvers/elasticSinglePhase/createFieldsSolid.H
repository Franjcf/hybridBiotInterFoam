#include "readMechanicalProperties.H"
#include "readThermalProperties.H"

Info<< "Reading field D\n" << endl;
volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Us\n" << endl;
volVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::ddt(D)
);

Info<< "Reading/calculating solid face flux field 'phis'\n" << endl;
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
      -phis*linearInterpolate(epss)
);

autoPtr<volScalarField> Tptr(nullptr);

if (thermalStress)
{
    Info<< "Reading field T\n" << endl;
    Tptr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );
}

Info<< "Reading field superficial displacement Dsup\n" << endl;
volVectorField Dsup
(
    IOobject
    (
        "Dsup",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    epss*D
);

Info<< "Reading superfcial velocity field UsSup\n" << endl;
volVectorField UsSup
(
    IOobject
    (
        "UsSup",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    epss*fvc::ddt(D)
);

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    epss*mu*twoSymm(fvc::grad(D)) + epss*lambda*(I*tr(fvc::grad(D)))
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);



if (compactNormalStress)
{
    divSigmaExp -= fvc::laplacian(epss*(2*mu + lambda), D, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div(epss*(2*mu + lambda)*fvc::grad(D), "div(sigmaD)");
}

mesh.setFluxRequired(D.name());

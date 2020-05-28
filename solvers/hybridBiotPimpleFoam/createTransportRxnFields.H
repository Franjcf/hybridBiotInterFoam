// Scalar Transport Fields //

Info<< "Reading diffusivity 'Df'\n" << endl;
dimensionedScalar Df
(
    "Df",
    transportProperties.lookup("Df")
);

Info<< "Reading diffusivity tourtuosity coefficient n\n" << endl;
dimensionedScalar n
(
    "n",
    transportProperties.lookup("n")
);

Info<< "Define if we want to solve Convection Eqn.\n" << endl;
dimensionedScalar SolveConvectionEqn
(
    "SolveConvectionEqn",
    transportProperties.lookup("SolveConvectionEqn")
);

Info<< "Reading Corrected Ion Diffusivity 'D'\n" << endl;
volScalarField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Df*pow((1-epss),(n-1))
);

Info<< "Reading Ionic Concentration 'C'\n" << endl;
volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading Superficial Ionic Concentration 'C'\n" << endl;
volScalarField Cepsf
(
    IOobject
    (
        "Cepsf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    C*epsf
);


// Rxn Fields //

Info<< "Reading equilibrium porosity\n" << endl;
dimensionedScalar epssEq
(
    "epssEq",
    transportProperties.lookup("epssEq")
);

Info<< "Reading reaction rate constant\n" << endl;
dimensionedScalar krxn
(
    "krxn",
    transportProperties.lookup("krxn")
);

Info<< "Reaction Rate Constant\n" << endl;
volScalarField Krxn
(
    IOobject
    (
        "Krxn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    krxn
);


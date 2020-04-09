//Creating fields for the swelling pressure model

Info<< "Reading avg. of clay particles per stack 'Ns'\n" << endl;
dimensionedScalar Ns
(
    "Ns",
    transportProperties.lookup("Ns")
);

Info<< "Reading irriducible void fraction 'E1'\n" << endl;
dimensionedScalar E1
(
    "E1",
    transportProperties.lookup("E1")
);

Info<< "Reading osmotic swelling transition void fraction 'E2'\n" << endl;
dimensionedScalar E2
(
    "E2",
    transportProperties.lookup("E2")
);

Info<< "Reading crystalline pressure constant 'Pcry0'\n" << endl;
dimensionedScalar Pcry0
(
    "Pcry0",
    transportProperties.lookup("Pcry0")
);

Info<< "Reading constant temperature 'Temp'\n" << endl;
dimensionedScalar Temp
(
    "Temp",
    transportProperties.lookup("Temp")
);

Info<< "Reading Clay Surface Charge Density 'SCharge'\n" << endl;
dimensionedScalar SCharge
(
    "SCharge",
    transportProperties.lookup("SCharge")
);

Info<< "Crystalline Swelling Decay length 'l'\n" << endl;
dimensionedScalar l
(
    "l",
    transportProperties.lookup("l")
);

Info<< "Reading Pore Pressure (Hydrostatic Pressure) 'Ppore'\n" << endl;
dimensionedScalar Ppore
(
    "Ppore",
    transportProperties.lookup("Ppore")
);

Info<< "Reading Contact Stresses (Total Pressure pressing on top of the pores) 'Pconf'\n" << endl;
dimensionedScalar Pconf
(
    "Pconf",
    transportProperties.lookup("Pconf")
);

/////////Creating fields in Mesh//////////////////////

Info<< "Reading Average distance between clay particles 'Havg'\n" << endl;
volScalarField Havg
(
    IOobject
    (
        "Havg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    ZeroLength //this is going to be overwritten first thing
);

Info<< "Reading distance between clay aggregates 'H'\n" << endl;
volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    ZeroLength //this is going to be overwritten first thing
);

Info<< "Reading swelling pressure 'Pswell'\n" << endl;
volScalarField Pswell
(
    IOobject
    (
        "Pswell",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    ZeroP //this is going to be overwritten first thing
);


Info<< "Reading crystaline swelling pressure 'Pccs'\n" << endl;
volScalarField Pccs
(
    IOobject
    (
        "Pccs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    ZeroP //this is going to be overwritten first thing
);

Info<< "Reading diffusive double layer swelling pressure 'Pddl'\n" << endl;
volScalarField Pddl
(
    IOobject
    (
        "Pddl",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    ZeroP //this is going to be overwritten first thing
);

Info<< "Reading van der Waals swelling pressure 'Pvdw'\n" << endl;
volScalarField Pvdw
(
    IOobject
    (
        "Pvdw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    ZeroP //this is going to be overwritten first thing
);

Info<< "Reading average thickness of clay stacks 'Dp'\n" << endl;
volScalarField Dp
(
    IOobject
    (
        "Dp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    ZeroLength //this is going to be overwritten first thing
);

Info<< "Reading inverse Debye Length 'Kdebye'\n" << endl;
volScalarField Kdebye
(
    IOobject
    (
        "Kdebye",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    Linv //this is going to be overwritten first thing
);

Info<< "Reading scaled surface potential of an isolated particle\n" << endl;
// electric potential (also called the electric field potential, potential drop or the electrostatic potential) is the amount of work needed to move a unit positive charge from a reference point to a specific point inside the field without producing any acceleration.
volScalarField y0inf
(
    IOobject
    (
        "y0inf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    ZeroValue //this is going to be overwritten first thing
);

///////////  DEFINING MODEL CONSTANTS ///////////////

//distance between unit layers (can calculate as weighted avg later)
dimensionedScalar dm("dm", dimensionSet(0,1,0,0,0,0,0), 1.18e-9); 

//valence of ions (can be calculated as weighted avg later)
dimensionedScalar v("v", dimensionSet(0,0,0,0,0,0,0), 1.06);

//Farady Constant
dimensionedScalar F("F", dimensionSet(0,0,1,0,-1,1,0), 9.64853e4);

//Permitivity of free vaccum
dimensionedScalar Perm("Perm", dimensionSet(-1,-3,4,0,0,2,0), 8.8542e-12);

//Ideal Gas Constant
dimensionedScalar R("R", dimensionSet(1,2,-2,-1,-1,0,0), 8.314);

//Hamaker's Constant
dimensionedScalar Ah("Ah", dimensionSet(1,2,-2,0,0,0,0), 2.2e-20);

//Relative Dialectric constant of water 
dimensionedScalar Dia("Dia", dimensionSet(0,0,0,0,0,0,0), 78.54);

//Minimum Avg distance between particles
dimensionedScalar HavgMin("HavgMin", dimensionSet(0,1,0,0,0,0,0), 0);


/////////////// CALCULATING VALUES OF CONSTANTS ////////

HavgMin = clayThickness*E1;

// dm = 1/CEC*(...)

// v = 1/CEC*(...)




//Reading the transportProperties Dictionary. This has basically all the user model inputs. 
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//Fluid Related Properties

Info<< "Reading viscosity of fluid 1\n" << endl;
dimensionedScalar mu1
(
    "mu1",
    transportProperties.lookup("mu1")
);

Info<< "Reading viscosity of fluid 2\n" << endl;
dimensionedScalar mu2
(
    "mu2",
    transportProperties.lookup("mu2")
);

Info<< "Reading diffusivity of ions in fluid 'Df'\n" << endl;
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

//Solid Properties

Info<< "Reading permeability coefficient\n" << endl;
dimensionedScalar k0
(
    "k0",
    transportProperties.lookup("k0")
);

Info<< "Reading clay sheet thickness 'clayThickness'\n" << endl;
dimensionedScalar clayThickness
(
    "clayThickness",
    transportProperties.lookup("clayThickness")
);

Info<< "Reading the swelling force Coeff. 'sig'\n" << endl;
dimensionedScalar sig
(
    "sig",
    transportProperties.lookup("sig")
);

Info<< "Reading surface tension sigma\n" << endl;
dimensionedScalar sigma
(
    "sigma",
    transportProperties.lookup("sigma")
);

Info<< "Reading density of clay 'rhoc'\n" << endl;
dimensionedScalar rhoc
(
    "rhoc",
    transportProperties.lookup("rhoc")
);

//Two-Phase Properties

Info<< "Reading Capillary Pressure Coeff 'm'\n" << endl;
dimensionedScalar m
(
    "m",
    transportProperties.lookup("m")
);

Info<< "Reading Capillary Pressure Coeff 'a'\n" << endl;
dimensionedScalar a
(
    "a",
    transportProperties.lookup("a")
);

Info<< "Reading alpha1IrrPerm'\n" << endl;
dimensionedScalar alpha1IrrPerm
(
    "alpha1IrrPerm",
    transportProperties.lookup("alpha1IrrPerm")
);

Info<< "Reading alpha2IrrPerm'\n" << endl;
dimensionedScalar alpha2IrrPerm
(
    "alpha2IrrPerm",
    transportProperties.lookup("alpha2IrrPerm")
);

Info<< "Reading alpha1IrrPres'\n" << endl;
dimensionedScalar alpha1IrrPres
(
    "alpha1IrrPres",
    transportProperties.lookup("alpha1IrrPres")
);

Info<< "Reading alpha1MaxPres'\n" << endl;
dimensionedScalar alpha1MaxPres
(
    "alpha1MaxPres",
    transportProperties.lookup("alpha1MaxPres")
);

Info<< "Reading capillary entry pressure\n" << endl;
dimensionedScalar PcEntry
(
    "PcEntry",
    transportProperties.lookup("PcEntry")
);

Info<< "Reading Solve With Ur Switch\n" << endl;
dimensionedScalar SolveWithUr
(
    "SolveWithUr",
    transportProperties.lookup("SolveWithUr")
);

Info<< "Reading VariablePc\n" << endl;
dimensionedScalar variablePc
(
    "variablePc",
    transportProperties.lookup("variablePc")
);

Info<< "Reading 2D extra permeability term\n" << endl;
dimensionedScalar twoDkinv
(
    "twoDkinv",
    transportProperties.lookup("twoDkinv")
);


//Initializing Some Dimensionalized Scalar and Vectors

dimensionedScalar ZeroValue("ZeroValue", dimensionSet(0,0,0,0,0,0,0), 0); 

dimensionedScalar ZeroVelocity("ZeroVelocity", dimensionSet(0,1,-1,0,0,0,0), 0);

dimensionedScalar ZeroLength("ZeroLength", dimensionSet(0,1,0,0,0,0,0), 1); 

dimensionedScalar ZeroP("ZeroP", dimensionSet(1,-1,-2,0,0,0,0), 0);

dimensionedScalar Linv("Linv", dimensionSet(0,-1,0,0,0,0,0), 0);

dimensionedScalar correction("correction", dimensionSet(0,0,0,0,0,0,0), 1e-20); 

/////////Creating Fields in Mesh//////////////////////

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "gradP_rgh\n" << endl;
volVectorField gradP_rgh
(
    IOobject
    (
        "gradP_rgh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(p_rgh)
);

Info<< "Reading porous media velocity field 'Up'\n" << endl;
volVectorField Up
(
    IOobject
    (
        "Up",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading fluid relative velocity field 'Ur'\n" << endl;
volVectorField Ur
(
    IOobject
    (
        "Ur",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    Up//overwriten before used
);


Info<< "Reading solid fraction field epss\n" << endl;
volScalarField epss
(
    IOobject
    (
        "epss",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading initial porous field 'epssInitial'\n" << endl;
volScalarField epssInitial
(
    IOobject
    (
        "epssInitial",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    epss
);

Info<< "Reading places where Solid Fraction has changed 'epssChange' \n" << endl;
volScalarField epssChange
(
    IOobject
    (
        "epssChange",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar ("epssChange", dimensionSet(0,0,0,0,0,0,0), 0)
);

Info<< "Reading field fluid fraction field epsf\n" << endl;
volScalarField epsf
(
    IOobject
    (
        "epsf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    (1-epss)
);

Info<< "Reading total solid pressure 'Ps'\n" << endl;
volScalarField Ps
(
    IOobject
    (
        "Ps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    ZeroP
);

Info<< "Reading indicator variable B\n" << endl;
volScalarField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("B", dimensionSet(0,0,0,0,0,0,0), 0) 
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"
#include "createSurfaceFluxesSolid.H"

///////Defining the Mixtures and Fluids

//setting the solid to a plastic rheology model
singlePhaseTransportModel fluid(Up, phip);

//settin the fluid to a two phase VOF model
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);
volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

//////

///////Defining Remaining fields

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
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

#include "createSurfaceFluxesTwoPhase.H"

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES flux from previous time-step
surfaceScalarField alphaPhi
(
    IOobject
    (
        "alphaPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    phi*fvc::interpolate(alpha1)
);


/////Permeability models fields

Info<< "Reading initial absolute permeability\n" << endl;
volScalarField k0real
(
    IOobject
    (
        "k0real",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    k0*((1-epss)*(1-epss)*(1-epss))/(max(epss*epss,SMALL))
);

Info<< "Reading k1rel\n" << endl;
volScalarField k1rel
(
    IOobject
    (
        "k1rel",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar ("k1rel", dimensionSet(0,0,0,0,0,0,0), 1) 
);

Info<< "Reading k2rel\n" << endl;
volScalarField k2rel
(
    IOobject
    (
        "k2rel",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar ("k2rel", dimensionSet(0,0,0,0,0,0,0), 1) 
);

Info<< "Reading inverse permeability field kinvMu\n" << endl;
volScalarField kinvMu
(
    IOobject
    (
        "kinvMu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (1/k0real)*(1/(k1rel/mu1+k2rel/mu2))
);

Info<< "Reading alpha1effPerm\n" << endl;
volScalarField alpha1effPerm
(
    IOobject
    (
        "alpha1effPerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha1
	
);

Info<< "Reading alpha1effPres\n" << endl;
volScalarField alpha1effPres
(
    IOobject
    (
        "alpha1effPres",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    alpha1
	
);			

/////Capillary Pressure Models Fields

Info<< "Reading Capillary Entry Pressure Pc0\n" << endl;
volScalarField Pc0
(
    IOobject
    (
        "Pc0",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    sigma/ZeroLength //will be changed first thing
);

Info<< "Reading Capillary Pressure Coefficient PcCoeff\n" << endl;
volScalarField PcCoeff
(
    IOobject
    (
        "PcCoeff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    Pc0 //will be changed first thing
);

Info<< "Reading gradPc\n" << endl;
volVectorField gradPc
(
    IOobject
    (
        "gradPc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    PcCoeff*fvc::grad(alpha1)/rhoc
);

Info<< "Reading gradP\n" << endl;
volVectorField gradP
(
    IOobject
    (
        "gradP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(p)
);

Info<< "Reading Transport Momentum Source Term\n" << endl;
volVectorField SourceTransport
(
    IOobject
    (
        "SourceTransport",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	epsf*kinvMu/rhoc*(U-Up*epss)

);

Info<< "Reading Swelling Momentum Source Term\n" << endl;
volVectorField SourceSwelling
(
    IOobject
    (
        "SourceSwelling",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    epss*(sig/rhoc)*fvc::grad(Ps)
);

volVectorField Bborder
(
    IOobject
    (
        "Bborder",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(mag(fvc::grad(B))) 
);

//Info<< "Solid Shear\n" << endl;
//volScalarField shear
//(
//    IOobject
//    (
//        "shear",
//        runTime.timeName(),
//        mesh,
//        IOobject::NO_READ,
//        IOobject::NO_WRITE
//    ),
//	1.4141*mag(symm(fvc::grad(Up)))
//);

	
// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;

#include "createMRF.H"

//updating indicator variable for capillary pressure
forAll(epss,celli) //P-pore only in clay
{
	if(epss[celli]>0.01) //0.001 Good and 0.0001 also good
	{	
	B[celli]=1; //in clay (use porous Pc models)
	}
	else
	{
	B[celli]=0; //in free pore (use large scale Pc model)
	}
}

//Reading Max Courant Number for Porous Media

scalar maxCoUp = runTime.controlDict().lookupOrDefault<scalar>("maxCoUp", 1.0);

maxCoUp=maxCoUp;//this prevents the "unused variable" warning


//including the swelling pressure model
#include "createSwellingPressureModelFields.H"

///////////////////////

//Creating Mass Conservation parameters (can erase this)

//Initial Values of Solid Fraction

//Epss of Places that gave changed
dimensionedScalar Epss0("Epss0", dimensionSet(0,0,0,0,0,0,0), 0); 
//Initial Epss of places that have changed
dimensionedScalar EpssIn("EpssIn", dimensionSet(0,0,0,0,0,0,0), 0);
dimensionedScalar localEpssChange("localEpssChange", dimensionSet(0,0,0,0,0,0,0), 0);

//for mass conservation indicators
dimensionedScalar initialEpss("initialEpss", dimensionSet(0,0,0,0,0,0,0), 0);
dimensionedScalar initialVol("initialEpss", dimensionSet(0,0,0,0,0,0,0), 0);
dimensionedScalar finalEpss("finalEpss", dimensionSet(0,0,0,0,0,0,0), 0); 
dimensionedScalar finalVol("finalEpss", dimensionSet(0,0,0,0,0,0,0), 0);

dimensionedScalar PercentEpssChange("PercentEpssChange", 0);
dimensionedScalar PercentVolChange("PercentVolChange", 0);

dimensionedScalar redCount("redCount", 0);




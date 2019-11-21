//Update Variables

//Calculating Changes in Permeability

Info<< "\nUpdating Relative Permeability \n" << endl;

#include "relativePerm.H"

//Pressure Models

Info<< "\nUpdating Capillary Pressure Models \n" << endl;

#include "PcModels.H"

//updating diffusivity and other variables

epsf=1-epss;

gradP = fvc::grad(p);

//shear = 1.4141*mag(symm(fvc::grad(Up)));

//D=Df*pow((1-epss),(n-1));





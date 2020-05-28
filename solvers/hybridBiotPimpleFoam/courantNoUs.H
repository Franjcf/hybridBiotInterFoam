//Determining the Courant Number of the Solid Phase Velocity

scalar CoNumUs = 0.0;
scalar meanCoNumUs = 0.0;

if (mesh.nInternalFaces())
{
    scalarField sumPhis
    (
        fvc::surfaceSum(mag(phis))().internalField()
    );

    CoNumUs = 0.5*gMax(sumPhis/mesh.V().field())*runTime.deltaTValue();

    meanCoNumUs =
        0.5*(gSum(sumPhis)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Us Courant Number mean: " << meanCoNumUs
    << " max: " << CoNumUs << endl;



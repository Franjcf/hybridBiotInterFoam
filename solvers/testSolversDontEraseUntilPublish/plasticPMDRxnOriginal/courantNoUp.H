//Determining the Courant Number of the Solid Phase Velocity

scalar CoNumUp = 0.0;
scalar meanCoNumUp = 0.0;

if (mesh.nInternalFaces())
{
    scalarField sumPhip
    (
        fvc::surfaceSum(mag(phip))().internalField()
    );

    CoNumUp = 0.5*gMax(sumPhip/mesh.V().field())*runTime.deltaTValue();

    meanCoNumUp =
        0.5*(gSum(sumPhip)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Up Courant Number mean: " << meanCoNumUp
    << " max: " << CoNumUp << endl;



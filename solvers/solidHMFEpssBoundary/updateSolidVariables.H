// Updating places that have porous media
forAll(epss,celli)
{
	if(epss[celli]>0.1) 
        {
        Solid[celli]=1;
        }
        else
        {
        Solid[celli]=0;
        }
}

Solid.correctBoundaryConditions();

// Updating boundary patches
forAll(mesh.boundary(), patchi)
{
        forAll(Solid.boundaryField()[patchi],facei)
        {
                if(epss.boundaryField()[patchi][facei]>0.1)
                {Solid.boundaryFieldRef()[patchi][facei]=1;}
                else
                {Solid.boundaryFieldRef()[patchi][facei]=0;}
        }
}

// Updating face-centered values of the Solid field
surfaceScalarField Solidf(fvc::interpolate(Solid));
forAll(Solidf,facei)
{
        if(Solidf[facei]<1) // if <1 ==0.
        {
        Solidf[facei]=0;
        }
}

// Updating Fluid Volume Fraction (porosity)
eps = (1-epss); 

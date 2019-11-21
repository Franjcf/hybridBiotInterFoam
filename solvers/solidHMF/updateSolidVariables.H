// Labeling places that have porous media
forAll(epss,celli)
{
	if(epss[celli]>0.01)
        {
        Solid[celli]=1;
        }
        else
        {
        Solid[celli]=0;
        }
}

//updating boundary patches
forAll(mesh.boundary(), patchi)
{
        forAll(Solid.boundaryField()[patchi],facei)
        {
                if(epss.boundaryField()[patchi][facei]>0.01)
                {Solid.boundaryFieldRef()[patchi][facei]=1;}
                else
                {Solid.boundaryFieldRef()[patchi][facei]=0;}
        }
}

surfaceScalarField Solidf(fvc::interpolate(Solid));

// Eliminating intermediate values of solid indicator function
forAll(Solidf,facei)
{
        if(Solidf[facei]<1) // if <1 ==0.
        {
        Solidf[facei]=0;
        }
}

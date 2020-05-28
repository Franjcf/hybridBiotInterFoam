// Updating Rxn Coefficients //
// defining reaction to just be at the interface

if(krxn.value()>0) //if growth
{

forAll(epss,celli)
{       
        if(epss[celli] >= epssEq.value())
        {
        Solid[celli]=1;
        }
        else
        {
        Solid[celli]=0;
        }
}

Solid.correctBoundaryConditions();

SolidBorder = fvc::grad(Solid);

forAll (Krxn, celli)
    {   
        if (mag(SolidBorder[celli]) != 0 && Solid[celli] ==0  )
           {Krxn[celli] = krxn.value();}
        else
           {Krxn[celli] = 0;}
    }
}

if(krxn.value()<0) //if decay
{

forAll(epss,celli)
{
        if(epss[celli] > 1e-20)
        {
        Solid[celli]=1;
        }
        else
        {
        Solid[celli]=0;
        }
}

Solid.correctBoundaryConditions();

SolidBorder = fvc::grad(Solid);

forAll (Krxn, celli)
    {
            if (mag(SolidBorder[celli]) != 0 && Solid[celli] ==1  )
            {Krxn[celli] = krxn.value();}
            else
            {Krxn[celli] = 0;}
     }
}



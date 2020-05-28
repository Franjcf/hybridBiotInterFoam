// Updating Relative Permeability Model

if (useKozenyCarman == 1)
{k0 = K*(1-epss)*(1-epss)*(1-epss)/max(epss*epss,SMALL);}

//Calculating Drag Coefficient
Drag = nuf/k0;

forAll(Drag,celli)
{
	if (Solid[celli] == 0)
		{Drag[celli]= 0;}
}

//Calculating A. Ferrari 2015 Coefficient (2.5 D term)
DragF = (1-Solid)*invFk0*nuf;

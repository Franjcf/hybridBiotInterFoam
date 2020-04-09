//UPDATING VARIABLES

kinv = k0inv*(epss*epss)/(max((1-epss)*(1-epss)*(1-epss),SMALL));

forAll (epss,celli)
{
	if (epss[celli]<1e-3)
	{
		kinv[celli]=1e-20;
	}
	else 
	{kinv[celli]=k0inv.value();}
}

d=Df*pow((1-epss),(n-1));

epsf=1-epss;






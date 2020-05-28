//UPDATING VARIABLES//

kinv = k0inv*(epss*epss)/(max((1-epss)*(1-epss)*(1-epss),SMALL));

forAll(kinv,celli)
{	
	if (epss[celli]<0.001)
	{kinv[celli]=0;}
}

D=Df*pow((1-epss),(n-1));
epsf=1-epss;






//Relative Permeability Model

//updating non-relative permeabilities

k0real =k0*((1-epss)*(1-epss)*(1-epss))/(max(epss*epss,SMALL));

//k0real.min(k0max);//k0max; 1e-4

//Calculating effective saturation
alpha1effPerm = (alpha1-alpha1IrrPerm)/(1-alpha1IrrPerm-alpha2IrrPerm);

alpha1effPerm.max(1e-3);
alpha1effPerm.min(1-1e-3);

////Brooks and Corey Model Relative Permeabilities
//k1rel=pow(alpha1effPerm,1);
//k2rel=pow((1-alpha1effPerm),1);

//VanGenutchen Model Relative Permeabilities
k1rel=pow(alpha1effPerm,0.5)*pow(1-pow(1-pow(alpha1effPerm,1/m),m),2);
k2rel=pow(1-alpha1effPerm,0.5)*pow(1-pow(alpha1effPerm,1/m),(2*m));

	//WITH DIFFERENT "M" (now "a") THAN P0 MODEL 
	//k1rel=pow(alpha1effPerm,0.5)*pow(1-pow(1-pow(alpha1effPerm,1/a),a),2);
	//k2rel=pow(1-alpha1effPerm,0.5)*pow(1-pow(alpha1effPerm,1/a),(2*a));

//calculating final comprehensive permeability 
kinvMu=(1/k0real)*(1/(k1rel/mu1+k2rel/mu2)); 

forAll(kinvMu,celli)
{
	if (B[celli] == 0) //if in free fluid. 
		{kinvMu[celli]= 0;}
}


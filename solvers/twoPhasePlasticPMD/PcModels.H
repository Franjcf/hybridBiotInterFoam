//Capillary Pressure Models

//Calculating effective saturation
alpha1effPres = (alpha1-alpha1IrrPres)/(alpha1MaxPres-alpha1IrrPres);

//This prevents us from dividing by zero. It basically means that, even though alpha does not reflect the "real" saturation values, are capilarry pressure model does.
alpha1effPres.max(1e-3);
alpha1effPres.min(1-1e-3);

//////Linear Model

//PcCoeff = -(PcMax-Pc0); 

//////Brooks-Corey

//PcCoeff = -(a*Pc0*pow(alpha1effPres,(-a)))/(alpha1MaxPres-alpha1IrrPres);

/////VanGenutchen Model

if (variablePc.value()==1)
{Pc0 = sigma/H;}
else
{Pc0 = PcEntry;}

PcCoeff =  -((1-m)/m)*
		(Pc0/(alpha1MaxPres-alpha1IrrPres))* 
             	pow(pow(alpha1effPres,(-1/m))-1,(-m))* 
	     	pow(alpha1effPres, (-(1+m)/m)  );



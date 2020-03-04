Info<< "\nCalculating Swelling Pressure" << endl;

//This is based on L. Liu Paper Prediction of Swelling Pressures of different types of bentonite in dilute solutions
//Also based on K. Norrish "The swelling of montmorillonite" 
// Also based on Aachanta, S "Nonequilibirum swelling and cappilarry pressure relatinos for colloidal systems"


///////////////CALCULATING CURRENT VALUES///////////////////////////////////////////

//Calculating the average distance between individual clay sheets

Havg = clayThickness*(1/(epss+correction)-1);

//Making sure Havg is not less than the physical minimum:

forAll(Havg,celli)
{
	if ( Havg[celli] < HavgMin.value() )
	{Havg[celli]=HavgMin.value();}
}

//Average distance between clay aggregates

H = (
     clayThickness*
     (Ns.value()-1)* 
     pow( ( (Havg/clayThickness)- E1), 2 )
     )
     /( 2*(E2- E1) )
     + clayThickness*(Havg/clayThickness);

//Calculating Dp, the avg. thickness of clay aggregates

Dp = Ns*clayThickness + (Ns-1)*dm;

//Crystaline Swelling Pressure Calculation

Pccs = Pcry0*exp(-H/l); 

//Van del Waals Swelling Pressure Calculation

Pvdw =  Ah/(6*3.14159)*
	(
	1/pow(H,3) - 2/pow((H+Dp),3) + 1/pow((H+2*Dp),3)
	);  

//Calculating the diffusive double layer contributions to the swelling pressure

Kdebye = sqrt(2*C*v*v*F*F/(Perm*Dia*R*Temp) );

y0inf = 2*asinh(v*F*SCharge/(2*Perm*Dia*Kdebye*R*Temp) );

Pddl =  2*C*R*Temp*
	(cosh
	     (
	     asinh
		  (			
		  2*sinh
		       (
		       4*atanh( tanh(y0inf/4)*exp(-Kdebye*H/2) )
		       )
		  +
		  4/(Kdebye*H)*sinh
	               (
    		       2*atanh( tanh(y0inf/4)*exp(-Kdebye*H) )
		       ) 
                  )
	     )
        -1);

//Calculating Total Swelling Pressure 

//Pswell = Pccs + (Pddl - Pvdw)*alpha1 + neg(-B*(1-B))*sigma/H; //last term is only active at the interface

//Pc0 = sigma/H;

//Pswell = Pccs + (Pddl - Pvdw)*alpha1 + (1-alpha1)*Pc0;
// + neg(-alpha1*(1-alpha1))*sigma/H; //last term is only active at the interface

Pswell = Pccs + Pddl - Pvdw;

//Applying A Murad and Terzaghi's principle:

Info<< "\nCalculating Momentum Source" << endl;

Ps = Pconf - Pswell - Ppore - p; 

//Setting Ps=0 in places with no clay
forAll(Ps,celli) 
{
	if (B[celli]==0) 
	{
		Ps[celli]=0;
	}
}

Ps.correctBoundaryConditions();

Info<< "\nCalculating Non-Newtonian Viscosity" << endl;

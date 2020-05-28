/////////////////////////Calculating Us////////////

Info<< "\nCalculating Us" << endl;

Terzaghi= -(sigma*rhof/rhos)*epss*fvc::grad(p);

fvVectorMatrix UsEqn
	(

	//Inertial Terms
	fvm::ddt(epss,Us) 

	+ fvm::div(phisEpss,Us)
	
	//Viscouus terms: div(mu grad(U*epss)) 
      - fvm::laplacian(solid.nu(),Us)

	//U Transpose viscocity term 
      - fvc::div(solid.nu()*dev2(T(fvc::grad(Us))))

	//SourceSwelling //same as terzaghi
      - Terzaghi 
   
	//Clay Transport Momentum Source
//      - fvc::Sp(epss*(rhof*nuf*kinv)/rhos, U) //basic
      - fvc::Sp(epss*(rhof*nuf*kinv)/rhos, U) 


      //Momentum Transfer to fluid do to relative velocity 
//      + fvm::Sp(epss*(rhof*nuf*kinv)/rhos, Us)  //basic
      + fvm::Sp(epsf*(rhof*nuf*kinv)/rhos, Us)  

);	

UsEqn.solve();
Us.correctBoundaryConditions();
					 
//creating phi equivalent flux in order to calculate divergence (phis,epss)
phis = fvc::interpolate(Us) & mesh.Sf();
phisEpss=fvc::interpolate(epss)*phis; 
phic= -phisEpss;


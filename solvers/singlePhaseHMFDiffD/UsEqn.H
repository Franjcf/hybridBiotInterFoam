/////////////////////////Calculating Up////////////

volVectorField Drag= rhof*eps/rhos*nuf*kinv*(U-Us*epss);
volVectorField Terzaghi= -epss*rhof/rhos*fvc::grad(p);

fvVectorMatrix UsEqn
	(

	//Inertial Terms
	fvm::ddt(epss,Us) 

	+ fvm::div(phisEpss,Us)
	
	//Viscouus terms: div(mu grad(U*epss)) 
      - fvm::laplacian(solid.nu(),Us)

	//U Transpose viscocity term 
      - fvc::div(solid.nu()*dev2(T(fvc::grad(Us))))

	//Bulk viscocity term //may add or not
         //- fvc::grad(zeta*epss*fvc::div(Up))
   
       //Effective Stress
      - Terzaghi

	//Gravity Term
      - epss*g
   
	//Momentum Transfer 1
      - fvc::Sp(epss*(rhof*nuf*kinv)/rhos, U) 
      //- fvc::Sp(eps*(rhof*nuf*kinv)/rhos, U) //standard one  

      //Momentum Trasfer 2
      //+ fvm::Sp(epss*(rhof*nuf*kinv)/rhos, Us)  
      + fvm::Sp(eps*epss*(rhof*nuf*kinv)/rhos, Us)  
);	

UsEqn.solve();
Us.correctBoundaryConditions();
					 
//creating phi equivalent flux in order to calculate divergence (phis,epss)
phis = fvc::interpolate(Us) & mesh.Sf();
phisEpss=fvc::interpolate(epss)*phis; 

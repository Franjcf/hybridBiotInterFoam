/////////////////////////Calculating Up////////////
volVectorField Terzaghi= -epss*rho/rho_s*fvc::grad(p);

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
//      - Terzaghi

	//Gravity Term
      - epss*g
   
	//Momentum Transfer 1
      - fvc::Sp(epss*(rho*Drag)/rho_s, U)   

      //Momentum Trasfer 2 
      + fvm::Sp(epss*(rho*Drag)/rho_s, Us)  
);	

UsEqn.solve();
Us.correctBoundaryConditions();
					 
phis = fvc::interpolate(Us) & mesh.Sf();
phisEpss=fvc::interpolate(epss)*phis; 

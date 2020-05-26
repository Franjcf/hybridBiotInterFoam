Info<< "\nCalculating Us" << endl;

//Calculate terzaghi pressure vector
volVectorField Terzaghi = -(epss*rho/rho_s)*fvc::grad(p); 

fvVectorMatrix UsEqn
	(

	//Inertial Terms
	fvm::ddt(epss,Us) 

      + fvm::div(phisEpss,Us)
	
	//Viscous terms:  
      - fvm::laplacian(solid.nu(),Us)
      - fvc::div(solid.nu()*dev2(T(fvc::grad(Us))))
      
        //Terzaghi Term
      - Terzaghi
   
	//Momentum Source from fluid to solid
      - fvc::Sp(eps*Drag*rho/rho_s, U) 
	
	//Momentum Transfer from solid to fluid 
      + fvm::Sp(epss*eps*Drag*rho/rho_s, Us) 
);	

UsEqn.solve();
Us.correctBoundaryConditions();
					 
phis = fvc::interpolate(Us) & mesh.Sf();
phisEpss=fvc::interpolate(epss)*phis; 


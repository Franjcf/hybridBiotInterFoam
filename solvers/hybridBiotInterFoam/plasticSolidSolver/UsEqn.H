Info<< "\nCalculating Us" << endl;

//Calculate terzaghi pressure vector
volVectorField Terzaghi = -(epss/rho_s)*fvc::grad(p);

//Calculating Explicit Capillary Interaction force
//volVectorField CapTerm (eps/rho_s*(mixture.sigmaK()*(1-Solid) + PcCoeff*Solid)*fvc::grad(alpha1)); 
volVectorField CapTerm ((Solid/rho_s)*
                         (
			   (eps*PcCoeff*dPcdS + epss*pc)*fvc::grad(alpha1)
                         +
                           eps*PcCoeff*(rho1-rho2)*g  
                         )
                       );
// it should be + pc in Cap Term 2 

fvVectorMatrix UsEqn
	(

	//Inertial Terms
	fvm::ddt(epss,Us) 

      + fvm::div(phisEpss,Us)
	
	//Viscouus terms: div(mu grad(U*epss)) 
      - fvm::laplacian(solid.nu(),Us)

	//U Transpose viscocity term 
      - fvc::div(solid.nu()*dev2(T(fvc::grad(Us))))
      
        //Terzaghi Term
      - Terzaghi

	//Gravity Term
      - epss*g 

	//Capillary presure source
      + CapTerm
   
	//Momentum Source from fluid to solid
      - fvc::Sp(eps*Drag/rho_s, U) 
	
	//Momentum Transfer from solid to fluid 
      + fvm::Sp(epss*eps*Drag/rho_s, Us) 
);	

UsEqn.solve();
Us.correctBoundaryConditions();
					 
phis = fvc::interpolate(Us) & mesh.Sf();
phisEpss=fvc::interpolate(epss)*phis; 


Info<< "\nCalculating Us" << endl;

//Calculate terzaghi pressure vector
Terzaghi = -(epss/rho_s)*fvc::grad(p);

//Calculating Explicit Capillary Interaction force
CapTerm = ((Solid/rho_s)*
                         (
			   (eps*PcCoeff*dPcdS - eps*pc)*fvc::grad(alpha1)
                         +
                           eps*PcCoeff*(rho1-rho2)*g  
                         )
                       ); 


SolidDrag =  - eps*Drag/rho_s*U + epss*eps*Drag/rho_s*Us;

StrainRateU = 1.414*mag(symm(fvc::grad(U)));
StrainRateUs = 1.414*mag(symm(fvc::grad(Us)));


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


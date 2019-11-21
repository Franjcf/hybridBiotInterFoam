    fvVectorMatrix UEqn
    (
	//Inertial Terms
        (1./eps)*(fvm::ddt(rho, U) + fvm::div(rhoPhiByEps, U) )

	//Viscous Terms
      //+ turbulence->divDevRhoReff((rho/eps), U)

      - (1/eps)*(fvm::laplacian((rho*mixture.nu()),U))
      - (1/eps)*(fvc::div((rho*mixture.nu())*dev2(T(fvc::grad(U)))))

	//Implicit part of the Darcy-Brinkman Drag Term
      + fvm::Sp(Drag,U) 

        //Explicit part of the Darcy-Brinkman Drag Term
      - fvc::Sp(epss*Drag,Us) 

        //2.5 D term from A. Ferrari 2015
      + fvm::Sp((twoDkinv*rho*mixture.nu()),U)

    );

    UEqn.relax();

    if (pimple.momentumPredictor()) 
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
		    mixture.surfaceTensionForce()*(1.-Solidf) //activate free fluid interface model
	          + fvc::interpolate(PcCoeff)*fvc::snGrad(alpha1)*Solidf //activate porous media model
                  - fvc::snGrad(p)
                ) * mesh.magSf()
		+fvc::interpolate(rho)*(g & mesh.Sf())
            )
        );
    }

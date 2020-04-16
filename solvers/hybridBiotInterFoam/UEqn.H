    fvVectorMatrix UEqn
    (
	//Inertial Terms
        (1/eps)*(fvm::ddt(rho, U) + fvm::div(rhoPhiByEps, U) )

	//Viscous Terms
      - (1/eps)*(fvm::laplacian((rho*mixture.nu()),U))
      - (1/eps)*(fvc::div((rho*mixture.nu())*dev2(T(fvc::grad(U)))))

	//Implicit part of the Darcy-Brinkman Drag Term
      + fvm::Sp(Drag,U) 

        //Explicit part of the Darcy-Brinkman Drag Term
      - fvc::Sp(epss*Drag,Us) 

        //2.5 D term from A. Ferrari 2015
      + fvm::Sp(DragF,U)

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
		    mixture.surfaceTensionForce()*(1.-Solidf) 
	          + fvc::interpolate(PcCoeff*dPcdS-pc)*fvc::snGrad(alpha1)*Solidf 
                  - fvc::snGrad(p)
                ) * mesh.magSf()
		+fvc::interpolate(rho)*(g & mesh.Sf())
            )
        );
    }

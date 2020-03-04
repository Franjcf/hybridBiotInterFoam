    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
	//Inertial Terms
        (1/epsf)*(fvm::ddt(rho, U) + fvm::div(rhoPhiByEpsf, U) )

	//Viscous Terms
      - (1/epsf)*(fvm::laplacian((rho*mixture.nu()),U))
      - (1/epsf)*(fvc::div((rho*mixture.nu())*dev2(T(fvc::grad(U)))))

	//Implicit Darcy-Brinkman Drag Term
      + fvm::Sp(kinvMu,U) 
	
	//Explicit Momentum Exchange Terms
      - fvc::Sp(epss*kinvMu,Up)

        //2.5 D term from A. Ferrari 2015
      + fvm::Sp((twoDkinv*rho*mixture.nu()),U)

    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

//if we use momentum predictor then we would have to change this part and volume average it correctly. 
    if (pimple.momentumPredictor()) 
    {

    surfaceScalarField Bsurf(fvc::interpolate(B));

    forAll (Bsurf,facei) //same correction we use in the alpha Eqn. 
    {
      if (Bsurf[facei]<1)
    	 {Bsurf[facei] = 0;}
    }

        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
		     mixture.surfaceTensionForce()*(1-Bsurf) //activate free fluid interface model
	          +  fvc::interpolate(PcCoeff)*fvc::snGrad(alpha1)*Bsurf //activate porous media model
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }

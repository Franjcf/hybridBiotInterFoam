    fvVectorMatrix UEqn
    (
	//Inertial Terms
        (1/eps)*(fvm::ddt(U) + fvm::div(phiByEps, U) )

	//Viscous Terms
      - (1/eps)*(fvm::laplacian(nuf,U))
      - (1/eps)*(fvc::div(nuf*dev2(T(fvc::grad(U)))))

	//Implicit part of the Darcy-Brinkman Drag Term
      + fvm::Sp(Drag,U) 

        //Explicit part of the Darcy-Brinkman Drag Term
      - fvc::Sp(epss*Drag,Us) 

        //2.5 D term from A. Ferrari 2015
      + fvm::Sp(DragF,U)

    );

    phiByEps=fvc::interpolate(1/max(eps,SMALL))*phi;

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
                  - fvc::snGrad(p)
                ) * mesh.magSf()
            )
        );
    }

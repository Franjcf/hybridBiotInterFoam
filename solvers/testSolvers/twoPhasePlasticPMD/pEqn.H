{
    volScalarField rAU("rAU", 1.0/UEqn.A()); //1/ap
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU)); //interpolated 1/a
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phiHbyA //Hbya with extra transient part
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );
    MRF.makeRelative(phiHbyA);
    adjustPhi(phiHbyA, U, p_rgh);

    surfaceScalarField Bsurf(fvc::interpolate(B));
  
    volScalarField PcCoeff2(PcCoeff);
    forAll (PcCoeff2, celli)
    {
    	if (mag(Bborder[celli]) > 0) 
    	  { PcCoeff2[celli] = 0;}
    }

    surfaceScalarField phig
    (
        (
	    //B=1 in clay, B=0 in free pore
            mixture.surfaceTensionForce()*(1-Bsurf) //activate free fluid interface model
	 +  fvc::interpolate(PcCoeff2)*fvc::snGrad(alpha1)*Bsurf //activate porous media model
	 -  ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
          
	    fvm::laplacian(rAUf, p_rgh) - fvc::ddt(epsf) == fvc::div(phiHbyA) //There should not be an alpha here     
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - p_rghEqn.flux();

            p_rgh.relax();

            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}

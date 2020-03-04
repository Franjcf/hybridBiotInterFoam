//Fluid Momentum equation
//MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn //setting the DB Equation
(
	//InertialTerms
        (1/epsf)*(fvm::ddt(U) + fvm::div(phiByEpsf, U))

	//Viscous Terms
      - (1/epsf)*(fvm::laplacian(nuf,U))
      - (1/epsf)*(fvc::div(nuf*dev2(T(fvc::grad(U)))))
	
	//Implicit Darcy-Brinkman Drag Term
      + fvm::Sp((nuf*kinv),U) 

	//Explicit Momentum Exchange Terms
	//Due to Relative Fluid Motion and Volume Averaging 
      - fvc::Sp(epss*nuf*kinv,Us) 

	//gravity
      - g

);

fvVectorMatrix& UEqn = tUEqn.ref();

phiByEpsf=fvc::interpolate(1/epsf)*phi;

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}

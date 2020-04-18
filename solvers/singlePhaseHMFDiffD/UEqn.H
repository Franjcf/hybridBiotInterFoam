//Fluid Momentum Equation
tmp<fvVectorMatrix> tUEqn
(
	//InertialTerms
        (1/eps)*(fvm::ddt(U) + fvm::div(phiByEps, U))

	//Viscous Terms
      - (1/eps)*(fvm::laplacian(nuf,U))
      - (1/eps)*(fvc::div(nuf*dev2(T(fvc::grad(U)))))
	
	//Implicit Darcy-Brinkman Drag Term
      + fvm::Sp((nuf*kinv),U) 

	//Explicit Momentum Exchange Terms
      - fvc::Sp(epss*nuf*kinv,Us) 

	//gravity
      - g

	//2.5 D term from A. Ferrari 2015
      + fvm::Sp((nuf*twoDkinv),U)

);

fvVectorMatrix& UEqn = tUEqn.ref();

phiByEps=fvc::interpolate(1/max(eps,SMALL))*phi;

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}

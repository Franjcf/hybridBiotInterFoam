//Fluid Momentum Equation

tmp<fvVectorMatrix> tUEqn 
(
	//InertialTerms
        (1/epsf)*(fvm::ddt(U) + fvm::div(phiByEpsf, U))

	//Viscous Terms
      - (1/epsf)*(fvm::laplacian(nuf,U))
      - (1/epsf)*(fvc::div(nuf*dev2(T(fvc::grad(U)))))
	
	//Implicit Drag Term
      + fvm::Sp((nuf*kinv),U) 

	//Explicit Drag Term 
      - fvc::Sp(epss*nuf*kinv,Us) //epss*epsf

	//2.5 D term from A. Ferrari 2015
      + fvm::Sp((nuf*Fkinv),U)

);

fvVectorMatrix& UEqn = tUEqn.ref();

phiByEpsf=fvc::interpolate(1/max(epsf,SMALL))*phi;

UEqn.relax();

//fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

//    fvOptions.correct(U);
}

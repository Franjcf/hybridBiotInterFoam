/////////////////////////Calculating Up////////////

Info<< "\nCalculating Up" << endl;

//Calculating Explicit Momentum Sources
SourceTransport= epsf/rhoc*kinvMu*(U-Up*epss);
SourceSwelling= -(sig/rhoc)*epss*fvc::grad(Ps);
gradPc  =  (mixture.sigmaK()*(1-B) + PcCoeff*B)*fvc::grad(alpha1)/rhoc; 

fvVectorMatrix UpEqn
	(

	//Inertial Terms
	fvm::ddt(epss,Up) 

	+ fvm::div(phipEpss,Up)
	
	//Viscouus terms: div(mu grad(U*epss)) 
      - fvm::laplacian(fluid.nu(),Up)

	//U Transpose viscocity term 
      - fvc::div(fluid.nu()*dev2(T(fvc::grad(Up))))

	//Bulk viscocity term if neccesary
      //- fvc::grad(zeta*epss*fvc::div(Up))

	//Gravity Term
     - epss*g //epss is multiplied 

	//SourceSwelling
      + SourceSwelling 

	//Capillary presure source
      + gradPc*epsf //This is definetily positive
   
	//Clay Transport Momentum Source
      - fvc::Sp(epsf*kinvMu/rhoc, U) 
	
	//Momentum Transfer to fluid do to relative velocity 
      + fvm::Sp(epss*epsf*kinvMu/rhoc, Up) 
);	

//To do capillary and Buckley-leverett solutions juts comment these out.
UpEqn.solve();
Up.correctBoundaryConditions();
					 
//creating phi equivalent flux in order to calculate divergence (phip,epss)
phip = fvc::interpolate(Up) & mesh.Sf();
phipEpss=fvc::interpolate(epss)*phip; 
phic= -phipEpss;


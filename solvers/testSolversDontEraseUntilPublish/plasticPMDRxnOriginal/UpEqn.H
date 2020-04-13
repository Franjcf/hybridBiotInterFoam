/////////////////////////Calculating Up////////////

Info<< "\nCalculating Up" << endl;

SourceTransport= rhof*epsf/rhoc*nuf*kinv*(U-Up*epss);
SourceSwelling= -(sigma/rhoc)*epss*fvc::grad(Ps);

fvVectorMatrix UpEqn
	(

	//Inertial Terms
	fvm::ddt(epss,Up) 

	+ fvm::div(phipEpss,Up)
	
	//Viscouus terms: div(mu grad(U*epss)) 
      - fvm::laplacian(fluid.nu(),Up)

	//U Transpose viscocity term 
      - fvc::div(fluid.nu()*dev2(T(fvc::grad(Up))))

	//Bulk viscocity term //may add or not
      //- fvc::grad(zeta*epss*fvc::div(Up))

	//Gravity Term
     - epss*g 

	//SourceSwelling
      + SourceSwelling 
   
	//Clay Transport Momentum Source
      - fvc::Sp(epss*(rhof*nuf*kinv)/rhoc, U) //epsf //maybe not even

      //Momentum Transfer to fluid do to relative velocity 
      + fvm::Sp(epss*(rhof*nuf*kinv)/rhoc, Up)  //epsf
);	

UpEqn.solve();
Up.correctBoundaryConditions();
					 
//creating phi equivalent flux in order to calculate divergence (phip,epss)

phip = fvc::interpolate(Up) & mesh.Sf();
phipEpss=fvc::interpolate(epss)*phip; 
phic= -phipEpss;


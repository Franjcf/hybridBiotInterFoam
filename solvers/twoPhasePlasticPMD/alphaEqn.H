{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    tmp<fv::ddtScheme<scalar>> ddtAlpha
    (
        fv::ddtScheme<scalar>::New
        (
            mesh,
            mesh.ddtScheme("ddt(alpha)")
        )
    );

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    if
    (
        isType<fv::EulerDdtScheme<scalar>>(ddtAlpha())
     || isType<fv::localEulerDdtScheme<scalar>>(ddtAlpha())
    )
    {
        ocCoeff = 0;
    }
    else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha()))
    {
        if (nAlphaSubCycles > 1)
        {
            FatalErrorInFunction
                << "Sub-cycling is not supported "
                   "with the CrankNicolson ddt scheme"
                << exit(FatalError);
        }

        ocCoeff =
            refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtAlpha())
           .ocCoeff();
    }
    else
    {
        FatalErrorInFunction
            << "Only Euler and CrankNicolson ddt schemes are supported"
            << exit(FatalError);
    }

    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf())); //defenition of compresion coeff. 

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

	Ur = -B*(k1rel/(max(alpha1,0.01)*mu1) - k2rel/(max((1-alpha1),0.01)*mu2) )*k0real*fvc::grad(p);

	//zeroing Ur corretion at the boundary (to prevent crazy instabilities)
	Bborder = fvc::grad(mag(fvc::grad(B)));
  
	Ur.correctBoundaryConditions();	

        surfaceScalarField phirPore(fvc::interpolate(Ur) & mesh.Sf());

        surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();

        surfaceScalarField::Boundary& phirPoreBf =
        phirPore.boundaryFieldRef(); 

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    forAll(phirPore.boundaryField(), patchi) 
    {
        fvsPatchScalarField& phirPorep = phirPoreBf[patchi];

        if (!phirPorep.coupled())
        {
            phirPorep == 0; 
        }
    }

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime();
    }

    if (MULESCorr)
    {
        fvScalarMatrix alpha1Eqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
              : fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phiCN,
                upwind<scalar>(mesh, phiCN)
            ).fvmDiv(phiCN, alpha1)
        );

        alpha1Eqn.solve();

        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        tmp<surfaceScalarField> talphaPhiUD(alpha1Eqn.flux());
        alphaPhi = talphaPhiUD();

        if (alphaApplyPrevCorr && talphaPhiCorr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            //MULES::correct(alpha1, alphaPhi, talphaPhiCorr0.ref(), 1, 0); 
	    MULES::correct(epsf ,alpha1, alphaPhi, talphaPhiCorr0.ref(), zeroField(),zeroField(), 1, 0); //**

            alphaPhi += talphaPhiCorr0();
        }

        // Cache the upwind-flux
        talphaPhiCorr0 = talphaPhiUD;

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }


    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
	surfaceScalarField Bsurf(fvc::interpolate(B)); //identifying solid space 
	//this loop is neccesary to completely segregate both models. 
	//the face at the interface has to be controled by the free fluid model
	//otherwise it will not work

        forAll (Bsurf, facei)
	{
		if (Bsurf[facei] < 1) //if free fluid
		{Bsurf[facei] = 0;}
	}
  
	surfaceScalarField phir(phic*mixture.nHatf()); //compressive velocity 

	if(SolveWithUr.value()==1) 
	{phir = (phirPore*Bsurf+phic*mixture.nHatf()*(1-Bsurf));}

        tmp<surfaceScalarField> talphaPhiUn
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, (alpha2*epsf), alpharScheme), 
		alpha1,
                alpharScheme
            )
        );

        // Calculate the Crank-Nicolson off-centred alpha flux
        if (ocCoeff > 0)
        {
            talphaPhiUn =
                cnCoeff*talphaPhiUn + (1.0 - cnCoeff)*alphaPhi.oldTime();
        }

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhiCorr(talphaPhiUn() - alphaPhi);
            volScalarField alpha10("alpha10", alpha1);

	    //Here we added epsf to the conservation eqn. 
	    //MULES::correct(alpha1, talphaPhiUn(), talphaPhiCorr.ref(),1, 0); 
	    MULES::correct(epsf ,alpha1, talphaPhiUn(), talphaPhiCorr.ref(), zeroField(),zeroField(), 1, 0); 

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi += talphaPhiCorr();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                alphaPhi += 0.5*talphaPhiCorr();
            }
        }
        else
        {
            alphaPhi = talphaPhiUn;

            //Here we added epsf to the conservation eqn
	    //MULES::explicitSolve(alpha1,phiCN,alphaPhi,1,0); 
	    MULES::explicitSolve(epsf,alpha1,phiCN,alphaPhi,zeroField(),zeroField(),1,0); 

        }

        alpha2 = 1.0 - alpha1;

        mixture.correct();
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        talphaPhiCorr0 = alphaPhi - talphaPhiCorr0;
    }

    if
    (
        word(mesh.ddtScheme("ddt(rho,U)"))
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        rhoPhi = alphaPhi*(rho1 - rho2) + phiCN*rho2;
    }
    else
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step alpha flux
            alphaPhi = (alphaPhi - (1.0 - cnCoeff)*alphaPhi.oldTime())/cnCoeff;
        }

        // Calculate the end-of-time-step mass flux
        rhoPhi = alphaPhi*(rho1 - rho2) + phi*rho2;
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}

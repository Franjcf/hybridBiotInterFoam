Info<< "\nCalculating Solid Displacement Field" << endl;

//Reading nCorrectors and stress type (thermal and/or plane or neither)
#include "readSolidDisplacementFoamControls.H"

//Initializing number of corrections and residuals
int iCorr = 0; 
scalar initialResidual = 0;

//Calculating Source Terms
volVectorField Terzaghi = -(epss/rho_s)*fvc::grad(p);
volVectorField CapTerm (eps/rho_s*PcCoeff*Solid*fvc::grad(alpha1)); 

//Calculating Displacement

        do
        {
            if (thermalStress)
            {
                volScalarField& T = Tptr();
                solve
                (
                    fvm::ddt(T) == fvm::laplacian(DT, T)
                );
            }

            { //This is always done
                fvVectorMatrix DEqn 
                (
		    //transient term, D is displacement vector
                    fvm::d2dt2(epss,D) 
                 ==
		    //implicit part of the stress tensor
                    fvm::laplacian((2*mu + lambda)*epss, D, "laplacian(DD,D)")

		    //Terzaghi Term
		    + Terzaghi

		    //Gravity Term
		    + epss*g

                    //Transport Momentum Source 
                    + fvc::Sp(eps*Drag/rho_s, U)
                    - fvc::Sp(epss*eps*Drag/rho_s, Us)

                    //Capillary Pressure Source
                    - CapTerm                   

                    //Explicit Part of Stress Tensor
                    + divSigmaExp 
                );


                if (thermalStress)
                {
                    const volScalarField& T = Tptr();
                    DEqn += fvc::grad(threeKalpha*T);
                }

                //DEqn.setComponentReference(1, 0, vector::X, 0);
                //DEqn.setComponentReference(1, 0, vector::Z, 0);

                initialResidual = DEqn.solve().max().initialResidual();

                if (!compactNormalStress)
                {
                    divSigmaExp = fvc::div(DEqn.flux()); 
                }
            }

            {//Calculating Explicit Stress Tensor

                volTensorField gradD(fvc::grad(D));
                sigmaD = (epss*mu*twoSymm(gradD) + epss*(lambda*I)*tr(gradD));

                if (compactNormalStress) 
                {
                    divSigmaExp = fvc::div
                    (
                        sigmaD - ((2*mu + lambda)*epss*gradD),
                        "div(sigmaD)"
                    );
                }
                else 
                {
                    divSigmaExp += fvc::div(sigmaD);
                }
            }

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

	//here we reconstruct the final stress tensor from everything to print
        #include "calculateStress.H"

// ************************************************************************* //

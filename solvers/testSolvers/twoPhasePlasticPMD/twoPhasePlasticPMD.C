/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

    For a two-fluid approach see twoPhaseEulerFoam.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "correctPhi.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //Counts total initial solid volume fraction
	//#include "totalEpssInitial.H"

    while (runTime.run())
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "courantNoUp.H"
            #include "alphaCourantNo.H" //just reports the number
            #include "setDeltaTUp.H" //modifies based on solid speed
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
///////////////////CLAY MOMEVMENT SECTION//////////////

if(runTime.value() > runTime.startTime().value() + runTime.deltaTValue())
//Only runs these after the first time step, this way we calculate
//Uf first. Uf is also calculated last based on the last movement 
//That is why Uf follows after this loop
{

	 #include "swellingPressureModel.H" //Swelling pressure model	

         fluid.correct(); //Upating Clay Non-Newtoninan Viscocity

         #include "UpEqn.H" //Porous Media Velocity

         #include "epssEqn.H" //Porous Media Conservation Equation

	 //#include "updateEpssVolConservation.H" //Check Mass conservation
}

/////////////////////FLUID VELOCITY SECTION//////////////////////////////////////////

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
	   
            #include "alphaControls.H"
            #include "alphaEqnSubCycle.H"

            mixture.correct(); //Updating VOF viscositites and densities
	
	    #include "updateVariables.H"
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        if (SolveConvectionEqn.value() ==1)
	{
	#include "convectionEqn.H" //Convection of Ions
	}

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;


       // #include "writeEpssVolConservation.H" //Write solid mass conservation

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

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
    pimpleFoam

Description
    Large time-step transient solver for incompressible, turbulent flow, using
    the PIMPLE (merged PISO-SIMPLE) algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createControlsSolid.H"
    #include "createFields.H"
    #include "createFieldsSolid.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

	//Counts total initial solid volume fraction
	#include "totalEpssInitial.H" 

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
        
///////////////////SOLID MECHANICS SECTION//////////////

if(runTime.value() > runTime.startTime().value() + runTime.deltaTValue())
//Only runs these after the first time step, this way we calculate
//Uf first. Uf is also calculated last based on the last movement 
//That is why Uf follows after this loop
{

	 #include "swellingPressureModel.H" //Swelling pressure model	
	 #include "solidMechanics.H"
	
	 #include "updateEpssVolConservation.H" //Check Mass conservation
	 #include "updateVariables.H" //Update epsf,k, ext...

}

/////////////////////FLUID VELOCITY SECTION//////////////////////////////////////////

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {   
    
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
                  
            if (pimple.turbCorr())
            {
          
            }
        }
             

	if (SolveConvectionEqn.value() ==1)
	{
	#include "convectionEqn.H" //Convection of Ions
	}
      
        runTime.write();

        Info<< "\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

         #include "writeEpssVolConservation.H" //Write solid mass conservation

    }

        Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

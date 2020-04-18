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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
	#include "courantNoUs.H"
        #include "setDeltaTUs.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        
///////////////////Solid Mechanics Section//////////////
if(runTime.value() > runTime.startTime().value() + runTime.deltaTValue())
{

         solid.correct(); //Non-Newtoninan Viscocity Update

         #include "UsEqn.H" //Porous Media Momentum Conservation Eqn. 

         #include "epssEqn.H" //Porous Media Continuity Eqn.

	 #include "updateVariables.H" //Update eps,k, ext...

}

/////////////////////Fluid Mechanics Section//////////////////////////////////////////

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {   
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
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
    }

        Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

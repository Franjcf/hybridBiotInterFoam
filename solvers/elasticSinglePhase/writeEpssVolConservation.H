//Writing Mass Conservation File

//Writing Epss change
	OFstream str(runTime.path()/"MassConservation.dat"); 
        str<< runTime.timeOutputValue()
        <<"\n\ninitialEpss: "   << initialEpss.value()
        <<"   finalEpss: "     << finalEpss.value()
        <<"   Total percentChange: " << PercentEpssChange.value()
        <<"   Local percentChange: " << localEpssChange.value()

//Writing Volume change
        <<"\n\ninitialVol: "     << initialVol.value() 
	<<"   finalVol: "      << finalVol.value()
	<<"   Total Vol percent Change: " << PercentVolChange.value()
        <<"\nredCount: " << redCount.value()

//Writing Time

	<< "\n\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

       


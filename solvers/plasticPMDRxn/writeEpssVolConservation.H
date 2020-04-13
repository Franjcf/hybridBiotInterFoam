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

//Writing Concentration change
        //<<"   C percent Change: " << Csum.value()/Csum0.value()
       // <<"   C*epsf percent Change: " << CepsfSum.value()/CepsfSum0.value()

//Writing Time

	<< "\n\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

       


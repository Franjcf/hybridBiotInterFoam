//CALCULATION OF INITIAL TOTAL POROSITY AND VOLUME

// initialEpss is the total amound of clay epss and initialVol is the volume fo the clay
initialEpss=0;
initialVol=0;

forAll(epssInitial,celli)
{
	initialEpss += epssInitial[celli];

	if (epss[celli]>0.001)
	{
		initialVol += 1;
	}

}


					 

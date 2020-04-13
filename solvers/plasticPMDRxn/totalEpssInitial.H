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

dimensionedScalar Csum0("Csum0", dimensionSet(0,0,0,0,0,0,0), 0);
dimensionedScalar CepsfSum0("CepsfSum0", dimensionSet(0,0,0,0,0,0,0), 0);

forAll(epss,celli)
{
        Csum0 += C[celli];
        CepsfSum0 += (epsf[celli]*C[celli]);

}



					 

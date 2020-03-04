//CONSERVATION OF EPSS AND VOLUME FILE

//Initializing the total changes
finalEpss=0;
finalVol=0;

//Initializing the local changes
Epss0=0;      //Epss0 is the current epss of the places that have expanded or contracted
EpssIn=0;     //EpssIn is the initial epss of clay in the places that have expanded or contracted

//////////////////////////////////////////

forAll(epss,celli) 
{

//labeling the places where porosity has changed
	if(mag(epss[celli]-epssInitial[celli])>0.005 && epssChange[celli]==0) 
		{
		epssChange[celli] = 1;
		}

//making sure values are less than one and more than zero
	if(epss[celli]>0.999)
		{
			epss[celli]=0.999;
			redCount+=1;
		}
	
//Just in case
	if (epss[celli]<0)
		{
			epss[celli]=1e-20;
		}

//ADDING UP THE TOTAL POROSITY AND VOLUME//////////////////////////////////////

	finalEpss += epss[celli]; //total epss in the total epss in clay

	if (epss[celli]>0.001)
	{
		finalVol  += 1; //finalVol is the total clay volume
	}
	                   
//CALCULATING THE LOCAL CHANGES IN POROSITY AND VOLUME///////////////////////
	
	//if porosity has changed, then add up the values here. 
	if(epssChange[celli]==1)
		{

		EpssIn += epssInitial[celli];
		Epss0 += epss[celli];		
		}			
}
		
//CALCULATING THE CHANGES/////////////////////////////////////////////////	

//Total Changes	 
PercentEpssChange = (finalEpss/(initialEpss+1e-20)-1)*100;
PercentVolChange = (finalVol/(initialVol+1e-20)-1)*100;

//Local changes 
localEpssChange=(Epss0/(EpssIn+1e-20)-1)*100;

//WRITING CHANGES IN THE TERMINAL////////////////////////////////////////

Info<< "\nTotal Percent Epss Change  = " << PercentEpssChange.value() << endl;
Info<< "Local Percent Epss Change  = " << localEpssChange.value() << endl;
Info<< "Total Percent Vol Change  = " << PercentVolChange.value() << endl;
Info<< "redCount  = " << redCount.value() << endl;

Info<< "\n" << endl;


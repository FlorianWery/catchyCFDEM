fileName resultsDir = mesh.time().path()/"postProcessing"/runTime.timeName();
if (!isDir(resultsDir))
{
	mkDir(resultsDir);
}

word resultsName=(word("circumAverage-")&word(std::to_string(planeZ)))&word(".csv");
OFstream resultsFile(resultsDir/resultsName);  
    
if (resultsFile.good())
{
	resultsFile <<  "# Circumferential averages at plane z = " << planeZ << endl;
	resultsFile <<  "radPos" 
                << tab << IOobject::groupName("magU",phaseName1)
                << tab << IOobject::groupName("radU",phaseName1)
                << tab << IOobject::groupName("tanU",phaseName1)
               /* << tab << IOobject::groupName("magU",phaseName2) 
                << tab << IOobject::groupName("radU",phaseName2) 
                << tab << IOobject::groupName("tanU",phaseName2)*/;
	forAll(fieldIDs, fi)
	{
		resultsFile << tab << fieldIDs[fi];
	}
	resultsFile << endl;
}
else
{
	Info<< "Cannot make results file" << endl;
};



fileName resultsDir = mesh.time().path()/"postProcessing"/"mixingCup"/regionName/runTime.timeName();
if (regionName == "region0") resultsDir = mesh.time().path()/"postProcessing"/"mixingCup"/runTime.timeName();
    
if (!isDir(resultsDir))
{
	mkDir(resultsDir);
}
 
OFstream resultsFile(resultsDir/"mixingCupResults.csv");  
    
if (resultsFile.good())
{
	Info << "Initialising results file" << endl;
	resultsFile <<  "distance" << "," << "weight";
	forAll(fieldIDs, fi)
	{
		resultsFile << "," << fieldIDs[fi];
	}
	resultsFile << endl;
}
else
{
	Info<< "Cannot make results file" << endl;
};



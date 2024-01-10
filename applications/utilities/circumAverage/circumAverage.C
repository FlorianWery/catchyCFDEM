/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Description
    Write the three components of the cell centres as volScalarFields so
    they can be used in postprocessing in thresholding.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of volScalarFields for averaging (required)"
    );
    argList::addOption
    (
		"nPoints",
		"number",
		"number of data intervals for averaging (default = 100)"
    );
    argList::addOption
    (
		"direction",
		"number",
		"specify the axis direction (x = 0, y = 1, z = 2, default: z)"
    );
    argList::addOption
    (
		"planeZ",
		"number",
		"specify the coordinate of the plane on which to average fields (required)"
    );
    argList::addBoolOption
    (
        "meanValue",
        "use mean velocity field value"
    );
    argList::addOption
    (
		"deltaZ",
		"number",
		"resolution for plane coordinate (default = 0.001)"
    );
    argList::addOption
    (
		"phaseName1",
		"word",
		"optional phase name for velocity field"
    );
    argList::addOption
    (
		"phaseName2",
		"word",
		"optional phase name for velocity field"
    );
    timeSelector::addOptions();

#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createNamedMesh.H"

    List<word> fieldIDs;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> fieldIDs;
    }
    if (fieldIDs.size())
    {
        Info<< "Calculating circumferential averages for fields " << fieldIDs << nl << endl;
    }
    else
    {
        FatalErrorIn("circumAverage::args")
            << " Please specify a list of fields" << exit(FatalError);
    }

    scalar nPoints;
    scalar readNPoints;
    if (args.optionReadIfPresent("nPoints", readNPoints)) {nPoints=readNPoints;}
    else {nPoints=100;}

    scalar planeZ;
    scalar readPlaneZ;
    if (args.optionReadIfPresent("planeZ", readPlaneZ)) {planeZ=readPlaneZ;}
    else 
    {
        FatalErrorIn("circumAverage::args")
            << " Please specify a the plane coordinate 'planeZ'" << exit(FatalError);
    }
   
    scalar deltaZ;
    scalar readDeltaZ;
    if (args.optionReadIfPresent("deltaZ", readDeltaZ)) {deltaZ=readDeltaZ;}
    else {deltaZ=0.001;}
    
    label axDir;
    label readDir;
    if (args.optionReadIfPresent("direction", readDir)) {axDir=readDir;}
    else {axDir=2;}
    
    word phaseName1;
    word readPhaseName1;
    if (args.optionReadIfPresent("phaseName1", readPhaseName1)) {phaseName1=readPhaseName1;}
    else {phaseName1=word::null;}
    
/*    word phaseName2;
    word readPhaseName2;
    if (args.optionReadIfPresent("phaseName2", readPhaseName2)) {phaseName2=readPhaseName2;}
    else {phaseName2=word::null;}*/

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        // Check for new mesh
        mesh.readUpdate();

		#include "createResultsFile.H"
		#include "createFields.H"

		Info << "Initializing averages" << endl;
		scalarList fieldAverages(fieldIDs.size()+3);
		scalar totWeight;
		forAll(fieldAverages, fa)
		{
			fieldAverages[fa]=0.0;
		}
		totWeight=0.0;	

		Info << "Calculating averages" << endl;
		scalar rStart=minR;
		label j =1;
        fieldAverages[0]=0.0+magU1[orderR[0]]*cv[orderR[0]];
        fieldAverages[1]=0.0+radU1[orderR[0]]*cv[orderR[0]];
        fieldAverages[2]=0.0+tanU1[orderR[0]]*cv[orderR[0]];
        /*fieldAverages[3]=0.0+magU2[orderR[0]]*cv[orderR[0]];
        fieldAverages[4]=0.0+radU2[orderR[0]]*cv[orderR[0]];
        fieldAverages[5]=0.0+tanU2[orderR[0]]*cv[orderR[0]];*/
		forAll(fieldIDs, fi)
		{
			fieldAverages[fi+3]=0.0+iFields[fi][orderR[0]]*cv[orderR[0]];
		}
		totWeight=0.0+cv[orderR[0]];
		forAll(orderR, i)
		{
			if (radPosZ[orderR[i]]<minR+j*deltaR)
			{
                fieldAverages[0]+=magU1[orderR[i]]*cv[orderR[i]];
                fieldAverages[1]+=radU1[orderR[i]]*cv[orderR[i]];
                fieldAverages[2]+=tanU1[orderR[i]]*cv[orderR[i]];
                /*fieldAverages[3]+=magU2[orderR[i]]*cv[orderR[i]];
                fieldAverages[4]+=radU2[orderR[i]]*cv[orderR[i]];
                fieldAverages[5]+=tanU2[orderR[i]]*cv[orderR[i]];*/
				forAll(fieldIDs, fi)
				{    
					fieldAverages[fi+3] += iFields[fi][orderR[i]]*cv[orderR[i]];
				}
				totWeight += cv[orderR[i]];
			}
			else
			{
				resultsFile << (radPosZ[orderR[i-1]]+rStart)/2.0;
				forAll(fieldAverages, fa)		
				{
					resultsFile << tab << fieldAverages[fa]/totWeight;
				}
				resultsFile << endl;

				rStart=radPosZ[orderR[i-1]];
				j=j+1;
                fieldAverages[0]=0.0+magU1[orderR[i]]*cv[orderR[i]];
                fieldAverages[1]=0.0+radU1[orderR[i]]*cv[orderR[i]];
                fieldAverages[2]=0.0+tanU1[orderR[i]]*cv[orderR[i]];
                /*fieldAverages[3]=0.0+magU2[orderR[i]]*cv[orderR[i]];
                fieldAverages[4]=0.0+radU2[orderR[i]]*cv[orderR[i]];
                fieldAverages[5]=0.0+tanU2[orderR[i]]*cv[orderR[i]];*/
				forAll(fieldIDs, fi)
				{    
					fieldAverages[fi+3]=0.0+iFields[fi][orderR[i]]*cv[orderR[i]];
				}
				totWeight=0.0+cv[orderR[i]];
			}
		}
		/*resultsFile << (maxR+rStart)/2.0 << tab << totWeight;
		forAll(fieldAverages, fa)		
		{
			resultsFile << tab << fieldAverages[fa]/totWeight;
		}
		resultsFile << endl;*/
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

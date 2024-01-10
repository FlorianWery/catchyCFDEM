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
#include "vectorIOField.H"
#include "volFields.H"
#include "fvCFD.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of volScalarFields for mixing-cup averaging (required)"
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
		"specify the axial direction (x = 0, y = 1, z = 2, default: z)"
    );
    argList::addOption
    (
		"inletPatch",
		"word",
		"name of inlet patch (default: inlet)"
    );
    argList::addOption
    (
		"outletPatch",
		"word",
		"name of outlet patch (default: outlet)"
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
        Info<< "Calculating mixing-cup averages for fields " << fieldIDs << nl << endl;
    }
    else
    {
        FatalErrorIn("calcMixingCup::args")
            << " Please specify a list of fields" << exit(FatalError);
    }

    word inName;
    word readInName;
    if (args.optionReadIfPresent("inletPatch", readInName)) {inName=readInName;}
    else {inName="inlet";}
	label inPatch=mesh.boundaryMesh().findPatchID(inName);

    word outName;
    word readOutName;
    if (args.optionReadIfPresent("outletPatch", readOutName)) {outName=readOutName;}
    else {outName="outlet";}
	label outPatch=mesh.boundaryMesh().findPatchID(outName);

    scalar nPoints;
    scalar readNPoints;
    if (args.optionReadIfPresent("nPoints", readNPoints)) {nPoints=readNPoints;}
    else {nPoints=100;}

    label axDir;
    label readDir;
    if (args.optionReadIfPresent("direction", readDir)) {axDir=readDir;}
    else {axDir=2;}

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        // Check for new mesh
        mesh.readUpdate();

		#include "createResultsFile.H"
		#include "createFields.H"

		Info << "Initializing averages" << endl;
		scalarList fieldAverages(fieldIDs.size());
		scalar totWeight;
		forAll(fieldIDs, fi)
		{
			fieldAverages[fi]=0.0;
		}
		totWeight=0.0;

		Info << "Calculating inlet boundary field averages" << endl;
		forAll(biZ, faceI)
		{
			forAll(fieldIDs, fi)
			{    
				fieldAverages[fi] += biFields[fi][faceI]*phi.boundaryField()[inPatch][faceI];
			}
			totWeight += phi.boundaryField()[inPatch][faceI];
		}
		resultsFile <<  inZ << "," << totWeight;
		forAll(fieldIDs, fi)		
		{
			resultsFile << "," << fieldAverages[fi]/totWeight;
		}
		resultsFile << endl;
	

		Info << "Calculating internal field averages" << endl;
		scalar zStart=inZ;
		label orderi = deltaZ>0 ? order[0] : order[order.size()];
		label j =1;
		forAll(fieldIDs, fi)
		{
			fieldAverages[fi]=0.0+iFields[fi][orderi]*phi.internalField()[orderi];
		}
		totWeight=0.0+phi.internalField()[orderi];
		for (label i=1; i<order.size(); i++)
		{
		    orderi = deltaZ>0 ? order[i] : order[order.size()-i];
			if ( deltaZ>0 ? cZ[orderi]<inZ+j*deltaZ : cZ[orderi]>inZ+j*deltaZ)
			{
				forAll(fieldIDs, fi)
				{    
					fieldAverages[fi] += iFields[fi][orderi]*phi.internalField()[orderi];
				}
				totWeight += phi.internalField()[orderi];
			}
			else
			{
				resultsFile << (cZ[orderi]+zStart)/2.0 << "," << totWeight;
				forAll(fieldIDs, fi)		
				{
					resultsFile << "," << fieldAverages[fi]/totWeight;
				}
				resultsFile << endl;

				zStart=cZ[orderi];
		        orderi = deltaZ>0 ? order[i] : order[order.size()-i];
				j=j+1;
				forAll(fieldIDs, fi)
				{    
					fieldAverages[fi]=0.0+iFields[fi][orderi]*phi.internalField()[orderi];
				}
				totWeight=0.0+phi.internalField()[orderi];
			}
		}
		resultsFile << (outZ+zStart)/2.0 << "," << totWeight;
		forAll(fieldIDs, fi)		
		{
			resultsFile << "," << fieldAverages[fi]/totWeight;
		}
		resultsFile << endl;

		Info << "Calculating outlet boundary field averages" << endl;
		forAll(boZ, faceI)
		{
			forAll(fieldIDs, fi)
			{    
				fieldAverages[fi] += boFields[fi][faceI]*phi.boundaryField()[outPatch][faceI];
			}
			totWeight += phi.boundaryField()[outPatch][faceI];
		}
		resultsFile <<  outZ << "," << totWeight;
		forAll(fieldIDs, fi)		
		{
			resultsFile << "," << fieldAverages[fi]/totWeight;
		}
		resultsFile << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

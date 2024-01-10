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
    


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
		"UName",
		"word",
		"optional name for velocity field"
    );
    argList::addOption
    (
		"UPrimeName",
		"word",
		"optional name for velocity variance tensorfield"
    );
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    
    label axDir;
    label readDir;
    if (args.optionReadIfPresent("direction", readDir)) {axDir=readDir;}
    else {axDir=2;}
    
    bool writeU = false;
    word Uname;
    word readUname;
    if (args.optionReadIfPresent("UName", readUname)) {Uname=readUname; writeU=true;}
    else {Uname=word::null;}
    
    bool writeUPrime = false;
    word UPrimename;
    word readUPrimename;
    if (args.optionReadIfPresent("UPrimeName", readUPrimename)) {UPrimename=readUPrimename; writeUPrime=true;}
    else {UPrimename=word::null;}

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "\nTime = " << runTime.timeName() << endl;

        // Check for new mesh
        mesh.readUpdate();
        #include "createFields.H"
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

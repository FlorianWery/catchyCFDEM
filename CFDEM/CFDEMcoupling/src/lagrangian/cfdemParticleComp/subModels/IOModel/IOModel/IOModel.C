/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "IOModel.H"
#include "particle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(IOModel, 0);

defineRunTimeSelectionTable(IOModel, dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

int IOModel::dumpDEMdata() const
{
    return -1;
}

bool IOModel::dumpNow() const
{
    //bool dmp(false);
    //if (time_.value()+SMALL > time_.endTime().value()-time_.deltaT().value() || time_.outputTime())
    //    dmp=true;

    return time_.outputTime();
}

fileName IOModel::createTimeDir(fileName path) const
{
    fileName timeDirPath(path/time_.timeName());
    mkDir(timeDirPath,0777);
    return timeDirPath;
}

fileName IOModel::createLagrangianDir(fileName path) const
{
    fileName lagrangianDirPath(path/"lagrangian");
    mkDir(lagrangianDirPath,0777);
    fileName cfdemCloudDirPath(lagrangianDirPath/"cfdemCloud1");
    mkDir(cfdemCloudDirPath,0777);
    return cfdemCloudDirPath;
}

fileName IOModel::buildFilePath(word dirName) const
{
    // create file structure
	fileName path("");
    if(parOutput_)
    {
    	path=fileName(particleCloud_.mesh().time().path()/particleCloud_.mesh().time().timeName()/dirName/"particleCloud");
    	mkDir(path,0777);
    } else
    {
		path=fileName("."/dirName);
    	mkDir(path,0777);
    	mkDir(fileName(path/"constant"),0777);
    	OFstream* stubFile = new OFstream(fileName(path/"particles.foam"));
    	delete stubFile;
        }
    return path;
}

void IOModel::streamDataToPath(fileName path, double** array,int nPProc,word name,word type,word className,word finaliser) const
{
    OFstream fileStream(fileName(path/name));

    fileStream
        << "/*--------------------------------*- C++ -*----------------------------------*\\" << nl
        << "  =========                 |" << nl
        << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox" << nl
        << "   \\\\    /   O peration     | Website:  https://openfoam.org" << nl
        << "    \\\\  /    A nd           | Version:  8" << nl
        << "     \\\\/     M anipulation  |" << nl
        << "\\*---------------------------------------------------------------------------*/" << nl
        << "FoamFile" << nl
        << "{" << nl
        << "    version     " << fileStream.version() << ";" << nl
        << "    format      " << fileStream.format() << ";" << nl
        << "    class       " << className << ";" << nl
        << "    location    " << path << ";" << nl
        << "    object      " << name << ";" << nl
        << "}" << nl
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << nl << nl;
    
    fileStream << nPProc <<"\n";

    if (type == "origId")
    {
        if (nPProc > 0) fileStream << "{0}" << "\n";
        else fileStream << "{}" << "\n";
        return;
    }

    if (type == "origProcId")
    {
        if (nPProc > 0) fileStream << "{0}" << "\n";
        else fileStream << "{}" << "\n";
        return;
    }

    fileStream << "(\n";

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        if (particleCloud_.cellIDs()[index][0] > -1) // particle Found
        {
            if (type=="scalar"){
                fileStream << array[index][0] << " \n";
            }
            else if (type == "position")
            {
                particle ofparticle
                (
                    particleCloud_.mesh(),
                    vector(array[index][0],array[index][1],array[index][2]),
                    particleCloud_.cellIDs()[index][0]
                );
                ofparticle.writePosition(fileStream);
                fileStream << nl;
            }
            else if (type == "vector")
            {
                fileStream << "( "<< array[index][0] << " " << array[index][1] << " " << array[index][2] << " ) " << " \n";
            }
        }
    }

    fileStream << ")\n";
}

void IOModel::streamDataToPath(fileName path, int** array,int nPProc,word name,word type,word className,word finaliser) const
{
    OFstream fileStream(fileName(path/name));

    fileStream
        << "/*--------------------------------*- C++ -*----------------------------------*\\" << nl
        << "  =========                 |" << nl
        << "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox" << nl
        << "   \\\\    /   O peration     | Website:  https://openfoam.org" << nl
        << "    \\\\  /    A nd           | Version:  8" << nl
        << "     \\\\/     M anipulation  |" << nl
        << "\\*---------------------------------------------------------------------------*/" << nl
        << "FoamFile" << nl
        << "{" << nl
        << "    version     " << fileStream.version() << ";" << nl
        << "    format      " << fileStream.format() << ";" << nl
        << "    class       " << className << ";" << nl
        << "    location    " << path << ";" << nl
        << "    object      " << name << ";" << nl
        << "}" << nl
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << nl << nl;
    
    fileStream << nPProc <<"\n";

    if (type == "origId")
    {
        if (nPProc > 0) fileStream << "{0}" << "\n";
        else fileStream << "{}" << "\n";
        return;
    }

    if (type == "origProcId")
    {
        if (nPProc > 0) fileStream << "{0}" << "\n";
        else fileStream << "{}" << "\n";
        return;
    }

    fileStream << "(\n";

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
        if (particleCloud_.cellIDs()[index][0] > -1) // particle Found
        {
            if (type=="scalar"){
                fileStream << array[index][0] << " \n";
            }
            else if (type == "position")
            {
                particle ofparticle
                (
                    particleCloud_.mesh(),
                    vector(array[index][0],array[index][1],array[index][2]),
                    particleCloud_.cellIDs()[index][0]
                );
                ofparticle.writePosition(fileStream);
                fileStream << nl;
            }
            else if (type == "vector")
            {
                fileStream << "( "<< array[index][0] << " " << array[index][1] << " " << array[index][2] << " ) " << " \n";
            }
        }
    }

    fileStream << ")\n";
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
IOModel::IOModel
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    dict_(dict),
    particleCloud_(sm),
    time_(sm.mesh().time()),
    parOutput_(true)
{
	if (
            particleCloud_.dataExchangeM().myType()=="oneWayVTK" ||
            dict_.found("serialOutput")
       )
    {
        parOutput_=false;
        Warning << "IO model is in serial write mode, only data on proc 0 is written" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

IOModel::~IOModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //


/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------

Application
    vtkDumpToLagrangian

Description
    Reads dump files from LIGGGHTS and converts it into OpenFOAM 
    lagrangian format   

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "IOstream.H"
#include "IFstream.H"
#include "vtkReader.H"
#include "vtkWriter.H"
#include "particle.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeHeader(fileName path, word name, word className)
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
}

void streamPointDataToPath(fvMesh& mesh, fileName path, pointField fld, word name, word type, word className)
{
    writeHeader(path, name, className);
    
    OFstream fileStream(fileName(path/name), IOstream::streamFormat::ASCII, IOstream::currentVersion, IOstream::compressionType::UNCOMPRESSED, true);
    fileStream << fld.size() <<"\n";
    
    if (type == "position")
    {
        fileStream << "(\n";
        for(int index = 0;index < fld.size() ; ++index)
        {
            particle ofparticle
            (
                mesh,
                fld[index],
                mesh.findCell(fld[index])
            );
            ofparticle.writePosition(fileStream);
            fileStream << nl;
        }
        fileStream << ")\n";
    }
}

void streamLabelDataToPath(fvMesh& mesh, fileName path, labelField fld, word name)
{
    writeHeader(path, name, "labelField");
    
    OFstream fileStream(fileName(path/name), IOstream::streamFormat::ASCII, IOstream::currentVersion, IOstream::compressionType::UNCOMPRESSED, true);

    if (name == "origProcId")
    {
        if (fld.size() > 0) fileStream << fld.size() << "{0}" << "\n";
        else fileStream << "{}" << "\n";
        return;
    }
    
    fileStream << fld << endl;
}

void streamScalarDataToPath(fvMesh& mesh, fileName path, scalarField fld, word name)
{
    writeHeader(path, name, "scalarField");
    
    OFstream fileStream(fileName(path/name), IOstream::streamFormat::ASCII, IOstream::currentVersion, IOstream::compressionType::UNCOMPRESSED, true);
    fileStream << fld << endl;
}

void streamVectorDataToPath(fvMesh& mesh, fileName path, vectorField fld, word name)
{
    writeHeader(path, name, "vectorField");
    
    OFstream fileStream(fileName(path/name));
    fileStream << fld << endl;
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    #include "readProperties.H"
    
    // loop over time steps
    forAll(timeDirs, i)
    {
        runTime.setTime(timeDirs[i], i);
        Info<< "Time = " << runTime.timeName() << endl;
        
        fileName cloudDir(runTime.path()/runTime.timeName()/"lagrangian"/cloudName_);
        if (!isDir(cloudDir))
        {
            mkDir(cloudDir);
        }

        std::string filename;

	label timeName_;

	if(selector=="TimeStep")
	{
	timeName_ = round(timeDirs[i].value()/demTimeStep);
	}
	else
	{
	timeName_ = demTimeStepNumber;
	}
        filename = dump_filename + std::to_string(timeName_) + ".vtk";
        Info << "Reading file "<< filename << endl;
        IFstream dstream(filename);
        Info << "Test 1" << endl;
        vtkReader reader(runTime, dstream);
        Info << "Test 2" << endl;
        pointField pts(reader.points());
        Info << "Test 3" << endl;
        labelField ids(reader.pointData().lookupObject<labelField>("id"));
        Info << "Test 4" << endl;
        
        streamPointDataToPath(mesh,cloudDir, pts,"positions","position","Cloud<passiveParticle>");

        Info << "Test 5" << endl;

        streamLabelDataToPath(mesh,cloudDir,ids,"origId");
        streamLabelDataToPath(mesh,cloudDir,ids,"origProcId");
        
        forAll(scalarFlds, fi)
        {
            scalarField sfld(reader.pointData().lookupObject<scalarField>(scalarFlds[fi]));
            streamScalarDataToPath(mesh,cloudDir,sfld,scalarFlds[fi]);
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "End\n" << endl;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    return 0;
}


// ************************************************************************* //

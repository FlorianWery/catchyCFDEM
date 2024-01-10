/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "vtkWriter.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "stringIOList.H"
#include "cellModeller.H"
#include "vectorIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkWriter::vtkWriter
(
    fileName& outFile,
    pointField pts,
    label n
)
:
    ostr_(outFile.c_str()),
    nattributes_(n)
{
    ostr_ << "# vtk DataFile Version 2.0" << std::endl;
    ostr_ << "Created by cfdemPost" << std::endl;
    ostr_ << "ASCII" << std::endl;
    ostr_ << "DATASET POLYDATA" << std::endl;
    writePoints(pts);
    writeVertices(pts.size());

    ostr_
        << "POINT_DATA " << pts.size() << std::endl
        << "FIELD attributes "<< nattributes_ << std::endl;
}


// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::vtkWriter::writePoints(pointField pts)
{
    ostr_<< "POINTS " << pts.size() << " float" << std::endl;

    forAll(pts, i)
    {
        ostr_  << pts[i].x() << ' '<< pts[i].y() << ' '<< pts[i].z() << std::endl;
    }
    ostr_ << std::endl;
}

void Foam::vtkWriter::writeVertices(label n)
{
    ostr_<< "VERTICES " << n << " " << 2*n << std::endl;

    for (label i=0; i<n; i++)
    {
        ostr_  << "1 " << i << std::endl;
    }
    ostr_ << std::endl;
}

void Foam::vtkWriter::writeVectorData(vectorField pdata, word name)
{
    ostr_<< name << " 3 " << pdata.size() << " float" << std::endl;
    forAll(pdata, i)
    {
        ostr_ << pdata[i].x() << ' '<< pdata[i].y() << ' '<< pdata[i].z() << std::endl;
    }
    ostr_ << std::endl;
}

void Foam::vtkWriter::writeScalarData(scalarField pdata, word name)
{
    ostr_<< name << " 1 " << pdata.size() << " float" << std::endl;
    forAll(pdata, i)
    {
        ostr_  << pdata[i];

        if (i > 0 && (i % 10) == 0)
        {
            ostr_  << std::endl;
        }
        else
        {
            ostr_  << ' ';
        }
    }
    ostr_ << std::endl;
}


void Foam::vtkWriter::writeLabelData(labelField pdata, word name)
{
    ostr_<< name << " 1 " << pdata.size() << " int" << std::endl;
    forAll(pdata, i)
    {
        ostr_  << pdata[i];

        if (i > 0 && (i % 10) == 0)
        {
            ostr_  << std::endl;
        }
        else
        {
            ostr_  << ' ';
        }
    }
    ostr_ << std::endl;
}



// ************************************************************************* //

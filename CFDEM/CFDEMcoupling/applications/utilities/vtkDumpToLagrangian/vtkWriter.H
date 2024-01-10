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

Class
    Foam::vtkWriter

Description


SourceFiles
    vtkWriter.C

\*---------------------------------------------------------------------------*/

#ifndef vtkWriter_H
#define vtkWriter_H

#include "OFstream.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class vtkWriter Declaration
\*---------------------------------------------------------------------------*/

class vtkWriter
{
    // Write stream
    std::ofstream ostr_;

    // Number of attributes to write
    label nattributes_;

public:

    //- Runtime type information
    ClassName("vtkWriter");

    // Constructors

        //- Construct from Istream, read all
        vtkWriter
        (
            fileName& outFile,
            pointField pts,
            label n
        );

    // Member Functions

        //- Write topology
        void writePoints(pointField pts);

        //- Write vertices        
        void writeVertices(label n);

        //- Write point data
        void writeVectorData(vectorField pdata, word name);
        void writeScalarData(scalarField pdata, word name);
        void writeLabelData(labelField pdata, word name);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

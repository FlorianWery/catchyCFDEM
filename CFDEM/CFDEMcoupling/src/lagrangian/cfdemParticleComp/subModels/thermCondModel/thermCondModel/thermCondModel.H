/*---------------------------------------------------------------------------*\
License

    This is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This code is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this code.  If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2015- Thomas Lichtenegger, JKU Linz, Austria


Description
    model for thermal conductivity of fluid phase in presence of particles

Class
    thermCondModel

SourceFiles
    thermCondModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermCondModel_H
#define thermCondModel_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class thermCondModel Declaration
\*---------------------------------------------------------------------------*/

class thermCondModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloudEnergy& particleCloud_;

public:

    //- Runtime type information
    TypeName("thermCondModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermCondModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloudEnergy& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        thermCondModel
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~thermCondModel();


    // Selector

        static autoPtr<thermCondModel> New
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Member Functions

        virtual tmp<volScalarField> thermCond() const = 0;

        virtual tmp<volScalarField> thermDiff() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

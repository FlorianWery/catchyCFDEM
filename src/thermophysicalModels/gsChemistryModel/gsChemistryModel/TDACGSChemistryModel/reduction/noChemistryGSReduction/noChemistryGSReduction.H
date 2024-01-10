/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::chemistryReductionMethods::none

Description

SourceFiles
    noChemistryGSReduction.C

\*---------------------------------------------------------------------------*/

#ifndef noChemistryGSReduction_H
#define noChemistryGSReduction_H

#include "chemistryGSReductionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryGSReductionMethods
{

/*---------------------------------------------------------------------------*\
                       Class none Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class GThermoType, class SThermoType>
class none
:
    public chemistryGSReductionMethod<CompType, GThermoType, SThermoType>
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        none
        (
            const IOdictionary& dict,
            TDACGSChemistryModel<CompType, GThermoType, SThermoType>& chemistry
        );


    //- Destructor
    virtual ~none();


    // Member Functions

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryGSReductionMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "noChemistryGSReduction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

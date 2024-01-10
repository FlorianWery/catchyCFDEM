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
    Foam::chemistryReductionMethod

Description
    An abstract class for methods of chemical mechanism reduction

SourceFiles
    chemistryReductionMethod.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryGSReductionMethod_H
#define chemistryGSReductionMethod_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class CompType, class GThermoType, class SThermoType>
class TDACGSChemistryModel;

/*---------------------------------------------------------------------------*\
                           Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class GThermoType, class SThermoType>
class chemistryGSReductionMethod
{

protected:

    const IOdictionary& dict_;

    //- Dictionary that store the algorithm data
    const dictionary coeffsDict_;

    //- Is mechanism reduction active?
    Switch active_;

    //- Switch to select performance logging
    Switch log_;

    TDACGSChemistryModel<CompType, GThermoType, SThermoType>& chemistry_;

    //- List of active species (active = true)
    List<bool> activeSpecies_;

    //- Number of active species
    label NsSimp_;
    label NsGasSimp_;
    label NsSolidSimp_;

    //- Number of species
    const label nSpecie_;
    const label nGasSpecie_;
    const label nSolidSpecie_;

    //- Tolerance for the mechanism reduction algorithm
    scalar tolerance_;


public:

    //- Runtime type information
    TypeName("chemistryGSReductionMethod");


    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        chemistryGSReductionMethod,
        dictionary,
        (
            const IOdictionary& dict,
            TDACGSChemistryModel<CompType, GThermoType, SThermoType>& chemistry
        ),
        (dict, chemistry)
    );


    // Constructors

        //- Construct from components
        chemistryGSReductionMethod
        (
            const IOdictionary& dict,
            TDACGSChemistryModel<CompType, GThermoType, SThermoType>& chemistry
        );


    // Selector

        static autoPtr<chemistryGSReductionMethod<CompType, GThermoType, SThermoType>> New
        (
            const IOdictionary& dict,
            TDACGSChemistryModel<CompType, GThermoType, SThermoType>& chemistry
        );


    //- Destructor
    virtual ~chemistryGSReductionMethod();


    // Member Functions

        //- Is mechanism reduction active?
        inline bool active() const;

        //- Is performance data logging enabled?
        inline bool log() const;

        //- Return the active species
        inline const List<bool>& activeSpecies() const;

        //- Return the number of active species
        inline label NsSimp();
        inline label NsGasSimp();
        inline label NsSolidSimp();

        //- Return the initial number of species
        inline label nSpecie();
        inline label nGasSpecie();
        inline label nSolidSpecie();

        //- Return the tolerance
        inline scalar tolerance() const;

        //- Reduce the mechanism
        virtual void reduceMechanism
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        ) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryGSReductionMethodI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "chemistryGSReductionMethod.C"
    #include "chemistryGSReductionMethodNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

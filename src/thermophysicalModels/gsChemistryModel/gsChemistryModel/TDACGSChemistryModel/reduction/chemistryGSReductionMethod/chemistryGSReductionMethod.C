/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "chemistryGSReductionMethod.H"
#include "TDACGSChemistryModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompType, class GThermoType, class SThermoType>
Foam::chemistryGSReductionMethod<CompType, GThermoType, SThermoType>::chemistryGSReductionMethod
(
    const Foam::IOdictionary& dict,
    Foam::TDACGSChemistryModel<CompType, GThermoType, SThermoType>& chemistry
)
:
    dict_(dict),
    coeffsDict_(dict.subDict("reduction")),
    active_(coeffsDict_.lookupOrDefault<Switch>("active", false)),
    log_(coeffsDict_.lookupOrDefault<Switch>("log", false)),
    chemistry_(chemistry),
    activeSpecies_(chemistry.nSpecie(), false),
    NsSimp_(chemistry.nSpecie()),
    NsGasSimp_(this->nGasSpecie_),
    NsSolidSimp_(this->nSolidSpecie_),
    nSpecie_(chemistry.nSpecie()),
    nGasSpecie_(chemistry.Y().size()),
    nSolidSpecie_(chemistry.Ys().size()),
    tolerance_(coeffsDict_.lookupOrDefault<scalar>("tolerance", 1e-4))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompType, class GThermoType, class SThermoType>
Foam::chemistryGSReductionMethod<CompType, GThermoType, SThermoType>::
~chemistryGSReductionMethod()
{}


// ************************************************************************* //

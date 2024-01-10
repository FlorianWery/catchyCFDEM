/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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

#include "PhaseThermophysicalTransportModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "makeThermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"

#include "fluidThermo.H"
#include "rhoThermo.H"

#include "laminarThermophysicalTransportModel.H"
#include "RASThermophysicalTransportModel.H"
#include "LESThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeThermophysicalTransportModelTypes
(
    PhaseThermophysicalTransportModel,
    fluidThermoPhaseCompressibleMomentumTransportModel,
    fluidThermo
);


makeThermophysicalTransportModels
(
    PhaseThermophysicalTransportModel,
    fluidThermoPhaseCompressibleMomentumTransportModel,
    fluidThermo
);


#define makeLaminarThermophysicalTransportModel(Type)                          \
    makeThermophysicalTransportModel                                           \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        fluidThermoPhaseCompressibleMomentumTransportModel,                    \
        fluidThermo,                                                           \
        laminar,                                                               \
        Type                                                                   \
    )

#define makeRASLESThermophysicalTransportModel(SType, Type)                    \
    makeTurbulenceThermophysicalTransportModel                                 \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        fluidThermoPhaseCompressibleMomentumTransportModel,                    \
        fluidThermo,                                                           \
        SType,                                                                 \
        Type                                                                   \
    )

#define makeRASThermophysicalTransportModel(Type)                              \
    makeThermophysicalTransportModel                                           \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        fluidThermoPhaseCompressibleMomentumTransportModel,                    \
        fluidThermo,                                                           \
        RAS,                                                                   \
        Type                                                                   \
    )

#define makeLESThermophysicalTransportModel(Type)                              \
    makeThermophysicalTransportModel                                           \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        fluidThermoPhaseCompressibleMomentumTransportModel,                    \
        fluidThermo,                                                           \
        LES,                                                                   \
        Type                                                                   \
    )


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

#include "Fourier.H"
makeLaminarThermophysicalTransportModel(Fourier);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "eddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(RAS, eddyDiffusivity);


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "eddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(LES, eddyDiffusivity);


// ************************************************************************* //

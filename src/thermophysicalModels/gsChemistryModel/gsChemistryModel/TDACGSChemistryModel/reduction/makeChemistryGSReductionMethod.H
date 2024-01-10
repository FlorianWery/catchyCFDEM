/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#ifndef makeChemistryGSReductionMethod_H
#define makeChemistryGSReductionMethod_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistryGSReductionMethod(ReactionThermo, SThermoPhysics, GThermoPhysics)          \
                                                                               \
    typedef chemistryGSReductionMethod<ReactionThermo, GThermoPhysics, SThermoPhysics>            \
        chemistryGSReductionMethod##ReactionThermo##GThermoPhysics##SThermoPhysics;               \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryGSReductionMethod##ReactionThermo##GThermoPhysics##SThermoPhysics,               \
        (                                                                      \
            word                                                               \
            (                                                                  \
                chemistryGSReductionMethod##ReactionThermo##GThermoPhysics##SThermoPhysics::      \
                typeName_()                                                    \
            )                                                                  \
          + "<" + word(ReactionThermo::typeName_())                            \
          + "," + GThermoPhysics::typeName()                                   \
          + "," + SThermoPhysics::typeName() + ">"                             \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        chemistryGSReductionMethod##ReactionThermo##GThermoPhysics##SThermoPhysics,               \
        dictionary                                                             \
    );                                                                         \

#define makeChemistryGSReductionMethod(Method, ReactionThermo, SThermoPhysics, GThermoPhysics)    \
                                                                               \
    typedef chemistryGSReductionMethods::Method<ReactionThermo, GThermoPhysics, SThermoPhysics>   \
        chemistryGSReductionMethod##Method##ReactionThermo##GThermoPhysics##SThermoPhysics;       \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryGSReductionMethod##Method##ReactionThermo##GThermoPhysics##SThermoPhysics,       \
        (                                                                      \
            word                                                               \
            (                                                                  \
                chemistryGSReductionMethod##Method##ReactionThermo##GThermoPhysics##SThermoPhysics::typeName_()                                                     \
            )                                                                  \
         + "<" + word(ReactionThermo::typeName_())                             \
         + "," + GThermoPhysics::typeName()                                    \
         + "," + SThermoPhysics::typeName() + ">"                               \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    chemistryGSReductionMethod<ReactionThermo, GThermoPhysics, SThermoPhysics>::                  \
        adddictionaryConstructorToTable                                        \
        <chemistryGSReductionMethod##Method##ReactionThermo##GThermoPhysics##SThermoPhysics>      \
        add##chemistryGSReductionMethod##Method##ReactionThermo##GThermoPhysics##SThermoPhysics##\
dictionaryConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

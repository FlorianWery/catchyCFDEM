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

#ifndef makeChemistryGSTabulationMethod_H
#define makeChemistryGSTabulationMethod_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistryGSTabulationMethod(ReactionThermo, SThermoPhysics, GThermoPhysics)         \
                                                                               \
    typedef chemistryGSTabulationMethod<ReactionThermo, GThermoPhysics, SThermoPhysics>           \
        chemistryGSTabulationMethod##ReactionThermo##GThermoPhysics##SThermoPhysics;              \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryGSTabulationMethod##ReactionThermo##GThermoPhysics##SThermoPhysics,              \
        (                                                                      \
            word                                                               \
            (                                                                  \
                chemistryGSTabulationMethod##ReactionThermo##GThermoPhysics##SThermoPhysics::     \
                typeName_()                                                    \
            )                                                                  \
          + "<"#ReactionThermo","                                              \
          + GThermoPhysics::typeName()                                         \
          + "," + SThermoPhysics::typeName() + ">"                             \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        chemistryGSTabulationMethod##ReactionThermo##GThermoPhysics##SThermoPhysics,              \
        dictionary                                                             \
    );                                                                         \

#define makeChemistryGSTabulationMethod(Method, ReactionThermo, SThermoPhysics, GThermoPhysics)   \
                                                                               \
    typedef chemistryGSTabulationMethods::Method<ReactionThermo, GThermoPhysics, SThermoPhysics>  \
        chemistryGSTabulationMethod##Method##ReactionThermo##GThermoPhysics##SThermoPhysics;      \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        chemistryGSTabulationMethod##Method##ReactionThermo##GThermoPhysics##SThermoPhysics,      \
        (                                                                      \
            word                                                               \
            (                                                                  \
                chemistryGSTabulationMethod##Method##ReactionThermo##\
GThermoPhysics##SThermoPhysics::typeName_()                                                     \
            )                                                                  \
         + "<"#ReactionThermo","                                               \
         + GThermoPhysics::typeName()                                          \
         + "," + SThermoPhysics::typeName() + ">"                              \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    chemistryGSTabulationMethod<ReactionThermo, GThermoPhysics, SThermoPhysics>::                 \
        adddictionaryConstructorToTable                                        \
        <chemistryGSTabulationMethod##Method##ReactionThermo##GThermoPhysics##SThermoPhysics>     \
        add##chemistryGSTabulationMethod##Method##ReactionThermo##GThermoPhysics##SThermoPhysics##\
dictionaryConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

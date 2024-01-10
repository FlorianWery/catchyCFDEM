/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::interstitialUniformInletVelocityFvPatchVectorField

Description
    Inlet velocity in which the actual interstitial velocity is calculated
    by dividing the specified inletVelocity field with the local phase-fraction.

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type              interstitialUniformInletVelocity;
        inletVelocity     uniform (0 0.2 0);// Non-interstitial inlet velocity
        alpha             alpha.particles;  // Name of the phase-fraction field
        value             uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    interstitialUniformInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef interstitialUniformInletVelocityFvPatchVectorField_H
#define interstitialUniformInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class interstitialUniformInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class interstitialUniformInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet velocity
        vectorField inletVelocity_;
        autoPtr<Function1<vector>> uniformValue_;

        //- Name of the flux transporting the field
        word alphaName_;


public:

   //- Runtime type information
   TypeName("interstitialUniformInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        interstitialUniformInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interstitialUniformInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  interstitialUniformInletVelocityFvPatchVectorField
        //  onto a new patch
        interstitialUniformInletVelocityFvPatchVectorField
        (
            const interstitialUniformInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        interstitialUniformInletVelocityFvPatchVectorField
        (
            const interstitialUniformInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new interstitialUniformInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interstitialUniformInletVelocityFvPatchVectorField
        (
            const interstitialUniformInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new interstitialUniformInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

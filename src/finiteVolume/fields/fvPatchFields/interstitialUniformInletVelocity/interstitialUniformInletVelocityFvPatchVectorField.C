/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "interstitialUniformInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interstitialUniformInletVelocityFvPatchVectorField::
interstitialUniformInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    inletVelocity_(p.size(), Zero),
    uniformValue_(),
    alphaName_("alpha")
{}


Foam::interstitialUniformInletVelocityFvPatchVectorField::
interstitialUniformInletVelocityFvPatchVectorField
(
    const interstitialUniformInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    inletVelocity_(mapper(ptf.inletVelocity_)),
    uniformValue_(ptf.uniformValue_, false),
    alphaName_(ptf.alphaName_)
{}


Foam::interstitialUniformInletVelocityFvPatchVectorField::
interstitialUniformInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    inletVelocity_(p.size(), Zero),
    uniformValue_(Function1<vector>::New("inletVelocity", dict)),
    alphaName_(dict.lookupOrDefault<word>("alpha", "alpha"))
{}


Foam::interstitialUniformInletVelocityFvPatchVectorField::
interstitialUniformInletVelocityFvPatchVectorField
(
    const interstitialUniformInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    inletVelocity_(ptf.inletVelocity_),
    uniformValue_(ptf.uniformValue_, false),
    alphaName_(ptf.alphaName_)
{}


Foam::interstitialUniformInletVelocityFvPatchVectorField::
interstitialUniformInletVelocityFvPatchVectorField
(
    const interstitialUniformInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    inletVelocity_(ptf.inletVelocity_),
    uniformValue_(ptf.uniformValue_, false),
    alphaName_(ptf.alphaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interstitialUniformInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<scalar>& alphap =
        patch().lookupPatchField<volScalarField, scalar>(alphaName_);

    const scalar t = this->db().time().timeOutputValue();

    forAll(inletVelocity_, facei)
    {
        inletVelocity_[facei] = uniformValue_->value(t);
    }

    operator==(inletVelocity_/alphap);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::interstitialUniformInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    writeEntryIfDifferent<word>(os, "alpha", "alpha", alphaName_);
    writeEntry(os, uniformValue_());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       interstitialUniformInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //

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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "noCatalyticChemistry.H"
#include "addToRunTimeSelectionTable.H"

#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noCatalyticChemistry, 0);

addToRunTimeSelectionTable(catalyticChemistryModel, noCatalyticChemistry, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
noCatalyticChemistry::noCatalyticChemistry
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    catalyticChemistryModel(dict,sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noCatalyticChemistry::~noCatalyticChemistry()
{}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void noCatalyticChemistry::execute()
{
    //do nothing
}

const volScalarField noCatalyticChemistry::RR(label i) const
{
    volScalarField RRi_
    (
        IOobject
        (
            "RR_",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
         particleCloud_.mesh(),
         dimensionedScalar("zero", dimMass/(dimVol*dimTime), 0.0)
    );

    return RRi_;
}

const volScalarField noCatalyticChemistry::Qdot() const
{
    volScalarField Qdot_
    (
        IOobject
        (
            "Qdot_",
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
         ),
         particleCloud_.mesh(),
         dimensionedScalar("zero", dimEnergy/(dimVol*dimTime), 0.0)
    );

    return Qdot_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

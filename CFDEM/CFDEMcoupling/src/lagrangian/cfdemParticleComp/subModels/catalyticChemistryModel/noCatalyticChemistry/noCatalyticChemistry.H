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

    Description
    Catalytic chemistry turned off

\*---------------------------------------------------------------------------*/

#ifndef noCatalyticChemistry_H
#define noCatalyticChemistry_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "catalyticChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noCatalyticChemistry Declaration
\*---------------------------------------------------------------------------*/


class noCatalyticChemistry
:
    public catalyticChemistryModel
{
public:

    //- Runtime type information
    TypeName("off");

    // Constructors

        //- Construct from components
        noCatalyticChemistry
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        ~noCatalyticChemistry();


    // Member Functions

        // - Calculate source terms
        virtual void execute();

        // - Return source terms for species i
        virtual const volScalarField RR(const label i) const;

        // - Return reaction heat
        virtual const volScalarField Qdot() const;

        virtual void postFlow() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

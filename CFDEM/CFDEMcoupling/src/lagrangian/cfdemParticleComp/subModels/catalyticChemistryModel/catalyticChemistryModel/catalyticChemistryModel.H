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

#ifndef catalyticChemistryModel_H
#define catalyticChemistryModel_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "probeModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class chemistryModel Declaration
\*---------------------------------------------------------------------------*/

class catalyticChemistryModel
{

protected:

    // Protected data
        const dictionary& dict_;

        cfdemCloudEnergy& particleCloud_;



public:

    //- Runtime type information
    TypeName("catalyticChemistryModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            catalyticChemistryModel,
            dictionary,
            (
                const dictionary& dict,
                cfdemCloudEnergy& sm
            ),
            (dict,sm)
        );


    // Constructors

        //- Construct from components
        catalyticChemistryModel
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~catalyticChemistryModel();


    // Selector

        static autoPtr<catalyticChemistryModel> New
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm,
            word catalyticChemistryModelType
        );


    // Member Functions

        // - Calculate source terms
        virtual void execute() = 0;

        // - Return source terms for species i
        virtual const volScalarField RR(const label i) const = 0;

        // - Return reaction heat
        virtual const volScalarField Qdot() const = 0;

        virtual void postFlow() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

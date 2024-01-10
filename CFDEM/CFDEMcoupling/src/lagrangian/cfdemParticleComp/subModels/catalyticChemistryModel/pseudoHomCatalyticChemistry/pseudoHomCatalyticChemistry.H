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
    Catalytic chemistry without mass transfer resistances

\*---------------------------------------------------------------------------*/

#ifndef pseudoHomCatalyticChemistry_H
#define pseudoHomCatalyticChemistry_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "catalyticChemistryModel.H"

#include "autoPtr.H"
#include "speciesTable.H"
#include "rhoReactionThermo.H"
#include "basicGSChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pseudoHomCatalyticChemistry Declaration
\*---------------------------------------------------------------------------*/


class pseudoHomCatalyticChemistry
:
    public catalyticChemistryModel
{
private:

    dictionary propsDict_;

    //- Name of the catalyst
    const word catalystName_;

    //- Pointer to thermo
    autoPtr<rhoReactionThermo> catThermoPtr_;
    rhoReactionThermo&  catThermo_;

    //- Pointer to chemistry model
    autoPtr<basicGSChemistryModel> catChemistryPtr_;

    //- Reference to the fluid temperature field
    const volScalarField& Tg_;

    //- Reference to the fluid density field
    const volScalarField& rhog_;

    //- Reference to the pressure field
    const volScalarField& p_;

    //- Reference to the field of fluid specie mass fractions
    const PtrList<volScalarField>& Yg_;

    //- Reference to the field of fluid specie mass fractions in pores
    PtrList<volScalarField>& Ygs_;

    //- Reference to the field of coverages
    PtrList<volScalarField>& Ys_;

    //- Number of gas species
    const label nGasSpecie_;

    //- Number of surface species
    const label nSolidSpecie_;

    //- Reference to surface species table
    const speciesTable& solidSpecies_;

    //- Switch to use fluid temeprature instead of particle temperature
    const bool useFluidTemperature_;

    //- Field of reaction rate per fluid specie [kg/m^3/s]
    PtrList<volScalarField> RR_;

    //- Field of reaction head [J/m^3/s = kg/m^2/s^3]
    volScalarField Qdot_;

    //- Particle temperature
    word partTempName_;
    mutable double **partTemp_;

    //- Temporary fields for calculations
    mutable scalarField dcdt_;
    mutable scalarField Yi_;
    mutable scalarField Ysi_;

    // Private Member Functions

        // - Allocate arrays
        virtual void allocateCatChemArrays();
        virtual void reallocateCatChemArrays();

public:

    //- Runtime type information
    TypeName("pseudoHomCatalyticChemistry");

    // Constructors

        //- Construct from components
        pseudoHomCatalyticChemistry
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        ~pseudoHomCatalyticChemistry();


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

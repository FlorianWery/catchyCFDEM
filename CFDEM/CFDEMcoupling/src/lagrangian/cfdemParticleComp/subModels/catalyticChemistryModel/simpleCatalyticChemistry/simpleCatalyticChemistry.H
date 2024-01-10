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

#ifndef simpleCatalyticChemistry_H
#define simpleCatalyticChemistry_H

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
                           Class simpleCatalyticChemistry Declaration
\*---------------------------------------------------------------------------*/


class simpleCatalyticChemistry
:
    public catalyticChemistryModel
{
private:

    dictionary propsDict_;

    //- Name of the catalyst
    const word catalystName_;

    //- Pointer to thermo
    autoPtr<rhoReactionThermo> catThermoPtr_;

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
    const speciesTable solidSpecies_;

    //- Reference to gas species table
    const speciesTable gasSpecies_;

    //- Switch to use fluid temperature instead of particle temperature
    const bool useFluidTemperature_;

    //- Field of reaction rate per fluid specie [kg/m^3/s]
    PtrList<volScalarField> RR_;

    //- Field of reaction heat [J/m^3/s = kg/m^2/s^3]
    volScalarField Qdot_;

    //- List of total reaction rate per fluid specie [kg/s] per particle
    mutable double **partRR_;
    mutable List<scalarList> ofPartRR_;

    //- List of surface species coverages on particles
    mutable double **partCoverages_;
    mutable List<scalarList> ofPartCoverages_;

    //- List of gas species concentrations in particles
    mutable double **partGas_;
    mutable List<scalarList> ofPartGas_;

    mutable double **partPressure_;

    //- Particle temperature
    word partTempName_;
    mutable double **partTemp_;

    //- Reaction heat per particle [J/s = kg/m^2/s^3]
    mutable double **partQdot_;
    word partHeatSourceName_;
    mutable double **partHeatSource_;

    //- Switch to release heat in particles
    const bool useParticleQdot_;

    //- Temporary fields for calculations
    mutable scalarField dcdt_;
    mutable scalarField Yi_;
    mutable scalarField Ysi_;

    //- Initialization flag
    bool ini_;

    word partPressName_;

    //- Flag for diffusive mass transfer
    bool massTransfer_;

    // Private Member Functions

        // - Allocate arrays
        virtual void allocateMyArrays() const;

public:

    //- Runtime type information
    TypeName("simpleCatalyticChemistry");

    // Constructors

        //- Construct from components
        simpleCatalyticChemistry
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        ~simpleCatalyticChemistry();


    // Member Functions

        inline double ** voidfractions();

        // - Calculate source terms
        virtual void execute();

        // - Return source terms for species i
        virtual const volScalarField RR(const label i) const;

        // - Return reaction heat
        virtual const volScalarField Qdot() const;

        virtual void postFlow();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

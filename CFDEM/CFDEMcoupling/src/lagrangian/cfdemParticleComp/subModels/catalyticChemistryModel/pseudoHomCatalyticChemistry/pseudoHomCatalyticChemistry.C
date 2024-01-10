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
#include "pseudoHomCatalyticChemistry.H"
#include "addToRunTimeSelectionTable.H"

#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pseudoHomCatalyticChemistry, 0);

addToRunTimeSelectionTable(catalyticChemistryModel, pseudoHomCatalyticChemistry, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pseudoHomCatalyticChemistry::pseudoHomCatalyticChemistry
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    catalyticChemistryModel(dict, sm),
    propsDict_(dict.subDict(typeName + "Props")),
    catalystName_(propsDict_.lookup("catalystName")),
    catThermoPtr_(rhoReactionThermo::New(particleCloud_.mesh(), catalystName_)),
    catThermo_(catThermoPtr_()),
    catChemistryPtr_(basicGSChemistryModel::New(catThermo_)),
    Tg_(sm.thermo().T()),
    rhog_(sm.thermo().rho()),
    p_(sm.thermo().p()),
    Yg_(dynamic_cast<const rhoReactionThermo&>(sm.thermo()).composition().Y()),
    Ygs_(catThermo_.composition().Y()),
    Ys_(catThermo_.composition().Ys()),
    nGasSpecie_(Yg_.size()),
    nSolidSpecie_(Ys_.size()),
    solidSpecies_(catThermo_.composition().solidSpecies()),
    useFluidTemperature_(propsDict_.lookupOrDefault("useFluidTemperature", false)),
    RR_(nGasSpecie_),
    Qdot_
    (
        IOobject
        (
            IOobject::groupName("Qdot", catalystName_),
            particleCloud_.mesh().time().timeName(),
            particleCloud_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        particleCloud_.mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
    ),
    partTempName_(propsDict_.lookupOrDefault<word>("partTempName", "Temp")),
    partTemp_(NULL),
    dcdt_(nGasSpecie_, 0.0),
    Yi_(nGasSpecie_, 0.0),
    Ysi_(nSolidSpecie_, 0.0)
{
    particleCloud_.checkCG(false);

    allocateCatChemArrays();

    forAll(RR_, fieldi)
    {
        RR_.set
        (
            fieldi,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "RR."+catThermo_.composition().Y()[fieldi].name(),
                        catalystName_
                    ),
                    particleCloud_.mesh().time().timeName(),
                    particleCloud_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                particleCloud_.mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0)
            )
        );
    }

    // initialize surface
    {
        label celli = 0.0;
        scalar pi = p_[celli];
        scalar Tgi = Tg_[celli];
        scalar Tsi = catThermo_.T()[celli];

        for (label i=0; i<nGasSpecie_; i++) Yi_[i] = Ygs_[i][celli];
        for (label i=0; i<nSolidSpecie_; i++) Ysi_[i] = Ys_[i][celli];

        Qdot_[celli] =
            catChemistryPtr_->getRatesQdotI
            (
                particleCloud_.mesh().time().deltaTValue(),
                celli,
                dcdt_,
                Yi_,
                Ysi_,
                rhog_[celli],
                Tgi,
                Tsi,
                pi,
                true
            );
    }

    Info<< "pseudoHomCatalyticChemistryModel: " << nl
        << indent << "Number of gas species = " << nGasSpecie_ << nl
        << indent << "Number of solid species = " << nSolidSpecie_ << nl << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pseudoHomCatalyticChemistry::~pseudoHomCatalyticChemistry()
{
    particleCloud_.dataExchangeM().destroy(partTemp_, 1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void pseudoHomCatalyticChemistry::allocateCatChemArrays()
{
    if (particleCloud_.dataExchangeM().maxNumberOfParticles() > 0)
    {
        if (!useFluidTemperature_) particleCloud_.dataExchangeM().allocateArray(partTemp_, 0.0, 1, "nparticles");
    }
}

void pseudoHomCatalyticChemistry::reallocateCatChemArrays()
{
    if (!useFluidTemperature_) particleCloud_.dataExchangeM().allocateArray(partTemp_, 0.0, 1, "nparticles");
}

void pseudoHomCatalyticChemistry::execute()
{
    reallocateCatChemArrays();

    forAll(Ygs_, i)
    {
        forAll(Ygs_[i], celli)
        {
            Ygs_[i][celli] = Yg_[i][celli];
        }
    }

    if (useFluidTemperature_)
    {
        catThermo_.T() = Tg_;
    }
    else
    {
        particleCloud_.dataExchangeM().getData(partTempName_, "scalar-atom", partTemp_);
        catThermo_.T().primitiveFieldRef() = 0.0;
        particleCloud_.averagingM().resetWeightFields();
        particleCloud_.averagingM().setScalarAverage
        (
            catThermo_.T(),
            partTemp_,
            particleCloud_.particleWeights(),
            particleCloud_.averagingM().UsWeightField(),
            NULL
        );
    }

    volScalarField alphap(1.0-particleCloud_.voidFractionM().voidFractionInterp());

    forAll(Tg_ ,celli)
    {
        scalar pi = p_[celli];
        scalar Tgi = Tg_[celli];
        scalar Tsi = catThermo_.T()[celli];

        for (label i=0; i<nGasSpecie_; i++) Yi_[i] = Ygs_[i][celli];
        for (label i=0; i<nSolidSpecie_; i++) Ysi_[i] = Ys_[i][celli];

        Qdot_[celli] =
            catChemistryPtr_->getRatesQdotI
            (
                particleCloud_.mesh().time().deltaTValue(),
                celli,
                dcdt_,
                Yi_,
                Ysi_,
                rhog_[celli],
                Tgi,
                Tsi,
                pi,
                false // initialize coverages
            )*alphap[celli];
        for (label i=0; i<nGasSpecie_; i++) RR_[i][celli] = dcdt_[i]*alphap[celli];
        for (label i=0; i<nSolidSpecie_; i++) Ys_[i][celli] = Ysi_[i];
    }
}

const volScalarField pseudoHomCatalyticChemistry::RR(label i) const
{
    return RR_[i];
}

const volScalarField pseudoHomCatalyticChemistry::Qdot() const
{
    return Qdot_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

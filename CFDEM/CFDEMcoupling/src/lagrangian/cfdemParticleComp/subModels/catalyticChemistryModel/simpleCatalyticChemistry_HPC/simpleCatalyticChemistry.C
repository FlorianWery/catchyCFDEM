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
#include "simpleCatalyticChemistry.H"
#include "addToRunTimeSelectionTable.H"

#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(simpleCatalyticChemistry, 0);

addToRunTimeSelectionTable(catalyticChemistryModel, simpleCatalyticChemistry, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
simpleCatalyticChemistry::simpleCatalyticChemistry
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    catalyticChemistryModel(dict, sm),
    propsDict_(dict.subDict(typeName + "Props")),
    catalystName_(propsDict_.lookup("catalystName")),
    catThermoPtr_(rhoReactionThermo::New(particleCloud_.mesh(), catalystName_)),
    catChemistryPtr_(basicGSChemistryModel::New(catThermoPtr_())),
    Tg_(sm.thermo().T()),
    rhog_(sm.thermo().rho()),
    p_(sm.thermo().p()),
    Yg_(dynamic_cast<const rhoReactionThermo&>(sm.thermo()).composition().Y()),
    Ygs_(catThermoPtr_->composition().Y()),
    Ys_(catThermoPtr_->composition().Ys()),
    nGasSpecie_(Yg_.size()),
    nSolidSpecie_(Ys_.size()),
    solidSpecies_(catThermoPtr_->composition().solidSpecies()),
    gasSpecies_(catThermoPtr_->composition().species()),
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
        dimensionedScalar(dimEnergy/dimVolume/dimTime, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    partRR_(NULL),
    ofPartRR_(nGasSpecie_),
    partCoverages_(NULL),
    ofPartCoverages_(nSolidSpecie_),
    partGas_(NULL),
    ofPartGas_(nGasSpecie_),
    partTempName_(propsDict_.lookupOrDefault<word>("partTempName", "Temp")),
    partTemp_(NULL),
    partQdot_(NULL),
    partHeatSourceName_(propsDict_.lookupOrDefault<word>("partHeatSourceName", "heatSource")),
    partHeatSource_(NULL),
    useParticleQdot_(propsDict_.lookupOrDefault("useParticleQdot", true)),
    dcdt_(nGasSpecie_, 0.0),
    Yi_(nGasSpecie_, 0.0),
    Ysi_(nSolidSpecie_, 0.0),
    ini_(propsDict_.lookupOrDefault("initialized", false)),
    massTransfer_(false)
{
    if(sm.nrMassTransferModels() != 0)
    {
        massTransfer_ = true;
    };

    allocateMyArrays();

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
                        "RR."+Yg_[fieldi].name(),
                        catalystName_
                    ),
                    particleCloud_.mesh().time().timeName(),
                    particleCloud_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                particleCloud_.mesh(),
                dimensionedScalar(dimMass/dimVolume/dimTime, 0),
                zeroGradientFvPatchField<scalar>::typeName
            )
        );
    }
    particleCloud_.checkCG(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

simpleCatalyticChemistry::~simpleCatalyticChemistry()
{
    if(massTransfer_)
    {
        delete partGas_;
    };
    delete partTemp_;
    delete partQdot_;
    delete partRR_;
    delete partCoverages_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void simpleCatalyticChemistry::allocateMyArrays() const
{
    double initVal = 0.0;
    particleCloud_.dataExchangeM().allocateArray(partTemp_, initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(partQdot_, initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(partHeatSource_, initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(partRR_, initVal, 1);
    particleCloud_.dataExchangeM().allocateArray(partCoverages_, initVal, 1);
    forAll(ofPartRR_, i) ofPartRR_[i].setSize(particleCloud_.numberOfParticles(), 0.0);
    forAll(ofPartCoverages_, i) ofPartCoverages_[i].setSize(particleCloud_.numberOfParticles(), Ys_[i][0]);

    if(massTransfer_)
    {
        forAll(ofPartGas_, i) ofPartGas_[i].setSize(particleCloud_.numberOfParticles(), Ygs_[i][0]);
        particleCloud_.dataExchangeM().allocateArray(partGas_, initVal, 1);
    };

    if (useParticleQdot_) particleCloud_.dataExchangeM().giveData(partHeatSourceName_, "scalar-atom", partHeatSource_);
}

void simpleCatalyticChemistry::execute()
{
    forAll(Ygs_, i)
    {
        forAll(Ygs_[i], celli)
        {
            Ygs_[i][celli] = Yg_[i][celli];
        }
    }

    // initialize surface
    if (!ini_)
    {
        label celli = 0;
        scalar pi = p_[celli];
        scalar Tgi = Tg_[celli];
        scalar Tsi = catThermoPtr_->T()[celli];

        forAll(Yi_, i) Yi_[i] = Ygs_[i][celli];
        forAll(Ysi_, i) Ysi_[i] = Ys_[i][celli];

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
                true, // initialize coverages
                !useParticleQdot_ // if particleQdot = true, use absolute energy (sensible energy = false)
            );

        allocateMyArrays();

        ini_ = true;
    }
    else
    {
        allocateMyArrays();

        forAll(ofPartCoverages_, i)
        {
            particleCloud_.dataExchangeM().getData(solidSpecies_[i], "scalar-atom", partCoverages_);
            forAll(ofPartCoverages_[i], index)
            {
                ofPartCoverages_[i][index] = partCoverages_[index][0];
            }
        }
        if(massTransfer_)
        {
            forAll(ofPartGas_, i)
            {
                particleCloud_.dataExchangeM().getData(gasSpecies_[i], "scalar-atom", partGas_);
                forAll(ofPartGas_[i], index)
                {
                    ofPartGas_[i][index] = partGas_[index][0];
                }
            } 
        };
    }

    if (!useFluidTemperature_)
    {
        particleCloud_.dataExchangeM().getData(partTempName_, "scalar-atom", partTemp_);
    }

    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        label cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
            if(!massTransfer_)
            {
                for (label i=0; i<nGasSpecie_; i++) Yi_[i] = Ygs_[i][cellI];
            }
            else
            {
                for (label i=0; i<nGasSpecie_; i++) Yi_[i] = ofPartGas_[i][index];
            };
            for (label i=0; i<nSolidSpecie_; i++) Ysi_[i] = ofPartCoverages_[i][index];

            if (useFluidTemperature_)
            {
                scalar pi = p_[cellI];
                scalar Tgi = Tg_[cellI];

                partQdot_[index][0] =
                    catChemistryPtr_->getRatesQdotI
                    (
                        particleCloud_.mesh().time().deltaTValue(),
                        cellI,
                        dcdt_,
                        Yi_,
                        Ysi_,
                        rhog_[cellI],
                        Tg_[cellI],
                        Tgi,
                        pi,
                        false // initialize coverages
                    );
                partQdot_[index][0] *= particleCloud_.particleVolume(index);
            }
            else
            {
                scalar pi = p_[cellI];

                partQdot_[index][0] =
                    catChemistryPtr_->getRatesQdotI
                    (
                        particleCloud_.mesh().time().deltaTValue(),
                        cellI,
                        dcdt_,
                        Yi_,
                        Ysi_,
                        rhog_[cellI],
                        Tg_[cellI],
                        partTemp_[index][0],
                        pi,
                        false, // initialize coverages
                        !useParticleQdot_ // if particleQdot = true, use absolute energy (sensible energy = false)
                    );
                partQdot_[index][0] *= particleCloud_.particleVolume(index);
            }
            for (label i=0; i<nGasSpecie_; i++) ofPartRR_[i][index] = dcdt_[i]*particleCloud_.particleVolume(index);
            for (label i=0; i<nSolidSpecie_; i++) ofPartCoverages_[i][index] = Ysi_[i];
        }
    }

    forAll(Ys_, i)
    {
        forAll(ofPartCoverages_[i], index)
        {
            partCoverages_[index][0] = ofPartCoverages_[i][index];
        }
        Ys_[i].primitiveFieldRef() = 0.0;
        particleCloud_.averagingM().resetWeightFields();
        particleCloud_.averagingM().setScalarAverage
        (
            Ys_[i],
            partCoverages_,
            particleCloud_.particleWeights(),
            particleCloud_.averagingM().UsWeightField(),
            NULL
        );
        Info << Ys_[i].name() << " max/min: " << max(Ys_[i]).value() << "/" << min(Ys_[i]).value() << endl;
    }

    Qdot_.primitiveFieldRef() = 0.0;
    if (!useParticleQdot_)
    {
        particleCloud_.averagingM().setScalarSum
        (
            Qdot_,
            partQdot_,
            particleCloud_.particleWeights(),
            NULL
        );
        Qdot_.primitiveFieldRef() /= Qdot_.mesh().V();
    }

    forAll(RR_, i)
    {
        forAll(ofPartRR_[i], index)
        {
            partRR_[index][0] = ofPartRR_[i][index];
        }
        RR_[i].primitiveFieldRef() = 0.0;
        particleCloud_.averagingM().setScalarSum
        (
            RR_[i],
            partRR_,
            particleCloud_.particleWeights(),
            NULL
        );
        RR_[i].primitiveFieldRef() /= RR_[i].mesh().V();
    }
}

const volScalarField simpleCatalyticChemistry::RR(label i) const
{
    return RR_[i];
}

const volScalarField simpleCatalyticChemistry::Qdot() const
{
    return Qdot_;
}

void simpleCatalyticChemistry::postFlow()
{
    forAll(ofPartCoverages_, i)
    {
        forAll(ofPartCoverages_[i], index)
        {
            partCoverages_[index][0] = particleCloud_.cellIDs()[index][0] >= 0 ? ofPartCoverages_[i][index] : 0.0;
        }
        particleCloud_.dataExchangeM().giveData(solidSpecies_[i], "scalar-atom", partCoverages_);
    }

    if(massTransfer_)
    {
        forAll(ofPartGas_, i)
        {
            forAll(ofPartGas_[i], index)
            {
                partGas_[index][0] = particleCloud_.cellIDs()[index][0] >= 0 ? ofPartGas_[i][index] : 0.0;
            }
            particleCloud_.dataExchangeM().giveData(gasSpecies_[i], "scalar-atom", partGas_);
        }
    };

    if (useParticleQdot_)
    {
        // energyModel::postFlow() has been called before this function (see cfdemCloudEnergy.C),
        // so the data in partHeatSource_ is the heat transfer from the current time step (e.g. via Gunn)
        // Here we add the reaction heat and send the field back to LIGGGHTS
        particleCloud_.dataExchangeM().getData(partHeatSourceName_,"scalar-atom", partHeatSource_);
        for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
        {
            label cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                partHeatSource_[index][0] = partHeatSource_[index][0] + partQdot_[index][0];
            }
            else
            {
                partHeatSource_[index][0] = 0.0;
            }
        }
        particleCloud_.dataExchangeM().giveData(partHeatSourceName_, "scalar-atom", partHeatSource_);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

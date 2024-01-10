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

    Copyright (C) 2015- Thomas Lichtenegger, JKU Linz, Austria

\*---------------------------------------------------------------------------*/

#include "cfdemCloudEnergy.H"
#include "energyModel.H"
#include "massTransferModel.H"
#include "thermCondModel.H"
#include "catalyticChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
cfdemCloudEnergy::cfdemCloudEnergy
(
    const fvMesh& mesh
)
:
    cfdemCloud(mesh),
    energyModels_(couplingProperties_.lookup("energyModels")),
    implicitEnergyModel_(false),
    massTransferModels_(couplingProperties_.lookup("massTransferModels")),
    catalyticChemistryModels_(couplingProperties_.lookup("catalyticChemistryModels")),
    energyModel_(nrEnergyModels()),
    massTransferModel_(nrMassTransferModels()),
    thermCondModel_
    (
        thermCondModel::New
        (
            couplingProperties_,
            *this
        )
    ),
    catalyticChemistryModel_(nrCatalyticChemistryModels())
{
    forAll(energyModels_, modeli)
    {
        energyModel_.set
        (
            modeli,
            energyModel::New
            (
                couplingProperties_,
                *this,
                energyModels_[modeli]
            )
        );
    }
    forAll(massTransferModels_, modeli)
    {
        massTransferModel_.set
        (
            modeli,
            massTransferModel::New
            (
                couplingProperties_,
                *this,
                massTransferModels_[modeli]
            )
        );
    }
        Info << "Creating new catchemmodel" << endl;
    forAll(catalyticChemistryModels_, modeli)
    {
        catalyticChemistryModel_.set
        (
            modeli,
            catalyticChemistryModel::New
            (
                couplingProperties_,
                *this,
                catalyticChemistryModels_[modeli]
            )
        );
        Info << "end of constructing cfdemCloudEnergy" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudEnergy::~cfdemCloudEnergy()
{

}

// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

void cfdemCloudEnergy::calcEnergyContributions()
{
    forAll(energyModel_, modeli)
    {
        energyModel_[modeli].calcEnergyContribution();
    }
}

void cfdemCloudEnergy::calcMassTransferContributions()
{
    forAll(massTransferModel_, modeli)
    {
        massTransferModel_[modeli].calcMassTransferContribution();
    }
}

void cfdemCloudEnergy::speciesExecute()
{
    forAll(catalyticChemistryModel_, modeli)
    {
        catalyticChemistryModel_[modeli].execute();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label cfdemCloudEnergy::nrEnergyModels() const
{
    return energyModels_.size();
}

label cfdemCloudEnergy::nrMassTransferModels() const
{
    return massTransferModels_.size();
}

int cfdemCloudEnergy::nrCatalyticChemistryModels()
{
    return catalyticChemistryModels_.size();
}

bool cfdemCloudEnergy::implicitEnergyModel() const
{
    return implicitEnergyModel_;
}

const energyModel& cfdemCloudEnergy::energyM(int i)
{
    return energyModel_[i];
}

const massTransferModel& cfdemCloudEnergy::massTransferM(int i)
{
    return massTransferModel_[i];
}

const catalyticChemistryModel& cfdemCloudEnergy::catalyticChemistryM(int i)
{
    return catalyticChemistryModel_[i];
}

const thermCondModel& cfdemCloudEnergy::thermCondM()
{
    return thermCondModel_;
}

void cfdemCloudEnergy::energyContributions(volScalarField& Qsource)
{
    Qsource.primitiveFieldRef()=0.0;
    Qsource.boundaryFieldRef()=0.0;
    forAll(energyModel_, modeli)
    {
        energyM(modeli).addEnergyContribution(Qsource);
    }
}

void cfdemCloudEnergy::energyCoefficients(volScalarField& Qcoeff)
{
    Qcoeff.primitiveFieldRef()=0.0;
    Qcoeff.boundaryFieldRef()=0.0;
    forAll(energyModel_, modeli)
    {
        energyM(modeli).addEnergyCoefficient(Qcoeff);
    }
}

void cfdemCloudEnergy::massTransferContributions(volScalarField& Msource, label i)
{
    Msource.primitiveFieldRef()=0.0;
    Msource.boundaryFieldRef()=0.0;
    forAll(massTransferModel_, modeli)
    {
        massTransferM(modeli).addMassTransferContribution(Msource, i);
    }
}

tmp<volScalarField> cfdemCloudEnergy::catalyticReactionRate(label i)
{
    tmp<volScalarField> tRRi
    (
        volScalarField::New
        (
            "RR.catalyticChemistry",
            mesh_,
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    );
    forAll(catalyticChemistryModel_, modeli)
    {
        tRRi = tRRi + catalyticChemistryM(modeli).RR(i);
    }
    return tRRi;
}

tmp<volScalarField> cfdemCloudEnergy::catalyticReactionHeat()
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot.catalyticChemistry",
            mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );
    forAll(catalyticChemistryModel_, modeli)
    {
        tQdot = tQdot + catalyticChemistryM(modeli).Qdot();
    }
    return tQdot;
}

bool cfdemCloudEnergy::evolve
(
    volScalarField& alpha,
    volVectorField& Us,
    volVectorField& U
)
{
    if (cfdemCloud::evolve(alpha, Us, U))
    {
        // calc energy contributions
        // position 26 was already defined as Flow in clockModels and RhoPimpleChem solver.
        clockM().start(27,"calcEnergyContributions");
        if(verbose_) Info << "- calcEnergyContributions" << endl;
        calcEnergyContributions();
        if(verbose_) Info << "calcEnergyContributions done." << endl;
        clockM().stop("calcEnergyContributions");

        // execute catalytic chemical model species
        clockM().start(28,"speciesExecute");
        if(verbose_) Info << "- speciesExecute()" << endl;
        speciesExecute();
        if(verbose_) Info << "speciesExecute done" << endl;
        clockM().stop("speciesExecute");

        clockM().start(29,"calcMassTransferContributions");
        if(verbose_) Info << "- calcMassTransferContributions" << endl;
        calcMassTransferContributions();
        if(verbose_) Info << "calcMassTransferContributions done." << endl;
        clockM().stop("calcMassTransferContributions");

        return true;
    }
    return false;
}

void cfdemCloudEnergy::postFlow()
{
    cfdemCloud::postFlow();
    forAll(energyModel_, modeli)
    {
        energyModel_[modeli].postFlow();
    }
    forAll(massTransferModel_, modeli)
    {
        massTransferModel_[modeli].postFlow();
    }
    forAll(catalyticChemistryModel_, modeli)
    {
        catalyticChemistryModel_[modeli].postFlow();
    }
}

void cfdemCloudEnergy::solve()
{
    forAll(energyModel_, modeli)
    {
        energyModel_[modeli].solve();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

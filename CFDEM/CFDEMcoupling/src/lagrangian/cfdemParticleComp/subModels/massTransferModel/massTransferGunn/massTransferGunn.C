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

#include "error.H"
#include "massTransferGunn.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(massTransferGunn, 0);

addToRunTimeSelectionTable(massTransferModel, massTransferGunn, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
massTransferGunn::massTransferGunn
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    massTransferModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    catalystName_(propsDict_.lookup("catalystName")),
    gasName_(propsDict_.lookup("gasName")),
    catThermoPtr_(rhoReactionThermo::New(particleCloud_.mesh(), catalystName_)),
    catChemistryPtr_(basicGSChemistryModel::New(catThermoPtr_())),
    Yg_(dynamic_cast<const rhoReactionThermo&>(sm.thermo()).composition().Y()),
    Tg_(sm.thermo().T()),
    nGasSpecie_(Yg_.size()),
    gasSpecies_(catThermoPtr_->composition().species()),
    ofPartGasFrac_(nGasSpecie_),
    partGasFrac_(NULL),
    massFlux_(nGasSpecie_),
    ReField_
    (   IOobject
        (
            "ReField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    ShField_
    (   IOobject
        (
            "ShField",
            sm.mesh().time().timeName(),
            sm.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sm.mesh(),
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),
    voidfractionFieldName_(propsDict_.lookupOrDefault<word>("voidfractionFieldName","voidfraction")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),
    velFieldName_(propsDict_.lookupOrDefault<word>("velFieldName","U")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
    densityFieldName_(propsDict_.lookupOrDefault<word>("densityFieldName","rho")),
    rho_(sm.mesh().lookupObject<volScalarField> (densityFieldName_)),
    pressureFieldName_(propsDict_.lookupOrDefault<word>("pressureFieldName","p")),
    p_(sm.mesh().lookupObject<volScalarField> (pressureFieldName_)),
    partTempName_(propsDict_.lookupOrDefault<word>("partTempName", "Temp")),
    partTemp_(NULL),
    partMassFlux_(NULL),
    partRe_(NULL),
    partSh_(NULL),
    ofPartMassFlux_(nGasSpecie_),
    partPressName_(propsDict_.lookupOrDefault<word>("partPressName", "p")),
    partPressure_(NULL),
    scaleDia_(1.),
    scaleNusselt_(propsDict_.lookupOrDefault<scalar>("scaleNusselt",1.0)),
    Le_(propsDict_.lookupOrDefault<scalar>("Le",1.0)),
    typeCG_(propsDict_.lookupOrDefault<scalarList>("coarseGrainingFactors",scalarList(1,1.0)))
{
    allocateMyArrays();

    forAll(massFlux_, fieldi)
    {
        massFlux_.set
        (
            fieldi,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "massFlux."+Yg_[fieldi].name(),
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

massTransferGunn::~massTransferGunn()
{
    particleCloud_.dataExchangeM().destroy(partGasFrac_,1);
    particleCloud_.dataExchangeM().destroy(partTemp_,1);
    particleCloud_.dataExchangeM().destroy(partPressure_,1);
    particleCloud_.dataExchangeM().destroy(partMassFlux_,1);
    particleCloud_.dataExchangeM().destroy(partRe_,1);
    particleCloud_.dataExchangeM().destroy(partSh_,1);
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //
void massTransferGunn::allocateMyArrays() const
{
    // get memory for 2d arrays
    double initVal=0.0;
    particleCloud_.dataExchangeM().allocateArray(partTemp_,initVal,1);  // field/initVal/width
    particleCloud_.dataExchangeM().allocateArray(partMassFlux_,initVal,1);
    particleCloud_.dataExchangeM().allocateArray(partPressure_,initVal,1);
    particleCloud_.dataExchangeM().allocateArray(partGasFrac_, initVal, 1);
    forAll(ofPartGasFrac_, i) ofPartGasFrac_[i].setSize(particleCloud_.numberOfParticles(), 0.0);
    forAll(ofPartMassFlux_, i) ofPartMassFlux_[i].setSize(particleCloud_.numberOfParticles(), 0.0);
}

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void massTransferGunn::calcMassTransferContribution()
{
   // realloc the arrays
    allocateMyArrays();

    forAll(ofPartGasFrac_, i)
    {
        particleCloud_.dataExchangeM().getData(gasSpecies_[i], "scalar-atom", partGasFrac_);
        forAll(ofPartGasFrac_[i], index)
        {
            ofPartGasFrac_[i][index] = partGasFrac_[index][0];
        }
    }
    particleCloud_.dataExchangeM().getData(partTempName_,"scalar-atom",partTemp_);
    particleCloud_.dataExchangeM().getData(partPressName_,"scalar-atom",partPressure_);

    if(particleCloud_.cg() > 1.)
    {
        scaleDia_ = particleCloud_.cg();
        Info << "Mass Transfer Gunn is using scale from liggghts cg = " << scaleDia_ << endl;
    }

    const volScalarField mufField = particleCloud_.turbulence().mu();
    const volScalarField CpField = particleCloud_.thermo().Cp();
    const volScalarField kappaField = particleCloud_.thermo().kappa();

    if (typeCG_.size()>1 || typeCG_[0] > 1)
    {
        Info << "massTransferGunn using scale = " << typeCG_ << endl;
    }
    else if (particleCloud_.cg() > 1)
    {
        scaleDia_=particleCloud_.cg();
        Info << "massTransferGunn using scale from liggghts cg = " << scaleDia_ << endl;
    }

    // calc mass flux
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    label cellI=0;
    vector Us(0,0,0);
    scalar ds(0);
    scalar muf(0);
    scalar rhof(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Pr(0);
    scalar Nup(0);
    scalar Shp(0);
    scalarField molefrac_(nGasSpecie_,0.0);
    scalar dt = particleCloud_.mesh().time().deltaTValue();

    for(int index = 0;index < particleCloud_.numberOfParticles(); ++index)
    {
            cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                voidfraction = voidfraction_[cellI];
                Ufluid = U_[cellI];

                if (voidfraction < 0.01)
                    voidfraction = 0.01;

                // calc relative velocity
                Us = particleCloud_.velocity(index);
                magUr = mag(Ufluid - Us);
                ds = 2.*particleCloud_.radius(index);
                muf = mufField[cellI];
                rhof = rho_[cellI];

                Rep = ds * magUr * voidfraction * rhof/ muf;
                Pr = max(SMALL, CpField[cellI] * muf / kappaField[cellI]);
                //scalar D0 = DOField_[celli]; //field implementation of D0 to be done
		//Sc = max(SMALL, muf/(rhof*D0));

                Nup = Nusselt(voidfraction, Rep, Pr);
                Shp = Nup*cbrt(Le_);

                scalar k_lump;
                scalar D_;
                scalar mu_;
                scalar Y_MW = 0;
                scalar mg_ptot0 = 0;
                scalar mg_ptot = 0;
                scalar MWav_ = 0;
                scalar ntot_ = 0;
                scalar R_ = 8.314462618;

                forAll(Yg_, i)
                {
                    Y_MW += ofPartGasFrac_[i][index]/catThermoPtr_->composition().Wi(i);
                }
                forAll(Yg_, i)
                {
                    molefrac_[i] = ofPartGasFrac_[i][index]/catThermoPtr_->composition().Wi(i)/Y_MW;
                    MWav_ += molefrac_[i]*catThermoPtr_->composition().Wi(i);
                }

                forAll(massFlux_, i)
                {
                    mu_ = catChemistryPtr_->mu(i, partTemp_[index][0], p_[cellI]);
                    D_ = mu_/rhof;
                    k_lump = D_ * Shp / ds * (ds * ds * M_PI * catChemistryPtr_->porosity()); //multiply with porosity since only through this fraction of surface mass is transferred?
                    mg_ptot0 = partPressure_[index][0]*catChemistryPtr_->porosity()*particleCloud_.particleVolume(index)/(R_*partTemp_[index][0])*MWav_*1e-3;
                    massFlux(index, i, k_lump, Yg_[i][cellI], ofPartGasFrac_[i][index], mg_ptot0, rhof);
                    ofPartGasFrac_[i][index] = max(ofPartGasFrac_[i][index]*mg_ptot0 + ofPartMassFlux_[i][index]*dt,0.0);
                    mg_ptot += ofPartGasFrac_[i][index];
                }

                forAll(massFlux_, i)
                {
                    ofPartGasFrac_[i][index] /= mg_ptot;
                }

                forAll(Yg_, i)
                {
                    ntot_ += ofPartGasFrac_[i][index]*mg_ptot/(catThermoPtr_->composition().Wi(i)*1e-3);
                }

                forAll(Yg_, i)
                {
                    partPressure_[index][0] = ntot_*R_*partTemp_[index][0]/(catChemistryPtr_->porosity()*particleCloud_.particleVolume(index));
                }
            }
            else
            {
                partPressure_[index][0] = 0;
            }
    }

    forAll(ofPartGasFrac_, i)
    {
        forAll(ofPartGasFrac_[i], index)
        {
            partGasFrac_[index][0] = particleCloud_.cellIDs()[index][0] >= 0 ? ofPartGasFrac_[i][index] : 0.0;
        }
        particleCloud_.dataExchangeM().giveData(gasSpecies_[i], "scalar-atom", partGasFrac_);
    }

    particleCloud_.dataExchangeM().giveData(partPressName_, "scalar-atom", partPressure_);

    forAll(massFlux_, i)
    {
        forAll(ofPartMassFlux_[i], index)
        {
            partMassFlux_[index][0] = ofPartMassFlux_[i][index];
        }
        massFlux_[i].primitiveFieldRef() = 0.0;
        particleCloud_.averagingM().setScalarSum
        (
            massFlux_[i],
            partMassFlux_,
            particleCloud_.particleWeights(),
            NULL
        );
        massFlux_[i].primitiveFieldRef() /= -massFlux_[i].mesh().V();
    }

    forAll(massFlux_, i)
    {
        massFlux_[i].correctBoundaryConditions();
    }
}

void massTransferGunn::addMassTransferContribution(volScalarField& Msource, label i) const
{
    Msource += massFlux_[i];
}

scalar massTransferGunn::Nusselt(scalar voidfraction, scalar Rep, scalar Pr) const
{
  return scaleNusselt_*((7 - 10 * voidfraction + 5 * voidfraction * voidfraction) *
                      (1 + 0.7 * Foam::pow(Rep,0.2) * Foam::pow(Pr,0.33)) +
                      (1.33 - 2.4 * voidfraction + 1.2 * voidfraction * voidfraction) *
                      Foam::pow(Rep,0.7) * Foam::pow(Pr,0.33));
}

void massTransferGunn::massFlux(label index, label i, scalar k, scalar Yg, scalar Ygs, scalar mg_ptot0, scalar rhof)
{
    ofPartMassFlux_[i][index] = - k * Ygs * mg_ptot0/(catChemistryPtr_->porosity()*particleCloud_.particleVolume(index));
    ofPartMassFlux_[i][index] += k * Yg * rhof;
}

void massTransferGunn::giveData()
{
    Info << "total convective particle-fluid mass flux [kg/s] (Eulerian) for first specie = " << gSum(massFlux_[0]*1.0*massFlux_[0].mesh().V()) << endl;
    Info << "net convective particle-fluid mass flux [kg/s] (Eulerian) for first specie = " << gSum((massFlux_[0]*1.0+massFlux_[0]*Yg_[0])*massFlux_[0].mesh().V()) << endl;
}

void massTransferGunn::postFlow()
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

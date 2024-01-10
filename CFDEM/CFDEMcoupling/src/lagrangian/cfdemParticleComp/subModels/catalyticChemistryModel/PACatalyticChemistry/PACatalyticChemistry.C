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
#include "PACatalyticChemistry.H"
#include "addToRunTimeSelectionTable.H"

#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PACatalyticChemistry, 0);

addToRunTimeSelectionTable(catalyticChemistryModel, PACatalyticChemistry, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
PACatalyticChemistry::PACatalyticChemistry
(
    const dictionary& dict,
    cfdemCloudEnergy& sm
)
:
    catalyticChemistryModel(dict, sm),
    propsDict_(dict.subDict(typeName + "Props")),
    catalystName_(propsDict_.lookup("catalystName")),
    PAtol_(propsDict_.lookupOrDefault<scalar>("tolerance", 1e-3)),
    PAtolp_(propsDict_.lookupOrDefault<scalar>("pressure_tolerance", 1e-3)),
    nBins_(propsDict_.lookup<scalar>("nBins")),
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
    ofPartRRav_(nGasSpecie_),
    partCoverages_(NULL),
    ofPartCoverages_(nSolidSpecie_),
    ofPartCoveragesav_(nSolidSpecie_),
    partGas_(NULL),
    ofPartGas_(nGasSpecie_),
    partPressure_(NULL),
    partTempName_(propsDict_.lookupOrDefault<word>("partTempName", "Temp")),
    partTemp_(NULL),
    partQdot_(NULL),
    partHeatSourceName_(propsDict_.lookupOrDefault<word>("partHeatSourceName", "heatSource")),
    partHeatSource_(NULL),
    useParticleQdot_(propsDict_.lookupOrDefault("useParticleQdot", true)),
    dcdt_(nGasSpecie_, 0.0),
    Yi_(nGasSpecie_, 0.0),
    Ysi_(nSolidSpecie_, 0.0),
    nVar_(nGasSpecie_ + nSolidSpecie_ + 5),
    ini_(propsDict_.lookupOrDefault("initialized", false)),
    partPressName_(propsDict_.lookupOrDefault<word>("partPressName", "p")),
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

PACatalyticChemistry::~PACatalyticChemistry()
{
	if(massTransfer_)
    {
        delete partGas_;
        delete partPressure_;
    };
    
    delete partTemp_;
    delete partQdot_;
    delete partRR_;
    delete partCoverages_;
}

// * * * * * * * * * * * * * * * private Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Member Fct  * * * * * * * * * * * * * * * //

void PACatalyticChemistry::allocateMyArrays() const
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
        particleCloud_.dataExchangeM().allocateArray(partPressure_, initVal, 1);
        forAll(ofPartGas_, i) ofPartGas_[i].setSize(particleCloud_.numberOfParticles(), Ygs_[i][0]);
        particleCloud_.dataExchangeM().allocateArray(partGas_, initVal, 1);
    };

    if (useParticleQdot_) particleCloud_.dataExchangeM().giveData(partHeatSourceName_, "scalar-atom", partHeatSource_);
}

void PACatalyticChemistry::execute()
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
		allocateMyArrays();

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

        for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
        {
            label cellI = particleCloud_.cellIDs()[index][0];
            if(cellI >= 0)
            {
                for (label i=0; i<nSolidSpecie_; i++) ofPartCoverages_[i][index] = Ysi_[i];
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
            particleCloud_.dataExchangeM().getData(partPressName_, "scalar-atom", partPressure_);
        };
        
        Info << "Solid species after QSSA: " << Ysi_ << endl;


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
            particleCloud_.dataExchangeM().getData(partPressName_, "scalar-atom", partPressure_);
        };

        if (!useFluidTemperature_)
        {
            particleCloud_.dataExchangeM().getData(partTempName_, "scalar-atom", partTemp_);
        }
    }

    if (!useFluidTemperature_)
    {
        particleCloud_.dataExchangeM().getData(partTempName_, "scalar-atom", partTemp_);
    }

    List<scalarList> phi_(particleCloud_.numberOfParticles());
    forAll(phi_, i) phi_[i].setSize(nVar_, 0.0);
    scalarField iPA_(particleCloud_.numberOfParticles(), 0.0);

    scalarField minPhi_(nVar_, 0.0);
    scalarField maxPhi_(nVar_, 2e-20);

    scalar si;
    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        label cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
			scalar pi;
            if(!massTransfer_)
            {
                for (label i=0; i<nGasSpecie_; i++) Yi_[i] = Ygs_[i][cellI];
                pi = p_[cellI];
            }
            else
            {
                for (label i=0; i<nGasSpecie_; i++) Yi_[i] = ofPartGas_[i][index];
                pi = partPressure_[index][0];
            }
            for (label i=0; i<nGasSpecie_; i++)
            {
                phi_[index][i] = Yi_[i];
                minPhi_[i] = max(min(phi_[index][i], minPhi_[i]),1e-20);
                maxPhi_[i] = max(phi_[index][i], maxPhi_[i]);
            }
            for (label i=0; i<nSolidSpecie_; i++)
            {
                si = i + nGasSpecie_;
                Ysi_[i] = ofPartCoverages_[i][index];
                phi_[index][si] = Ysi_[i];
                minPhi_[si] = max(min(phi_[index][si], minPhi_[si]),1e-20);
                maxPhi_[si] = max(phi_[index][si], maxPhi_[si]);
            }
            phi_[index][nVar_ - 5 ] = particleCloud_.particleVolume(index);
            minPhi_[nVar_ - 5 ] = max(min(phi_[index][nVar_ - 5 ], minPhi_[nVar_ - 5 ]),1e-20);
            maxPhi_[nVar_ - 5 ] = max(phi_[index][nVar_ - 5 ], maxPhi_[nVar_ - 5 ]);
            scalar Tgi = Tg_[cellI];
            phi_[index][nVar_ - 4 ] = Tgi;
            minPhi_[nVar_ - 4 ] = max(min(phi_[index][nVar_ - 4 ], minPhi_[nVar_ - 4 ]),1e-20);
            maxPhi_[nVar_ - 4 ] = max(phi_[index][nVar_ - 4 ], maxPhi_[nVar_ - 4 ]);
            if(useFluidTemperature_)
            {
                phi_[index][nVar_ - 3 ] = Tgi;
                minPhi_[nVar_ - 3 ] = max(min(phi_[index][nVar_ - 3 ], minPhi_[nVar_ - 3 ]),1e-20);
                maxPhi_[nVar_ - 3 ] = max(phi_[index][nVar_ - 3 ], maxPhi_[nVar_ - 3 ]);
            }
            else
            {
                phi_[index][nVar_ - 3 ] = partTemp_[index][0];
                minPhi_[nVar_ - 3 ] = max(min(phi_[index][nVar_ - 3 ], minPhi_[nVar_ - 3 ]),1e-20);
                maxPhi_[nVar_ - 3 ] = max(phi_[index][nVar_ - 3 ], maxPhi_[nVar_ - 3 ]);
            }
			phi_[index][nVar_ - 2 ] = pi;
			minPhi_[nVar_ - 2 ] = max(min(phi_[index][nVar_ - 2 ], minPhi_[nVar_ - 2 ]),1e-20);
            maxPhi_[nVar_ - 2 ] = max(phi_[index][nVar_ - 2 ], maxPhi_[nVar_ - 2 ]);
            phi_[index][nVar_ - 1 ] = rhog_[cellI];
            minPhi_[nVar_ - 1 ] = max(min(phi_[index][nVar_ - 1 ], minPhi_[nVar_ - 1 ]),1e-20);
            maxPhi_[nVar_ - 1 ] = max(phi_[index][nVar_ - 1 ], maxPhi_[nVar_ - 1 ]);
        }
    }

    HashTable<scalarList,label> phiTable_;
    phiTable_.clear();
    HashTable<label,label> partIndexTable_;
    partIndexTable_.clear();
    HashTable<label,label> partiPATable_;
    partiPATable_.clear();

    scalar iPAmin_ = GREAT;
    scalar iPAmax_ = SMALL;
    
    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        label cellI = particleCloud_.cellIDs()[index][0];
        if(cellI >= 0)
        {
            for (label i=0; i < nVar_ ;i++)
            {
				if(massTransfer_ && i == nVar_ -2)
				{
					iPA_[index] += (i+1)*std::log(nVar_/PAtolp_*(std::log(max(phi_[index][i],2e-20))-std::log(minPhi_[i]))/(std::log(maxPhi_[i])-std::log(minPhi_[i]))); //hash mapping function
				}
                else
                {
					iPA_[index] += (i+1)*std::log(nVar_/PAtol_*(std::log(max(phi_[index][i],2e-20))-std::log(minPhi_[i]))/(std::log(maxPhi_[i])-std::log(minPhi_[i]))); //hash mapping function
				}
            }
            iPAmin_ = floor(min(iPA_[index], iPAmin_));
            iPAmax_ = ceil(max(iPA_[index], iPAmax_));
            phiTable_.insert(iPA_[index],phi_[index]);
            partIndexTable_.insert(iPA_[index],index);
            partiPATable_.insert(index,iPA_[index]);
        }
    }
    
    scalarField binCount_(nBins_,0.0);
    scalarField phiBin_(nVar_,0.0);

    HashTable<label,label> partBinTable_;
    partBinTable_.clear();

    List<scalarField> Yiav_(nBins_);
    List<scalarField> Ysiav_(nBins_);
    List<scalarField> YsiavOld_(nBins_);
    scalarField partVolav_(nBins_,0.0);
    scalarField Tgav_(nBins_,0.0);
    scalarField Tsav_(nBins_,0.0);
    scalarField pav_(nBins_,0.0);
    scalarField rhogav_(nBins_,0.0);

    forAll(Yiav_, i) Yiav_[i].setSize(nGasSpecie_, 0.0);
    forAll(Ysiav_, i) Ysiav_[i].setSize(nSolidSpecie_, 0.0);
    forAll(YsiavOld_, i) YsiavOld_[i].setSize(nSolidSpecie_, 0.0);
    forAll(ofPartRRav_, i) ofPartRRav_[i].setSize(nBins_, 0.0);
    forAll(ofPartCoveragesav_, i) ofPartCoveragesav_[i].setSize(nBins_, 0.0);

	scalar binNumber_ = 0;
	
    for(int index = 0; index < particleCloud_.numberOfParticles(); ++index)
    {
        if(particleCloud_.cellIDs()[index][0] >= 0)
        {
            scalar partiPA_ = partiPATable_.find(index)();
            binNumber_ = floor(nBins_*(partiPA_-iPAmin_)/(max(iPAmax_-iPAmin_,SMALL)));
            phiBin_ = phiTable_.find(partiPA_)();
            partBinTable_.insert(index,binNumber_);
            binCount_[binNumber_]++;
            for(label k=0; k < nGasSpecie_; k++)
            {
                Yiav_[binNumber_][k] += phiBin_[k];
            }
            for(label k=0; k < nSolidSpecie_; k++)
            {
                Ysiav_[binNumber_][k] += phiBin_[k + nGasSpecie_];
            }
            partVolav_[binNumber_] += phiBin_[nVar_ - 5];
            Tgav_[binNumber_] += phiBin_[nVar_ - 4];
            Tsav_[binNumber_] += phiBin_[nVar_ - 3];
            pav_[binNumber_] += phiBin_[nVar_ - 2];
            rhogav_[binNumber_] += phiBin_[nVar_ - 1];
        }
    }

    for(int i = 0; i < nBins_; i++)
    {
        if(binCount_[i] > 0)
        {
            Yiav_[i] /= binCount_[i];
            Ysiav_[i] /= binCount_[i];
            partVolav_[i] /= binCount_[i];
            Tgav_[i] /= binCount_[i];
            Tsav_[i] /= binCount_[i];
            pav_[i] /= binCount_[i];
            rhogav_[i] /= binCount_[i];
        }
    }
    
    YsiavOld_ = Ysiav_;
    scalarField pavOld_(nBins_,0.0);
    List<scalarField> YiavOld_(nBins_);
    
    if(massTransfer_)
    {
		forAll(YiavOld_, i) YiavOld_[i].setSize(nSolidSpecie_, 0.0);
		YiavOld_ = Yiav_;
		pavOld_ = pav_;
	}

    scalarField partQdotav_(nBins_,0.0);
    for(int index = 0; index < nBins_; index++)
    {
        if(binCount_[index] > 0)
        {
            partQdotav_[index] =
                catChemistryPtr_->getRatesQdotI
                (
                    particleCloud_.mesh().time().deltaTValue(),
                    1,
                    dcdt_,
                    Yiav_[index],
                    Ysiav_[index],
                    rhogav_[index],
                    Tgav_[index],
                    Tsav_[index],
                    pav_[index],
                    false, // initialize coverages
                    !useParticleQdot_ // if particleQdot = true, use absolute energy (sensible energy = false)
              );
            partQdotav_[index] *= partVolav_[index];
            
            for (label i=0; i<nGasSpecie_; i++) ofPartRRav_[i][index] = dcdt_[i]*partVolav_[index]; //OK FOR UNIFORM PARTICLE DIAMETER
            for (label i=0; i<nSolidSpecie_; i++) ofPartCoveragesav_[i][index] = Ysiav_[index][i];
        }
    }
    
    if(massTransfer_)
    {
		scalar dt = particleCloud_.mesh().time().deltaTValue();
		
		for(int index = 0; index < nBins_; index++)
        {
			if(binCount_[index] > 0)
			{
				scalar RRnet_ = 0.0;
				scalar R_ = 8.314462618;
				scalar T_;
				scalar Y_MW = 0;
				if(useFluidTemperature_)
				{
					T_ = Tgav_[index];
				}
				else
				{
					T_ = Tsav_[index];
				}
            
				forAll(Yg_, i)
				{
					Y_MW += Yiav_[index][i]/catThermoPtr_->composition().Wi(i);
				}
            
				forAll(Yg_, i)
				{
					Yiav_[index][i] = Yiav_[index][i]/catThermoPtr_->composition().Wi(i)/Y_MW;
				}

				forAll(Yg_, i)
				{
					RRnet_  += ofPartRRav_[i][index]*dt/(catThermoPtr_->composition().Wi(i)*1e-3); //in mol/s
                
					Yiav_[index][i] = max(Yiav_[index][i]*pav_[index]/(R_*T_)*catChemistryPtr_->porosity()*partVolav_[index] + ofPartRRav_[i][index]*dt/(catThermoPtr_->composition().Wi(i)*1e-3),0.0);
               	}
            
				scalar Mt_ = 0.0;
				scalar nt_ = 0.0;
				forAll(Yg_, i)
				{
					Mt_ += Yiav_[index][i]*catThermoPtr_->composition().Wi(i);
					nt_ += Yiav_[index][i];
				}
            
				forAll(Yg_, i)
				{
					Yiav_[index][i] = Yiav_[index][i]*catThermoPtr_->composition().Wi(i)/Mt_;
				}
				
				pav_[index] += RRnet_*R_*T_/(catChemistryPtr_->porosity()*partVolav_[index]);
            }
		}
		
		scalar Ygtot;
        for(int index = 0; index < particleCloud_.numberOfParticles(); index++)
        {
			Ygtot = 0.0;
			if(particleCloud_.cellIDs()[index][0] >= 0)
			{
				
				scalar binNumberFound_ = partBinTable_.find(index)();

				for(label k=0; k < nGasSpecie_; k++)
				{
					ofPartGas_[k][index] = max(ofPartGas_[k][index] +  Yiav_[binNumberFound_][k] - YiavOld_[binNumberFound_][k],0.0);
					Ygtot += ofPartGas_[k][index];
				}
				for(label k=0; k < nGasSpecie_; k++)
				{
					ofPartGas_[k][index] /= Ygtot;
				}
				partPressure_[index][0] = max(partPressure_[index][0] + pav_[binNumberFound_] - pavOld_[binNumberFound_],0.0);
				if(index==0) Info << "partPressure_: " << partPressure_[index][0] << endl;
			}
			else
            {
                partPressure_[index][0] = 0;
            }
		}
		
		forAll(ofPartGas_, i)
        {
            forAll(ofPartGas_[i], index)
            {
                partGas_[index][0] = particleCloud_.cellIDs()[index][0] >= 0 ? ofPartGas_[i][index] : 0.0;
            }
            particleCloud_.dataExchangeM().giveData(gasSpecies_[i], "scalar-atom", partGas_);
        }
        particleCloud_.dataExchangeM().giveData(partPressName_, "scalar-atom", partPressure_);
	}

    scalar Ystot;
    for(int index = 0; index < particleCloud_.numberOfParticles(); index++)
    {
		Ystot = 0.0;
        if(particleCloud_.cellIDs()[index][0] >= 0)
        {
            scalar binNumberFound_ = partBinTable_.find(index)();
            for(label k=0; k < nGasSpecie_; k++)
            {
                ofPartRR_[k][index] = ofPartRRav_[k][binNumberFound_];
            }
            for(label k=0; k < nSolidSpecie_; k++)
            {
                ofPartCoverages_[k][index] = max(ofPartCoverages_[k][index] +  ofPartCoveragesav_[k][binNumberFound_] - YsiavOld_[binNumberFound_][k],0.0);
                if(!massTransfer_) Ystot += ofPartCoverages_[k][index];
            }
            Ystot = max(SMALL, Ystot);
            if(!massTransfer_)
            {
				for(label k=0; k < nSolidSpecie_; k++)
				{
					ofPartCoverages_[k][index] /= Ystot;
				}
			}
            partQdot_[index][0] = partQdotav_[binNumberFound_];
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

    if(!massTransfer_)
    {
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
    else
    {
        forAll(RR_, i)
        {
            RR_[i].primitiveFieldRef() = 0.0;
        }
    }
}

const volScalarField PACatalyticChemistry::RR(label i) const
{
    return RR_[i];
}

const volScalarField PACatalyticChemistry::Qdot() const
{
    return Qdot_;
}

void PACatalyticChemistry::postFlow()
{
    forAll(ofPartCoverages_, i)
    {
        forAll(ofPartCoverages_[i], index)
        {
            partCoverages_[index][0] = particleCloud_.cellIDs()[index][0] >= 0 ? ofPartCoverages_[i][index] : 0.0;
        }
        particleCloud_.dataExchangeM().giveData(solidSpecies_[i], "scalar-atom", partCoverages_);
    }

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

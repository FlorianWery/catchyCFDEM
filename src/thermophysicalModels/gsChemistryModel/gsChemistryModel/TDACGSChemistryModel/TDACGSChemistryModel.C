/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "TDACGSChemistryModel.H"
#include "multiComponentGSMixture.H"
#include "UniformField.H"
#include "localEulerDdtScheme.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::TDACGSChemistryModel
(
    rhoReactionThermo& thermo
)
:
    StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>(thermo),
    variableTimeStep_
    (
        this->mesh().time().controlDict().lookupOrDefault
        (
            "adjustTimeStep",
            false
        )
     || fv::localEulerDdt::enabled(this->mesh())
    ),
    timeSteps_(0),
    Treact_
    (
        ReactionThermo::template lookupOrDefault<scalar>
        (
            "Treact",
            0
        )
    ),
    NsDAC_(this->nSpecie_),
    completeC_(this->nSpecie_, 0),
    gasReactionsDisabled_(this->nGasReaction_, false),
    solidReactionsDisabled_(this->nSolidReaction_, false),
    specieComp_(this->nSpecie_),
    completeToSimplifiedIndex_(this->nSpecie_, -1),
    simplifiedToCompleteIndex_(this->nSpecie_),
    tabulationResults_
    (
        IOobject
        (
            thermo.phasePropertyName("TabulationResults"),
            this->time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        scalar(0)
    )
{
    const basicSpecieMixture& composition = this->thermo().composition();

    const HashTable<List<specieElement>>& specCompG =
        dynamicCast<const multiComponentGSMixture<SThermoType, GThermoType>&>(this->thermo())
       .specieComposition();

    const HashTable<List<specieElement>>& specCompS =
        dynamicCast<const multiComponentGSMixture<SThermoType, GThermoType>&>(this->thermo())
       .solidSpecieComposition();

    forAll(specieComp_, i)
    {
        if (i < this->nGasSpecie_)
        {
            specieComp_[i] = specCompG[this->Y()[i].member()];
        }
        else
        {
            specieComp_[i] = specCompS[this->Ys()[i-this->nGasSpecie_].member()];
        }

    }

    mechRed_ = chemistryGSReductionMethod<ReactionThermo, GThermoType, SThermoType>::New
    (
        *this,
        *this
    );

    // When the mechanism reduction method is used, the 'active' flag for every
    // species should be initialized (by default 'active' is true)
    if (mechRed_->active())
    {
        forAll(this->Y(), i)
        {
            IOobject header
            (
                this->Y()[i].name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ
            );

            // Check if the species file is provided, if not set inactive
            // and NO_WRITE
            if (!header.typeHeaderOk<volScalarField>(true))
            {
                composition.setInactive(i);
            }
        }

        forAll(this->Ys(), i)
        {
            IOobject header
            (
                this->Ys()[i].name(),
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ
            );

            // Check if the species file is provided, if not set inactive
            // and NO_WRITE
            if (!header.typeHeaderOk<volScalarField>(true))
            {
                composition.setInactive(i+this->nGasSpecie_);
            }
        }
    }

    tabulation_ = chemistryGSTabulationMethod<ReactionThermo, GThermoType, SThermoType>::New
    (
        *this,
        *this
    );

    if (mechRed_->log())
    {
        cpuReduceFile_ = logFile("cpu_reduce.out");
        nActiveSpeciesFile_ = logFile("nActiveSpecies.out");
    }

    if (tabulation_->log())
    {
        cpuAddFile_ = logFile("cpu_add.out");
        cpuGrowFile_ = logFile("cpu_grow.out");
        cpuRetrieveFile_ = logFile("cpu_retrieve.out");
    }

    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_ = logFile("cpu_solve.out");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::~TDACGSChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::omega
(
    const scalar p,
    const scalar T,
    const scalarField& c, // Contains all species even when mechRed is active
    const label li,
    scalarField& dcdt
) const
{
    const bool reduced = mechRed_->active();

    dcdt = Zero;

    forAll(this->solidReactions_, i)
    {
      if (!solidReactionsDisabled_[i])
      {
            const SurfaceReaction<SThermoType>& R = this->solidReactions_[i];
            R.omega(p, T, c, li, dcdt, reduced, completeToSimplifiedIndex_);
      }
    }

    for(label i = 0; i < this->nGasSpecie_; i++)
    {
        dcdt[i] *= this->catarea_*this->multiplyAV_;
    }

    if (this->por_*this->multiplyPor_ > 0.0)
    {
        scalarField dcdtg(this->nGasSpecie_+2, 0.0);

        forAll(this->gasReactions_, i)
        {
            if (!gasReactionsDisabled_[i])
            {
                const Reaction<GThermoType>& R = this->gasReactions_[i];

                R.omega(p, T, c, li, dcdtg, reduced, completeToSimplifiedIndex_);
            }
        }

        for(label i = 0; i < this->nGasSpecie_; i++)
        {
            dcdt[i] += dcdtg[i]*this->por_*this->multiplyPor_;
        }
    }

    /*forAll(this->gasReactions_, i)
    {
        if (!reactionsDisabled_[i])
        {
            const Reaction<GThermoType>& R = this->gasReactions_[i];

            scalar omegai = R.omega
            (
                 p, T, c, li, pf, cf, lRef, pr, cr, rRef
            );

            forAll(R.lhs(), s)
            {
                label si = R.lhs()[s].index;
                if (reduced)
                {
                    si = completeToSimplifiedIndex_[si];
                }

                const scalar sl = R.lhs()[s].stoichCoeff;
                dcdt[si] -= sl*omegai;
            }
            forAll(R.rhs(), s)
            {
                label si = R.rhs()[s].index;
                if (reduced)
                {
                    si = completeToSimplifiedIndex_[si];
                }

                const scalar sr = R.rhs()[s].stoichCoeff;
                dcdt[si] += sr*omegai;
            }
        }
    }

    forAll(this->solidReactions_, i)
    {
        if (!reactionsDisabled_[i])
        {
            const SurfaceReaction<SThermoType>& R = this->solidReactions_[i];

            scalar omegai = R.omega
            (
                 p, T, c, li, pf, cf, lRef, pr, cr, rRef
            );

            forAll(R.lhs(), s)
            {
                label si = R.lhs()[s].index;
                if (reduced)
                {
                    si = completeToSimplifiedIndex_[si];
                }

                const scalar sl = R.lhs()[s].stoichCoeff;
                dcdt[si] -= sl*omegai;
            }
            forAll(R.rhs(), s)
            {
                label si = R.rhs()[s].index;
                if (reduced)
                {
                    si = completeToSimplifiedIndex_[si];
                }

                const scalar sr = R.rhs()[s].stoichCoeff;
                dcdt[si] += sr*omegai;
            }
        }
    }*/
}

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::omega
(
    const Reaction<GThermoType>& R,
    const scalar p,
    const scalar T,
    const scalarField& c, // Contains all species even when mechRed is active
    const label li,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const scalar kf = R.kf(p, T, c, li);
    const scalar kr = R.kr(kf, p, T, c, li);

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s=1; s<Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(c[lRef], 0), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(c[si], 0), exp);
        }
    }
    cf = max(c[lRef], 0);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1)
        {
            if (cf > small)
            {
                pf *= pow(cf, exp - 1);
            }
            else
            {
                pf = 0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // Find the matrix element and element position for the rhs
    pr = kr;
    for (label s=1; s<Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(c[rRef], 0), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(c[si], 0), exp);
        }
    }
    cr = max(c[rRef], 0);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1)
        {
            if (cr > small)
            {
                pr *= pow(cr, exp - 1);
            }
            else
            {
                pr = 0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1);
        }
    }

    return (pf*cf - pr*cr)*this->por_*this->multiplyPor_;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::omega
(
    const SurfaceReaction<SThermoType>& R,
    const scalar p,
    const scalar T,
    const scalarField& c, // Contains all species even when mechRed is active
    const label li,
    scalar& pf,
    scalar& cf,
    label& lRef,
    scalar& pr,
    scalar& cr,
    label& rRef
) const
{
    const scalar kf = R.kf(p, T, c, li);
    const scalar kr = R.kr(kf, p, T, c, li);

    const label Nl = R.lhs().size();
    const label Nr = R.rhs().size();

    label slRef = 0;
    lRef = R.lhs()[slRef].index;

    pf = kf;
    for (label s=1; s<Nl; s++)
    {
        const label si = R.lhs()[s].index;

        if (c[si] < c[lRef])
        {
            const scalar exp = R.lhs()[slRef].exponent;
            pf *= pow(max(c[lRef], 0), exp);
            lRef = si;
            slRef = s;
        }
        else
        {
            const scalar exp = R.lhs()[s].exponent;
            pf *= pow(max(c[si], 0), exp);
        }
    }
    cf = max(c[lRef], 0);

    {
        const scalar exp = R.lhs()[slRef].exponent;
        if (exp < 1)
        {
            if (cf > small)
            {
                pf *= pow(cf, exp - 1);
            }
            else
            {
                pf = 0;
            }
        }
        else
        {
            pf *= pow(cf, exp - 1);
        }
    }

    label srRef = 0;
    rRef = R.rhs()[srRef].index;

    // Find the matrix element and element position for the rhs
    pr = kr;
    for (label s=1; s<Nr; s++)
    {
        const label si = R.rhs()[s].index;
        if (c[si] < c[rRef])
        {
            const scalar exp = R.rhs()[srRef].exponent;
            pr *= pow(max(c[rRef], 0), exp);
            rRef = si;
            srRef = s;
        }
        else
        {
            const scalar exp = R.rhs()[s].exponent;
            pr *= pow(max(c[si], 0), exp);
        }
    }
    cr = max(c[rRef], 0);

    {
        const scalar exp = R.rhs()[srRef].exponent;
        if (exp < 1)
        {
            if (cr > small)
            {
                pr *= pow(cr, exp - 1);
            }
            else
            {
                pr = 0;
            }
        }
        else
        {
            pr *= pow(cr, exp - 1);
        }
    }

    return (pf*cf - pr*cr)*this->catarea_*this->multiplyAV_;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::derivatives
(
    const scalar time,
    const scalarField& c,
    const label li,
    scalarField& dcdt
) const
{
    const bool reduced = mechRed_->active();

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    if (reduced)
    {
        // When using DAC, the ODE solver submit a reduced set of species the
        // complete set is used and only the species in the simplified mechanism
        // are updated
        this->c_ = completeC_;

        // Update the concentration of the species in the simplified mechanism
        // the other species remain the same and are used only for third-body
        // efficiencies
        for (label i=0; i<NsDAC_; i++)
        {
            this->c_[simplifiedToCompleteIndex_[i]] = max(c[i], 0);
        }
    }
    else
    {
        for (label i=0; i<this->nSpecie(); i++)
        {
            this->c_[i] = max(c[i], 0);
        }
    }

    dcdt = Zero;
    omega(p, T, this->c_, li, dcdt);

    // Constant pressure
    // dT/dt = ...
    /*scalar rho = 0;
    for (label i=0; i<this->c_.size(); i++)
    {
        const scalar W = this->specieThermo_[i].W();
        rho += W*this->c_[i];
    }

    scalar cp = 0;
    for (label i=0; i<this->c_.size(); i++)
    {
        // cp function returns [J/kmol/K]
        cp += this->c_[i]*this->specieThermo_[i].cp(p, T);
    }
    cp /= rho;

    // When mechanism reduction is active
    // dT is computed on the reduced set since dcdt is null
    // for species not involved in the simplified mechanism
    scalar dT = 0;
    for (label i=0; i<this->nSpecie_; i++)
    {
        label si;
        if (reduced)
        {
            si = simplifiedToCompleteIndex_[i];
        }
        else
        {
            si = i;
        }

        // ha function returns [J/kmol]
        const scalar hi = this->specieThermo_[si].ha(p, T);
        dT += hi*dcdt[i];
    }
    dT /= rho*cp;

    dcdt[this->nSpecie_] = -dT;

    // dp/dt = ...
    dcdt[this->nSpecie_ + 1] = 0;*/
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::jacobian
(
    const scalar t,
    const scalarField& c,
    const label li,
    scalarField& dcdt,
    scalarSquareMatrix& J
) const
{
    const bool reduced = mechRed_->active();

    // If the mechanism reduction is active, the computed Jacobian
    // is compact (size of the reduced set of species)
    // but according to the information of the complete set
    // (i.e. for the third-body efficiencies)

    const scalar T = c[this->nSpecie_];
    const scalar p = c[this->nSpecie_ + 1];

    if (reduced)
    {
        this->c_ = completeC_;

        for (label i=0; i<NsDAC_; i++)
        {
            this->c_[simplifiedToCompleteIndex_[i]] = max(c[i], 0);
        }
    }
    else
    {
        forAll(this->c_, i)
        {
            this->c_[i] = max(c[i], 0);
        }
    }

    //scalarField hi(this->c_.size());
    //scalarField cpi(this->c_.size());
    /*forAll(hi, i)
    {
        hi[i] = this->specieThermo_[i].ha(p, T);
        cpi[i] = this->specieThermo_[i].cp(p, T);
    }*/

    J = Zero;
    dcdt = Zero;

    scalar omegaI = 0;
    //List<label> dummy;

    forAll(this->solidReactions_, sri)
    {
        if (!solidReactionsDisabled_[sri])
        {
            const SurfaceReaction<SThermoType>& R = this->solidReactions_[sri];
            scalar kfwd, kbwd;
            R.dwdc(p, T, this->c_, li, J, dcdt, omegaI, kfwd, kbwd, reduced, completeToSimplifiedIndex_);
        }
    }

    for(label j = 0; j < this->nGasSpecie_; j++)
    {
        for(label i = 0; i < this->nSpecie_+2; i++)
        {
            J(j,i) *= this->catarea_*this->multiplyAV_;
        }
        dcdt[j] *= this->catarea_*this->multiplyAV_;
    }

    if (this->por_*this->multiplyPor_>0.0)
    {
        scalarSquareMatrix Jg(this->nGasSpecie_+2, 0.0);
        scalarField dcdtg(this->nGasSpecie_+2, 0.0);

        forAll(this->gasReactions_, gri)
        {
          if (!gasReactionsDisabled_[gri])
          {
              const Reaction<GThermoType>& R = this->gasReactions_[gri];
              scalar kfwd, kbwd;
              R.dwdc(p, T, this->c_, li, Jg, dcdtg, omegaI, kfwd, kbwd, reduced, completeToSimplifiedIndex_);
          }
        }

        for(label j = 0; j < this->nGasSpecie_; j++)
        {
            for(label i = 0; i < this->nGasSpecie_; i++)
            {
                J(j,i) += Jg(j,i)*this->por_*this->multiplyPor_;
            }
            dcdt[j] += dcdtg[j]*this->por_*this->multiplyPor_;

            J(j,this->nSpecie_) += Jg(j,this->nGasSpecie_)*this->por_*this->multiplyPor_;
        }
    }

    // The species derivatives of the temperature term are partially computed
    // while computing dwdc, they are completed hereunder:
    /*scalar cpMean = 0;
    scalar dcpdTMean = 0;
    forAll(this->c_, i)
    {
        cpMean += this->c_[i]*cpi[i]; // J/(m^3 K)
        // Already multiplied by rho
        dcpdTMean += this->c_[i]*this->specieThermo_[i].dcpdT(p, T);
    }

    scalar dTdt = 0;
    forAll(hi, i)
    {
        if (reduced)
        {
            const label si = completeToSimplifiedIndex_[i];
            if (si != -1)
            {
                dTdt += hi[i]*dcdt[si]; // J/(m^3 s)
            }
        }
        else
        {
            dTdt += hi[i]*dcdt[i]; // J/(m^3 s)
        }
    }
    dTdt /= -cpMean; // K/s
    dcdt[this->nSpecie_] = dTdt;*/

    /*for (label i = 0; i < this->nSpecie_; i++)
    {
        J(this->nSpecie_, i) = 0;
        for (label j = 0; j < this->nSpecie_; j++)
        {
            const label sj = reduced ? simplifiedToCompleteIndex_[j] : j;
            J(this->nSpecie_, i) += hi[sj]*J(j, i);
        }
        const label si = reduced ? simplifiedToCompleteIndex_[i] : i;
        J(this->nSpecie_, i) += cpi[si]*dTdt; // J/(mol s)
        J(this->nSpecie_, i) /= -cpMean;    // K/s / (mol/m^3)
    }*/

    // ddT of dTdt
    /*J(this->nSpecie_, this->nSpecie_) = 0;
    for (label i = 0; i < this->nSpecie_; i++)
    {
        const label si = reduced ? simplifiedToCompleteIndex_[i] : i;
        J(this->nSpecie_, this->nSpecie_) +=
            cpi[si]*dcdt[i]
          + hi[si]*J(i, this->nSpecie_);
    }
    J(this->nSpecie_, this->nSpecie_) += dTdt*dcpdTMean;
    J(this->nSpecie_, this->nSpecie_) /= -cpMean;
    J(this->nSpecie_, this->nSpecie_) += dTdt/T;*/
}


template<class ReactionThermo, class GThermoType, class SThermoType>
template<class DeltaTType>
Foam::scalar Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
(
    const DeltaTType& deltaT
)
{
	Info << "solve function called" << endl;
    // Increment counter of time-step
    timeSteps_++;

    const bool reduced = mechRed_->active();

    label nAdditionalEqn = (tabulation_->variableTimeStep() ? 3 : 2);

    const basicSpecieMixture& composition = this->thermo().composition();

    // CPU time analysis
    const clockTime clockTime_= clockTime();
    clockTime_.timeIncrement();
    scalar reduceMechCpuTime_ = 0;
    scalar addNewLeafCpuTime_ = 0;
    scalar growCpuTime_ = 0;
    scalar solveChemistryCpuTime_ = 0;
    scalar searchISATCpuTime_ = 0;

    this->resetTabulationResults();

    // Average number of active species
    scalar nActiveSpecies = 0;
    scalar nAvg = 0;

    ReactionThermo::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    if (!this->initialized_)
    {
        this->initializeSurface();
    }

    Info << "surface initialized" << endl;

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    scalarField c(this->nSpecie_);
    scalarField c0(this->nSpecie_);

    // Composition vector (Yi, T, p)
    scalarField phiq(this->nEqns() + nAdditionalEqn);

    scalarField Rphiq(this->nEqns() + nAdditionalEqn);

    forAll(T, celli)
    {
        scalar pi = p[celli];
        scalar Ti = T[celli];
        const scalar rhoi = dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
            (
                this->thermo()
            ).gasCellMixture(celli).rho(pi, Ti);

        for (label i=0; i<this->nGasSpecieOrig_; i++)
        {
                const scalar Yi = this->Y()[i][celli];
                c[i] = rhoi*Yi/this->gasSpecieThermo_[i].W();
                c0[i] = c[i];
                phiq[i] = Yi;
        }

        for (label i=0; i<this->nSolidSpecieOrig_; i++)
        {
                label si = i + this->nGasSpecieOrig_;
                const scalar Ysi = this->Ys()[i][celli];
                c[si] = Ysi*this->solidSpecieThermo_[i].sden()/this->solidSpecieThermo_[i].size();
                c0[si] = c[si];
                phiq[si] = Ysi;
        }
        label nSpecieOrig_ = this->nGasSpecieOrig_ + this->nSolidSpecieOrig_;
        phiq[nSpecieOrig_]=Ti;
        phiq[nSpecieOrig_ + 1]=pi;
        phiq[nSpecieOrig_ + 2]=this->multiplyAV_*this->catarea_;
        phiq[nSpecieOrig_ + 3]=this->multiplyPor_*this->por_;
        if (tabulation_->variableTimeStep())
        {
            phiq[nSpecieOrig_ + 4] = deltaT[celli];
        }


        // Initialise time progress
        scalar timeLeft = deltaT[celli];

        // Not sure if this is necessary
        Rphiq = Zero;

        clockTime_.timeIncrement();

        // When tabulation is active (short-circuit evaluation for retrieve)
        // It first tries to retrieve the solution of the system with the
        // information stored through the tabulation method

        if (tabulation_->active() && tabulation_->retrieve(phiq, Rphiq))
        {
            // Retrieved solution stored in Rphiq
            for (label i=0; i<this->nGasSpecieOrig_; i++)
            {
                c[i] = rhoi*Rphiq[i]/this->gasSpecieThermo_[i].W();
            }

            for (label i=0; i<this->nSolidSpecieOrig_; i++)
            {
                label si = i + this->nGasSpecieOrig_;
                c[si] = Rphiq[si]*this->solidSpecieThermo_[i].sden()/this->solidSpecieThermo_[i].size();
            }

            searchISATCpuTime_ += clockTime_.timeIncrement();
        }
        // This position is reached when tabulation is not used OR
        // if the solution is not retrieved.
        // In the latter case, it adds the information to the tabulation
        // (it will either expand the current data or add a new stored point).
        else
        {

			Info << "entering else loop" << endl;

            // Reset the time
            clockTime_.timeIncrement();

            if (reduced)
            {
                // Reduce mechanism change the number of species (only active)
                mechRed_->reduceMechanism(pi, Ti, c, celli);

                nActiveSpecies += mechRed_->NsSimp();
                nAvg++;
                reduceMechCpuTime_ += clockTime_.timeIncrement();
            }

            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                if (reduced)
                {
                    // completeC_ used in the overridden ODE methods
                    // to update only the active species
                    completeC_ = c;

                    // Solve the reduced set of ODE

                    this->solve
                    (
                        pi,
                        Ti,
                        simplifiedC_,
                        celli,
                        dt,
                        this->deltaTChem_[celli]
                    );

                    for (label i=0; i<NsDAC_; i++)
                    {
                        c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                    }
                }
                else
                {
                    this->solve(pi, Ti, c, celli, dt, this->deltaTChem_[celli]);
                }
                timeLeft -= dt;
            }
                solveChemistryCpuTime_ += clockTime_.timeIncrement();

            // If tabulation is used, we add the information computed here to
            // the stored points (either expand or add)
            if (tabulation_->active())
            {
                for(label i=0; i < this->Y_.size();i++)
                {
                    Rphiq[i] = c[i]/rhoi*this->gasSpecieThermo_[i].W();
                }
                for(label i=0; i < this->Ys_.size();i++)
                {
                    label si = i + this->nGasSpecieOrig_;
                    Rphiq[si] = c[si]/this->solidSpecieThermo_[i].sden()*this->solidSpecieThermo_[i].size();
                }
                if (tabulation_->variableTimeStep())
                {
                  Rphiq[Rphiq.size()-5] = Ti;
                  Rphiq[Rphiq.size()-4] = pi;
                  Rphiq[Rphiq.size()-3] = this->multiplyAV_*this->catarea_;
                  Rphiq[Rphiq.size()-2] = this->multiplyPor_*this->por_;
                  Rphiq[Rphiq.size()-1] = deltaT[celli];
                }
                else
                {
                  Rphiq[Rphiq.size()-4] = Ti;
                  Rphiq[Rphiq.size()-3] = pi;
                  Rphiq[Rphiq.size()-2] = this->multiplyAV_*this->catarea_;
                  Rphiq[Rphiq.size()-1] = this->multiplyPor_*this->por_;
                }
                label growOrAdd =
                    tabulation_->add(phiq, Rphiq, celli, rhoi, deltaT[celli]);
                if (growOrAdd)
                {
                    this->setTabulationResultsAdd(celli);
                    addNewLeafCpuTime_ += clockTime_.timeIncrement();
                }
                else
                {
                    this->setTabulationResultsGrow(celli);
                    growCpuTime_ += clockTime_.timeIncrement();
                }
            }
            // When operations are done and if mechanism reduction is active,
            // the number of species (which also affects nEqns) is set back
            // to the total number of species (stored in the mechRed object)
            if (reduced)
            {
              this->nSpecie_ = mechRed_->nSpecie();
              this->nGasSpecie_ = mechRed_->nGasSpecie();
              this->nSolidSpecie_ = mechRed_->nSolidSpecie();
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);
        }
        // Set the RR vector (used in the solver)
        for (label i=0; i<this->Y_.size(); i++)
        {
            this->RRg_[i][celli] =
                (c[i] - c0[i])*this->gasSpecieThermo_[i].W()/deltaT[celli];
        }

        for (label i=0; i<this->Ys_.size(); i++)
        {
            label si = i + this->nGasSpecieOrig_;
            //if (!this->thermo().composition().active(i+this->nGasSpecieOrig_))
            //{
            this->Ys_[i][celli] =
                c[si]/this->solidSpecieThermo_[i].sden()*this->solidSpecieThermo_[i].size();
            //}
            this->RRs_[i][celli] =
               (c[si] - c0[si])/deltaT[celli]
               /this->solidSpecieThermo_[i].sden()*this->solidSpecieThermo_[i].size();
        }
    }

    if (mechRed_->log() || tabulation_->log())
    {
        cpuSolveFile_()
            << this->time().timeOutputValue()
            << "    " << solveChemistryCpuTime_ << endl;
    }

    if (mechRed_->log())
    {
        cpuReduceFile_()
            << this->time().timeOutputValue()
            << "    " << reduceMechCpuTime_ << endl;
    }

    if (tabulation_->active())
    {
        // Every time-step, look if the tabulation should be updated
        tabulation_->update();

        // Write the performance of the tabulation
        tabulation_->writePerformance();

        if (tabulation_->log())
        {
            cpuRetrieveFile_()
                << this->time().timeOutputValue()
                << "    " << searchISATCpuTime_ << endl;

            cpuGrowFile_()
                << this->time().timeOutputValue()
                << "    " << growCpuTime_ << endl;

            cpuAddFile_()
                << this->time().timeOutputValue()
                << "    " << addNewLeafCpuTime_ << endl;
        }
    }

    if (reduced && nAvg && mechRed_->log())
    {
        // Write average number of species
        nActiveSpeciesFile_()
            << this->time().timeOutputValue()
            << "    " << nActiveSpecies/nAvg << endl;
    }

    if (Pstream::parRun())
    {
        List<bool> active(composition.active());
        Pstream::listCombineGather(active, orEqOp<bool>());
        Pstream::listCombineScatter(active);

        forAll(active, i)
        {
            if (active[i])
            {
                composition.setActive(i);
            }
        }
    }

    return deltaTMin;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setTabulationResultsAdd
(
    const label celli
)
{
    tabulationResults_[celli] = 0;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setTabulationResultsGrow(const label celli)
{
    tabulationResults_[celli] = 1;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setTabulationResultsRetrieve
(
    const label celli
)
{
    tabulationResults_[celli] = 2;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
correctCatalyticWallFluxes
(
    PtrList<volScalarField>& Yg,
    const volScalarField& rhog,
    volScalarField& heg,
    volScalarField& Tg,
    const volScalarField& pg,
    const wordList& patches
)
{
    if (!this->chemistry_)
    {
        Info << "Chemistry is off" << endl;
        return;
    }
    Info << "Correcting catalytic wall fluxes" << endl;

    // Increment counter of time-step
    timeSteps_++;

    label nAdditionalEqn = (tabulation_->variableTimeStep() ? 3 : 2);

    // CPU time analysis
    const clockTime clockTime_= clockTime();
    clockTime_.timeIncrement();
    scalar reduceMechCpuTime_ = 0;
    scalar addNewLeafCpuTime_ = 0;
    scalar growCpuTime_ = 0;
    scalar solveChemistryCpuTime_ = 0;
    scalar searchISATCpuTime_ = 0;

    this->resetTabulationResults();

    scalarField c(this->nSpecie_, 0.0);
    scalarField c0(this->nSpecie_, 0.0);

    scalar deltaT = this->mesh().time().deltaTValue();

    // Composition vector (Yi, T, p)
    scalarField phiq(this->nEqns() + nAdditionalEqn);
    scalarField Rphiq(this->nEqns() + nAdditionalEqn);

    forAll(patches, patchii)
    {
        label patchi = this->mesh().boundary().findPatchID(patches[patchii]);
        const fvPatchScalarField& pT = Tg.boundaryField()[patchi];
        forAll(pT, facei)
        {
            label celli = pT.patch().faceCells()[facei];
           // if(this->setYSolution_)
                this->multiplyAV_ = pT.patch().magSf()[facei]/this->mesh().V()[celli];
            this->multiplyPor_ = 0.0;

            scalar Ti = Tg.primitiveField()[celli];
            scalar pi = pg.primitiveField()[celli];
            scalar rhoi = rhog.primitiveField()[celli];
            this->thermo().T().boundaryFieldRef()[patchi][facei] = Ti;
            this->thermo().p().boundaryFieldRef()[patchi][facei] = pi;

            forAll(Yg, i)
            {
                this->RRgwall_[i][patchi][facei] = 0.0;
                scalar Ygi = Yg[i].primitiveField()[celli];
                this->Y_[i].primitiveFieldRef()[celli] = Ygi;
                c[i] = rhoi*Ygi/this->gasSpecieThermo_[i].W();// kmol/m3
                c0[i] = c[i];
                phiq[i] = Ygi;
            }

            for (label i=0; i<this->nSolidSpecie_; i++)
            {
                scalar Ysi = this->Ys_[i].primitiveField()[celli];
                label si = i + this->nGasSpecie_;
                c[si] = Ysi*this->solidSpecieThermo_[i].sden()/this->solidSpecieThermo_[i].size(); // kmol/m2
                c0[si] = c[si];
                phiq[si] = c[si];
            }

            phiq[this->nSpecie_]=Ti;
            phiq[this->nSpecie_ + 1]=pi;
            phiq[this->nSpecie_ + 2]=this->multiplyAV_*this->catarea_;
            phiq[this->nSpecie_ + 3]=this->multiplyPor_*this->por_;
            if (tabulation_->variableTimeStep())
            {
                phiq[this->nSpecie_ + 4] = deltaT;
            }

            // Initialise time progress
            scalar timeLeft = deltaT;

            // Not sure if this is necessary
            Rphiq = Zero;

            clockTime_.timeIncrement();

            // When tabulation is active (short-circuit evaluation for retrieve)
            // It first tries to retrieve the solution of the system with the
            // information stored through the tabulation method
            if (tabulation_->active() && tabulation_->retrieve(phiq, Rphiq))
            {
                // Retrieved solution stored in Rphiq
                for (label i=0; i<this->Y_.size(); i++)
                {
                    c[i] = rhoi*Rphiq[i]/this->gasSpecieThermo_[i].W();
                }
                for (label i=0; i<this->Ys_.size(); i++)
                {
                    label si = i + this->nGasSpecieOrig_;
                    c[si] = Rphiq[si]*this->solidSpecieThermo_[i].sden()
                            /this->solidSpecieThermo_[i].size();
                }

                searchISATCpuTime_ += clockTime_.timeIncrement();
            }
            // This position is reached when tabulation is not used OR
            // if the solution is not retrieved.
            // In the latter case, it adds the information to the tabulation
            // (it will either expand the current data or add a new stored point).
            else
            {
                // Store total time waiting to attribute to add or grow
                scalar timeTmp = clockTime_.timeIncrement();

                // calculate the chemical source terms
                while (timeLeft > small)
                {
                    scalar dt = timeLeft;
                    this->solve(pi, Ti, c, celli, dt, this->deltaTChem_[celli]);
                    timeLeft -= dt;
                }

                {
                    scalar timeIncr = clockTime_.timeIncrement();
                    solveChemistryCpuTime_ += timeIncr;
                    timeTmp += timeIncr;
                }

                // If tabulation is used, we add the information computed here to
                // the stored points (either expand or add)
                if (tabulation_->active())
                {
                    for (label i=0; i<this->Y_.size(); i++)
                    {
                        Rphiq[i] = c[i]/rhoi*this->gasSpecieThermo_[i].W();
                    }
                    for (label i=0; i<this->Ys_.size(); i++)
                    {
                        label si = i+this->nGasSpecieOrig_;
                        Rphiq[si] = c[si]/this->solidSpecieThermo_[i].sden()
                                *this->solidSpecieThermo_[i].size();
                    }
                    if (tabulation_->variableTimeStep())
                    {
                        Rphiq[Rphiq.size()-5] = Ti;
                        Rphiq[Rphiq.size()-4] = pi;
                        Rphiq[Rphiq.size()-3] = this->multiplyAV_*this->catarea_;
                        Rphiq[Rphiq.size()-2] = this->multiplyPor_*this->por_;
                        Rphiq[Rphiq.size()-1] = deltaT;
                    }
                    else
                    {
                        Rphiq[Rphiq.size()-4] = Ti;
                        Rphiq[Rphiq.size()-3] = pi;
                        Rphiq[Rphiq.size()-2] = this->multiplyAV_*this->catarea_;
                        Rphiq[Rphiq.size()-1] = this->multiplyPor_*this->por_;
                    }
                    label growOrAdd =
                        tabulation_->add(phiq, Rphiq, celli, rhoi, deltaT);
                    if (growOrAdd)
                    {
                        this->setTabulationResultsAdd(celli);
                        addNewLeafCpuTime_ += clockTime_.timeIncrement() + timeTmp;
                    }
                    else
                    {
                        this->setTabulationResultsGrow(celli);
                        growCpuTime_ += clockTime_.timeIncrement() + timeTmp;
                    }
                }

                this->deltaTChem_[celli] =
                    min(this->deltaTChem_[celli], this->deltaTChemMax_);
            }

            for (label i=0; i<this->Y_.size(); i++)
            {
              //  if(this->setYSolution_)
              //  {
                    Yg[i].primitiveFieldRef()[celli] = c[i]*this->gasSpecieThermo_[i].W()/rhoi;
                    Yg[i].boundaryFieldRef()[patchi][facei] = c[i]*this->gasSpecieThermo_[i].W()/rhoi;
                    this->Y_[i].primitiveFieldRef()[celli] = c[i]*this->gasSpecieThermo_[i].W()/rhoi;
                    this->Y_[i].boundaryFieldRef()[patchi][facei] = c[i]*this->gasSpecieThermo_[i].W()/rhoi;
               // }
                this->RRgwall_[i][patchi][facei] =
                    (c[i] - c0[i])*this->gasSpecieThermo_[i].W()/deltaT/this->multiplyAV_;
            }

            scalar sumCS = 0.0;
            for (label i=1; i<this->nSolidSpecie_; i++)
            {
                label si = i + this->nGasSpecieOrig_;
                this->Ys_[i].boundaryFieldRef()[patchi][facei] = max(c[si]
                    /this->solidSpecieThermo_[i].sden()*this->solidSpecieThermo_[i].size(),0.0);
                this->Ys_[i].primitiveFieldRef()[celli] = max(c[si]
                    /this->solidSpecieThermo_[i].sden()*this->solidSpecieThermo_[i].size(),0.0);
                sumCS += this->Ys_[i].boundaryFieldRef()[patchi][facei];
            }
            this->Ys_[0].primitiveFieldRef()[celli] = 1.0 - sumCS;
            this->Ys_[0].boundaryFieldRef()[patchi][facei] = 1.0 - sumCS;
            for (label i=0; i<this->Ys_.size(); i++)
            {
                label si = i+this->nGasSpecieOrig_;
                this->RRswall_[i][patchi][facei] =
                    (c[si] - c0[si])/deltaT;
            }
        }
    }

    if (tabulation_->log())
    {
        cpuSolveFile_()
            << this->time().timeOutputValue()
            << "    " << solveChemistryCpuTime_ << endl;
    }

    if (tabulation_->active())
    {
        // Every time-step, look if the tabulation should be updated
        tabulation_->update();

        // Write the performance of the tabulation
        tabulation_->writePerformance();

        if (tabulation_->log())
        {
            cpuRetrieveFile_()
                << this->time().timeOutputValue()
                << "    " << searchISATCpuTime_ << endl;

            cpuGrowFile_()
                << this->time().timeOutputValue()
                << "    " << growCpuTime_ << endl;

            cpuAddFile_()
                << this->time().timeOutputValue()
                << "    " << addNewLeafCpuTime_ << endl;
        }
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
correctRates
(
    PtrList<volScalarField>& Yg,
    const volScalarField& rhog,
    const volScalarField& Tg,
    const volScalarField& pg
)
{
    if (!this->chemistry_)
    {
        Info << "Chemistry is off" << endl;
        return;
    }
    this->setGasMassFractions(Yg);
    scalarField& T = this->thermo().T().primitiveFieldRef();
    scalarField& p = this->thermo().p().primitiveFieldRef();

    this->multiplyAV_ = 1.0;
    this->multiplyPor_ = 1.0;

    forAll(T, celli)
    {
        T[celli] = Tg[celli];
        p[celli] = pg[celli];
    }

    this->solve(this->mesh().time().deltaTValue());

    /*if(this->setYSolution_)
    {
        forAll(T, celli)
        {
            forAll(this->Y_, i)
            {
                Yg[i].primitiveFieldRef()[celli] = this->Y_[i][celli];
            }
        }
    }*/
}

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
getRatesQdotI
(
    const scalar deltaT,
    const label celli,
    scalarField& RR,
    const scalarField& Yg,
    scalarField& Ys,
    const scalar rhog,
    const scalar Tg,
    scalar Ts,
    scalar p,
    const bool ini,
    const bool sensibleHeat
)
{
    if (!this->chemistry_)
    {
        return 0.0;
    }

    // Checks
    if (deltaT < 0.0)
      FatalErrorInFunction << "Negative deltaT" << abort(FatalError);

    if (celli < 0 || celli>=this->mesh().V().size())
      FatalErrorInFunction << "Invalid celli" << abort(FatalError);

    if (RR.size() != this->nGasSpecieOrig_)
      FatalErrorInFunction << "Wrong RR size" << abort(FatalError);

    if (Yg.size() != this->nGasSpecieOrig_)
      FatalErrorInFunction << "Wrong Yg size" << abort(FatalError);

    if (Ys.size() != this->nSolidSpecieOrig_)
      FatalErrorInFunction << "Wrong Ys size" << abort(FatalError);

    if (Tg < 10.0 || Ts < 10.0)
      FatalErrorInFunction << "Too low temperature: \n Tg = "<< Tg << nl << "Ts = " << Ts << abort(FatalError);

    if (rhog < 1e-8)
      FatalErrorInFunction << "Too low density" << abort(FatalError);

    if (p < 10)
      FatalErrorInFunction << "Too low pressure" << abort(FatalError);

    if (ini)
    {
        this->initializeSurface();
        for (label i=0; i<this->nSolidSpecie_; i++)
        {
            Ys[i] = this->Ys_[i][0];
            Info << "Ys[i]: " << Ys[i] << endl;
        }
        return 0.0;
    }

    if(Ts > Treact_)
    {
        // Increment counter of time-step
        timeSteps_++;

        const bool reduced = mechRed_->active();

        const basicSpecieMixture& composition = this->thermo().composition();

        label nAdditionalEqn = (tabulation_->variableTimeStep() ? 3 : 2);

        //const basicSpecieMixture& composition = this->thermo().composition();

        // CPU time analysis
        const clockTime clockTime_= clockTime();
        clockTime_.timeIncrement();
        scalar reduceMechCpuTime_ = 0;
        scalar addNewLeafCpuTime_ = 0;
        scalar growCpuTime_ = 0;
        scalar solveChemistryCpuTime_ = 0;
        scalar searchISATCpuTime_ = 0;

        this->resetTabulationResults();

        // Average number of active species
        scalar nActiveSpecies = 0;
        scalar nAvg = 0;

        ReactionThermo::correct();

        scalar deltaTMin = great;

        scalarField c(this->nSpecie_);
        scalarField c0(this->nSpecie_);

        // Composition vector (Yi, T, p)
        scalarField phiq(this->nEqns() + nAdditionalEqn);
        scalarField Rphiq(this->nEqns() + nAdditionalEqn);

        scalar rhogs_ = rhog*Tg/Ts;

        for (label i=0; i<this->nGasSpecieOrig_; i++)
        {
            c[i] = rhogs_*Yg[i]/this->gasSpecieThermo_[i].W();
            c0[i] = c[i];
            phiq[i] = Yg[i];
        }

        for (label i=0; i<this->nSolidSpecieOrig_; i++)
        {
            label si = i + this->nGasSpecieOrig_;
            c[si] = Ys[i]*this->solidSpecieThermo_[i].sden()/this->solidSpecieThermo_[i].size();
            c0[si] = c[si];
            phiq[si] = Ys[i];
        }

        label nSpecieOrig_ = this->nGasSpecieOrig_ + this->nSolidSpecieOrig_;
        phiq[nSpecieOrig_]=Ts;
        phiq[nSpecieOrig_ + 1]=p;
        phiq[nSpecieOrig_ + 2]=this->multiplyAV_*this->catarea_;
        phiq[nSpecieOrig_ + 3]=this->multiplyPor_*this->por_;
        if (tabulation_->variableTimeStep())
        {
            phiq[nSpecieOrig_ + 4] = deltaT;
        }

        // Initialise time progress
        scalar timeLeft = deltaT;

        // Not sure if this is necessary
        Rphiq = Zero;

        clockTime_.timeIncrement();

        // When tabulation is active (short-circuit evaluation for retrieve)
        // It first tries to retrieve the solution of the system with the
        // information stored through the tabulation method

        if (tabulation_->active() && tabulation_->retrieve(phiq, Rphiq))
        {
            // Retrieved solution stored in Rphiq
            for (label i=0; i<this->nGasSpecieOrig_; i++)
            {
                c[i] = rhogs_*Rphiq[i]/this->gasSpecieThermo_[i].W();
            }

            for (label i=0; i<this->nSolidSpecieOrig_; i++)
            {
                label si = i + this->nGasSpecieOrig_;
                c[si] = Rphiq[si]*this->solidSpecieThermo_[i].sden()/this->solidSpecieThermo_[i].size();
            }

            searchISATCpuTime_ += clockTime_.timeIncrement();
        }

        // This position is reached when tabulation is not used OR
        // if the solution is not retrieved.
        // In the latter case, it adds the information to the tabulation
        // (it will either expand the current data or add a new stored point).
        else
        {

            // Reset the time
            clockTime_.timeIncrement();

            if (reduced)
            {
                // Reduce mechanism change the number of species (only active)
                Info << "Calling reduceMech" << endl;
                mechRed_->reduceMechanism(p, Ts, c, celli);

                nActiveSpecies += mechRed_->NsSimp();
                nAvg++;
                reduceMechCpuTime_ += clockTime_.timeIncrement();
            }
            // Calculate the chemical source terms
            while (timeLeft > small)
            {
                scalar dt = timeLeft;
                if (reduced)
                {
                    // completeC_ used in the overridden ODE methods
                    // to update only the active species
                    completeC_ = c;

                    // Solve the reduced set of ODE

                    this->solve
                    (
                        p,
                        Ts,
                        simplifiedC_,
                        celli,
                        dt,
                        this->deltaTChem_[celli]
                    );

                    for (label i=0; i<NsDAC_; i++)
                    {
                        c[simplifiedToCompleteIndex_[i]] = simplifiedC_[i];
                    }
                }
                else
                {
                    this->solve(p, Ts, c, celli, dt, this->deltaTChem_[celli]);
                }
                timeLeft -= dt;
            }
                //solveChemistryCpuTime_ += clockTime_.timeIncrement();

            // If tabulation is used, we add the information computed here to
            // the stored points (either expand or add)
            if (tabulation_->active())
            {
                for(label i=0; i < Yg.size();i++)
                {
                    Rphiq[i] = c[i]/rhogs_*this->gasSpecieThermo_[i].W();
                }
                for(label i=0; i < Ys.size();i++)
                {
                    label si = i + this->nGasSpecieOrig_;
                    Rphiq[si] = c[si]/this->solidSpecieThermo_[i].sden()*this->solidSpecieThermo_[i].size();
                }
                if (tabulation_->variableTimeStep())
                {
                  Rphiq[Rphiq.size()-5] = Ts;
                  Rphiq[Rphiq.size()-4] = p;
                  Rphiq[Rphiq.size()-3] = this->multiplyAV_*this->catarea_;
                  Rphiq[Rphiq.size()-2] = this->multiplyPor_*this->por_;
                  Rphiq[Rphiq.size()-1] = deltaT;
                }
                else
                {
                  Rphiq[Rphiq.size()-4] = Ts;
                  Rphiq[Rphiq.size()-3] = p;
                  Rphiq[Rphiq.size()-2] = this->multiplyAV_*this->catarea_;
                  Rphiq[Rphiq.size()-1] = this->multiplyPor_*this->por_;
                }
                label growOrAdd =
                    tabulation_->add(phiq, Rphiq, celli, rhogs_, deltaT);
                if (growOrAdd)
                {
                    this->setTabulationResultsAdd(celli);
                    addNewLeafCpuTime_ += clockTime_.timeIncrement();
                }
                else
                {
                    this->setTabulationResultsGrow(celli);
                    growCpuTime_ += clockTime_.timeIncrement();
                }
            }
            // When operations are done and if mechanism reduction is active,
            // the number of species (which also affects nEqns) is set back
            // to the total number of species (stored in the mechRed object)
            if (reduced)
            {
              this->nSpecie_ = mechRed_->nSpecie();
              this->nGasSpecie_ = mechRed_->nGasSpecie();
              this->nSolidSpecie_ = mechRed_->nSolidSpecie();
            }

            deltaTMin = min(this->deltaTChem_[celli], deltaTMin);

            this->deltaTChem_[celli] =
                min(this->deltaTChem_[celli], this->deltaTChemMax_);
        }
        // Set the RR vector (used in the solver)
        for (label i=0; i<Yg.size(); i++)
        {
            RR[i] = (c[i] - c0[i])*this->gasSpecieThermo_[i].W()/deltaT;
        }

        for (label i=0; i<this->Ys_.size(); i++)
        {
            //label si = i + this->nGasSpecieOrig_;
            //if (!this->thermo().composition().active(i+this->nGasSpecieOrig_))
            //{
            Ys[i] = c[i + this->nGasSpecieOrig_]/this->solidSpecieThermo_[i].sden()*this->solidSpecieThermo_[i].size();
            //}
            //this->RRs_[i][celli] =
            //   (c[si] - c0[si])/deltaT[celli]
            //   /this->solidSpecieThermo_[i].sden()*this->solidSpecieThermo_[i].size();
        }

        if (mechRed_->log() || tabulation_->log())
        {
            cpuSolveFile_()
                << this->time().timeOutputValue()
                << "    " << solveChemistryCpuTime_ << endl;
        }

        if (mechRed_->log())
        {
            cpuReduceFile_()
                << this->time().timeOutputValue()
                << "    " << reduceMechCpuTime_ << endl;
        }

        if (tabulation_->active())
        {
            // Every time-step, look if the tabulation should be updated
            tabulation_->update();

            // Write the performance of the tabulation
            tabulation_->writePerformance();

            if (tabulation_->log())
            {
                cpuRetrieveFile_()
                    << this->time().timeOutputValue()
                    << "    " << searchISATCpuTime_ << endl;

                cpuGrowFile_()
                    << this->time().timeOutputValue()
                    << "    " << growCpuTime_ << endl;

                cpuAddFile_()
                    << this->time().timeOutputValue()
                    << "    " << addNewLeafCpuTime_ << endl;
            }
        }

        if (reduced && nAvg && mechRed_->log())
        {
            // Write average number of species
            nActiveSpeciesFile_()
                << this->time().timeOutputValue()
                << "    " << nActiveSpecies/nAvg << endl;
        }

        /*if (Pstream::parRun())
        {
            List<bool> active(composition.active());
            Pstream::listCombineGather(active, orEqOp<bool>());
            Pstream::listCombineScatter(active);

            forAll(active, i)
            {
                if (active[i])
                {
                    composition.setActive(i);
                }
            }
        }*/
    }
    else
    {
        for (label i=0; i<this->nGasSpecieOrig_; i++)
        {
            RR[i] = 0.0;
        }
    }

    // Calculate and return reaction heat
    scalar Qdot = 0.0;
    forAll(Yg, i)
    {
        const scalar hi = sensibleHeat? this->gasSpecieThermo_[i].Hf() : this->gasSpecieThermo_[i].Ha(p, Ts);
        Qdot -= hi*RR[i];
    }
    //forAll(Ys, i)
    //{
    //    const scalar hi = solidSpecieThermo_[i].Hf();
    //    Qdot -= hi*RR[i+nGasSpecie_];
    //}
    return Qdot;

}

// ************************************************************************* //

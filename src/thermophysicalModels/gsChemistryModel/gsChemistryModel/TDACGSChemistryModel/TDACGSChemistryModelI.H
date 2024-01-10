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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
inline bool
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::variableTimeStep() const
{
    return variableTimeStep_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::label
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::timeSteps() const
{
    return timeSteps_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::autoPtr<Foam::OFstream>
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
logFile(const word& name) const
{
    mkDir(this->ReactionThermo::mesh().time().path()/"TDACGS"/this->group());
    return autoPtr<OFstream>
    (
        new OFstream
        (
            this->ReactionThermo::mesh().time().path()/"TDACGS"/this->group()/name
        )
    );
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline const Foam::PtrList<Foam::volScalarField>&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::Y()
{
    return this->Y_;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
inline const Foam::PtrList<Foam::volScalarField>&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::Ys()
{
    return this->Ys_;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
inline
Foam::autoPtr<Foam::chemistryGSReductionMethod<ReactionThermo, GThermoType, SThermoType>>&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::mechRed()
{
    return mechRed_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::setActive
(
    const label i
)
{
    this->thermo().composition().setActive(i);
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline bool Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::active
(
    const label i
) const
{
    return this->thermo().composition().active(i);
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setNsDAC(const label newNsDAC)
{
    NsDAC_ = newNsDAC;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setNSpecie(const label newNs)
{
    this->nSpecie_ = newNs;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
inline void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setNGasSpecie(const label newNsGas)
{
    this->nGasSpecie_ = newNsGas;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
inline void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setNSolidSpecie(const label newNsSolid)
{
    this->nSolidSpecie_ = newNsSolid;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::scalar
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::Treact() const
{
    return Treact_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::scalar&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::Treact()
{
    return Treact_;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::DynamicList<Foam::label>&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
simplifiedToCompleteIndex()
{
    return simplifiedToCompleteIndex_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::Field<Foam::label>&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
completeToSimplifiedIndex()
{
    return completeToSimplifiedIndex_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline const Foam::Field<Foam::label>&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
completeToSimplifiedIndex() const
{
    return completeToSimplifiedIndex_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::Field<bool>&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::gasReactionsDisabled()
{
    return gasReactionsDisabled_;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::Field<bool>&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solidReactionsDisabled()
{
    return solidReactionsDisabled_;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::scalarField&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::completeC()
{
    return completeC_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::scalarField&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::simplifiedC()
{
    return simplifiedC_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
inline Foam::List<Foam::List<Foam::specieElement>>&
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::specieComp()
{
    return specieComp_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
resetTabulationResults()
{
    forAll(tabulationResults_, tabi)
    {
        tabulationResults_[tabi] = 2;
    }
}

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::mu(label i, scalar T, scalar p)
{
    return this->gasSpecieThermo_[i].mu(T, p);
}

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar
Foam::TDACGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::porosity()
{
    return this->por_;
}
// ************************************************************************* //

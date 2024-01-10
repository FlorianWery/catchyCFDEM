/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "ANNChemistryModel.H"
#include "multiComponentGSMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"


// * * * * * * * * * * *  Private member functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
void Foam::ANNChemistryModel<ReactionThermo, ThermoType>::relu
(
    scalarField& y
)
{
    forAll(y, i)
    {
        y[i] = max(y[i], 0.0);
    }
}


template<class ReactionThermo, class ThermoType>
void Foam::ANNChemistryModel<ReactionThermo, ThermoType>::tanhAct
(
    scalarField& y
)
{
    forAll(y, i)
    {
        y[i] = (exp(y[i])-exp(-y[i]))/(exp(y[i])+exp(-y[i]));
    }
}

template<class ReactionThermo, class ThermoType>
void Foam::ANNChemistryModel<ReactionThermo, ThermoType>::swish
(
    scalarField& y
)
{
    forAll(y, i)
    {
        y[i] = y[i]/(1.0+exp(-y[i]));
    }
}

template<class ReactionThermo, class ThermoType>
void Foam::ANNChemistryModel<ReactionThermo, ThermoType>::sigmoid
(
    scalarField& y
)
{
    forAll(y, i)
    {
        y[i] = 1.0/(1.0+exp(-y[i]));
    }
}

template<class ReactionThermo, class ThermoType>
void Foam::ANNChemistryModel<ReactionThermo, ThermoType>::hardSigmoid
(
    scalarField& y
)
{
    forAll(y, i)
    {
        y[i] = max(min(0.2*y[i] + 0.5, 1.0), 0.0);
    }
}


template<class ReactionThermo, class ThermoType>
void Foam::ANNChemistryModel<ReactionThermo, ThermoType>::matrixMultiplySum
(
    scalarField& y,
    const scalarField& x,
    const scalarField& M,
    const scalarField& b
)
{
    y = b;
    label n(y.size());
    forAll(x, j)
    {
        forAll(y, i)
        {
            y[i] += x[j]*M[i+j*n];
        }
    }
}

template<class ReactionThermo, class ThermoType>
void Foam::ANNChemistryModel<ReactionThermo, ThermoType>::calculateANNRates
(
    const scalar dt,
    const scalarField& c,
    const scalar T,
    const scalar P,
    scalarField& dcdt
)
{

  //for 0 to 1
    /*input_[0] = max(min(1.0, (dt-dtmin_)/(dtmax_-dtmin_)), 0.0);
    input_[1] = max(min(1.0, (T-Tmin_)/(Tmax_-Tmin_)), 0.0);
    input_[2] = max(min(1.0, (P-Pmin_)/(Pmax_-Pmin_)), 0.0);

    for(label i=0; i<nGasSpecie_; i++)
    {
        input_[i+3] = max(min(1.0, (c[i]-Cmin_[i])/(Cmax_[i]-Cmin_[i])), 0.0);
    }

    matrixMultiplySum(l1_, input_, w1_, b1_);
    relu(l1_);
    matrixMultiplySum(l2_, l1_, w2_, b2_);
    relu(l2_);
    //matrixMultiplySum(l3_, l2_, w3_, b3_);
    //relu(l3_);
    matrixMultiplySum(lout_, l2_, wout_, bout_);

    for(label i=0; i<nGasSpecie_; i++)
    {
        scalar ci = lout_[i]*(Cmax_[i]-Cmin_[i]) + Cmin_[i];
        dcdt[i] = (ci - c[i])/dt;
    }*/

  //for -1 to 1
    input_[0] = max(min(1.0, ((dt-dtmin_)/(dtmax_-dtmin_)-0.5)*2), -1.0);
    input_[1] = max(min(1.0, ((T-Tmin_)/(Tmax_-Tmin_)-0.5)*2), -1.0);
    input_[2] = max(min(1.0, ((P-Pmin_)/(Pmax_-Pmin_)-0.5)*2), -1.0);

    for(label i=0; i<this->nSpecie_; i++)
    {
        input_[i+3] = max(min(1.0, ((c[i]-Cmin_[i])/(Cmax_[i]-Cmin_[i])-0.5)*2), -1.0);
    }

    matrixMultiplySum(l1_, input_, w1_, b1_);
    sigmoid(l1_);
    matrixMultiplySum(l2_, l1_, w2_, b2_);
    sigmoid(l2_);
    //matrixMultiplySum(l3_, l2_, w3_, b3_);
    //sigmoid(l3_);
    matrixMultiplySum(lout_, l2_, wout_, bout_);

    for(label i=0; i<this->nSpecie_; i++)
    {
        //dcdt[i] = (max(lout_[i],0.0) - c[i])/dt;
        scalar ci = (max(min(lout_[i],1.0),-1.0)/2+0.5)*(Cmax_[i]-Cmin_[i]) + Cmin_[i];
        dcdt[i] = (ci - c[i])/dt;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::ANNChemistryModel<ReactionThermo, ThermoType>::ANNChemistryModel
(
    const ReactionThermo& thermo
)
:
    StandardChemistryModel<ReactionThermo, ThermoType>(thermo),
    b1_(IFstream(this->subDict("ANN").lookup("b1"))()),
    b2_(IFstream(this->subDict("ANN").lookup("b2"))()),
    //b3_(IFstream(this->subDict("ANN").lookup("b3"))()),
    bout_(IFstream(this->subDict("ANN").lookup("bout"))()),
    w1_(IFstream(this->subDict("ANN").lookup("w1"))()),
    w2_(IFstream(this->subDict("ANN").lookup("w2"))()),
    //w3_(IFstream(this->subDict("ANN").lookup("w3"))()),
    wout_(IFstream(this->subDict("ANN").lookup("wout"))()),
    Trhominmax_(IFstream(this->subDict("ANN").lookup("minmaxTPrho"))()),
    Tmin_(Trhominmax_[0]),
    Tmax_(Trhominmax_[1]),
    Pmin_(Trhominmax_[2]),
    Pmax_(Trhominmax_[3]),
    rhomin_(Trhominmax_[4]),
    rhomax_(Trhominmax_[5]),
    dtmin_(Trhominmax_[6]),
    dtmax_(Trhominmax_[7]),
    Cmin_(IFstream(this->subDict("ANN").lookup("minC"))()),
    Cmax_(IFstream(this->subDict("ANN").lookup("maxC"))()),
    input_(this->nSpecie_+3, 0.0),
    l1_(b1_.size(), 0.0),
    l2_(b2_.size(), 0.0),
    //l3_(b3_.size(), 0.0),
    lout_(bout_.size(), 0.0)
{
    Info<< "Number of species in ANNChemistryModel: " << this->nSpecie_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::ANNChemistryModel<ReactionThermo, ThermoType>::
~ANNChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::ANNChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();
    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(T, celli)
    {
        if (!this->thermo().cellReacting(celli)) continue;

        scalar Ti = T[celli];

        if (Ti > this->Treact_)
        {
            scalar pi = p[celli];
            scalar rhoi = rho[celli];

            for (label i=0; i<this->nSpecie_; i++)
            {
                this->c_[i] = rhoi*this->Y_[i][celli]/this->specieThermos_[i].W();
            }

            // Calculate rates via ANN
            calculateANNRates(deltaT[celli], this->c_, Ti, pi, this->dcdt_);

            for (label i=0; i<this->nSpecie_; i++)
            {
                this->RR_[i][celli] = this->dcdt_[i]*this->specieThermos_[i].W();
            }
        }
        else
        {
            for (label i=0; i<this->nSpecie_; i++)
            {
                this->RR_[i][celli] = 0;
            }
        }
    }

    return deltaTMin;
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::ANNChemistryModel<ReactionThermo, ThermoType>::solve
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


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::ANNChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}

// ************************************************************************* //

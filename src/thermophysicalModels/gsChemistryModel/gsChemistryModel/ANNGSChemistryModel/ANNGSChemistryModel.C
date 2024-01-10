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

#include "ANNGSChemistryModel.H"
#include "multiComponentGSMixture.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"


// * * * * * * * * * * *  Private member functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::relu
(
    scalarField& y
)
{
    forAll(y, i)
    {
        y[i] = max(y[i], 0.0);
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::tanhAct
(
    scalarField& y
)
{
    forAll(y, i)
    {
        y[i] = (exp(y[i])-exp(-y[i]))/(exp(y[i])+exp(-y[i]));
    }
}

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::swish
(
    scalarField& y
)
{
    forAll(y, i)
    {
        y[i] = y[i]/(1.0+exp(-y[i]));
    }
}

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::sigmoid
(
    scalarField& y
)
{
    forAll(y, i)
    {
        y[i] = 1.0/(1.0+exp(-y[i]));
    }
}

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::hardSigmoid
(
    scalarField& y
)
{
    forAll(y, i)
    {
        y[i] = max(min(0.2*y[i] + 0.5, 1.0), 0.0);
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::matrixMultiplySum
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

template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::calculateANNRates
(
    const scalar dt,
    scalarField& c,
    scalarField& Ys,
    const scalar Tg,
    const scalar Ts,
    const scalar P,
    scalarField& dcdt
)
{
    //Assumes input is scaled from -1 to 1

    /*input_[0] = max(min(1.0, ((dt-dtmin_)/(dtmax_-dtmin_)-0.5)*2), -1.0);
    input_[1] = max(min(1.0, ((Tg-Tmin_)/(Tmax_-Tmin_)-0.5)*2), -1.0);
    //input_[2] = max(min(1.0, ((Ts-Tsmin_)/(Tsmax_-Tsmin_)-0.5)*2), -1.0);
    input_[2] = max(min(1.0, ((P-Pmin_)/(Pmax_-Pmin_)-0.5)*2), -1.0);*/

    input_[0] = max(min(1.0, ((std::log10(dt)-dtmin_)-0.5)*2), -1.0);
    input_[1] = max(min(1.0, ((std::log10(Tg)-Tmin_)-0.5)*2), -1.0);
    //input_[2] = max(min(1.0, ((Ts-Tsmin_)/(Tsmax_-Tsmin_)-0.5)*2), -1.0);
    input_[2] = max(min(1.0, ((std::log10(P)-Pmin_)-0.5)*2), -1.0);

    for(label i=0; i<this->nGasSpecie_; i++)
    {
        input_[i+3] = max(min(1.0, ((c[i]-Cmingas_[i])/(Cmaxgas_[i]-Cmingas_[i])-0.5)*2), -1.0);
    }

    for(label i=0; i<this->nSolidSpecie_; i++)
    {
        scalar si = i + 3 + this->nGasSpecie_;
        //input_[si] = max(min(1.0, ((Ys[i]-Cminsurf_[i])/(Cmaxsurf_[i]-Cminsurf_[i])-0.5)*2), -1.0);
        input_[si] = max(min(1.0, ((std::log10(max(Ys[i],1e-30))-Cminsurf_[i])/(Cmaxsurf_[i]-Cminsurf_[i])-0.5)*2), -1.0);
    }

    //Define architecture IF NO SPLIT BETWEEN GAS AND SURFAC
    matrixMultiplySum(l1_, input_, w1_, b1_);
    sigmoid(l1_);
    matrixMultiplySum(l2_, l1_, w2_, b2_);
    sigmoid(l2_);
    matrixMultiplySum(l3_, l2_, w3_, b3_);
    sigmoid(l3_);
    matrixMultiplySum(l4_, l3_, w4_, b4_);
    sigmoid(l4_);
    matrixMultiplySum(loutgas_, l4_, woutgas_, boutgas_);
    matrixMultiplySum(loutsurf_, l4_, woutsurf_, boutsurf_);

    //Define architecture IF SPLIT BETWEEN GAS AND SURFAC
    /*matrixMultiplySum(l1_, input_, w1_, b1_);
    sigmoid(l1_);
    matrixMultiplySum(l2gas_, l1_, w2gas_, b2gas_);
    matrixMultiplySum(l2surf_, l1_, w2surf_, b2surf_);
    sigmoid(l2gas_);
    sigmoid(l2surf_);
    matrixMultiplySum(l3gas_, l2gas_, w3gas_, b3gas_);
    matrixMultiplySum(l3surf_, l2surf_, w3surf_, b3surf_);
    sigmoid(l3gas_);
    sigmoid(l3surf_);
    matrixMultiplySum(loutgas_, l3gas_, woutgas_, boutgas_);
    matrixMultiplySum(loutsurf_, l3surf_, woutsurf_, boutsurf_);*/

    for(label i=0; i<this->nGasSpecie_; i++)
    {
        scalar ci = (max(min(loutgas_[i],1.0),-1.0)/2+0.5)*(Cmaxgas_[i]-Cmingas_[i]) + Cmingas_[i];
        dcdt[i] = (ci - c[i])/dt;
    }

    scalar Yst = 0;
    for(label i=0; i<this->nSolidSpecie_; i++)
    {
        //Ys[i] = (max(min(loutsurf_[i],1.0),-1.0)/2+0.5)*(Cmaxsurf_[i]-Cminsurf_[i]) + Cminsurf_[i];
        Ys[i] = pow(10,(max(min(loutsurf_[i],1.0),-1.0)/2+0.5)*(Cmaxsurf_[i]-Cminsurf_[i]) + Cminsurf_[i]);
        Yst += Ys[i];
    }

    for(label i=0; i<this->nSolidSpecie_; i++)
    {
        Ys[i] = Ys[i]/Yst;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::ANNGSChemistryModel
(
    rhoReactionThermo& thermo
)
:
    StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>(thermo),
    b1_(IFstream(this->subDict("ANN").lookup("b1"))()),
    b2_(IFstream(this->subDict("ANN").lookup("b2"))()),
    //b2gas_(IFstream(this->subDict("ANN").lookup("b2gas"))()),
    //b2surf_(IFstream(this->subDict("ANN").lookup("b2surf"))()),
    b3_(IFstream(this->subDict("ANN").lookup("b3"))()),
    b4_(IFstream(this->subDict("ANN").lookup("b4"))()),
    //b3gas_(IFstream(this->subDict("ANN").lookup("b3gas"))()),
    //b3surf_(IFstream(this->subDict("ANN").lookup("b3surf"))()),
    boutgas_(IFstream(this->subDict("ANN").lookup("boutgas"))()),
    boutsurf_(IFstream(this->subDict("ANN").lookup("boutsurf"))()),
    w1_(IFstream(this->subDict("ANN").lookup("w1"))()),
    w2_(IFstream(this->subDict("ANN").lookup("w2"))()),
    //w2gas_(IFstream(this->subDict("ANN").lookup("w2gas"))()),
    //w2surf_(IFstream(this->subDict("ANN").lookup("w2surf"))()),
    w3_(IFstream(this->subDict("ANN").lookup("w3"))()),
    w4_(IFstream(this->subDict("ANN").lookup("w4"))()),
    //w3gas_(IFstream(this->subDict("ANN").lookup("w3gas"))()),
    //w3surf_(IFstream(this->subDict("ANN").lookup("w3surf"))()),
    woutgas_(IFstream(this->subDict("ANN").lookup("woutgas"))()),
    woutsurf_(IFstream(this->subDict("ANN").lookup("woutsurf"))()),
    Trhominmax_(IFstream(this->subDict("ANN").lookup("minmaxTPrho"))()),
    Tmin_(Trhominmax_[0]),
    Tmax_(Trhominmax_[1]),
    //Tsmin_(Trhominmax_[2]),
    //Tsmax_(Trhominmax_[3]),
    Pmin_(Trhominmax_[2]),
    Pmax_(Trhominmax_[3]),
    rhomin_(Trhominmax_[4]),
    rhomax_(Trhominmax_[5]),
    dtmin_(Trhominmax_[6]),
    dtmax_(Trhominmax_[7]),
    Cmingas_(IFstream(this->subDict("ANN").lookup("minCgas"))()),
    Cmaxgas_(IFstream(this->subDict("ANN").lookup("maxCgas"))()),
    Cminsurf_(IFstream(this->subDict("ANN").lookup("minCsurf"))()),
    Cmaxsurf_(IFstream(this->subDict("ANN").lookup("maxCsurf"))()),
    input_(this->nSpecie_+4, 0.0),
    //input_(this->nGasSpecie_+4, 0.0),
    l1_(b1_.size(), 0.0),
    l2_(b2_.size(), 0.0),
    //l2gas_(b2gas_.size(), 0.0),
    //l2surf_(b2surf_.size(), 0.0),
    l3_(b3_.size(), 0.0),
    l4_(b4_.size(), 0.0),
    //l3gas_(b3gas_.size(), 0.0),
    //l3surf_(b3surf_.size(), 0.0),
    loutgas_(boutgas_.size(), 0.0),
    loutsurf_(boutsurf_.size(), 0.0),
    rhogs_(0.0)
{
    Info<< "ANNGSChemistryModel: " << nl
        << indent << "Number of gas species = " << this->nGasSpecie_ << nl
        << indent << "Catalyst porosity = " << this->por_ << nl
        << indent << "Catalyst surface area to volume = " << this->catarea_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
~ANNGSChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::tmp<Foam::volScalarField>
Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::Qdot() const
{
    tmp<volScalarField> tQdot
    (
        volScalarField::New
        (
            "Qdot",
            this->mesh_,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
        )
    );

    if (this->chemistry_)
    {
        scalarField& Qdot = tQdot.ref();

        forAll(this->Y_, i)
        {
            forAll(Qdot, celli)
            {
                const scalar hi = this->gasSpecieThermo_[i].Hf();
                Qdot[celli] -= hi*this->RRg_[i][celli];
            }
        }
    }

    return tQdot;
}

template<class ReactionThermo, class GThermoType, class SThermoType>
template<class DeltaTType>
Foam::scalar Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    //IGNORE THIS FUNCTION, GETRATESQDOTI DOES THE WORK
    ReactionThermo::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(T, celli)
    {
        if (!this->thermo().cellReacting(celli)) continue;

        scalar Ti = T[celli];

        if (Ti > this->Treact_)
        {
            scalar pi = p[celli];
            const scalar rhoi = dynamic_cast<const multiComponentGSMixture<SThermoType, GThermoType>&>
            (
                this->thermo()
            ).gasCellMixture(celli).rho(pi, Ti);

            for (label i=0; i<this->nGasSpecie_; i++)
            {
                this->c_[i] = rhoi*this->Y_[i][celli]/this->gasSpecieThermo_[i].W();
            }

            // Calculate rates via ANN
            //calculateANNRates(deltaT[celli], this->c_, Ti, Ti, pi, this->dcdt_);

            for (label i=0; i<this->nGasSpecie_; i++)
            {
                this->RRg_[i][celli] = this->dcdt_[i]*this->gasSpecieThermo_[i].W();
            }
        }
        else
        {
            for (label i=0; i<this->nGasSpecie_; i++)
            {
                this->RRg_[i][celli] = 0;
            }
        }
    }

    return deltaTMin;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
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
Foam::scalar Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setGasMassFractions
(
    const PtrList<volScalarField>& Yg
)
{
    forAll(this->Y_, i)
    {
        forAll(this->Y_[i], celli)
        {
            this->Y_[i][celli] = Yg[i][celli];
        }
        this->Y_[i].correctBoundaryConditions();
        Info << "min/max "<< this->Y_[i].name() << ": " << gMin(this->Y_[i]) << " / " << gMax(this->Y_[i]) << endl;
    }
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
disableSurfaceSpeciesEqn()
{
    Info << "Disabling solution of surface species equation" << nl
        << indent << "NOTE: this is obvious when using ANN chemistry model" << endl;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
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
    Info << "Correcting catalytic zone rates" << endl;

    setGasMassFractions(Yg);
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
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
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
    NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalarField
Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
catalyticWallSpeciesFlux
(
    const word speciei,
    const label patchi
) const
{
    NotImplemented;
    return this->c_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalarField
Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
catalyticWallHeatFlux
(
    const label patchi
) const
{
    NotImplemented;
    return this->c_;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setAlpha
(
    const volScalarField& alpha
)
{
    NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
void Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
setTemperature
(
    const volScalarField& T
)
{
    NotImplemented;
}


template<class ReactionThermo, class GThermoType, class SThermoType>
Foam::scalar Foam::ANNGSChemistryModel<ReactionThermo, GThermoType, SThermoType>::
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

    if (RR.size() != this->nGasSpecie_)
      FatalErrorInFunction << "Wrong RR size" << abort(FatalError);

    if (Yg.size() != this->nGasSpecie_)
      FatalErrorInFunction << "Wrong Yg size" << abort(FatalError);

    if (Ys.size() != this->nSolidSpecie_)
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
        return 0.0;
    }
    // Calculate reaction rates
    scalarField c0(this->nSpecie_);
    if (Ts > this->Treact_)
    {
        rhogs_ = rhog*Tg/Ts;
        for (label i=0; i<this->nGasSpecie_; i++)
        {
            this->c_[i] = rhogs_*Yg[i]/this->gasSpecieThermo_[i].W();
            c0[i] = this->c_[i];
        }
        /*for (label i=0; i<this->nSolidSpecie_; i++)
        {
            //Coverages are predicted by ANN, not surface concentration

            this->c_[i + this->nGasSpecie_] = Ys[i]; //*this->solidSpecieThermo_[i].sden()/this->solidSpecieThermo_[i].size();
            c0[i + this->nGasSpecie_] = this->c_[i + this->nGasSpecie_];
        }*/

        // Calculate the chemical source terms
        calculateANNRates(deltaT, this->c_, Ys, Tg, Ts, p, this->dcdt_);

        for (label i=0; i<this->nGasSpecie_; i++)
        {
            RR[i] = this->dcdt_[i]*this->gasSpecieThermo_[i].W();
        }

        //Coverages are predicted by ANN and assigned to Ys directly in calculateANNRates
        /*for (label i=0; i<this->nSolidSpecie_; i++)
        {
            //RR[i+nGasSpecie_] = (c_[i+nGasSpecie_] - c0[i+nGasSpecie_])/deltaT
            //        /solidSpecieThermo_[i].sden()*solidSpecieThermo_[i].size();
            Ys[i] = this->c_[i+this->nGasSpecie_]; //  /this->solidSpecieThermo_[i].sden()*this->solidSpecieThermo_[i].size();
        }*/
    }
    else
    {
        for (label i=0; i<this->nGasSpecie_; i++)
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

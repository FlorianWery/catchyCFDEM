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

Class
    Foam::ANNGSChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    ANNGSChemistryModelI.H
    ANNGSChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef ANNGSChemistryModel_H
#define ANNGSChemistryModel_H

#include "StandardGSChemistryModel.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ANNGSChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class GThermoType, class SThermoType>
class ANNGSChemistryModel
:
    public StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>
{
    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

        //- Gas density in pores [kg/m3]
        scalar rhogs_;

        //- Relu activation
        void relu
        (
            scalarField& y
        );

        //- Tanh activation
        void tanhAct
        (
            scalarField& y
        );

        //- Swish activation
        void swish
        (
            scalarField& y
        );

        //- Sigmoid activation
        void sigmoid
        (
            scalarField& y
        );

        //- Hard sigmoid activation
        void hardSigmoid
        (
            scalarField& y
        );

        //- Matrix multiplication
        void matrixMultiplySum
        (
            scalarField& y,
            const scalarField& x,
            const scalarField& M,
            const scalarField& b
        );

        //- Calculate ANN-based rates
        void calculateANNRates
        (
            const scalar dt,
            scalarField& c,
            scalarField& Ys,
            const scalar Tg,
            const scalar Ts,
            const scalar P,
            scalarField& dcdt
        );


protected:

    typedef SThermoType thermoType;


    // Protected data

        //- Bias vectors
        const scalarField b1_;
        const scalarField b2_;
        const scalarField b2gas_;
        const scalarField b2surf_;
        const scalarField b3_;
        const scalarField b4_;
        const scalarField b3gas_;
        const scalarField b3surf_;
        const scalarField boutgas_;
        const scalarField boutsurf_;

        //- Weights
        const scalarField w1_;
        const scalarField w2_;
        const scalarField w2gas_;
        const scalarField w2surf_;
        const scalarField w3_;
        const scalarField w4_;
        const scalarField w3gas_;
        const scalarField w3surf_;
        const scalarField woutgas_;
        const scalarField woutsurf_;

        //- Minima and maxima
        const scalarField Trhominmax_;
        const scalar Tmin_;
        const scalar Tmax_;
        //const scalar Tsmin_;
        //const scalar Tsmax_;
        const scalar Pmin_;
        const scalar Pmax_;
        const scalar rhomin_;
        const scalar rhomax_;
        const scalar dtmin_;
        const scalar dtmax_;
        const scalarField Cmingas_;
        const scalarField Cmaxgas_;
        const scalarField Cminsurf_;
        const scalarField Cmaxsurf_;

        //- Temporary matrices for coverages calculation
        scalarField input_;
        scalarField l1_;
        scalarField l2_;
        scalarField l2gas_;
        scalarField l2surf_;
        scalarField l3_;
        scalarField l4_;
        scalarField l3gas_;
        scalarField l3surf_;
        scalarField loutgas_;
        scalarField loutsurf_;

public:

    //- Runtime type information
    TypeName("gsANN");


    // Constructors

        //- Construct from thermo
        ANNGSChemistryModel(rhoReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        ANNGSChemistryModel(ANNGSChemistryModel&);


    //- Destructor
    virtual ~ANNGSChemistryModel();


    // Member Functions

        // Chemistry model functions (overriding abstract functions in
        // StandardGSChemistryModel.H)
            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volScalarField> Qdot() const;

        // Other functions

            //- Correct catalytic reaction rates
            virtual void correctRates
            (
                PtrList<volScalarField>& Yg,
                const volScalarField& rhog,
                const volScalarField& Tg,
                const volScalarField& pg
            );

            //- Set gas phase mass fractions (e.g. from other phase)
            virtual void setGasMassFractions(const PtrList<volScalarField>& Yg);

            //- Disable surface species equations
            virtual void disableSurfaceSpeciesEqn();

            //- Correct catalytic wall fluxes
            virtual void correctCatalyticWallFluxes
            (
                PtrList<volScalarField>& Yg,
                const volScalarField& rhog,
                volScalarField& heg,
                volScalarField& Tg,
                const volScalarField& pg,
                const wordList& patches
            );

            //- Return catalytic wall species flux on patchi
            virtual scalarField catalyticWallSpeciesFlux
            (
                const word speciei,
                const label patchi
            ) const;

            //- Return catalytic wall heat flux on patchi
            virtual scalarField catalyticWallHeatFlux(const label patchi) const;

            //- Set volume fraction alpha
            virtual void setAlpha(const volScalarField& alpha);

            //- Set temperatures (e.g. from other phase)
            virtual void setTemperature(const volScalarField& Tg);

            //- Get rate in one cell or for one particle,
            //  given the conditions (Y, Ys, rho, T, p).
            //  This function returns the reaction heat
            //  Optionally, initialize surface coverages
            virtual scalar getRatesQdotI
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
                const bool initialize = false,
                const bool sensibleHeat = true
            );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ANNGSChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ANNGSChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Class
    Foam::TDACChemistryModel

Description
    Extends StandardChemistryModel by adding the TDAC method.

    References:
    \verbatim
        Contino, F., Jeanmart, H., Lucchini, T., & D’Errico, G. (2011).
        Coupling of in situ adaptive tabulation and dynamic adaptive chemistry:
        An effective method for solving combustion in engine simulations.
        Proceedings of the Combustion Institute, 33(2), 3057-3064.

        Contino, F., Lucchini, T., D'Errico, G., Duynslaegher, C.,
        Dias, V., & Jeanmart, H. (2012).
        Simulations of advanced combustion modes using detailed chemistry
        combined with tabulation and mechanism reduction techniques.
        SAE International Journal of Engines,
        5(2012-01-0145), 185-196.

        Contino, F., Foucher, F., Dagaut, P., Lucchini, T., D’Errico, G., &
        Mounaïm-Rousselle, C. (2013).
        Experimental and numerical analysis of nitric oxide effect on the
        ignition of iso-octane in a single cylinder HCCI engine.
        Combustion and Flame, 160(8), 1476-1483.

        Contino, F., Masurier, J. B., Foucher, F., Lucchini, T., D’Errico, G., &
        Dagaut, P. (2014).
        CFD simulations using the TDAC method to model iso-octane combustion
        for a large range of ozone seeding and temperature conditions
        in a single cylinder HCCI engine.
        Fuel, 137, 179-184.
    \endverbatim

SourceFiles
    TDACChemistryModelI.H
    TDACChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef TDACGSChemistryModel_H
#define TDACGSChemistryModel_H

#include "StandardGSChemistryModel.H"
#include "chemistryGSReductionMethod.H"
#include "chemistryGSTabulationMethod.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class TDACChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class GThermoType, class SThermoType>
class TDACGSChemistryModel
:
    public StandardGSChemistryModel<ReactionThermo, GThermoType, SThermoType>
{
    // Private member data

        bool variableTimeStep_;

        label timeSteps_;

        //- Temperature below which the reaction rates are assumed 0
        scalar Treact_;

        // Mechanism reduction
        label NsDAC_;
        scalarField completeC_;
        scalarField simplifiedC_;
        Field<bool> gasReactionsDisabled_;
        Field<bool> solidReactionsDisabled_;
        List<List<specieElement>> specieComp_;
        Field<label> completeToSimplifiedIndex_;
        DynamicList<label> simplifiedToCompleteIndex_;
        autoPtr<chemistryGSReductionMethod<ReactionThermo, GThermoType, SThermoType>>
            mechRed_;

        // Tabulation
        autoPtr<chemistryGSTabulationMethod<ReactionThermo, GThermoType, SThermoType>>
            tabulation_;

        // Log file for the average time spent reducing the chemistry
        autoPtr<OFstream> cpuReduceFile_;

        // Write average number of species
        autoPtr<OFstream> nActiveSpeciesFile_;

        //- Log file for the average time spent adding tabulated data
        autoPtr<OFstream> cpuAddFile_;

        //- Log file for the average time spent growing tabulated data
        autoPtr<OFstream> cpuGrowFile_;

        //- Log file for the average time spent retrieving tabulated data
        autoPtr<OFstream> cpuRetrieveFile_;

        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;

        // Field containing information about tabulation:
        // 0 -> add (direct integration)
        // 1 -> grow
        // 2 -> retrieve
        volScalarField tabulationResults_;


    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


public:

    //- Runtime type information
    TypeName("TDACGS");


    // Constructors

        //- Construct from thermo
        TDACGSChemistryModel(rhoReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        TDACGSChemistryModel(TDACGSChemistryModel&);


    //- Destructor
    virtual ~TDACGSChemistryModel();


    // Member Functions

        //- Return true if the time-step is variable and/or non-uniform
        inline bool variableTimeStep() const;

        //- Return the number of chemistry evaluations (used by ISAT)
        inline label timeSteps() const;

        //- Create and return a TDAC log file of the given name
        inline autoPtr<OFstream> logFile(const word& name) const;

        inline const PtrList<volScalarField>& Y();

        inline const PtrList<volScalarField>& Ys();

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual void omega
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& dcdt
        ) const;


        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction<GThermoType>& r,
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;


        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const SurfaceReaction<SThermoType>& r,
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

        // Chemistry model functions (overriding functions in
        // StandardGSChemistryModel to use the private solve function)

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);


        // ODE functions (overriding functions in StandardChemistryModel to take
        // into account the variable number of species)

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt,
                scalarSquareMatrix& J
            ) const;

            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& c,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;


        // Mechanism reduction access functions

            inline void setNsDAC(const label newNsDAC);

            inline void setNSpecie(const label newNs);

            inline void setNGasSpecie(const label newNsGas);

            inline void setNSolidSpecie(const label newNsSolid);

            inline scalarField& completeC();

            inline scalarField& simplifiedC();

            inline Field<bool>& gasReactionsDisabled();

            inline Field<bool>& solidReactionsDisabled();

            inline bool active(const label i) const;

            inline void setActive(const label i);

            inline DynamicList<label>& simplifiedToCompleteIndex();

            inline Field<label>& completeToSimplifiedIndex();

            //- Temperature below which the reaction rates are assumed 0
            inline scalar Treact() const;

            //- Temperature below which the reaction rates are assumed 0
            inline scalar& Treact();

            inline const Field<label>& completeToSimplifiedIndex() const;

            inline List<List<specieElement>>& specieComp();

            inline
                autoPtr<chemistryGSReductionMethod<ReactionThermo, GThermoType, SThermoType>>&
                mechRed();

            tmp<volScalarField> tabulationResults() const
            {
                return tabulationResults_;
            }

            void setTabulationResultsAdd(const label celli);

            void setTabulationResultsGrow(const label celli);

            void setTabulationResultsRetrieve(const label celli);

            inline void resetTabulationResults();

	    //- Return species mu fields
	    virtual scalar mu(label i, scalar T, scalar p);

      //- Catalyst porosity
      virtual scalar porosity();

            //- Correct catalytic wall fluxes
            void correctCatalyticWallFluxes
            (
              PtrList<volScalarField>& Yg,
              const volScalarField& rhog,
              volScalarField& heg,
              volScalarField& Tg,
              const volScalarField& pg,
              const wordList& patches
            );

            //- Correct catalytic reaction rates
            virtual void correctRates
            (
              PtrList<volScalarField>& Yg,
              const volScalarField& rhog,
              const volScalarField& Tg,
              const volScalarField& pg
            );

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
        void operator=(const TDACGSChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TDACGSChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TDACGSChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

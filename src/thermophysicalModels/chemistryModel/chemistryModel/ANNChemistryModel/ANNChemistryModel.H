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
    Foam::ANNChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    ANNChemistryModelI.H
    ANNChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef ANNChemistryModel_H
#define ANNChemistryModel_H

#include "StandardChemistryModel.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class ANNChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class ANNChemistryModel
:
    public StandardChemistryModel<ReactionThermo, ThermoType>
{
    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

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
            const scalarField& c,
            const scalar T,
            const scalar P,
            scalarField& dcdt
        );


protected:

    // Protected data

        //- Bias vectors
        const scalarField b1_;
        const scalarField b2_;
        const scalarField b3_;
        const scalarField bout_;

        //- Weights
        const scalarField w1_;
        const scalarField w2_;
        const scalarField w3_;
        const scalarField wout_;

        //- Minima and maxima
        const scalarField Trhominmax_;
        const scalar Tmin_;
        const scalar Tmax_;
        const scalar Pmin_;
        const scalar Pmax_;
        const scalar rhomin_;
        const scalar rhomax_;
        const scalar dtmin_;
        const scalar dtmax_;
        const scalarField Cmin_;
        const scalarField Cmax_;

        //- Temporary matrices for coverages calculation
        scalarField input_;
        scalarField l1_;
        scalarField l2_;
        scalarField l3_;
        scalarField lout_;


public:

    //- Runtime type information
    TypeName("ANN");


    // Constructors

        //- Construct from thermo
        ANNChemistryModel(const ReactionThermo& thermo);

        //- Disallow default bitwise copy construction
        ANNChemistryModel(const ANNChemistryModel&);


    //- Destructor
    virtual ~ANNChemistryModel();


    // Member Functions
        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ANNChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ANNChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

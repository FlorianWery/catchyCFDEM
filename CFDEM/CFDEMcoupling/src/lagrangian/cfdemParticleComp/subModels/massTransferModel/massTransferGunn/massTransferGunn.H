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

    Description
    Correlation for Nusselt number according to
    Gunn, D. J. International Journal of Heat and Mass Transfer 21.4 (1978)

\*---------------------------------------------------------------------------*/

#ifndef massTransferGunn_H
#define massTransferGunn_H

#include "fvCFD.H"
#include "cfdemCloudEnergy.H"
#include "massTransferModel.H"

#include "autoPtr.H"
#include "speciesTable.H"
#include "rhoReactionThermo.H"

#include "fluidThermo.H"
#include "solidThermo.H"

#include "basicGSChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class massTransferGunn Declaration
\*---------------------------------------------------------------------------*/

class massTransferGunn
:
    public massTransferModel
{
protected:

    dictionary propsDict_;

    const word catalystName_;

    const word gasName_;

    //- Pointer to thermo
    autoPtr<rhoReactionThermo> catThermoPtr_;

    autoPtr<basicGSChemistryModel> catChemistryPtr_;

    //- Reference to the field of fluid specie mass fractions
    const PtrList<volScalarField>& Yg_;

    const volScalarField& Tg_;

    //- Number of gas species
    const label nGasSpecie_;

    const speciesTable gasSpecies_;

    mutable List<scalarList> ofPartGasFrac_;

    mutable double **partGasFrac_;

    PtrList<volScalarField> massFlux_;

    volScalarField ReField_;

    volScalarField ShField_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;        // ref to voidfraction field

    word velFieldName_;

    const volVectorField& U_;

    word densityFieldName_;

    const volScalarField& rho_;

    word pressureFieldName_;

    const volScalarField& p_;

    word partTempName_;

    mutable double **partTemp_;

    mutable double **partMassFlux_;

    mutable double **partRe_;

    mutable double **partSh_;

    mutable List<scalarList> ofPartMassFlux_;

    word partPressName_;

    mutable double **partPressure_;

    mutable scalar scaleDia_;

    mutable scalar scaleNusselt_;

    mutable scalar Le_;

    scalarList typeCG_;

    void allocateMyArrays() const;

    scalar Nusselt(scalar, scalar, scalar) const;

    virtual void giveData();

    virtual void massFlux(label, label, scalar, scalar, scalar, scalar, scalar);

public:

    //- Runtime type information
    TypeName("massTransferGunn");

    // Constructors

        //- Construct from components
        massTransferGunn
        (
            const dictionary& dict,
            cfdemCloudEnergy& sm
        );


    // Destructor

        virtual ~massTransferGunn();

    // Member Functions

        void addMassTransferContribution(volScalarField&, label) const;

        void calcMassTransferContribution();

        void postFlow();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

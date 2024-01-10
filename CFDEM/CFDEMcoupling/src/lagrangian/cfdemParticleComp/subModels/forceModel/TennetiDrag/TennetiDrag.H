/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

    Description
    Drag model of Tenneti et al. for flow past fixed assemblies of monodisperse 
    spheres with direct numerical simulations with the Immersed boundary method.
    Range of applicability: 0 < alphaParticles < 0.5, Re <= 300

    Reference:
        Tenneti, S., R. Garg, and S. Subramaniam. (2011) 
	Drag law for monodisperse gas–solid systems using particle-resolved 
	direct numerical simulation of flow past fixed assemblies of spheres
        International journal of multiphase flow, 37(9), 1072-1092. 


Class
    TennetiDrag

SourceFiles
    TennetiDrag.C

\*---------------------------------------------------------------------------*/

#ifndef TennetiDrag_H
#define TennetiDrag_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TennetiDrag Declaration
\*---------------------------------------------------------------------------*/

class TennetiDrag
:
    public forceModel
{
private:
    dictionary propsDict_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    const scalar minVoidfraction_;

    word UsFieldName_;

    const volVectorField& UsField_;  // the average particle velocity field (for implicit/expliti force split)

    mutable scalar scaleDia_;

    mutable scalar scaleDrag_;

    double F(double, double) const;

public:

    //- Runtime type information
    TypeName("TennetiDrag");


    // Constructors

        //- Construct from components
        TennetiDrag
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~TennetiDrag();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

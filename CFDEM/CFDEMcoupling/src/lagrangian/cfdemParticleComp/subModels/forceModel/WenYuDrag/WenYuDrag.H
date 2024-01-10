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
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

		Enwald, H., Peirano, E., & Almstedt, A. E. (1996).
        Eulerian two-phase flow theory applied to fluidization.
        International Journal of Multiphase Flow, 22, 21-66.
        Eq. 86-87, p. 40.
Class
    WenYuDrag

SourceFiles
    WenYuDrag.C

\*---------------------------------------------------------------------------*/

#ifndef WenYuDrag_H
#define WenYuDrag_H

#include "forceModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class WenYuDrag Declaration
\*---------------------------------------------------------------------------*/

class WenYuDrag
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
    TypeName("WenYuDrag");


    // Constructors

        //- Construct from components
        WenYuDrag
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~WenYuDrag();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

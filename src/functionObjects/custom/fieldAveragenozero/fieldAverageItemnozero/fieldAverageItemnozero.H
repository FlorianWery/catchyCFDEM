/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::functionObjects::fieldAverageItem

Description
    Helper class to describe what form of averaging to apply.  A set will be
    applied to each base field in Foam::fieldAverage, of the form:

    \verbatim
    {
        mean            on;   // (default = on)
        prime2Mean      on;   // (default = off)
        base            time; // time or iteration (default = time)
        window          200;  // optional averaging window
        windowName      w1;   // optional window name (default = "")
    }
    \endverbatim

    The averaging window corresponds to the averaging interval (iters or time)
    If not specified, the averaging is over 'all iters/time'

SourceFiles
    fieldAverageItem.C
    fieldAverageItemIO.C

\*---------------------------------------------------------------------------*/

#ifndef fieldAverageItemnozero_H
#define fieldAverageItemnozero_H

#include "NamedEnum.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;

namespace functionObjects
{

class fieldAveragenozero;

/*---------------------------------------------------------------------------*\
                      Class fieldAverageItem Declaration
\*---------------------------------------------------------------------------*/

class fieldAverageItemnozero
{
    // Private Data

        //- Field name
        word fieldName_;

        //- Compute mean flag
        Switch mean_;

        //- Name of mean field
        word meanFieldName_;

        //- Compute prime-squared mean flag
        Switch prime2Mean_;

        //- Name of prime-squared mean field
        word prime2MeanFieldName_;


public:

    // Static data

        // File and field name extensions

            //- Mean average
            static const word meanExt;

            //- Prime-squared average
            static const word prime2MeanExt;


    // Constructors

        //- Construct null
        fieldAverageItemnozero();

        //- Construct from Istream for given fieldAverage
        fieldAverageItemnozero(const fieldAveragenozero& fa, Istream& is);

        //- Copy constructor
        fieldAverageItemnozero(const fieldAverageItemnozero&) = default;

        //- Clone function
        autoPtr<fieldAverageItemnozero> clone() const
        {
            return autoPtr<fieldAverageItemnozero>(new fieldAverageItemnozero(*this));
        }

        //- Class used for the read-construction of
        //  PtrLists of fieldAverageItem
        class iNew
        {
            const fieldAveragenozero& fa_;

        public:

            iNew(const fieldAveragenozero& fa)
            :
                fa_(fa)
            {}

            autoPtr<fieldAverageItemnozero> operator()(Istream& is) const
            {
                return autoPtr<fieldAverageItemnozero>(new fieldAverageItemnozero(fa_, is));
            }
        };


    //- Destructor
    ~fieldAverageItemnozero();


    // Member Functions

        // Access

            //- Return const access to the field name
            const word& fieldName() const
            {
                return fieldName_;
            }

            //- Return const access to the mean flag
            const Switch& mean() const
            {
                return mean_;
            }

            //- Return non-const access to the mean flag
            Switch& mean()
            {
                return mean_;
            }

            //- Return const access to the mean field name
            const word& meanFieldName() const
            {
                return meanFieldName_;
            }

            //- Return const access to the prime-squared mean flag
            const Switch& prime2Mean() const
            {
                return prime2Mean_;
            }

            //- Return non-const access to the prime-squared mean flag
            Switch& prime2Mean()
            {
                return prime2Mean_;
            }

            //- Return const access to the prime-squared mean field name
            const word& prime2MeanFieldName() const
            {
                return prime2MeanFieldName_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

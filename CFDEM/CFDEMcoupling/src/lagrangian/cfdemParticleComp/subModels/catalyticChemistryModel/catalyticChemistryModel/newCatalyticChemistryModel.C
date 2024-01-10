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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "catalyticChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

autoPtr<catalyticChemistryModel> catalyticChemistryModel::New
(
    const dictionary& dict,
    cfdemCloudEnergy& sm,
    word catalyticChemistryModelType
)
{
    Info<< "Selecting catalyticChemistryModel "
         << catalyticChemistryModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(catalyticChemistryModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "catalyticChemistryModel::New(const dictionary&, cfdemCloudEnergy&): "
            << endl
            << "    unknown catalyticChemistryModelType type "
            << catalyticChemistryModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid chemistryModel types are :"
            << endl;

        Info<< dictionaryConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    return autoPtr<catalyticChemistryModel>(cstrIter()(dict,sm));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

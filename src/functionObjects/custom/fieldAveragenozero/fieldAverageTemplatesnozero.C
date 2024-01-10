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

\*---------------------------------------------------------------------------*/

#include "fieldAverageItemnozero.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldAveragenozero::addMeanFieldType(const label fieldi)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();

    Log << "    Reading/initialising field " << meanFieldName << endl;

    if (obr_.foundObject<Type>(meanFieldName))
    {}
    else if (obr_.found(meanFieldName))
    {
        Log << "    Cannot allocate average field " << meanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        faItems_[fieldi].mean() = false;
    }
    else
    {
        const Type& baseField = obr_.lookupObject<Type>(fieldName);

        // Store on registry
        obr_.store
        (
            new Type
            (
                IOobject
                (
                    meanFieldName,
                    obr_.time().timeName(obr_.time().startTime().value()),
                    obr_,
                    restartOnOutput_
                  ? IOobject::NO_READ
                  : IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                1*baseField
            )
        );
    }
}


template<class Type>
void Foam::functionObjects::fieldAveragenozero::addMeanField(const label fieldi)
{
    if (faItems_[fieldi].mean())
    {
        typedef GeometricField<Type, fvPatchField, volMesh>
            VolFieldType;

        typedef GeometricField<Type, fvsPatchField, surfaceMesh>
            SurfaceFieldType;

        const word& fieldName = faItems_[fieldi].fieldName();

        if (obr_.foundObject<VolFieldType>(fieldName))
        {
            addMeanFieldType<VolFieldType>(fieldi);
        }
        else if (obr_.foundObject<SurfaceFieldType>(fieldName))
        {
            addMeanFieldType<SurfaceFieldType>(fieldi);
        }
    }
}

template<class Type>
void Foam::functionObjects::fieldAveragenozero::fieldSize(const label fieldi)
{

    typedef GeometricField<Type, fvPatchField, volMesh>
            VolFieldType;

    typedef GeometricField<Type, fvsPatchField, surfaceMesh>
            SurfaceFieldType;

    const word& fieldName = faItems_[fieldi].fieldName();

    if (obr_.foundObject<VolFieldType>(fieldName))
    {
        const VolFieldType& baseField = obr_.lookupObject<VolFieldType>(fieldName);
        fieldSize_ = baseField.size();
    }

    if (obr_.foundObject<SurfaceFieldType>(fieldName))
    {
        const SurfaceFieldType& baseField = obr_.lookupObject<SurfaceFieldType>(fieldName);
        fieldSize_ = baseField.size();
    }

}

template<class Type>
void Foam::functionObjects::fieldAveragenozero::modifyTotalIterTime
(
    const label fieldi
)
{

    const word& fieldName = faItems_[fieldi].fieldName();

    typedef GeometricField<Type, fvPatchField, volMesh>
            VolFieldType;

    typedef GeometricField<Type, fvsPatchField, surfaceMesh>
            SurfaceFieldType;

    if (obr_.foundObject<VolFieldType>(fieldName))
    {
        const VolFieldType& baseField = obr_.lookupObject<VolFieldType>(fieldName);
	for(int i = 0; i < baseField.size(); i++)
        {
            if(mag(baseField[i]) != 0)
	        {
                    totalTime_[fieldi][i] += obr_.time().deltaTValue();
		    if (iterBase())
       		    {
			totalIter_[fieldi][i]++;
       		    }
                }
        }
    }

    if (obr_.foundObject<SurfaceFieldType>(fieldName))
    {
        const SurfaceFieldType& baseField = obr_.lookupObject<SurfaceFieldType>(fieldName);
	for(int i = 0; i < baseField.size(); i++)
        {
            if(mag(baseField[i]) != 0)
	        {
                    totalTime_[fieldi][i] += obr_.time().deltaTValue();
		    if (iterBase())
       		    {
			totalIter_[fieldi][i]++;
       		    }
                }
        }
    }
}

template<class Type1, class Type2>
void Foam::functionObjects::fieldAveragenozero::addPrime2MeanFieldType
(
    const label fieldi
)
{
    const word& fieldName = faItems_[fieldi].fieldName();
    const word& meanFieldName = faItems_[fieldi].meanFieldName();
    const word& prime2MeanFieldName = faItems_[fieldi].prime2MeanFieldName();

    Log << "    Reading/initialising field " << prime2MeanFieldName << nl;

    if (obr_.foundObject<Type2>(prime2MeanFieldName))
    {}
    else if (obr_.found(prime2MeanFieldName))
    {
        Log << "    Cannot allocate average field " << prime2MeanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << nl;

        faItems_[fieldi].prime2Mean() = false;
    }
    else
    {
        const Type1& baseField = obr_.lookupObject<Type1>(fieldName);
        const Type1& meanField = obr_.lookupObject<Type1>(meanFieldName);

        // Store on registry
        obr_.store
        (
            new Type2
            (
                IOobject
                (
                    prime2MeanFieldName,
                    obr_.time().timeName(obr_.time().startTime().value()),
                    obr_,
                    restartOnOutput_
                  ? IOobject::NO_READ
                  : IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                sqr(baseField) - sqr(meanField)
            )
        );
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAveragenozero::addPrime2MeanField(const label fieldi)
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;

    if (faItems_[fieldi].prime2Mean())
    {
        const word& fieldName = faItems_[fieldi].fieldName();

        if (!faItems_[fieldi].mean())
        {
            FatalErrorInFunction
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << fieldName << nl << exit(FatalError);
        }

        if (obr_.foundObject<VolFieldType1>(fieldName))
        {
            addPrime2MeanFieldType<VolFieldType1, VolFieldType2>(fieldi);
        }
        else if (obr_.foundObject<SurfaceFieldType1>(fieldName))
        {
            addPrime2MeanFieldType<SurfaceFieldType1, SurfaceFieldType2>
            (
                fieldi
            );
        }
    }
}

template<class Type>
void Foam::functionObjects::fieldAveragenozero::calculateMeanFieldType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (obr_.foundObject<Type>(fieldName))
    {
        const Type& baseField = obr_.lookupObject<Type>(fieldName);

        Type& meanField =
            obr_.lookupObjectRef<Type>(faItems_[fieldi].meanFieldName());

        scalar dt = obr_.time().deltaTValue();

	for(int i = 0; i < meanField.size(); i++)
	{
		if(mag(baseField[i]) != 0)
		{
		scalar Dt = totalTime_[fieldi][i];
    		if (iterBase())
        	{
            		dt = 1;
            		Dt = scalar(totalIter_[fieldi][i]);
        	}

        	scalar beta = dt/Dt;


        	if (window() > 0)
        	{
            		const scalar w = window();

            		if (Dt - dt >= w)
            		{
                		beta = dt/w;
            		}
        	}

        	meanField[i] = (1 - beta)*meanField[i] + beta*baseField[i];
		}
	}
    }
}


template<class Type>
void Foam::functionObjects::fieldAveragenozero::calculateMeanFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            calculateMeanFieldType<VolFieldType>(i);
            calculateMeanFieldType<SurfaceFieldType>(i);
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAveragenozero::calculatePrime2MeanFieldType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (obr_.foundObject<Type1>(fieldName))
    {
        const Type1& baseField = obr_.lookupObject<Type1>(fieldName);
        const Type1& meanField =
            obr_.lookupObject<Type1>(faItems_[fieldi].meanFieldName());

        Type2& prime2MeanField =
            obr_.lookupObjectRef<Type2>(faItems_[fieldi].prime2MeanFieldName());

        scalar dt = obr_.time().deltaTValue();

	for(int i = 0; i < prime2MeanField.size(); i++)
	{
		if(mag(baseField[i]) !=0)
		{
		scalar Dt = totalTime_[fieldi][i];
    		if (iterBase())
        	{
            		dt = 1;
            		Dt = scalar(totalIter_[fieldi][i]);
        	}

        	scalar beta = dt/Dt;


        	if (window() > 0)
        	{
            		const scalar w = window();

            		if (Dt - dt >= w)
            		{
                		beta = dt/w;
            		}
        	}
		prime2MeanField[i] = (1 - beta)*prime2MeanField[i] + beta*sqr(baseField[i]) - sqr(meanField[i]);
		}
	}
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAveragenozero::calculatePrime2MeanFields() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;

    forAll(faItems_, i)
    {
        if (faItems_[i].prime2Mean())
        {
            calculatePrime2MeanFieldType<VolFieldType1, VolFieldType2>(i);
            calculatePrime2MeanFieldType<SurfaceFieldType1, SurfaceFieldType2>
            (
                i
            );
        }
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAveragenozero::addMeanSqrToPrime2MeanType
(
    const label fieldi
) const
{
    const word& fieldName = faItems_[fieldi].fieldName();

    if (obr_.foundObject<Type1>(fieldName))
    {
        const Type1& baseField = obr_.lookupObject<Type1>(fieldName);
        const Type1& meanField =
            obr_.lookupObject<Type1>(faItems_[fieldi].meanFieldName());

        Type2& prime2MeanField =
            obr_.lookupObjectRef<Type2>(faItems_[fieldi].prime2MeanFieldName());

	for(int i = 0; i < prime2MeanField.size(); i++)
	{
		if(mag(baseField[i]) !=0)
		{
        		prime2MeanField[i] += sqr(meanField[i]);
		}
	}
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAveragenozero::addMeanSqrToPrime2Mean() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;

    forAll(faItems_, i)
    {
        if (faItems_[i].prime2Mean())
        {
            addMeanSqrToPrime2MeanType<VolFieldType1, VolFieldType2>(i);
            addMeanSqrToPrime2MeanType<SurfaceFieldType1, SurfaceFieldType2>(i);
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAveragenozero::writeFieldType
(
    const word& fieldName
) const
{
    if (obr_.foundObject<Type>(fieldName))
    {
        const Type& f = obr_.lookupObject<Type>(fieldName);
        f.write();
    }
}


template<class Type>
void Foam::functionObjects::fieldAveragenozero::writeFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    forAll(faItems_, i)
    {
        if (faItems_[i].mean())
        {
            const word& fieldName = faItems_[i].meanFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
        }
        if (faItems_[i].prime2Mean())
        {
            const word& fieldName = faItems_[i].prime2MeanFieldName();
            writeFieldType<VolFieldType>(fieldName);
            writeFieldType<SurfaceFieldType>(fieldName);
        }
    }
}


// ************************************************************************* //

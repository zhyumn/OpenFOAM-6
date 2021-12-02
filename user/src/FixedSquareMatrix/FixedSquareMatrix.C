/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "FixedSquareMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::FixedSquareMatrix::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class Type, unsigned Size>
Foam::FixedSquareMatrix<Type, Size>::FixedSquareMatrix()
{}

/*
Foam::FixedSquareMatrix::FixedSquareMatrix(const dataType& data)
:
    baseClassName(),
    data_(data)
{}
*/
template<class Type, unsigned Size>
Foam::FixedSquareMatrix<Type, Size>::FixedSquareMatrix(const FixedSquareMatrix& M)
{
    for (label i = 0; i < Size * Size; i++)
    {
        v_[i] = M.v_[i];
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::FixedSquareMatrix>
Foam::FixedSquareMatrix::New()
{
    return autoPtr<FixedSquareMatrix>(new FixedSquareMatrix);
}*/


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
template<class Type, unsigned Size>
Foam::FixedSquareMatrix<Type, Size>::~FixedSquareMatrix()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //
/*
void Foam::FixedSquareMatrix::operator=(const FixedSquareMatrix& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}
*/
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //




// ************************************************************************* //

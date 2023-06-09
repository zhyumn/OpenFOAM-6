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

#include "FixedLUscalarMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::FixedLUscalarMatrix::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
/*
Foam::FixedLUscalarMatrix::FixedLUscalarMatrix()
:
    baseClassName(),
    data_()
{}


Foam::FixedLUscalarMatrix::FixedLUscalarMatrix(const dataType& data)
:
    baseClassName(),
    data_(data)
{}


Foam::FixedLUscalarMatrix::FixedLUscalarMatrix(const FixedLUscalarMatrix&)
:
    baseClassName(),
    data_()
{}
*/
template<unsigned Size>
Foam::FixedLUscalarMatrix<Size>::FixedLUscalarMatrix(const FixedSquareMatrix<scalar, Size>& matrix)
    :
    FixedSquareMatrix<scalar, Size>(matrix)
{

}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::FixedLUscalarMatrix>
Foam::FixedLUscalarMatrix::New()
{
    return autoPtr<FixedLUscalarMatrix>(new FixedLUscalarMatrix);
}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
/*
Foam::FixedLUscalarMatrix::~FixedLUscalarMatrix()
{}
*/

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //
/*
void Foam::FixedLUscalarMatrix::operator=(const FixedLUscalarMatrix& rhs)
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

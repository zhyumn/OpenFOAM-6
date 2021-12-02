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

Class
    Foam::FixedSquareMatrix

Description

SourceFiles
    FixedSquareMatrixI.H
    FixedSquareMatrix.C
    FixedSquareMatrixIO.C

\*---------------------------------------------------------------------------*/

#ifndef FixedSquareMatrix_H
#define FixedSquareMatrix_H
#include "bool.H"
#include "label.H"
#include "uLabel.H"
#include <type_traits>
#include <initializer_list>
#include "scalar.H"
//#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of classes
    //class Istream;
    //class Ostream;

    // Forward declaration of friend functions and operators
    //class FixedSquareMatrix;
    //Istream& operator>>(Istream&, FixedSquareMatrix&);
    //Ostream& operator<<(Ostream&, const FixedSquareMatrix&);


    /*---------------------------------------------------------------------------*\
                             Class FixedSquareMatrix Declaration
    \*---------------------------------------------------------------------------*/

    template<class Type, unsigned Size>
    class FixedSquareMatrix
    {
        static_assert
            (
                Size&& Size <= INT_MAX,
                "Size must be positive (non-zero) and also fit as a signed value"
                );
        // Private data

            //- Description of data_
        Type v_[Size * Size];


        // Private Member Functions

            //- Disallow default bitwise copy construct
        //FixedSquareMatrix(const FixedSquareMatrix&);

        //- Disallow default bitwise assignment
        //void operator=(const FixedSquareMatrix&);


    public:

        // Static data members

            //- Static data staticData
            //static const dataType staticData;


        // Constructors

            //- Construct null
        FixedSquareMatrix();

        //- Construct from components
        //FixedSquareMatrix(const dataType& data);

        //- Construct from Istream
        //FixedSquareMatrix(Istream&);

        //- Construct as copy
        FixedSquareMatrix(const FixedSquareMatrix&);


    // Selectors

        //- Select null constructed
        //static autoPtr<FixedSquareMatrix> New();


    //- Destructor
        ~FixedSquareMatrix();


        // Member Functions

            // Access

            // Check
        inline void check(const label i) const;
        // Edit

        // Write


    // Member Operators

        //void operator=(const FixedSquareMatrix&);

        inline Type& operator()(const label i, const label j);

        inline const Type& operator()(const label i, const label j) const;

        // Friend Functions

        // Friend Operators

        // IOstream Operators

        //friend Istream& operator>>(Istream&, FixedSquareMatrix&);
        //friend Ostream& operator<<(Ostream&, const FixedSquareMatrix&);
    };

    // * * * * * * * * * * * * * Explicitly Instantiate * * * * * * * * * * * * //

    template class Foam::FixedSquareMatrix<Foam::scalar, 10>;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FixedSquareMatrixI.H"

/*
#ifdef NoRepository
#include "FixedSquareMatrix.C"
#endif
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

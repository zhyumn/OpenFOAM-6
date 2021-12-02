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
    Foam::FixedLUscalarMatrix

Description

SourceFiles
    FixedLUscalarMatrixI.H
    FixedLUscalarMatrix.C
    FixedLUscalarMatrixIO.C

\*---------------------------------------------------------------------------*/

#ifndef FixedLUscalarMatrix_H
#define FixedLUscalarMatrix_H

#include "FixedSquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Forward declaration of classes
    //class Istream;
    //class Ostream;

    // Forward declaration of friend functions and operators
    //class FixedLUscalarMatrix;
    //Istream& operator>>(Istream&, FixedLUscalarMatrix&);
    //Ostream& operator<<(Ostream&, const FixedLUscalarMatrix&);


    /*---------------------------------------------------------------------------*\
                             Class FixedLUscalarMatrix Declaration
    \*---------------------------------------------------------------------------*/
    template<unsigned Size>
    class FixedLUscalarMatrix
        :
        public FixedSquareMatrix<scalar, Size>
    {
        // Private data

            //- Description of data_
        //dataType data_;


        // Private Member Functions

            //- Disallow default bitwise copy construct
        //FixedLUscalarMatrix(const FixedLUscalarMatrix&);

        //- Disallow default bitwise assignment
        //void operator=(const FixedLUscalarMatrix&);


    public:

        // Static data members

            //- Static data staticData
        //static const dataType staticData;


        // Constructors

            //- Construct null
        //FixedLUscalarMatrix();
        FixedLUscalarMatrix(const FixedSquareMatrix<scalar, Size>& M);

        //- Construct from components
        //FixedLUscalarMatrix(const dataType& data);

        //- Construct from Istream
        //FixedLUscalarMatrix(Istream&);

        //- Construct as copy
        //FixedLUscalarMatrix(const FixedLUscalarMatrix&);


        // Selectors

            //- Select null constructed
        //static autoPtr<FixedLUscalarMatrix> New();


        //- Destructor
        ~FixedLUscalarMatrix();


        // Member Functions

            // Access

            // Check

            // Edit

            // Write


        // Member Operators

        //void operator=(const FixedLUscalarMatrix&);


        // Friend Functions

        // Friend Operators

        // IOstream Operators

        //friend Istream& operator>>(Istream&, FixedLUscalarMatrix&);
        //friend Ostream& operator<<(Ostream&, const FixedLUscalarMatrix&);
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "FixedLUscalarMatrixI.H"
/*
#ifdef NoRepository
    #include "FixedLUscalarMatrix.C"
#endif
*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

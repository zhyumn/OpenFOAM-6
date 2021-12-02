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
    Foam::parallelISAT_1D

Description

SourceFiles
    parallelISAT_1DI.H
    parallelISAT_1D.C
    parallelISAT_1DIO.C

\*---------------------------------------------------------------------------*/

#ifndef parallelISAT_1D_H
#define parallelISAT_1D_H

#include "parallelISAT.H"
#include <iostream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


    class parallelISAT_1D;
    /*---------------------------------------------------------------------------*\
                             Class parallelISAT_1D Declaration
    \*---------------------------------------------------------------------------*/
    struct ISAT_1D
    {
        typedef scalar inputType;
        typedef scalar outputType;
        typedef scalar gradientType;
        static parallelISAT_1D* pISAT;
        class nodeData;
        struct leafData
        {
        private:
            inputType v;
            outputType Rv;
            scalar A;//gradient
            scalar EOA;
        public:
            //void set(const inputType& x) { v = x; Rv = x * x; gradient = 2 * x; EOA = pISAT->tolerance_ / max(fabs(gradient), 1); }
            void set(const inputType& x);
            bool checkSolution(const inputType& x, const outputType& y);
            void set(const inputType& x, const outputType& y);
            template<typename ...Args>
            static void computeA(const inputType& x, const outputType& y, gradientType& A, Args ...args);
            bool grow(const inputType& x) { EOA = fabs(x - v); return true; }
            static outputType func(const inputType& x);
            static gradientType gradFunc(const inputType& x);
            bool collide(const inputType& x) { return x == v; }
            void retrieve(const inputType& x, outputType& y) { y = Rv + A * (x - v); }
            bool inEOA(const inputType& x) { return fabs(x - v) < EOA; }
            friend nodeData;
            friend std::ostream& operator<<(std::ostream& out, ISAT_1D::leafData& A);
        };
        struct nodeData
        {
        private:
            scalar v;
        public:
            void set() { v = 0; }
            void set(const inputType& x, const inputType& y) { v = (x + y) / 2; }
            void set(const leafData& x, const leafData& y) { v = (x.v + y.v) / 2; }
            bool goLeft(inputType& x) { return x < v; }
            bool goLeft(leafData& x) { return x.v < v; }
            friend leafData;
            friend std::ostream& operator<<(std::ostream& out, ISAT_1D::nodeData& A);
        };
    };
    std::ostream& operator<<(std::ostream& out, ISAT_1D::leafData& A);
    std::ostream& operator<<(std::ostream& out, ISAT_1D::nodeData& A);
    struct emptyClass
    {};
    class parallelISAT_1D
        :
        public parallelISAT<ISAT_1D, emptyClass>
    {
        // Private data

            //- Description of data_
        //dataType data_;


        // Private Member Functions

            //- Disallow default bitwise copy construct

        scalar tolerance_;
        //- Disallow default bitwise assignment



    public:

        // Static data members

            //- Static data staticData
        //static const dataType staticData;

        typedef ISAT_1D DataType;
        // Constructors

            //- Construct null


        //- Construct from components
        parallelISAT_1D();

        //- Construct from Istream
        //parallelISAT_1D(Istream&);

        //- Construct as copy
        //parallelISAT_1D(const parallelISAT_1D&);


        // Selectors

            //- Select null constructed
        //static autoPtr<parallelISAT_1D> New();


        //- Destructor
        ~parallelISAT_1D();


        // Member Functions

            // Access
        const scalar& tolerance() { return tolerance_; }
        // Check

        // Edit
    //void computeA(scalarSquareMatrix& A, const calarField& Rphiq, const scalar rhoi, const scalar dt);
        // Write


    // Member Operators

    //void operator=(const parallelISAT_1D&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

    //friend Istream& operator>>(Istream&, parallelISAT_1D&);
    //friend Ostream& operator<<(Ostream&, const parallelISAT_1D&);
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#include "parallelISAT_1D.C"
//#endif
//#include "parallelISAT_1DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
